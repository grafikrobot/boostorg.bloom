/* Common base for all boost::bloom::filter instantiations.
 *
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_FILTER_DETAIL_CORE_HPP
#define BOOST_BLOOM_FILTER_DETAIL_CORE_HPP

#include <boost/bloom/detail/mulx.hpp>
#include <boost/bloom/detail/sse2.hpp>
#include <boost/config.hpp>
#include <boost/core/bit.hpp>
#include <boost/core/empty_value.hpp>
#include <boost/core/allocator_traits.hpp>
#include <cstring>
#include <memory>
#include <type_traits>

/* We use BOOST_BLOOM_PREFETCH[_WRITE] macros rather than proper
 * functions because of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=109985
 */

#if defined(BOOST_GCC)||defined(BOOST_CLANG)
#define BOOST_BLOOM_PREFETCH(p) __builtin_prefetch((const char*)(p))
#define BOOST_BLOOM_PREFETCH_WRITE(p) __builtin_prefetch((const char*)(p),1)
#elif defined(BOOST_BLOOM_SSE2)
#define BOOST_BLOOM_PREFETCH(p) _mm_prefetch((const char*)(p),_MM_HINT_T0)
#if defined(_MM_HINT_ET0)
#define BOOST_BLOOM_PREFETCH_WRITE(p) \
_mm_prefetch((const char*)(p),_MM_HINT_ET0)
#else
#define BOOST_BLOOM_PREFETCH_WRITE(p) \
_mm_prefetch((const char*)(p),_MM_HINT_T0)
#endif
#else
#define BOOST_BLOOM_PREFETCH(p) ((void)(p))
#define BOOST_BLOOM_PREFETCH_WRITE(p) ((void)(p))
#endif

namespace boost{
namespace bloom{
namespace detail{

#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable:4714) /* marked as __forceinline not inlined */
#endif

/*  mcg_and_fastrange produces (pos,hash') from hash, where
 *   - x=mulx(hash,range), mulx denotes extended multiplication
 *   - pos=high(x)
 *   - hash'=low(x)
 *  pos is uniformly distributed in [0,range) (see
 *  https://arxiv.org/pdf/1805.10941), whereas hash'<-hash is a multiplicative
 *  congruential generator of the form hash'<-hash*rng mod 2^n, where n is the
 *  size in bits of std::size_t. This MCG generates long cycles when the
 *  initial value of hash is odd and rng = +-3 (mod 8), which is why we adjust
 *  hash and rng as seen below. As a result, the low bits of hash' are of poor
 *  quality, and the least significant bit in particular is always one.
 */

struct mcg_and_fastrange
{
  mcg_and_fastrange(std::size_t m)noexcept:
    rng{
      m+(
        (m%8<=3)?3-(m%8):
        (m%8<=5)?5-(m%8):
                 8-(m%8)+3)
    }
    {}

  inline std::size_t range()const noexcept{return rng;}

  inline void prepare_hash(std::size_t& hash)const noexcept
  {
    hash|=1u;
  }

  inline std::size_t next_position(std::size_t& hash)const noexcept
  {
    std::size_t hi;
    hash=mulx(hash,rng,hi);
    return hi;
  }

  std::size_t rng;
};

/* used_block_size<Subfilter>::value is Subfilter::used_value_size if it
 * exists, or sizeof(Subfilter::value_type) otherwise. This covers the
 * case where a subfilter only operates on the first bytes of its entire
 * value_type (e.g. fast_multiblock32<K> with K<8).
 */

template<typename Subfilter,typename=void>
struct used_block_size
{
  static constexpr std::size_t value=sizeof(typename Subfilter::value_type);
};

template<typename Subfilter>
struct used_block_size<
  Subfilter,
  typename std::enable_if<Subfilter::used_value_size!=0>::type
>
{
  static constexpr std::size_t value=Subfilter::used_value_size;
};

/* GCD with x,p > 1, p a power of two */

inline constexpr std::size_t gcd_pow2(std::size_t x,std::size_t p)
{
  /* x&-x: maximum power of two dividing x */
  return (x&(0-x))<p?(x&(0-x)):p;
}

template<
  std::size_t K,typename Subfilter,std::size_t BucketSize,typename Allocator
>
class filter_core:empty_value<Allocator,0>
{
  static_assert(K>0,"K must be >= 1");
  static_assert(
    std::is_same<allocator_value_type_t<Allocator>,unsigned char>::value,
    "Allocator value_type must be unsigned char");

public:
  static constexpr std::size_t k=K;
  using subfilter=Subfilter;

private:
  using block_type=typename subfilter::value_type;
  static constexpr std::size_t used_block_size=
    detail::used_block_size<subfilter>::value;
  static constexpr std::size_t bucket_size=
    BucketSize?BucketSize:used_block_size;
  static constexpr std::size_t tail_size=sizeof(block_type)-bucket_size;
  static constexpr bool are_blocks_aligned=
    (bucket_size%alignof(block_type)==0);
  static constexpr std::size_t cacheline=64; /* unknown at compile time */
  static constexpr std::size_t initial_alignment=
    are_blocks_aligned?
      alignof(block_type)>cacheline?alignof(block_type):cacheline:
      1;
  static constexpr std::size_t prefetched_cachelines=
    1+(used_block_size+cacheline-1-
       detail::gcd_pow2(bucket_size,cacheline))/cacheline;
  using hash_strategy=detail::mcg_and_fastrange;

protected:
  using allocator_type=Allocator;
  using size_type=std::size_t;
  using difference_type=std::ptrdiff_t;
  using pointer=allocator_pointer_t<allocator_type>;
  using const_pointer=allocator_const_pointer_t<allocator_type>;

  filter_core(std::size_t m,const allocator_type& al_):
    allocator_base{empty_init,al_},
    hs{((m+CHAR_BIT-1)/CHAR_BIT+bucket_size-1)/bucket_size}
  {
    std::size_t spc=space_for(hs.range());
    data_=allocator_allocate(al(),spc);
    std::memset(data_,0,spc);
    buckets=buckets_for(data_);
  }

  ~filter_core()
  {
    allocator_deallocate(al(),data_,space_for(hs.range()));
  }

  std::size_t capacity()const noexcept
  {
    return hs.range()*CHAR_BIT;
  }

  BOOST_FORCEINLINE void insert(std::size_t hash)
  {
    hs.prepare_hash(hash);
    for(auto n=k;n--;){
      auto p=next_element(hash); /* modifies hash */
      set(p,hash);
    }
  }

  BOOST_FORCEINLINE bool may_contain(std::size_t hash)const
  {
    hs.prepare_hash(hash);
#if 1
    auto p0=next_element(hash);
    for(std::size_t n=k-1;n--;){
      auto p=p0;
      auto hash0=hash;
      p0=next_element(hash);
      if(!get(p,hash0))return false;
    }
    if(!get(p0,hash))return false;
    return true;
#else
    for(auto n=k;n--;){
      auto p=next_element(hash); /* modifies hash */
      if(!get(p,hash))return false;
    }
    return true;
#endif
  }

private:
  using allocator_base=empty_value<Allocator,0>;

  Allocator& al(){return allocator_base::get();}

  static constexpr std::size_t space_for(std::size_t rng)
  {
    return (initial_alignment-1)+rng*bucket_size+tail_size;
  }

  static unsigned char* buckets_for(unsigned char* p)
  {
    return p+
      (boost::uintptr_t(initial_alignment)-
       boost::uintptr_t(p))%initial_alignment;
  }

  BOOST_FORCEINLINE bool get(const unsigned char* p,std::size_t hash)const
  {
    return get(p,hash,std::integral_constant<bool,are_blocks_aligned>{});
  }

  BOOST_FORCEINLINE bool get(
    const unsigned char* p,std::size_t hash,
    std::true_type /* blocks aligned */)const
  {
    return subfilter::check(*reinterpret_cast<const block_type*>(p),hash);
  }

  BOOST_FORCEINLINE bool get(
    const unsigned char* p,std::size_t hash,
    std::false_type /* blocks not aligned */)const
  {
    block_type x;
    std::memcpy(&x,p,used_block_size);
    return subfilter::check(x,hash);
  }

  BOOST_FORCEINLINE void set(unsigned char* p,std::size_t hash)
  {
    return set(p,hash,std::integral_constant<bool,are_blocks_aligned>{});
  }

  BOOST_FORCEINLINE void set(
    unsigned char* p,std::size_t hash,
    std::true_type /* blocks aligned */)const
  {
    subfilter::mark(*reinterpret_cast<block_type*>(p),hash);
  }

  BOOST_FORCEINLINE void set(
    unsigned char* p,std::size_t hash,
    std::false_type /* blocks not aligned */)const
  {
    block_type x;
    std::memcpy(&x,p,used_block_size);
    subfilter::mark(x,hash);
    std::memcpy(p,&x,used_block_size);
  }

  BOOST_FORCEINLINE unsigned char* next_element(std::size_t& hash)
  {
    auto p=buckets+hs.next_position(hash)*bucket_size;
    for(std::size_t i=0;i<prefetched_cachelines;++i){
      BOOST_BLOOM_PREFETCH_WRITE((unsigned char*)p+i*cacheline);
    }
    return p;
  }

  BOOST_FORCEINLINE const unsigned char* next_element(std::size_t& hash)const
  {
    auto p=buckets+hs.next_position(hash)*bucket_size;
    for(std::size_t i=0;i<prefetched_cachelines;++i){
      BOOST_BLOOM_PREFETCH((unsigned char*)p+i*cacheline);
    }
    return p;
  }

  hash_strategy hs;
  pointer       data_,buckets;
};

#if defined(BOOST_MSVC)
#pragma warning(pop) /* C4714 */
#endif

} /* namespace detail */
} /* namespace bloom */
} /* namespace boost */
#endif
