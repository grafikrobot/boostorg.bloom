/* Configurable Bloom filter.
 * 
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_FILTER_HPP
#define BOOST_BLOOM_FILTER_HPP

#include <boost/bloom/block.hpp>
#include <boost/bloom/detail/core.hpp>
#include <boost/config.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/core/empty_value.hpp>
#include <boost/cstdint.hpp>
#include <boost/unordered/hash_traits.hpp> // TODO: internalize?

namespace boost{
namespace bloom{
namespace detail{

/* Mixing policies: no_mix_policy is the identity function, and
 * mulx_mix_policy uses the mulx_mix function from
 * <boost/bloom/detail/mulx.hpp>.
 *
 * filter mixes hash results with mulx64_mix if the hash is not marked as
 * avalanching, i.e. it's not of good quality (see
 * <boost/unordered/hash_traits.hpp>), or if std::size_t is less than 64 bits
 * (mixing policies promote to boost::uint64_t).
 */

struct no_mix_policy
{
  template<typename Hash,typename T>
  static inline boost::uint64_t mix(const Hash& h,const T& x)
  {
    return (boost::uint64_t)h(x);
  }
};

struct mulx64_mix_policy
{
  template<typename Hash,typename T>
  static inline boost::uint64_t mix(const Hash& h,const T& x)
  {
    return mulx64_mix((boost::uint64_t)h(x));
  }
};

} /* namespace detail */

#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable:4714) /* marked as __forceinline not inlined */
#endif

template<
  typename T,typename Hash,std::size_t K,
  typename Subfilter=block<unsigned char,1>,std::size_t BucketSize=0,
  typename Allocator=std::allocator<unsigned char>
>
class

#if defined(_MSC_VER)&&_MSC_FULL_VER>=190023918
__declspec(empty_bases) /* activate EBO with multiple inheritance */
#endif

filter:
  public detail::filter_core<K,Subfilter,BucketSize,Allocator>, //TODO: revert to private
  empty_value<Hash,0>
{
  using super=detail::filter_core<K,Subfilter,BucketSize,Allocator>;
  using mix_policy=typename std::conditional<
    unordered::hash_is_avalanching<Hash>::value&&
    sizeof(std::size_t)>=sizeof(boost::uint64_t),
    detail::no_mix_policy,
    detail::mulx64_mix_policy
  >::type;

public:
  using value_type=T;
  using super::k;
  using subfilter=typename super::subfilter;
  using allocator_type=typename super::allocator_type;
  using size_type=typename super::size_type;
  using difference_type=typename super::difference_type;
  using reference=value_type&;
  using const_reference=const value_type&;
  using pointer=typename super::pointer;
  using const_pointer=typename super::const_pointer;

  filter(std::size_t m,const allocator_type& al={}):super{m,al}{}

  using super::capacity;

  BOOST_FORCEINLINE void insert(const T& x)
  {
    super::insert(hash_for(x));
  }

  BOOST_FORCEINLINE bool may_contain(const T& x)const
  {
    return super::may_contain(hash_for(x));
  }

private:
  using hash_base=empty_value<Hash,0>;

  const Hash& h()const{return hash_base::get();}

  inline boost::uint64_t hash_for(const T& x)const
  {
    return mix_policy::mix(h(),x);
  }
};

#if defined(BOOST_MSVC)
#pragma warning(pop) /* C4714 */
#endif

} /* namespace bloom */
} /* namespace boost */
#endif
