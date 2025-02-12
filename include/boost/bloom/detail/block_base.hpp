/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_BLOCK_BASE_HPP
#define BOOST_BLOOM_BLOCK_BASE_HPP

#include <boost/bloom/detail/mulx64.hpp>
#include <boost/core/bit.hpp>
#include <boost/cstdint.hpp>
#include <cstddef>

namespace boost{
namespace bloom{
namespace detail{

// TODO: describe

template<typename Block,std::size_t K>
struct block_base
{
  static constexpr std::size_t hash_width=sizeof(boost::uint64_t)*CHAR_BIT;
  static constexpr std::size_t block_width=sizeof(Block)*CHAR_BIT;
  static constexpr std::size_t mask=block_width-1;
  static constexpr std::size_t shift=boost::core::bit_width(mask);
  static constexpr std::size_t rehash_k=(hash_width-shift)/shift;

  static inline void next(
    std::size_t i,boost::uint64_t& h,boost::uint64_t& hash)
  {
    if(i&&(i%rehash_k)==0){
      h=hash=detail::mulx64_mix(hash);
    }
    else{
      h>>=shift;
    }
  }
};

} /* namespace detail */
} /* namespace bloom */
} /* namespace boost */
#endif
