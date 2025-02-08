/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_MULTIBLOCK_HPP
#define BOOST_BLOOM_MULTIBLOCK_HPP

#include <boost/bloom/detail/block_base.hpp>
#include <cstddef>

namespace boost{
namespace bloom{

template<typename Block,std::size_t K>
struct multiblock:private detail::block_base<Block,K>
{
  static constexpr std::size_t k=K;
  using value_type=Block[k];

  static inline void mark(value_type& x,std::size_t hash)
  {
    auto h=hash;
    for(std::size_t i=0;i<k;++i){
      next(i,h,hash);
      x[i]|=Block(1)<<(h&mask);
    }
  }

  static inline bool check(const value_type& x,std::size_t hash)
  {
    Block res=1;
    auto h=hash;
    for(std::size_t i=0;i<k;++i){
      next(i,h,hash);
      res&=(x[i]>>(h&mask));
    }
    return res;
  }

private:
  using super=detail::block_base<Block,K>;
  using super::mask;
  using super::next;
};

} /* namespace bloom */
} /* namespace boost */
#endif
