/* Copyright 2022 Peter Dimov.
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_DETAIL_MULX_HPP
#define BOOST_BLOOM_DETAIL_MULX_HPP

#include <boost/bloom/detail/64bit_arch.hpp>
#include <boost/cstdint.hpp>
#include <climits>
#include <cstddef>

#if defined(_MSC_VER)&&!defined(__clang__)
#include <intrin.h>
#endif

namespace boost{
namespace bloom{
namespace detail{

#if defined(_MSC_VER)&&defined(_M_X64)&&!defined(__clang__)

__forceinline boost::uint64_t mulx64(
  boost::uint64_t x,boost::uint64_t y,boost::uint64_t& hi)
{
  return _umul128(x,y,&hi);
}

#elif defined(_MSC_VER)&&defined(_M_ARM64)&&!defined(__clang__)

__forceinline boost::uint64_t mulx64(
  boost::uint64_t x,boost::uint64_t y,boost::uint64_t& hi)
{
  hi=__umulh(x,y);
  return x*y;
}

#elif defined(__SIZEOF_INT128__)

inline boost::uint64_t mulx64(
  boost::uint64_t x,boost::uint64_t y,boost::uint64_t& hi)
{
  __uint128_t r=(__uint128_t)x*y;
  hi=(boost::uint64_t)(r>>64);
  return (boost::uint64_t)r;
}

#else

inline boost::uint64_t mulx64(
  boost::uint64_t x,boost::uint64_t y,boost::uint64_t& hi)
{
  boost::uint64_t x1=(boost::uint32_t)x;
  boost::uint64_t x2=x >> 32;

  boost::uint64_t y1=(boost::uint32_t)y;
  boost::uint64_t y2=y >> 32;

  boost::uint64_t r3=x2*y2;

  boost::uint64_t r2a=x1*y2;

  r3+=r2a>>32;

  boost::uint64_t r2b=x2*y1;

  r3+=r2b>>32;

  boost::uint64_t r1=x1*y1;

  boost::uint64_t r2=(r1>>32)+(boost::uint32_t)r2a+(boost::uint32_t)r2b;

  r1=(r2<<32)+(boost::uint32_t)r1;
  r3+=r2>>32;

  hi=r3;
  return r1;
}

#endif

inline boost::uint32_t mulx32(
  boost::uint32_t x,boost::uint32_t y,boost::uint32_t& hi)
{
  boost::uint64_t r=(boost::uint64_t)x*y;
  hi=(boost::uint32_t)(r>>32);
#if defined(__MSVC_RUNTIME_CHECKS)
  return (boost::uint32_t)(r&UINT32_MAX);
#else
  return (boost::uint32_t)r;
#endif
}

inline std::size_t mulx(std::size_t x,std::size_t y,std::size_t& hi)noexcept
{
#if defined(BOOST_BLOOM_64B_ARCHITECTURE)
  boost::uint64_t hi_;
  boost::uint64_t r=mulx64((boost::uint64_t)x,(boost::uint64_t)y,hi_);
  hi=(std::size_t)hi_;
  return (std::size_t)r;
#else /* 32 bits assumed */
  return mulx32(x,y,hi);
#endif
}

inline std::size_t mulx_mix(std::size_t x)noexcept
{
#if defined(BOOST_BLOOM_64B_ARCHITECTURE)
  /* multiplier is phi */
  boost::uint64_t hi;
  boost::uint64_t lo=mulx64((boost::uint64_t)x,0x9E3779B97F4A7C15ull,hi);
  return (std::size_t)(hi^lo);
#else /* 32 bits assumed */
  /* multiplier from https://arxiv.org/abs/2001.05304 */
  boost::uint32_t hi;
  boost::uint32_t lo=mulx32(x,0xE817FB2Du,hi);
  return (std::size_t)(hi^lo);
#endif
}

} /* namespace detail */
} /* namespace bloom */
} /* namespace boost */
#endif
