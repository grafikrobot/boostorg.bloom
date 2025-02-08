/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_DETAIL_64BIT_ARCH_HPP
#define BOOST_BLOOM_DETAIL_64BIT_ARCH_HPP

#include <climits>

#if defined(SIZE_MAX)
#if ((((SIZE_MAX >> 16) >> 16) >> 16) >> 15) != 0
#define BOOST_BLOOM_64B_ARCHITECTURE /* >64 bits assumed as 64 bits */
#endif
#elif defined(UINTPTR_MAX) /* used as proxy for std::size_t */
#if ((((UINTPTR_MAX >> 16) >> 16) >> 16) >> 15) != 0
#define BOOST_BLOOM_64B_ARCHITECTURE
#endif
#endif

#endif
