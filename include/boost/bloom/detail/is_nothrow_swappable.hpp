/* Copyright 2022-2023 Christian Mazakas.
 * Copyright 2024 Braden Ganetsky.
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_DETAIL_IS_NOTHROW_SWAPPABLE_HPP
#define BOOST_BLOOM_DETAIL_IS_NOTHROW_SWAPPABLE_HPP

#include <boost/type_traits/make_void.hpp>
#include <type_traits>
#include <utility>

namespace boost{
namespace bloom{
namespace detail{
namespace is_nothrow_swappable_helper_detail{

using std::swap;

template<typename T,typename=void>
struct is_nothrow_swappable_helper
{
  constexpr static bool value=false;
};

template <typename T>
struct is_nothrow_swappable_helper<
  T,
  boost::void_t<decltype(swap(std::declval<T&>(),std::declval<T&>()))>
>
{
  constexpr static bool value=
    noexcept(swap(std::declval<T&>(),std::declval<T&>()));
};

} /* namespace is_nothrow_swappable_helper_detail */

template <class T>
struct is_nothrow_swappable:std::integral_constant<
  bool,
  is_nothrow_swappable_helper_detail::is_nothrow_swappable_helper<T>::value
>{};

#define BOOST_BLOOM_STATIC_ASSERT_IS_NOTHROW_SWAPPABLE(T) \
static_assert(                                            \
  boost::bloom::detail::is_nothrow_swappable<T>::value,   \
  #T " must be nothrow swappable")

} /* namespace detail */
} /* namespace bloom */
} /* namespace boost */

#endif
