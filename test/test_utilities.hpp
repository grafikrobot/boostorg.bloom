/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_TEST_TEST_UTILITIES_HPP
#define BOOST_BLOOM_TEST_TEST_UTILITIES_HPP

#include <string>

namespace test_utilities{

template<typename T>
struct value_factory
{
  T operator()(){return n++;}
  T n=0;
};

template<>
struct value_factory<std::string>
{
  std::string operator()()
  {
    return std::to_string(n++);
  }

  int n=0;
};

} /* namespace test_utilities */
#endif
