/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_TEST_TEST_UTILITIES_HPP
#define BOOST_BLOOM_TEST_TEST_UTILITIES_HPP

#include <boost/bloom/filter.hpp>
#include <boost/core/allocator_traits.hpp>
#include <boost/core/lightweight_test.hpp>
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

template<typename Filter,typename U>
struct revalue_filter_impl;

template<
  typename T,typename H,std::size_t K,typename S,std::size_t B,typename A,
  typename U
>
struct revalue_filter_impl<boost::bloom::filter<T,H,K,S,B,A>,U>
{
  using type=boost::bloom::filter<U,H,K,S,B,boost::allocator_rebind_t<A,U>>;
};

template<typename Filter,typename U>
using revalue_filter=typename revalue_filter_impl<Filter,U>::type;

template<typename Filter,typename Hash>
struct rehash_filter_impl;

template<
  typename T,typename H,std::size_t K,typename S,std::size_t B,typename A,
  typename Hash
>
struct rehash_filter_impl<boost::bloom::filter<T,H,K,S,B,A>,Hash>
{
  using type=boost::bloom::filter<T,Hash,K,S,B,A>;
};

template<typename Filter,typename Hash>
using rehash_filter=typename rehash_filter_impl<Filter,Hash>::type;

template<typename Filter,typename Allocator>
struct realloc_filter_impl;

template<
  typename T,typename H,std::size_t K,typename S,std::size_t B,typename A,
  typename Allocator
>
struct realloc_filter_impl<boost::bloom::filter<T,H,K,S,B,A>,Allocator>
{
  using type=boost::bloom::filter<T,H,K,S,B,Allocator>;
};

template<typename Filter,typename Allocator>
using realloc_filter=typename realloc_filter_impl<Filter,Allocator>::type;

template<typename Filter,typename Input>
void check_may_contain(const Filter& f,const Input& input)
{
  std::size_t res=0;
  for(const auto& x:input)res+=f.may_contain(x);
  BOOST_TEST_EQ(res,input.size());
}

} /* namespace test_utilities */
#endif
