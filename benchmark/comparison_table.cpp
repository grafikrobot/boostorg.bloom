/* Comparison table for several configurations of boost::bloom::filter.
 * 
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <numeric>

std::chrono::high_resolution_clock::time_point measure_start,measure_pause;

template<typename F>
double measure(F f)
{
  using namespace std::chrono;

  static const int              num_trials=10;
  static const milliseconds     min_time_per_trial(10);
  std::array<double,num_trials> trials;

  for(int i=0;i<num_trials;++i){
    int                               runs=0;
    high_resolution_clock::time_point t2;
    volatile decltype(f())            res; /* to avoid optimizing f() away */

    measure_start=high_resolution_clock::now();
    do{
      res=f();
      ++runs;
      t2=high_resolution_clock::now();
    }while(t2-measure_start<min_time_per_trial);
    trials[i]=duration_cast<duration<double>>(t2-measure_start).count()/runs;
  }

  std::sort(trials.begin(),trials.end());
  return std::accumulate(
    trials.begin()+2,trials.end()-2,0.0)/(trials.size()-4);
}

void pause_timing()
{
  measure_pause=std::chrono::high_resolution_clock::now();
}

void resume_timing()
{
  measure_start+=std::chrono::high_resolution_clock::now()-measure_pause;
}

#include <boost/bloom/block.hpp>
#include <boost/bloom/filter.hpp>
#include <boost/bloom/multiblock.hpp>
#include <boost/core/detail/splitmix64.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

template<typename T>
struct unordered_flat_set_filter
{
  using value_type=T;

  unordered_flat_set_filter(std::size_t){}
  void insert(const T& x){s.insert(x);}
  bool may_contain(const T& x){return s.contains(x);}

  boost::unordered_flat_set<T> s;
};

struct test_results
{
  double fpr;                      /* % */
  double insertion_time;           /* ns per element */
  double successful_lookup_time;   /* ns per element */
  double unsuccessful_lookup_time; /* ns per element */
};

template<typename Filter>
test_results test(std::size_t c)
{
  using value_type=typename Filter::value_type;
  static constexpr std::size_t N=10000000;

  std::vector<value_type> data_in,data_out;
  {
    boost::detail::splitmix64             rng;
    boost::unordered_flat_set<value_type> unique;
    for(std::size_t i=0;i<N;++i){
      for(;;){
        auto x=value_type(rng());
        if(unique.insert(x).second){
          data_in.push_back(x);
          break;
        }
      }
    }
    for(std::size_t i=0;i<N;++i){
      for(;;){
        auto x=value_type(rng());
        if(!unique.contains(x)){
          data_out.push_back(x);
          break;
        }
      }
    }
  }

  double fpr=0.0;
  {
    std::size_t res=0;
    Filter f(c*N);
    for(const auto& x:data_in)f.insert(x);
    for(const auto& x:data_out)res+=f.may_contain(x);
    fpr=(double)res*100/N;
  }

  double insertion_time=0.0;
  {
    double t=measure([&]{
      pause_timing();
      {
        Filter f(c*N);
        resume_timing();
        for(const auto& x:data_in)f.insert(x);
        pause_timing();
      }
      resume_timing();
      return 0;
    });
    insertion_time=t/N*1E9;
  }

  double successful_lookup_time=0.0;
  double unsuccessful_lookup_time=0.0;
  {
    Filter f(c*N);
    for(const auto& x:data_in)f.insert(x);
    double t=measure([&]{
      std::size_t res=0;
      for(const auto& x:data_in)res+=f.may_contain(x);
      return 0;
    });
    successful_lookup_time=t/N*1E9;
    t=measure([&]{
      std::size_t res=0;
      for(const auto& x:data_out)res+=f.may_contain(x);
      return 0;
    });
    unsuccessful_lookup_time=t/N*1E9;
  }

  return {fpr,insertion_time,successful_lookup_time,unsuccessful_lookup_time};
}

struct print_double
{
  print_double(double x_,int precision_=2):x{x_},precision{precision_}{}

  friend std::ostream& operator<<(std::ostream& os,const print_double& pd)
  {
    const auto default_precision{std::cout.precision()};
    os<<std::fixed<<std::setprecision(pd.precision)<<pd.x;
    std::cout.unsetf(std::ios::fixed);
    os<<std::setprecision(default_precision);
    return os;
  }

  double x;
  int    precision;
};

template<std::size_t K> void row(std::size_t c)
{
  using namespace boost::bloom;
  using filters=boost::mp11::mp_list<
    boost::bloom::filter<
      int,boost::hash<int>,K
    >,
    boost::bloom::filter<
      int,boost::hash<int>,1,
      block<boost::uint64_t,K>
    >,
    boost::bloom::filter<
      int,boost::hash<int>,1,
      multiblock<boost::uint64_t,K>
    >
  >;

  std::cout<<
    "  <tr>\n"
    "    <td align=\"center\">"<<c<<"</td>\n"
    "    <td align=\"center\">"<<K<<"</td>\n";

  boost::mp11::mp_for_each<
    boost::mp11::mp_transform<boost::mp11::mp_identity,filters>
  >([&](auto i){
    using filter=typename decltype(i)::type;
    auto res=test<filter>(c);
    std::cout<<
      "    <td align=\"right\">"<<print_double(res.fpr,4)<<"</td>\n"
      "    <td align=\"right\">"<<print_double(res.insertion_time)<<"</td>\n"
      "    <td align=\"right\">"<<print_double(res.successful_lookup_time)<<"</td>\n"
      "    <td align=\"right\">"<<print_double(res.unsuccessful_lookup_time)<<"</td>\n";
  });

  std::cout<<
    "  </tr>\n";
}

int main()
{
  /* reference table: boost::unordered_flat_set */

  auto res=test<unordered_flat_set_filter<int>>(0);
  std::cout<<
    "<table>\n"
    "  <tr><th colspan=\"3\"><code>boost::unordered_flat_set</code></tr>\n"
    "  <tr>\n"
    "    <th>insertion</th>\n"
    "    <th>successful</br>lookup</th>\n"
    "    <th>unsuccessful</br>lookup</th>\n"
    "  </tr>\n"
    "  <tr>\n"
    "    <td align=\"right\">"<<print_double(res.insertion_time)<<"</td>\n"
    "    <td align=\"right\">"<<print_double(res.successful_lookup_time)<<"</td>\n"
    "    <td align=\"right\">"<<print_double(res.unsuccessful_lookup_time)<<"</td>\n"
    "  </tr>\n"
    "</table>\n";

  /* filter table */

  auto subheader=
    "    <th>FPR [%]</th>\n"
    "    <th>ins.</th>\n"
    "    <th>succ.</br>lookup</th>\n"
    "    <th>unsucc.</br>lookup</th>\n";

  std::cout<<
    "<table>\n"
    "  <tr>\n"
    "    <th colspan=\"2\"></th>\n"
    "    <th colspan=\"4\"><code>filter&lt;K></code></th>\n"
    "    <th colspan=\"4\"><code>filter&lt;1, block&lt;uint64_t, K>></code></th>\n"
    "    <th colspan=\"4\"><code>filter&lt;1, multiblock&lt;uint64_t, K>></code></th>\n"
    "  </tr>\n"
    "  <tr>\n"
    "    <th>c</th>\n"
    "    <th>K</th>\n";
  std::cout<<
    subheader<<
    subheader<<
    subheader<<
    "  </tr>\n";

  row<6>(8);
  row<9>(12);
  row<11>(16);
  row<14>(20);

  std::cout<<"</table>\n";
}
