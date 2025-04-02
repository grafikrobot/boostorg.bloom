# Candidate Boost Bloom Library

[![Branch](https://img.shields.io/badge/branch-develop-brightgreen.svg)](https://github.com/joaquintides/bloom/tree/develop) [![CI](https://github.com/joaquintides/bloom/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/joaquintides/bloom/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/joaquintides/bloom/develop?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/joaquintides/bloom) </br>
[![BSL 1.0](https://img.shields.io/badge/license-BSL_1.0-blue.svg)](https://www.boost.org/users/license.html) <img alt="C++11 required" src="https://img.shields.io/badge/standard-C%2b%2b11-blue.svg"> <img alt="Header-only library" src="https://img.shields.io/badge/build-header--only-blue.svg">

(Candidate) Boost.Bloom provides the class template `boost::bloom::filter` that
can be configured to implement a classical [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter)
as well as variations discussed in the literature such as _blocked_ filters, _split block_/_multi-block_
filters, and more.

## Example

```cpp
#include <boost/bloom/filter.hpp>
#include <cassert>
#include <string>

int main()
{
  // Bloom filter of strings with 5 bits set per insertion
  using filter = boost::bloom::filter<std::string, 5>;

  // create filter with a capacity of 1'000'000 **bits**
  filter f(1'000'000);

  // insert elements (they can't be erased, Bloom filters are insert-only)
  f.insert("hello");
  f.insert("Boost");
  //...

  // elements inserted are always correctly checked as such
  assert(f.may_contain("hello") == true);

  // elements not inserted may incorrectly be identified as such with a
  // false positive rate (FPR) which is a function of the array capacity,
  // the number of bits set per element and generally how the boost::bloom::filter
  // was specified
  if(f.may_contain("bye")) { // likely false
    //...
  }
}
```

## Filter definition

A `boost::bloom::filter` can be regarded as an array of _buckets_ selected pseudo-randomly
(based on a hash function) upon insertion: each of the buckets is passed to a _subfilter_
that marks one or more of its bits according to some associated strategy.

```cpp
template<
  typename T, std::size_t K,
  typename Subfilter = block<unsigned char, 1>, std::size_t BucketSize = 0,
  typename Hash = boost::hash<T>, typename Allocator = std::allocator<T>  
>
class filter;
```

* `T`: type of the elements inserted.
* `K` number of buckets marked per insertion.
* `Subfilter`: type of subfilter used (more on this later).
* `BucketSize`: the number of buckets is just the capacity of the underlying
array (in bytes), divided by `BucketSize`. When `BucketSize` is specified as zero,
the value `sizeof(Subfilter::value_type)` is used.

The default configuration with `block<unsigned char,1>` corresponds to a
classical Bloom filter setting `K` bits per elements uniformly distributed across
the array.

### Overlapping buckets

`BucketSize` can be any value other (and typically less) than `sizeof(Subfilter::value_type)`.
When this is the case, subfilters in adjacent buckets act on overlapping byte ranges.
Far from being a problem, this improves the resulting _false positive rate_ (FPR) of
the filter. The downside is that buckets won't be in general properly aligned in memory,
which may result in more cache misses.

## Provided subfilters

### `block<Block, K'>`

Sets `K'` bits in an underlying value of the unsigned integral type `Block`
(e.g. `unsigned char`, `uint32_t`, `uint64_t`). So,
a `filter<T, K, block<Block, K'>>` will set `K*K'` bits per element.
The tradeoff here is that insertion/lookup will be (much) faster than
with `filter<T, K*K'>` while the FPR will be worse (larger).
FPR is better the wider `Block` is.

### `multiblock<Block, K'>`

Instead of setting `K'` bits in a `Block` value, this subfilter sets
one bit on each of the elements of a `Block[K']` subarray. This improves FPR
but impacts performance with respect to `block<Block, K'>`, among other
things because cacheline boundaries can be crossed when accessing the subarray.

### `fast_multiblock32<K'>`

Statistically equivalent to `multiblock<uint32_t, K'>`, but uses
faster SIMD-based algorithms when SSE2, AVX2 or Neon are available. The approach is
similar to that of [Apache Kudu](https://kudu.apache.org/)
[`BlockBloomFilter`](https://github.com/apache/kudu/blob/master/src/kudu/util/block_bloom_filter_avx2.cc),
but that implementation is fixed to `K'` = 8 whereas we accept any value.

### `fast_multiblock64<K'>`

Statistically equivalent to `multiblock<uint64_t, K'>`, but uses a
faster SIMD-based algorithm when AVX2 is available.

## Estimating FPR

For a classical Bloom filter, the theoretical false positive rate, under some simplifying assumptions,
is given by

$$FPR(n,m,k)=\left(1 - \left(1 - \frac{1}{m}\right)^{kn}\right)^k \approx \left(1 - e^{-kn/m}\right)^k \text{   for large } m,$$

where $$n$$ is the number of elements inserted in the filter, $$m$$ its capacity in bits and $$k$$ the
number of bits set per insertion (see a [derivation](https://en.wikipedia.org/wiki/Bloom_filter#Probability_of_false_positives)
of this formula). For a given inverse load factor $$c=m/n$$, the optimum $$k$$ is
the integer closest to:

$$k_{\text{opt}}=c\cdot\ln2,$$

yielding a minimum attainable FPR of $$1/2^{k_{\text{opt}}} \approx 0.6185^{c}$$.

In the case of a Boost.Bloom block filter of the form `filter<T, K, block<Block, K'>>`, we can extend
the approach from [Putze et al.](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=f376ff09a64b388bfcde2f5353e9ddb44033aac8)
to derive the (approximate but very precise) formula:

$$FPR_{\text{block}}(n,m,b,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}(i,nbk/m) \cdot FPR(i,b,k')\right)^{k},$$

where

$$\text{Pois}(i,\lambda)=\frac{\lambda^i e^{-\lambda}}{i!}$$

is the probability mass function of a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution)
with mean $$\lambda$$, and $$b$$ is the size of `Block` in bits. If we're using `multiblock<Block,K'>`, we have

$$FPR_\text{multiblock}(n,m,b,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}(i,nbkk'/m) \cdot FPR(i,b,1)^{k'}\right)^{k}.$$

As we have commented before, in general 

$$FPR_\text{block}(n,m,b,k,k') \geq FPR_\text{multiblock}(n,m,b,k,k') \geq FPR(n,m,kk'),$$

that is, block and multi-block filters have worse FPR than the classical filter for the same number of bits
set per insertion, but they will be much faster. We have the particular case

$$FPR_{\text{block}}(n,m,b,k,1)=FPR_{\text{multiblock}}(n,m,b,k,1)=FPR(n,m,k),$$

which follows simply from the observation that using `{block|multiblock}<Block, 1>` behaves exactly as
a classical Bloom filter.

We don't know of any closed, simple formula for the FPR of block and multiblock filters when
`Bucketsize` is not its "natural" size (`sizeof(Block)` for `block<Block, K'>`,
`sizeof(Block)*K'` for `multiblock<Block, K'>`), that is, when subfilter values overlap.
We can use the following approximations ($$s$$ = `BucketSize` in bits):

$$FPR_{\text{block}}(n,m,b,s,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}\left(i,\frac{n(2b-s)k}{m}\right) \cdot FPR(i,2b-s,k')\right)^{k},$$
$$FPR_\text{multiblock}(n,m,b,s,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}\left(i,\frac{n(2bk'-s)k}{m}\right) \cdot FPR\left(i,\frac{2bk'-s}{k'},1\right)^{k'}\right)^{k},$$

where the replacement of $$b$$ with $$2b-s$$ (or $$bk'$$ with $$2bk'-s$$ for multiblock filters) accounts
for the fact that the window of hashing positions affecting a particular bit spreads due to
overlapping. Note that the formulas reduce to the non-ovelapping case when $$s$$ takes its
default value ($$b$$ for block, $$bk´$$ for multiblock). These approximations are acceptable for
low values of $$k'$$ but tend to underestimate the actual FPR as $$k'$$ grows.
In general, the use of overlapping improves (decreases) FPR by a factor ranging from
0.6 to 0.9 for typical filter configurations.

## Experimental results

Provided in a [dedicated repo](https://github.com/joaquintides/boost_bloom_benchmarks).
</table>

