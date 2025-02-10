# Candidate Boost Bloom Library

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
  using filter = boost::bloom::filter<std::string, boost::hash<std::string>, 5>;

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
  typename T, typename Hash, std::size_t K,
  typename Subfilter = block<unsigned char, 1>, std::size_t BucketSize = 0,
  typename Allocator = std::allocator<unsigned char>  
>
class filter;
```

* `T`: type of the elements inserted.
* `Hash`: a hash function for `T`.
* `K` number of buckets marked per insertion.
* `Subfilter`: type of subfilter used (more on this later).
* `BucketSize`: the number of buckets is just the capacity of the underlying
array (in bytes), divided by `BucketSize`. When `BucketSize` is specified as zero,
the value `sizeof(Subfilter::value_type)` is used.

The default configuration with `block<unsigned char,1>` corresponds to a
classical Bloom filter setting K bits per elements uniformly distributed across
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
a `filter<T, Hash, K, block<Block, K'>>` will set `K*K'` bits per element.
The tradeoff here is that insertion/lookup will be (much) faster than
with `filter<T, Hash, K*K'>` while the FPR will be worse (larger).
FPR is better the wider `Block` is.

### `multiblock<Block, K'>`

Instead of setting `K'` bits in a `Block` value, this subfilter sets
one bit on each of the elements of a `Block[K']` subarray. This improves FPR
without noticeably impacting performance with respect to `block<Block, K'>`,
although cacheline boundaries can be crossed when accessing the subarray.

### `fast_multiblock32<K'>`

Statistically equivalent to `multiblock<uint32_t, K'>`, but uses a much
faster SIMD-based algorithm when AVX2 is available. This algorithm is
pretty much the same as that of [Apache Kudu](https://kudu.apache.org/)
[`BlockBloomFilter`](https://github.com/apache/kudu/blob/master/src/kudu/util/block_bloom_filter_avx2.cc),
but that implementation is fixed to `K'` = 8 whereas we accept
any value between 1 and 8.

## Estimating FPR

For a classical Bloom filter, the theoretical false positive rate, under some simplifying assumptions,
is given by

$$FPR(n,m,k)=\left(1 - \left(1 - \frac{1}{m}\right)^{kn}\right)^k \approx \left(1 - e^{-kn/m}\right)^k \text{   for large } m,$$

where $$n$$ is the number of elements inserted in the filter, $$m$$ its capacity in bits and $$k$$ the
number of bits set per insertion (see a [derivation](https://en.wikipedia.org/wiki/Bloom_filter#Probability_of_false_positives)
of this formula). For a given inverse load factor $$c=m/n$$, the optimum $$k$$ is
the integer closest to:

$$k_{\text{opt}}=c\cdot\ln2,$$

yielding a minimum attainable FPR of $$0.5^{c \cdot \ln 2} \approx 1/1.6168^{c}$$.

In the case of a Boost.Bloom block filter of the form `filter<T, Hash, K, block<Block, K'>>`, we can extend
the approach from [Putze et al.](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=f376ff09a64b388bfcde2f5353e9ddb44033aac8)
to derive the (approximate but very precise) formula:

$$FPR_{\text{block}}(n,m,b,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}(i,kbn/m) \cdot FPR(i,b,k')\right)^{k},$$

where

$$\text{Pois}(i,\lambda)=\frac{\lambda^i e^{-\lambda}}{i!}$$

is the probability mass function of a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution)
with mean $$\lambda$$, and $$b$$ is the size of `Block` in bits. If we're using `multiblock<Block,K'>`, we have

$$FPR_\text{multiblock}(n,m,b,k,k')=\left(\sum_{i=0}^{\infty} \text{Pois}(i,kk'bn/m) \cdot FPR(i,b,1)^{k'}\right)^{k}.$$

As we have commented before, in general 

$$FPR_\text{block}(n,m,b,k,k'), FPR_\text{multiblock}(n,m,b,k,k') \geq FPR(n,m,kk'),$$

that is, block and multi-block filters have worse FPR than the classical filter for the same number of bits
set per insertion, but they will be much faster. We have the particular case

$$FPR_{\text{block}}(n,m,b,k,1)=FPR_{\text{multiblock}}(n,m,b,k,1)=FPR(n,m,k),$$

which follows simply from the fact that using `{block|multiblock}<Block, 1>` behaves exactly as
a classical Bloom filter.

We don't know of any closed, simple formula for the FPR of block and multiblock filters when
`Bucketsize` is not its "natural" size (`sizeof(Block)` for `block<Block, K'>`,
`K'*sizeof(Block)` for `multiblock<Block, K'>`), that is, when subfilter values overlap,
but empirical calculations show that $$FPR/FPR_\text{baseline}$$ improves (reduces) with smaller values of `BucketSize`
and larger values of $$k$$, $$k'$$ and $$c=m/n$$. Some examples:

* `filter<T, Hash, 1, multiblock<unsigned char, 9>, BucketSize>`:
$$\frac{FPR(c=12,\texttt{BucketSize}=1)}{FPR(c=12,\texttt{BucketSize}=0)}=0.65$$
* `filter<T, Hash, 1, multiblock<uint32_t, 9>, BucketSize>`: 
$$\frac{FPR(c=12,\texttt{BucketSize}=1)}{FPR(c=12,\texttt{BucketSize}=0)}=0.80$$
* `filter<T, Hash, 1, multiblock<uint64_t, 9>, BucketSize>`: 
$$\frac{FPR(c=12,\texttt{BucketSize}=1)}{FPR(c=12,\texttt{BucketSize}=0)}=0.87$$

(Remember that `BucketSize` = 0 selects the non-overlapping case.)
