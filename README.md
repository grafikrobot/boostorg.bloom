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
  typename Subfilter=block<unsigned char, 1>, std::size_t BucketSize=0,
  typename Allocator=std::allocator<unsigned char>  
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
a `bloom::filter<T, Hash, K, block<Block, K'>>` will set `K*K'` bits per element.
The tradeoff here is that insertion/lookup will be (much) faster than
with `bloom::filter<T, Hash, K*K'>` while the FPR will be worse (larger).
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
