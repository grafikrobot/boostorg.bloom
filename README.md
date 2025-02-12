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
but impacts performance with respect to `block<Block, K'>`, among other
things because cacheline boundaries can be crossed when accessing the subarray.

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

yielding a minimum attainable FPR of $$1/2^{k_{\text{opt}}} \approx 0.6185^{c}$$.

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

$$FPR_\text{block}(n,m,b,k,k') \geq FPR_\text{multiblock}(n,m,b,k,k') \geq FPR(n,m,kk'),$$

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

## Experimental results

All benchmarks run on a Windows 10 machine with 8GB RAM and an Intel Core i5-8265U CPU @1.60GHz
(base frequency) using Clang 13.0.1 for Visual Studio (clang-cl) in 64-bit (x64), `/arch:AVX2`,
release mode (see [benchmarking code](benchmark/comparison_table.cpp)).

The tables show the FPR and execution times in nanoseconds per operatoin 
for three different configurationsof `boost::bloom::filter<int, boost::hash<int>, ...>`.
The number of inserted elements is `N` = 10 million. Filters are constructed with a capacity
`c*N` (bits), so `c` is the number of bits used per element. For each combination of `c` and
a given filter configuration, the optimum value of `K` (that yielding the minimum FPR)
has been used.

For reference, we provide also insertion, successful lookup and unsuccessful lookup times
for a `boost::unordered_flat_set<int>` with `N` = 10 million elements.
This container uses 63.75 bits per element.

### 64-bit mode

<table>
  <tr><th colspan="3"><code>boost::unordered_flat_set</code></tr>
  <tr>
    <th>insertion</th>
    <th>successful<br/>lookup</th>
    <th>unsuccessful<br/>lookup</th>
  </tr>
  <tr>
    <td align="right">81.96</td>
    <td align="right">27.52</td>
    <td align="right">16.41</td>
  </tr>
</table>

<table>
  <tr>
    <th></th>
    <th colspan="5"><code>filter&lt;K></code></th>
    <th colspan="5"><code>filter&lt;1, block&lt;uint64_t, K>></code></th>
    <th colspan="5"><code>filter&lt;1, multiblock&lt;uint64_t, K>></code></th>
  </tr>
  <tr>
    <th>c</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
  </tr>
  <tr>
    <td align="center">8</td>
    <td align="center">6</td>
    <td align="right">2.1566</td>
    <td align="right">75.06</td>
    <td align="right">46.92</td>
    <td align="right">35.14</td>
    <td align="center">4</td>
    <td align="right">3.3462</td>
    <td align="right">18.14</td>
    <td align="right">12.25</td>
    <td align="right">12.23</td>
    <td align="center">5</td>
    <td align="right">2.4515</td>
    <td align="right">20.74</td>
    <td align="right">19.51</td>
    <td align="right">17.57</td>
  </tr>
  <tr>
    <td align="center">12</td>
    <td align="center">9</td>
    <td align="right">0.3146</td>
    <td align="right">136.26</td>
    <td align="right">76.64</td>
    <td align="right">41.55</td>
    <td align="center">5</td>
    <td align="right">1.0310</td>
    <td align="right">20.79</td>
    <td align="right">14.11</td>
    <td align="right">14.86</td>
    <td align="center">8</td>
    <td align="right">0.4244</td>
    <td align="right">20.24</td>
    <td align="right">16.80</td>
    <td align="right">17.59</td>
  </tr>
  <tr>
    <td align="center">16</td>
    <td align="center">11</td>
    <td align="right">0.0456</td>
    <td align="right">184.34</td>
    <td align="right">98.16</td>
    <td align="right">42.74</td>
    <td align="center">6</td>
    <td align="right">0.4035</td>
    <td align="right">21.57</td>
    <td align="right">15.88</td>
    <td align="right">15.55</td>
    <td align="center">11</td>
    <td align="right">0.0770</td>
    <td align="right">89.57</td>
    <td align="right">82.17</td>
    <td align="right">80.68</td>
  </tr>
  <tr>
    <td align="center">20</td>
    <td align="center">14</td>
    <td align="right">0.0066</td>
    <td align="right">235.31</td>
    <td align="right">130.24</td>
    <td align="right">44.70</td>
    <td align="center">7</td>
    <td align="right">0.1879</td>
    <td align="right">23.97</td>
    <td align="right">16.85</td>
    <td align="right">17.09</td>
    <td align="center">14</td>
    <td align="right">0.0153</td>
    <td align="right">98.09</td>
    <td align="right">88.73</td>
    <td align="right">83.27</td>
  </tr>
  <tr>
    <th></th>
    <th colspan="5"><code>filter&lt;K, fast_multiblock32&ltK>></code></th>
    <th colspan="5"><code>filter&lt;1, block&lt;uint64_t, K>, 1></code></th>
    <th colspan="5"><code>filter&lt;1, multiblock&lt;uint64_t, K>, 1></code></th>
  </tr>
  <tr>
    <th>c</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
  </tr>
  <tr>
    <td align="center">8</td>
    <td align="center">5</td>
    <td align="right">2.7367</td>
    <td align="right">16.71</td>
    <td align="right">12.90</td>
    <td align="right">12.77</td>
    <td align="center">5</td>
    <td align="right">3.0580</td>
    <td align="right">15.98</td>
    <td align="right">12.07</td>
    <td align="right">11.96</td>
    <td align="center">5</td>
    <td align="right">2.3256</td>
    <td align="right">16.49</td>
    <td align="right">13.73</td>
    <td align="right">13.59</td>
  </tr>
  <tr>
    <td align="center">12</td>
    <td align="center">8</td>
    <td align="right">0.5399</td>
    <td align="right">15.40</td>
    <td align="right">11.49</td>
    <td align="right">11.41</td>
    <td align="center">6</td>
    <td align="right">0.8253</td>
    <td align="right">25.64</td>
    <td align="right">15.20</td>
    <td align="right">15.30</td>
    <td align="center">8</td>
    <td align="right">0.3777</td>
    <td align="right">23.02</td>
    <td align="right">16.70</td>
    <td align="right">16.18</td>
  </tr>
  <tr>
    <td align="center">16</td>
    <td align="center">8</td>
    <td align="right">0.1314</td>
    <td align="right">16.91</td>
    <td align="right">12.97</td>
    <td align="right">13.22</td>
    <td align="center">7</td>
    <td align="right">0.2891</td>
    <td align="right">22.21</td>
    <td align="right">17.15</td>
    <td align="right">17.24</td>
    <td align="center">11</td>
    <td align="right">0.0654</td>
    <td align="right">86.78</td>
    <td align="right">53.44</td>
    <td align="right">52.94</td>
  </tr>
  <tr>
    <td align="center">20</td>
    <td align="center">8</td>
    <td align="right">0.0421</td>
    <td align="right">17.29</td>
    <td align="right">12.34</td>
    <td align="right">12.38</td>
    <td align="center">8</td>
    <td align="right">0.1213</td>
    <td align="right">24.80</td>
    <td align="right">16.85</td>
    <td align="right">16.62</td>
    <td align="center">14</td>
    <td align="right">0.0111</td>
    <td align="right">99.33</td>
    <td align="right">98.62</td>
    <td align="right">98.65</td>
  </tr>
</table>

### 32-bit mode

<table>
  <tr><th colspan="3"><code>boost::unordered_flat_set</code></tr>
  <tr>
    <th>insertion</th>
    <th>successful<br/>lookup</th>
    <th>unsuccessful<br/>lookup</th>
  </tr>
  <tr>
    <td align="right">84.74</td>
    <td align="right">31.87</td>
    <td align="right">16.97</td>
  </tr>
</table>
<table>
  <tr>
    <th></th>
    <th colspan="5"><code>filter&lt;K></code></th>
    <th colspan="5"><code>filter&lt;1, block&lt;uint64_t, K>></code></th>
    <th colspan="5"><code>filter&lt;1, multiblock&lt;uint64_t, K>></code></th>
  </tr>
  <tr>
    <th>c</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
  </tr>
  <tr>
    <td align="center">8</td>
    <td align="center">6</td>
    <td align="right">3.1793</td>
    <td align="right">68.78</td>
    <td align="right">58.73</td>
    <td align="right">37.94</td>
    <td align="center">4</td>
    <td align="right">3.7080</td>
    <td align="right">30.30</td>
    <td align="right">15.12</td>
    <td align="right">15.34</td>
    <td align="center">5</td>
    <td align="right">2.7967</td>
    <td align="right">25.68</td>
    <td align="right">21.18</td>
    <td align="right">23.07</td>
  </tr>
  <tr>
    <td align="center">12</td>
    <td align="center">9</td>
    <td align="right">1.6736</td>
    <td align="right">135.91</td>
    <td align="right">92.26</td>
    <td align="right">49.24</td>
    <td align="center">5</td>
    <td align="right">1.4633</td>
    <td align="right">42.49</td>
    <td align="right">23.55</td>
    <td align="right">23.69</td>
    <td align="center">8</td>
    <td align="right">0.8633</td>
    <td align="right">25.04</td>
    <td align="right">20.94</td>
    <td align="right">21.09</td>
  </tr>
  <tr>
    <td align="center">16</td>
    <td align="center">11</td>
    <td align="right">1.4012</td>
    <td align="right">166.22</td>
    <td align="right">113.25</td>
    <td align="right">49.71</td>
    <td align="center">6</td>
    <td align="right">0.8558</td>
    <td align="right">47.21</td>
    <td align="right">30.13</td>
    <td align="right">30.51</td>
    <td align="center">11</td>
    <td align="right">0.5380</td>
    <td align="right">97.28</td>
    <td align="right">96.97</td>
    <td align="right">95.76</td>
  </tr>
  <tr>
    <td align="center">20</td>
    <td align="center">14</td>
    <td align="right">1.4417</td>
    <td align="right">240.92</td>
    <td align="right">160.94</td>
    <td align="right">53.82</td>
    <td align="center">7</td>
    <td align="right">0.6439</td>
    <td align="right">47.01</td>
    <td align="right">29.86</td>
    <td align="right">29.87</td>
    <td align="center">14</td>
    <td align="right">0.4777</td>
    <td align="right">101.01</td>
    <td align="right">99.93</td>
    <td align="right">100.89</td>
  </tr>
  <tr>
    <th></th>
    <th colspan="5"><code>filter&lt;K, fast_multiblock32&ltK>></code></th>
    <th colspan="5"><code>filter&lt;1, block&lt;uint64_t, K>, 1></code></th>
    <th colspan="5"><code>filter&lt;1, multiblock&lt;uint64_t, K>, 1></code></th>
  </tr>
  <tr>
    <th>c</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
    <th>K</th>
    <th>FPR<br/>[%]</th>
    <th>ins.</th>
    <th>succ.<br/>lkp.</th>
    <th>uns.<br/>lkp.</th>
  </tr>
  <tr>
    <td align="center">8</td>
    <td align="center">5</td>
    <td align="right">3.0370</td>
    <td align="right">17.99</td>
    <td align="right">13.53</td>
    <td align="right">13.94</td>
    <td align="center">5</td>
    <td align="right">3.4309</td>
    <td align="right">33.84</td>
    <td align="right">17.59</td>
    <td align="right">17.27</td>
    <td align="center">5</td>
    <td align="right">2.7102</td>
    <td align="right">20.21</td>
    <td align="right">23.60</td>
    <td align="right">26.72</td>
  </tr>
  <tr>
    <td align="center">12</td>
    <td align="center">8</td>
    <td align="right">1.0461</td>
    <td align="right">20.19</td>
    <td align="right">12.64</td>
    <td align="right">12.30</td>
    <td align="center">6</td>
    <td align="right">1.2554</td>
    <td align="right">43.27</td>
    <td align="right">23.88</td>
    <td align="right">24.42</td>
    <td align="center">8</td>
    <td align="right">0.8259</td>
    <td align="right">27.83</td>
    <td align="right">22.00</td>
    <td align="right">21.87</td>
  </tr>
  <tr>
    <td align="center">16</td>
    <td align="center">8</td>
    <td align="right">0.5877</td>
    <td align="right">17.22</td>
    <td align="right">12.84</td>
    <td align="right">12.79</td>
    <td align="center">7</td>
    <td align="right">0.7493</td>
    <td align="right">45.95</td>
    <td align="right">26.00</td>
    <td align="right">24.96</td>
    <td align="center">11</td>
    <td align="right">0.5273</td>
    <td align="right">97.09</td>
    <td align="right">106.32</td>
    <td align="right">106.63</td>
  </tr>
  <tr>
    <td align="center">20</td>
    <td align="center">8</td>
    <td align="right">0.5039</td>
    <td align="right">17.81</td>
    <td align="right">13.19</td>
    <td align="right">13.44</td>
    <td align="center">8</td>
    <td align="right">0.5794</td>
    <td align="right">49.02</td>
    <td align="right">24.39</td>
    <td align="right">24.73</td>
    <td align="center">14</td>
    <td align="right">0.4747</td>
    <td align="right">108.56</td>
    <td align="right">120.47</td>
    <td align="right">120.42</td>
  </tr>
</table>
