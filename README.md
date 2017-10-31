# goldbach-partition

Goldbach Strong Conjecture (GSC) states that every even integer greater than 2 can be expressed as the sum of two primes.

_2n_ = _p_ + _q_

where:

  * _n_ >= 2, 
  * both _p_ and _q_ are prime numbers 

So far (by 2017) GSC has not been proven yet, although has been computationally verified up to _4 × 10^18_ in: 
Tomás Oliveira e Silva, Siegfried Herzog, and Silvio Pardi, _Empirical verification of the even Goldbach conjecture 
and computation of prime gaps up to 4·10^18_, Mathematics of Computation, vol. 83, no. 288, pp. 2033-2060, July 2014 
(published electronically on November 18, 2013).
Expression of even integer as a sum of two primes is called Goldbach Partition (GP).

This framework is used to find the fastest computational confirmation method of GSC.
  
## Methods based on primality check of components of 2n

_goldbach-fast_conf-prim_check.py_ is examining 6 different approach which are based on iteration over even integers, decomposition to two integers and final primality check of both components.
Search for GP is successful if both components are primes. _goldbach-fast_conf-prim_check-distr.py_ is a distributed version of this approach.

## Methods based on sum building from already known prime numbers

_goldbach-fast_conf-sum_build.py_ is examining how fast it is possible to build even integers from a set of already known prime numbers.

## Auxiliary script - reduction of primes

_goldbach-reduction.py_ verifies if two primes from at least one GP can be lowered by 2^x (x >= 1), still remaining prime numbers. It is examining hypothesis that:

_2n_ = _p_ + _q_ + _2^a_ + _2^b_

where:

  * _n_ >= 2
  * _a_ >= 1
  * _b_ >= 1
  * both _p_ and _q_ are primes

## Auxiliary script - GSC statistics

_goldbach_stats.py_ is calculating all possible sum of primes for consequtive even integers first, then is presenting various interesting statistical properties of GPs like:

  * number of elements in _R(n)_ (so called Goldbach comet)
  * min/max/avg difference in _R(n)_
  * trens of max difference in _R(n)_
  * minimal prime in _R(n)_
  * histogram of minimal primes
  * participation of max/min/avg prime from _R(n)_ in _n_

where:

  * _R(n)_ is a set of all possible GPs for _n_
  * _n_ is an even integer

## Auxiliary script - symmetrical primes

_goldbach-prime_gap.py_ is checking relations between GP and so called symmetrical primes. We can rewrite GSC to the following form: 

 "All positive integers > 1 can be expressed as a half of sum of two primes."

This means that every positive integer _n_ > 1 can be treated as a symmetry point for two primes _p_ and _q_:

_n_ = (_p_ + _q_) / 2
 
This also means that for every integer _n_ > 1 there exists yet another integer _i_ >= 0 that:

 * _p_ = _n_ - _i_
 * _q_ = _n_ + _i_

where _p_ and _q_ are primes (_p_ <= _q_).
  
## Dependencies

Framework depends on https://github.com/mbarylsk/primes which is supporting all required operations on prime numbers.