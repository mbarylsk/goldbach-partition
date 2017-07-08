# goldbach-partition

Goldbach Strong Conjecture (GSC) states that every even integer greater than 2 can be expressed as the sum of two primes.

_2n_ = _p_ + _q_

where:

  * _n_ >=2, 
  * _p_ and _q_ are prime numbers 

So far (by 2017) GSC has not been proven yet, although has been computationally verified up to _4 × 10^18_.
Expression of even integer as a sum of two primes is called Goldbach Partition (GP).

This framework is used to find the fastest computational confirmation method of GSC.
  
## Methods based on primality check of components of 2n

_goldbach-fast_conf-prim_check.py_ is examining 6 different approach which are based on iteration over even integers, decomposition to two integers and final primality check of both components.
Search for GP is successful if both components are primes.

## Methods based on possible sum building from already known prime numbers

_goldbach-fast_conf-sum_build.py_ is examining how fast it is possible to build even integers from a set of already known prime numbers.

## Dependencies

Framework depends on https://github.com/mbarylsk/primes which is supporting all required operations on prime numbers.