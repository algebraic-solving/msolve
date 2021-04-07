- contains usolve => real root isolation and refinement for univariate
  polynomials with integer coefficients

Here is our (prelim.) todo list: 

* clean usolve and allow it to isolate real roots of rational parametrizations 

* let bsolve take more than 2 polynomials as input ; work parallelism there

* in bsolve, fix the bug when one input polynomial has degree <= 2

* in fglm, parallelism to be implemented in the linear algebra part 

* in fglm, memory management should be better (urgent)

* in fglm, finish degenerate cases DONE => but needs to be tested more.

* in fglm, bug when a variable is a leading monomial of the ideal (for grevlex ordering).
[FIXED]

* in fglm, implement a version with 16-bit primes

* in neogb parser with mpz coefficients [DONE]

* multi-modular rational reconstruction [DONE]

* in neogb, integrate (and rework) gbla 

* in neogb, elimination and weighted orderings. 

* in neogb, make a version with unsigned 32-bit primes [DONE]

* hilbert stuff to be implemented (to get the dimension, degree, etc.)

* try to reconstruct polynomials from the first m steps of the GB algorithm
  (for generic systems those may have smaller coefficients). then use this
  information to apply the tracer for further primes only from step m+1 onwards.

* Smaller examples + collect timings from other computer algebra systems
  during summer. [DONE]

* investigate the right functions to use in computer algebra systems [DONE]

* Finish to reconstruct the parametrizations + start using the tracer [DONE]

* connection with usolve todo

* about the paper: start by writing the tracer part to figure out what we
  actually need for the description of F4.

* prepare the conics example. Collect the timing of the tracer run on the
  conics.

* modular computation using the tracer and the prob linear algebra are ready. 
needs to store the parametrizations in appropriate variables and write them in 
files. [DONE]

* generate_lucky_primes

* normalization of parametrizations with derivative of the eliminating polynomial [DONE]

* non-radical case => manage the case where there are lin. forms in the GB

* check shape lemma position

* fix behaviour when input gb contains linear forms [DONE]

* fix bug when the input is a linear system (!) [DONE]

* output format when the GB is [1] (!)

* calls to compute_parametrizations_non_radical_case to be done 
  in nmod_fglm_compute trace functions

* multi-threading for multi-mod computations (duplicte_data function to finish to write). 
