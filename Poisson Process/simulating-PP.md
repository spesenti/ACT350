# How to simulate a Poisson process

Recall that a Poisson process with rate *λ* &gt; 0 is a counting process
{*N*(*t*), *t* ≥ 0} that has the following properties

1.  *starts at 0*: *N*(0) = 0 with probability 1  
2.  *is increasing*: *N*(*s*) ≤ *N*(*t*) for all *s* ≤ *t*  
3.  *has Poisson distributed increments*:
    *P*(*N*(*s* + *t*) − *N*(*s*) = *n*)∼ Pois(*λ**t*) for all
    *s*, *t* ≥ 0.

However, this definition does not help if we want to simulate a Poisson
process. Thus, to simulate a Poisson process we use its alternative
definition. That is, for independent random variables
*X*<sub>1</sub>, *X*<sub>2</sub>, … that are Exponentially distributed
with rate *λ* &gt; 0, a Poisson process with rate *λ* &gt; 0 can be
written as

Thus, for simulating a Poisson process we first need to simulate
exponential random variables, add them up and check whether they are
smaller or equal to *t*. However, we will run into problems when we want
to store the value of *N*(*t*) for every 0 ≤ *t* &lt; ∞. The way out is
to implement *N*(*t*) as a function of *t* in the following way. First
generate the jump points of the Poisson process. Second, we know that a
Poisson process can only jump one unit at a time. Thus, we can define
the Poisson process as a step function that jumps exactly one unit at
every jump point.

## 1. Generating jump times

This is an algorithm for generating the jump times of a Poisson process

1.  Set *τ*<sub>0</sub> = 0.  
2.  For *k* ≥ 1 do
    1.  generate *X*∼Exp(*λ*)  
    2.  set *τ*<sub>*k*</sub> = *τ*<sub>*k* − 1</sub> + *X*

The *τ*<sub>*k*</sub>, *k* = 0, 1, … are the time points, where the
Poisson process jumps.

**Numerical implementation: R code**

    # set the seed for generating random variables (for reproducability) 
    # set.seed(2019)

    # set the rate of the Poisson process/Exponential distribution
    lambda <- 2
    # number of jumps 
    n <- 50
    # initialise the jump times
    jump_times <- rep(0, lenthought = n)

    # note that the first value of the jump_times is already 0
    for(i in 2:n){
      jump_times[i] <- jump_times[i - 1] + rexp(1, lambda)
    }
    jump_times

    ##  [1]  0.00000000  0.02448321  0.24809148  0.35258581  1.39493350  1.57780882
    ##  [7]  1.69025705  1.72996782  2.29834041  2.36717230  2.71170119  2.85881459
    ## [13]  3.07426387  3.88015562  5.62763056  6.66078917  6.75034561  6.96574035
    ## [19]  7.18283283  8.08959753  8.94817806  8.99011209  9.10927751  9.19064317
    ## [25]  9.67051369 11.55664340 11.57695735 12.21075934 13.44668427 13.94012428
    ## [31] 14.53353901 15.53672306 15.66627549 15.92759279 16.74850066 17.01516505
    ## [37] 17.35273957 17.95037395 18.64770376 19.41962784 19.58733482 20.14885744
    ## [43] 20.28945536 21.07511818 21.29603414 21.31948105 21.36684063 22.05026880
    ## [49] 22.25317902 23.24447851

A more efficient implementation, avoiding the for loop, is

    jump_times_2 = c(0, cumsum(rexp(n-1, lambda)))
    # note that the two generated sample_paths are different, 
    # since the Exponential random variables are different.
    jump_times_2

    ##  [1]  0.0000000  0.3761636  0.5205277  1.2355153  1.8795389  2.8855586
    ##  [7]  3.1941576  4.8642651  5.5442555  5.5924231  5.6143261  6.5785516
    ## [13]  7.0027804  7.0173085  8.4408377  8.7897951  9.3788216 10.0003501
    ## [19] 11.1333459 12.6937935 12.7761392 13.2042240 13.3196488 13.5066792
    ## [25] 14.0419907 14.1945179 14.2587654 14.9511660 15.1625553 15.1936642
    ## [31] 16.7366934 16.7603917 16.9517960 17.0956475 18.2244326 18.2816330
    ## [37] 18.7197244 19.8831578 20.5477862 20.9298677 21.3082175 21.4944692
    ## [43] 21.5955008 22.5190349 22.8327660 23.3417422 24.6036850 25.0562683
    ## [49] 25.2067298 25.2116927

## 2. Defining the Poisson process as a function

Next, we can define the Poisson process as a step function that jumps
exactly one unit at every value in `jump_times`. Note that the
implemented Poisson process is now a *function* of the time *t*.

    # define the Poisson process as a function
    poisson_process <- stepfun(jump_times[2 : 50], seq(0, n - 1, by = 1))
    # the first argument of `stepfun` is where the step function jumps
    # the second argument is the hight of the step function

    # Let us check some values.
    # The starting value is 0.
    poisson_process(0)

    ## [1] 0

    # What is the value of N(3)? (Note that we use jump_times and not jump_times_2)
    poisson_process(3)

    ## [1] 11

    poisson_process(10)

    ## [1] 24

    poisson_process(15)

    ## [1] 30

The first jump time is 0.0244832. Let us check whether the Poisson
process jumps at that time.

    # The first jump of the Poisson process is 
    first_jump <- jump_times[2]
    poisson_process(first_jump - 0.001)

    ## [1] 0

    poisson_process(first_jump + 0.001)

    ## [1] 1

Note you can evaluate the Poisson process for any time point *t*,
however, keep in mind that we only simulated 50 jumps of the Poisson
process. Thus, evaluating `N(t)` for *t* large does not make sense.

## 3. Plotting a sample path

To plot a sample path we have to plot the implemented step function
`poisson_process`.

    plot(poisson_process, xlab = "time", ylab = "Value of Poisson process", main = NULL,
         verticals = FALSE, do.points = FALSE, xlim = c(0,jump_times[n]))

<img src="simulating-PP_files/figure-markdown_strict/simulated sample path-1.png" style="display: block; margin: auto;" />

**Question:** Rerun the code, what do you observe?  
**Question:** Rerun the code with a different *λ* or change the number
of jumps. How does the sample path change?

## Tasks:

**a)** Plot multiple sample paths in one plot and provide an
interpretation.

**b)** Plot a histogram of the Poisson process at a time *t*. Verify
that it is indeed a Poisson random variable with parameter *λ*. E.g.,
calculate the sample mean, sample variance, compare the histograms.

**b)** Plot a histogram of an increment of a Poisson process and verify
that it is Poisson distributed.

**c)** How can you show that disjoint increments are indeed independent?

# Renewal Theory

A renewal process is defined as follows: Consider a sequence
*T*<sub>1</sub>, *T*<sub>2</sub>, … of independent identically
distributed (i.i.d.) non-negative random variables. Define the
stochastic process
*X*<sub>0</sub> = 0 ,  *X*<sub>*n*</sub> = *T*<sub>1</sub> + … + *T*<sub>*n*</sub> ,  *n* = 1, 2, … ,
and set
$$N(t) = \\max\\Big\\{ n \\in \\mathbb{N} ~\\Big|~ \\sum\_{i = 1}^n X\_i \\leq t \\Big\\}, \\quad \\text{for all } t \\geq 0. $$
Then the process {*N*(*t*) | *t* ≥ 0} is called a renewal process.

If the *T*<sub>*i*</sub> , *i* = 1, 2, … are i.i.d. Exponential with
parameter *λ* &gt; 0, then the renewal process {*N*(*t*) | *t* ≥ 0} is a
Poisson process.

## Tasks:

**a)** Generate a code to simulate from a renewal process, where the
interarrival times *T*<sub>*i*</sub> , *i* = 1, 2, … are i.i.d.
LogNormal distributed with mean 2 and standard deviation 0.5. *Hint: use
the function `rlnorm`*.

**b)** Plot a sample path of the renewal process define in **a)**.

**c)** What do you observe compared to the Poisson process?

**d)** Repeat tasks **a)** to **b)** with interarrival times being Gamma
distributed with the same mean and standard deviation as in **a)**.
*Hint: use the function `rgamma`*.

**e)** Compare the Poisson process with the renewal process from **a)**
and the renewal process from **d)** and give interpretations.
