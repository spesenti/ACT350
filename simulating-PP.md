















How to simulate a Poisson process
=================================

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

1. Generating jump times
------------------------

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

    # note that the first value of the jump_times is aleady 0
    for(i in 2:n){
      jump_times[i] <- jump_times[i - 1] + rexp(1, lambda)
    }
    jump_times

    ##  [1]  0.0000000  0.5364204  0.6797592  1.1976821  1.4547997  1.6618271
    ##  [7]  2.0279133  2.1273068  3.4025758  3.8068813  4.0546129  4.5056082
    ## [13]  5.0440875  5.6096541  5.7136987  6.4815797  6.9479823  7.2108230
    ## [19]  9.7099614  9.7107387  9.9232506 10.1332872 10.3478895 10.6781867
    ## [25] 11.1893327 11.3628134 11.7966425 12.4920872 12.7023133 12.8938252
    ## [31] 13.4776079 13.4793536 13.6448219 14.0062459 14.6350073 14.9966568
    ## [37] 16.0899974 16.6408139 17.2828964 17.5200645 17.7641976 18.7120167
    ## [43] 19.2365888 20.2029024 20.6987375 20.7324573 21.2953410 21.8781391
    ## [49] 22.2128972 23.0060851

A more efficient implementation, avoiding the for loop, is

    jump_times_2 = cumsum(rexp(n-1, lambda))
    # note that the two generated sample_paths are different, 
    # since the Exponential random variables are different.
    jump_times_2

    ##  [1]  0.6544514  1.7731351  2.5005603  2.7452029  4.1155929  4.1822179
    ##  [7]  4.8691621  5.0903334  5.3616160  5.9824905  6.2469722  6.4557381
    ## [13]  7.2963775  7.7890900  8.3599224  8.8945182 10.0275949 10.6340718
    ## [19] 10.6483474 11.1267670 11.8834403 12.5927588 13.5410685 13.8921572
    ## [25] 14.2962352 14.4068175 14.9627015 15.8998191 15.9262117 17.0412885
    ## [31] 17.2213172 17.8005522 18.1404292 18.4289045 18.6027566 20.4972512
    ## [37] 20.5012016 21.0312639 21.2505673 22.7531266 23.2049816 23.2783735
    ## [43] 23.4024019 23.6441928 24.3085241 25.7634094 26.3752509 26.7056330
    ## [49] 26.7160101

2. Defining the Poisson process as a function
---------------------------------------------

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

    ## [1] 7

    poisson_process(10)

    ## [1] 20

    poisson_process(15)

    ## [1] 35

The first jump time is 0.5364204. Let us check whether the Poisson
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

3. Plotting a sample path
-------------------------

To plot a sample path we have to plot the implemented step function
`poisson_process`.

    plot(poisson_process, xlab = "time", ylab = "Value of Poisson process", main = NULL,
         verticals = FALSE, do.points = FALSE, xlim = c(0,jump_times[n]))

![](simulating-PP_files/figure-markdown_strict/simulated%20sample%20path-1.png)

**Question:** Rerun the code, what do you observe?  
**Question:** Rerun the code with a different *λ* or change the number
of jumps. How does the sample path change?

Renewal Theory
==============

A renewal process is defined as follows: Consider a sequence
*T*<sub>1</sub>, *T*<sub>2</sub>, … of independent identically
distributed (i.i.d.) non-negative random variables. Define the
stochastic process
*X*<sub>0</sub> = 0 , *X*<sub>*n*</sub> = *T*<sub>1</sub> + … + *T*<sub>*n*</sub> , *n* = 1, 2, … ,
and set
$$N(t) = \\max\\Big\\{ n \\in \\mathbb{N} ~\\Big|~ \\sum\_{i = 1}^n X\_i \\leq t \\Big\\}, \\quad \\text{for all } t \\geq 0. $$
Then the process {*N*(*t*) | *t* ≥ 0} is called a renewal process.

If the *T*<sub>*i*</sub> , *i* = 1, 2, … are i.i.d. Exponential with
parameter *λ* &gt; 0, then the renewal process {*N*(*t*) | *t* ≥ 0} is a
Poisson process.

**Task:**
=========

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
