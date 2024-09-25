Learning goals: 1. simulation and intuition of Poisson processes, 2.
simulation and intuition of compound Poisson processes, 3. introduction
to rMarkdown, 4. introduction to tidyverse/dplyr, 5. introduction to
ggplot.

# How to simulate a Poisson process

Recall that a Poisson process with rate *λ* &gt; 0 is a counting process
{*N*(*t*) : *t* ≥ 0} that has the following properties

1.  *starts at 0*: *N*(0) = 0 with probability 1
2.  *has independent increments*
3.  *is increasing*: *N*(*s*) ≤ *N*(*t*) for all *s* ≤ *t*
4.  *has Poisson distributed increments*:
    *P*(*N*(*s*+*t*)−*N*(*s*)=*n*)∼ Pois(*λ**t*) for all *s*, *t* ≥ 0.

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

We demonstrate how to generate sample jump times in two ways:

    # set the seed for generating random variables (for reproducibility) 
    set.seed(0310)

    # set the rate of the Poisson process/Exponential distribution
    lambda <- 2

    # number of jumps 
    n <- 50

    # initialize the jump times 
    # let's do in matrix - we will change to tibble/dataframe later
    jumps <- matrix(ncol = 3, nrow = n)

    head(jumps)

    ##      [,1] [,2] [,3]
    ## [1,]   NA   NA   NA
    ## [2,]   NA   NA   NA
    ## [3,]   NA   NA   NA
    ## [4,]   NA   NA   NA
    ## [5,]   NA   NA   NA
    ## [6,]   NA   NA   NA

    # we have do step 1., setting $\tau_0 = 0$
    jumps[1,] = c(1,0,1)

    # now we can generate the rest of the sample
    for(i in 2:n){
      jumps[i, 2] <- jumps[i - 1, 2] + rexp(1, lambda)
      jumps[i,3] = as.factor(1)
      jumps[i, 1] = i
    }

    head(jumps)

    ##      [,1]      [,2] [,3]
    ## [1,]    1 0.0000000    1
    ## [2,]    2 0.6936561    1
    ## [3,]    3 1.0369326    1
    ## [4,]    4 1.3594916    1
    ## [5,]    5 1.5650738    1
    ## [6,]    6 1.8360996    1

and

    # directly concatenate to our previous simulation
    jumps = rbind(cbind(c(1:n), c(0, cumsum(rexp(n-1, lambda))), rep(2,n)), jumps)

    head(jumps)

    ##      [,1]     [,2] [,3]
    ## [1,]    1 0.000000    2
    ## [2,]    2 2.229272    2
    ## [3,]    3 2.994915    2
    ## [4,]    4 4.458804    2
    ## [5,]    5 5.044406    2
    ## [6,]    6 5.565207    2

Which one is more efficient? Why?

We can easily add more sample jumps:

    jumps = rbind(cbind(c(1:n), c(0, cumsum(rexp(n-1, lambda))), rep(3,n)), jumps)

    jumps = rbind(cbind(c(1:n), c(0, cumsum(rexp(n-1, lambda))), rep(4,n)), jumps)

    head(jumps)

    ##      [,1]      [,2] [,3]
    ## [1,]    1 0.0000000    4
    ## [2,]    2 0.1073899    4
    ## [3,]    3 0.3720459    4
    ## [4,]    4 2.6008797    4
    ## [5,]    5 2.6462263    4
    ## [6,]    6 2.7696518    4

    jumps = data.frame(jumps)

    # name the columns
    colnames(jumps) = c("index", "sim", "path")

    # you can also make sample in a wider matrix that would look like this
    # how would you implement that?
    jumps_w = pivot_wider(jumps, names_from = path, values_from = sim)

    head(jumps_w)

    ## # A tibble: 6 x 5
    ##   index   `4`    `3`   `2`   `1`
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     1 0     0       0    0    
    ## 2     2 0.107 0.0904  2.23 0.694
    ## 3     3 0.372 1.06    2.99 1.04 
    ## 4     4 2.60  1.13    4.46 1.36 
    ## 5     5 2.65  1.77    5.04 1.57 
    ## 6     6 2.77  2.66    5.57 1.84

Now, let’s plot our sample jumps:

    # coerce our path label into factor type so we can use ggplot
    jumps$path = as.factor(jumps$path)

    # note that the four generated sample_paths are different, 
    # since the Exponential random variables are different samples.
    ggplot(data = jumps, aes(x = index, y = sim)) + 
      geom_step(aes(color = path, linetype = path)) +
      xlab("Index") + ylab("Jump Time") +
      ggtitle("Realizations of Poisson Process Jump times")

![](ACT350-simulation_PP_files/figure-markdown_strict/unnamed-chunk-4-1.png)

## 2. Defining the Poisson process as a function

Once the jump times have been simulated, we can go on to generate the
Poisson process.

Next, we can define the Poisson process as a step function that jumps
exactly one unit at every value in `jump_times`. Note that the
implemented Poisson process is now a *function* of the time *t*.

    # need to massage the data some more, R complains about numeric names
    names(jumps_w) = make.names(names(jumps_w))
    # the default X# names are fine for our purposes

    head(jumps_w)

    ## # A tibble: 6 x 5
    ##   index    X4     X3    X2    X1
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     1 0     0       0    0    
    ## 2     2 0.107 0.0904  2.23 0.694
    ## 3     3 0.372 1.06    2.99 1.04 
    ## 4     4 2.60  1.13    4.46 1.36 
    ## 5     5 2.65  1.77    5.04 1.57 
    ## 6     6 2.77  2.66    5.57 1.84

    # define a Poisson process as a function taking inputs of jump times
    # it outputs the value of a sample PP
    poisson_process = function(jumpVector){
      stepfun(jumpVector[2 : n], seq(0, n - 1, by = 1))(0:n)
    }

    # the first argument of `stepfun` is where the step function jumps
    # the second argument is the height of the step function

    # Let us check some values.
    # The starting value is 0.
    poisson_process(jumpVector = jumps_w$X3)[1]

    ## [1] 0

    # What is the value of N(3)?
    poisson_process(jumpVector = jumps_w$X3)[3]

    ## [1] 4

    poisson_process(jumpVector = jumps_w$X3)[10]

    ## [1] 17

    poisson_process(jumpVector = jumps_w$X3)[15]

    ## [1] 32

Note you can evaluate the Poisson process for any time point *t*,
however, keep in mind that we only simulated 50 jumps of the Poisson
process. Thus, evaluating *N*(*t*) for large *t* does not make sense.

## 3. Plotting a sample path

We have to apply our ‘poisson\_process’ function to each of our sample
paths. So let’s create a new dataframe that will contain this
information.

    # apply the poisson process function to each of our sample jump times
    PP = jumps_w %>% 
      reframe(across(2:5, poisson_process)) %>% 
      cbind(time = 0:n) 

    head(PP)

    ##   X4 X3 X2 X1 time
    ## 1  0  0  0  0    0
    ## 2  2  1  0  1    1
    ## 3  2  4  0  5    2
    ## 4  5  5  2  8    3
    ## 5 10  6  2 10    4
    ## 6 14  8  3 13    5

    # make a long version of the table with renaming (so we can use ggplot2)
    PP_l = PP %>% 
      pivot_longer(cols = 1:4) %>% 
      rename(pPath = name, pVal = value)

    head(PP_l)

    ## # A tibble: 6 x 3
    ##    time pPath  pVal
    ##   <int> <chr> <dbl>
    ## 1     0 X4        0
    ## 2     0 X3        0
    ## 3     0 X2        0
    ## 4     0 X1        0
    ## 5     1 X4        2
    ## 6     1 X3        1

To plot a sample path we have to plot the implemented step function
`poisson_process`.

    ggplot(data = PP_l) + 
      geom_step(aes(x = time, y = pVal, colour = pPath, linetype = pPath)) + 
      xlim(c(0,25)) +
      xlab("Time") + 
      ylab("Value of Poisson Process") + 
      ggtitle("Simulated Sample Paths of a Poisson Process") +
      scale_colour_discrete(name  ="Sample",
                              breaks=c("X1", "X2", "X3", "X4"),
                              labels=c("1", "2", "3", "4")) +
      scale_linetype_discrete(name  ="Sample",
                              breaks=c("X1", "X2", "X3", "X4"),
                              labels=c("1", "2", "3", "4"))

![](ACT350-simulation_PP_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    # what happens if I don't manually rename the legend?

## Confirming that samples are Poisson distributed

To (heuristically) check that our simulated process indeed follows a
Poisson distribution, can do a few “eyeball tests”.

We have to generate samples from a Poisson process at a fixed time *t*.

    # set a fixed time
    t = 8

    # set a reasonable sample size
    N = 1000

    # intialize a vector to collect our samples at time t, 
    # and a vector to temporarily store jump times
    tPois = rep(NA, N)
    tempJumps = rep(NA, n)

    # generate the samples
    for (i in 1:N){
      tempJumps = c(0, cumsum(rexp(n-1, lambda))) 
      # we overwrite this vector at every loop
      tPois[i] = poisson_process(jumpVector = tempJumps)[t + 1]
      # need to add 1 to account for index differences: R is 1 indexed
    }

First, we can check the mean and variance.

    mean(tPois)

    ## [1] 15.934

    var(tPois)

    ## [1] 15.9516

This is consistent with the theory, noting that the variance and
expectation of Poisson distributed random variables are equal. In a
numerical setting, we can expect the values to get closer as the number
of samples increases (recall: law of large numbers)

Also, recall for a Poisson process {*N*(*t*) : *t* ≥ 0}:
𝔼\[*N*(*t*)\] = *λ**t* = Var(*N*(*t*))
In particular, we have that *λ**t* = 2 ⋅ 8 = 16

Now we have the sample. What the next most obvious thing to look at if
we want to confirm the distribution of a sample? Histogram!

    ggplot() + 
      geom_histogram(data = tibble(tPois), 
                     aes(x = tPois, y = after_stat(density)), 
                     binwidth = 1, fill = "#56B4E9", alpha = 0.7)  + 
      geom_density() +
      stat_function(data = tibble(0:t*lambda*2), 
                    geom = "point", n = t*lambda*2 + 1, fun = dpois, 
                    args = list(lambda = t*lambda), aes(colour="#E69F00"), 
                    show.legend = FALSE) +
      xlim(c(0,t*lambda*2)) +
      xlab(" ")

![](ACT350-simulation_PP_files/figure-markdown_strict/unnamed-chunk-10-1.png)
The orange dots are the true values of a Poisson density. Why are they
dots? Are you convince that our sample is Poisson distributed? Why or
why not? How could we improve these results?

## Independent Increments

Let’s intuitively confirm that the increments are independent. Recall
that for a Poisson Process {*N*(*t*) : *t* ≥ 0}, the increments
*N*(*t*+*s*) − *N*(*s*) is independent of *N*(*s*). To test this
numerically, I arbitrarily select *t* = 4 and *s* = 3
i.e. *t* + *s* = 7.

Let’s simulate a sample of *N*(*s*) = *N*(3) and *N*(*t*+*s*) = *N*(7).

    # set s 
    s = 3
    r = 7 # denote r := t + s

    #initialize our vectors
    sPois = rep(NA, N)
    rPois = rep(NA, N)

    # generate the samples
    for (i in 1:N){
      tempJumps = c(0, cumsum(rexp(n-1, lambda))) #using the same vec as before
      sPois[i] = poisson_process(jumpVector = tempJumps)[s + 1]
      rPois[i] = poisson_process(jumpVector = tempJumps)[r + 1]
    }

How do we show independence between two samples? We can always eyeball a
scatterplot.

    ggplot() + geom_bin2d(aes(x = rPois - sPois, y = sPois), bins = 50) +
      xlab("N(t - s) - N(s)") +
      ylab("N(s)") + 
      ggtitle("Scatterplot of Poisson process Increments")

![](ACT350-simulation_PP_files/figure-markdown_strict/unnamed-chunk-12-1.png)

Do you think these samples are independent based off of the scatter
plot? Why or why not?

We can also check some independence tests… Consider the null hypothesis
that our two samples are statistically independent. We perform a Pearson
Chi-squared test (does our experiment satisfy the assumptions of this
test?).

    chisq.test(x = rPois - sPois, y = sPois, correct = TRUE)

    ## Warning in chisq.test(x = rPois - sPois, y = sPois, correct = TRUE): Chi-squared
    ## approximation may be incorrect

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  rPois - sPois and sPois
    ## X-squared = 245.83, df = 270, p-value = 0.8518

We fail to reject the null hypothesis. This may or may not be convincing
to you. Become confident in this fact by proving that Poisson increments
are independent (HW).

# How to simulate a compound Poisson process

In order to construct a compound Poisson process, we need a (regular)
Poisson process and a *jump size distribution*. Like Poisson processes,
a compound Poisson process is also a jump process, meaning that it is a
piecewise, increasing process. We define the compound Poisson process
constructively:

Let {*N*(*t*) : *t* ≥ 0} be a Poisson process with rate *λ* &gt; 0, let
*F* be a jump size distribution, and let random variables
*X*<sub>*i*</sub>, *i* = 1, 2, … follow the distribution *F* identically
and independently. Then, the compound Poisson distribution
{*W*(*t*) : *t* ≥ 0} is defined as
$$
W(t) := \\sum\_{i = 1}^{N(t)} X\_i
$$
for every *t* ≥ 0.

We already know how to generate a Poisson process. Let’s make a new
sample:

    # we can use the same sample sizes and parameters as earlier
    # so we're not going to re-declare the variables

    # generate the Poisson jumps
    jumps_cPP = c(0, cumsum(rexp(n-1, lambda)))

    # generate the Poisson Process
    PP_cPP = poisson_process(jumpVector = jumps_cPP)

Now, we need to sample the jump sizes, sometimes this is called the
*severity* of the compound Poisson process - especially in the context
of actuarial science (recall the severity of insurance claims). We will
consider the Weibull distribution with shape and parameters to be 0.75
and 1 respectively to be the severity. Note that the Weibull
distribution with shape between 0 and 1 is a heavy tailed distribution.

    # generate samples from the jump size distribution
    # let's use the Weibull(shape = 3/4, scale = 1)
    shape = 3/4
    scale = 1
    severity = rweibull(500, shape, scale)
    # we don't know how many jump size samples we'll need
    # so I choose 500 to make sure we have enough

Now we have all the components necessary to generate one sample path of
our compound Poisson process.

    # find the values of the compound Poisson process
    value_cPP = c(0, cumsum(severity[1 : n - 1]))

And we plot this

    ggplot() + geom_step(aes(x = jumps_cPP, y = value_cPP)) +
      xlab("Time") +
      ylab("Value") +
      ggtitle("Compound Poisson Process")

![](ACT350-simulation_PP_files/figure-markdown_strict/unnamed-chunk-17-1.png)

## Tasks:

1.  Plot multiple sample paths of a compound Poisson Process in one plot
    and provide an interpretation.

2.  Show that a compound Poisson process has the properties of a Poisson
    process, or provide a counterexample (in R and/or on paper),
    i.e. that it i) has value 0 at time 0, ii) is non-decreasing in
    time, iii) is Poisson distributed with rate *t**λ* at any time
    *t* ≥ 0, and iv) that disjoint increments are independent.

3.  An insurance company models a total claim portfolio via a compound
    Poisson process
    $$C(t) = \\sum\_{i = 1}^{N(t)} X\_i$$
    where *C*(*t*) is the total value of the claims at time *t* (with
    weeks as unit), *X*<sub>1</sub>, *X*<sub>2</sub>, … are
    independently and identically distributed log-normally with a mean
    claim value of $1500 and standard deviation of $2000. The frequency
    of the claims is modelled by an (independent) Poisson process
    {*N*(*t*:*t*≥0)} with an average of 4 claims per week.

    1.  Implement this process in ‘R’ and plot 4 sample paths (in one
        plot) for a time period of at least 1 year. Discuss the results.
    2.  What is the expected total claim to the insurance company after
        6 months? What about 1 year? Solve on paper and approximate this
        in ’R”.
    3.  Now, say that the intensity of the claims is not constant,
        i.e. consider *λ* to be a function of time *t*. In particular,
        use $\\lambda(t) = 3\\left(\\cos(\\frac{t}{2\\pi}) + 1\\right)$.
        Why would an insurance company want to consider such a scenario?
        Repeat parts a) and b), and provide an interpretation.
