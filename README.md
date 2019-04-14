trialr - Clinical Trial Designs in `RStan`
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# trialr

[![cran
version](http://www.r-pkg.org/badges/version/trialr)](https://cran.r-project.org/web/packages/trialr)
<https://cranlogs.r-pkg.org/badges/trialr>
<https://cranlogs.r-pkg.org/badges/grand-total/trialr>

`trialr` is a collection of Bayesian clinical trial designs implemented
in Stan and R. The documentation is available at
<https://brockk.github.io/trialr/>

Many notable Bayesian designs for clinical trials have been published.
However, one of the factors that has constrained their adoption is
availability of software. We present here some of the most popular,
implemented and demonstrated in a consistent style, leveraging the
powerful Stan environment.

It should be stressed that Bayesian trialists are not generally without
code. Often authors make available code with their design publication.
There are also some fantastic packages that aid the use of certain
designs. However, challenges to use still persist. The disparate methods
are naturally presented in a style that appeals to the particular
author. Features implemented in one package for one design may be
missing in another. Sometimes the technology chosen may only be
available on one particular operating system, or the chosen technology
may have fallen into disuse.

`trialr` seeks to address these problems. Models are specified in
[Stan](http://mc-stan.org/), a state-of-the-art environment for Bayesian
analysis. It uses Hamiltonian Monte Carlo to take samples from the
posterior distributions. This method is more efficient than Gibbs
sampling, for instance, and reliable inference can be performed on a few
thousand posterior samples. R, Stan and `trialr` are each available on
Mac, Linux, and Windows, so all of the examples presented here should
work on each operating system. Furthermore, Stan offers a very simple
method to split the sampling across *n* cores, taking full advantage of
the modern multicore processor in your computer (probably).

The designs implemented in `trialr` are introduced briefly below, and
developed more fully in vignettes. We focus on real-life usage,
including:

  - fitting models to observed data using your prior;
  - processing posterior samples to produce useful inferences;
  - and visualising inferences using modern `ggplot` graphics.

# Examples

In all examples, we will need to load `trialr`

``` r
library(trialr)
```

## CRM

The Continual Reassessment Method (CRM) was first published by
@OQuigley1990. It assumes a smooth mathematical form for the
dose-toxicity curve to conduct a dose-finding trial seeking a maximum
tolerable dose. There are many variations to suit different clinical
scenarios and the design has enjoyed *relatively* common use (although
nowhere near as common as the ubiquitous and inferior 3+3 design).

We will demonstrate the method using a notional trial example. In a
scenario of five potential doses, let us assume that we seek the dose
with probability of toxicity closest to 25% where our prior guesses of
the rates of toxicity can be represented:

``` r
target <- 0.25
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
```

Let us assume that we have already treated 2 patients each at doses 2, 3
and 4, having only seen toxicity at dose-level 4. What dose should we
give to the next patient or cohort? We can fit the data to the popular
empiric model

``` r
mod1 <- stan_crm(outcome_str = '2NN 3NN 4TT', skeleton = skeleton, 
                 target = target, model = 'empiric', beta_sd = sqrt(1.34), 
                 seed = 123)
```

The fitted model contains lots of useful of information:

``` r
mod1
#>   Patient Dose Toxicity Weight
#> 1       1    2        0      1
#> 2       2    2        0      1
#> 3       3    3        0      1
#> 4       4    3        0      1
#> 5       5    4        1      1
#> 6       6    4        1      1
#> 
#>   DoseLevel Skeleton N Tox ProbTox MedianProbTox ProbMTD
#> 1         1     0.05 0   0   0.108        0.0726  0.2140
#> 2         2     0.15 2   0   0.216        0.1900  0.2717
#> 3         3     0.25 2   0   0.310        0.2972  0.2657
#> 4         4     0.40 2   2   0.444        0.4484  0.2090
#> 5         5     0.60 0   0   0.624        0.6395  0.0395
#> 
#> The model targets a toxicity level of 0.25.
#> The dose with estimated toxicity probability closest to target is 2.
#> The dose most likely to be the MTD is 2.
```

``` r
library(ggplot2)
plot_df = data.frame(DoseLevel = 1:length(skeleton),
                     ProbTox = mod1$prob_tox)
ggplot(plot_df, aes(x = DoseLevel, y = ProbTox)) +
  geom_point() + geom_line() + ylim(0, 1) + 
  geom_hline(yintercept = target, col = 'orange', linetype = 'dashed') +
  labs(title = 'Posterior dose-toxicity curve under empiric CRM model')
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

Several variants of the CRM [are implemented in
‘trialr’](https://brockk.github.io/trialr/articles/CRM.html).
Further visualisation techniques are demonstrated in the [Visualisation
in CRM](https://brockk.github.io/trialr/articles/CRM-visualisation.html)
vignette.

## EffTox

EffTox by @Thall2004 is a dose-finding design that uses binary efficacy
and toxicity outcomes to select a dose with a high utility score. We
present it briefly here but there is a much more thorough examination in
the [EffTox
vignette](https://brockk.github.io/trialr/articles/EffTox.html).

For demonstration, We fit the model parameterisation introduced by
@Thall2014 to the following notional outcomes:

| Patient | Dose-level | Toxicity | Efficacy |
| :-----: | :--------: | :------: | :------: |
|    1    |     1      |    0     |    0     |
|    2    |     1      |    0     |    0     |
|    3    |     1      |    0     |    1     |
|    4    |     2      |    0     |    1     |
|    5    |     2      |    0     |    1     |
|    6    |     2      |    1     |    1     |

``` r
outcomes <- '1NNE 2EEB'
mod <- stan_efftox_demo(outcomes, seed = 123)
```

``` r
mod
#>   Patient Dose Toxicity Efficacy
#> 1       1    1        0        0
#> 2       2    1        0        0
#> 3       3    1        0        1
#> 4       4    2        0        1
#> 5       5    2        0        1
#> 6       6    2        1        1
#> 
#>   DoseLevel   ProbEff    ProbTox ProbAccEff ProbAccTox    Utility
#> 1         1 0.4045039 0.08990953    0.33175    0.92700 -0.3397885
#> 2         2 0.7917219 0.09875146    0.94575    0.92250  0.4237935
#> 3         3 0.9313427 0.21522248    0.98475    0.72900  0.5249445
#> 4         4 0.9572788 0.30606939    0.98475    0.62925  0.4380717
#> 5         5 0.9657038 0.36255571    0.98350    0.57725  0.3685257
#>   Acceptable
#> 1       TRUE
#> 2       TRUE
#> 3       TRUE
#> 4      FALSE
#> 5      FALSE
#> 
#> The model recommends selecting dose-level 3.
```

In this instance, after evaluation of our six patients, the dose
advocated for the next group is dose-level 3. This is contained in the
fitted object:

``` r
mod$recommended_dose
#> [1] 3
```

This is not surprising because dose 3 has the highest utility score:

``` r
mod$utility
#> [1] -0.3397885  0.4237935  0.5249445  0.4380717  0.3685257
```

Sometimes, doses other than the maximal-utility dose will be recommended
because of the dose-admissibility rules. See the papers for details.

Functions are provided to create useful plots. For instance, it is
illuminating to plot the posterior means of the probabilities of
efficacy and toxicity at each of the doses on the trade-off contours.
The five doses are shown in red. Doses closer to the lower-right corner
have higher
utility.

``` r
efftox_contour_plot(mod$dat, prob_eff = mod$prob_eff, prob_tox = mod$prob_tox)
title('EffTox utility contours')
```

![](man/figures/README-unnamed-chunk-11-1.png)<!-- -->

This example continues in the [EffTox
vignette](https://brockk.github.io/trialr/articles/EffTox.html).

There are many publications related to EffTox but the two most important
are @Thall2004 and @Thall2014.

## Hierachical analysis of response in related cohorts

Sticking with Peter Thall’s huge contribution to Bayesian clinical
trials, @Thall2003 described a method for analysing treatment effects of
a single intervention in several sub-types of a single disease.

We demonstrate the method for partially-pooling response rates to a
single drug in various subtypes of sarcoma. The following convenience
function returns the necessary data:

Fitting the data to the model:

``` r
mod0 <- stan_hierarchical_response_thall(
  group_responses = c(0, 0, 1, 3, 5, 0, 1, 2, 0, 0), 
  group_sizes = c(0, 2 ,1, 7, 5, 0, 2, 3, 1, 0), 
  mu_mean = -1.3863,
  mu_sd = sqrt(1 / 0.1),
  tau_alpha = 2,
  tau_beta = 20)
```

The returned object is the same type as the fits returned by rstan:

``` r
mod0
#> Inference for Stan model: ThallHierarchicalBinary.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                     mean se_mean    sd   2.5%    25%    50%    75%  97.5%
#> mu                 -0.05    0.03  1.38  -2.85  -0.89  -0.03   0.82   2.63
#> sigma2             11.96    0.33 10.35   3.35   6.21   9.15  14.15  37.99
#> rho[1]             -0.02    0.07  3.84  -7.80  -2.25   0.03   2.24   7.56
#> rho[2]             -2.96    0.06  2.49  -8.93  -4.20  -2.52  -1.26   0.59
#> rho[3]              2.28    0.05  2.67  -1.88   0.44   1.96   3.59   8.53
#> rho[4]             -0.30    0.01  0.78  -1.84  -0.82  -0.29   0.22   1.22
#> rho[5]              3.62    0.05  2.24   0.58   2.06   3.23   4.68   9.11
#> rho[6]             -0.12    0.06  3.67  -7.74  -2.31  -0.04   2.09   7.25
#> rho[7]             -0.02    0.02  1.49  -2.99  -0.96  -0.05   0.90   3.02
#> rho[8]              0.75    0.02  1.26  -1.58  -0.08   0.68   1.51   3.43
#> rho[9]             -2.32    0.06  2.74  -8.93  -3.68  -1.95  -0.48   1.88
#> rho[10]            -0.02    0.06  3.59  -7.27  -2.19  -0.02   2.13   7.01
#> sigma               3.26    0.04  1.14   1.83   2.49   3.02   3.76   6.16
#> prob_response[1]    0.50    0.01  0.38   0.00   0.09   0.51   0.90   1.00
#> prob_response[2]    0.15    0.00  0.18   0.00   0.01   0.07   0.22   0.64
#> prob_response[3]    0.77    0.00  0.26   0.13   0.61   0.88   0.97   1.00
#> prob_response[4]    0.43    0.00  0.17   0.14   0.31   0.43   0.56   0.77
#> prob_response[5]    0.92    0.00  0.10   0.64   0.89   0.96   0.99   1.00
#> prob_response[6]    0.49    0.01  0.38   0.00   0.09   0.49   0.89   1.00
#> prob_response[7]    0.49    0.00  0.26   0.05   0.28   0.49   0.71   0.95
#> prob_response[8]    0.64    0.00  0.22   0.17   0.48   0.66   0.82   0.97
#> prob_response[9]    0.23    0.00  0.26   0.00   0.02   0.13   0.38   0.87
#> prob_response[10]   0.50    0.01  0.38   0.00   0.10   0.50   0.89   1.00
#> lp__              -34.05    0.12  3.53 -42.58 -35.93 -33.62 -31.49 -28.70
#>                   n_eff Rhat
#> mu                 2222    1
#> sigma2              997    1
#> rho[1]             3188    1
#> rho[2]             1721    1
#> rho[3]             2446    1
#> rho[4]             3753    1
#> rho[5]             1963    1
#> rho[6]             3218    1
#> rho[7]             4347    1
#> rho[8]             4196    1
#> rho[9]             2109    1
#> rho[10]            3619    1
#> sigma              1018    1
#> prob_response[1]   4008    1
#> prob_response[2]   4496    1
#> prob_response[3]   4859    1
#> prob_response[4]   3822    1
#> prob_response[5]   4683    1
#> prob_response[6]   3848    1
#> prob_response[7]   4207    1
#> prob_response[8]   4644    1
#> prob_response[9]   3702    1
#> prob_response[10]  3870    1
#> lp__                866    1
#> 
#> Samples were drawn using NUTS(diag_e) at Sun Apr  7 19:34:36 2019.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

``` r
library(rstan)
library(ggplot2)

plot(mod0, pars = 'prob_response') + 
  geom_vline(xintercept = 0.3, col = 'orange', linetype = 'dashed') +
  labs(title = 'Partially-pooled analysis of response rate in 10 sarcoma subtypes')
```

![](man/figures/README-unnamed-chunk-14-1.png)<!-- -->

The hierarchical model for binary responses is developed in [its own
vignette](https://brockk.github.io/trialr/articles/HierarchicalBayesianResponse.html).

## BEBOP in PePS2

@Thall2008 introduced an extension of EffTox that allows dose-finding by
efficacy and toxicity outcomes and adjusts for covariate information.
Brock, et al. simplified the method by removing the dose-finding
components to leave a design that studies associated co-primary and
toxicity outcomes in an arbitrary number of cohorts determined by the
basline covariates. They refered to the simplifed design as BEBOP, for
*Bayesian Evaluation of Bivariate binary Outcomes with Predictive
variables*.

The investigators implement the design is a phase II trial of
pembrolizumab in non-small-cell lung cancer. A distinct feature of the
trial is the availability of predictive baseline covariates, the most
notwworthy of which is the PD-L1 tumour proportion score, shown by
@Garon2015 to be a predictive biomarker.

This example is demonstrated in the [BEBOP
vignette](https://brockk.github.io/trialr/articles/BEBOP.html).

## Installation

You can install trialr from github with:

``` r
# install.packages("devtools")
devtools::install_github("brockk/trialr")
```

If the latest CRAN build is what you seek then instead run:

``` r
install.packages("trialr")
```

It should go without saying that the CRAN release will be older than the
github version.

## Extending trialr and getting in touch

If there is a published Bayesian design you want implemented in Stan,
get in touch. Contact @brockk on github.

## References

Garon, Edward B, Naiyer a Rizvi, Rina Hui, Natasha Leighl, Ani S
Balmanoukian, Joseph Paul Eder, Amita Patnaik, et al. 2015.
“Pembrolizumab for the treatment of non-small-cell lung cancer.” The
New England Journal of Medicine 372 (21): 2018–28.
<doi:10.1056/NEJMoa1501824>.

O’Quigley, J, M Pepe, and L Fisher. 1990. “Continual reassessment
method: a practical design for phase 1 clinical trials in cancer.”
Biometrics 46 (1): 33–48. <doi:10.2307/2531628>.

Thall, Peter F., Hoang Q. Nguyen, and Elihu H. Estey. 2008.
“Patient-specific dose finding based on bivariate outcomes and
covariates.” Biometrics 64 (4): 1126–36.
<doi:10.1111/j.1541-0420.2008.01009.x>.

Thall, Peter F., J. Kyle Wathen, B. Nebiyou Bekele, Richard E. Champlin,
Laurence H. Baker, and Robert S. Benjamin. 2003. “Hierarchical Bayesian
approaches to phase II trials in diseases with multiple subtypes.”
Statistics in Medicine 22 (5): 763–80. <doi:10.1002/sim.1399>.

Thall, PF, and JD Cook. 2004. “Dose-Finding Based on Efficacy-Toxicity
Trade-Offs.” Biometrics 60 (3): 684–93.

Thall, PF, RC Herrick, HQ Nguyen, JJ Venier, and JC Norris. 2014.
“Effective sample size for computing prior hyperparameters in Bayesian
phase I-II dose-finding.” Clinical Trials 11 (6): 657–66.
<doi:10.1177/1740774514547397>.
