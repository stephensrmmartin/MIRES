
# MIRES

<!-- badges: start -->
<!-- badges: end -->

MIRES is a package for *M*easurement *I*nvariance assessment using *R*andom *E*ffects and *S*hrinkage.
This is the official implementation of the approach taken in [Measurement Invariance Assessment with Bayesian Hierarchical Inclusion Modeling](https://osf.io/preprints/psyarxiv/qbdjt_v1).

Unlike most other MI assessment methods, MIRES assumes all measurement parameters (loadings, intercepts, and residual SDs) can randomly vary across groups.
Therefore, it is a fully mixed effects CFA model.
If invariance is tenable, then the random effect SDs should all be effectively zero.

Unlike most mixed effects models, MIRES further regularizes the random effect SDs (RE-SDs) via an inclusion model.
In sum, it includes an adaptive, dependent regularization model that shares regularization intensity across several RE-SDs.
This allows the RE-SDs to rapidly shrink toward zero if they need to (i.e., if invariance is met), and to be largely untouched if they need to be non-zero (i.e., if invariance is not met).
Consequently, MIRES can gain evidence in favor of invariance quicker than a mixed effects model otherwise could, especially with few groups.

Altogether, MIRES is a mixed effects CFA model, where the RE-SDs are also hierarchically modeled with a dependent regularizing prior (See Details for more information).
When invariance is plausible, it will collapse into a fixed effects CFA model; when it is not, then a partially invariant or fully non-invariant model will be produced.


## Installation

You can install the development version of MIRES from [github](https://github.com/stephenSRMMartin/MIRES) with:

``` r
remotes::install_github("stephenSRMMartin/MIRES")
```

or from CRAN with:
``` r
install.packages("MIRES")
```

## Example

Simulated data are provided in the package.

This data has non-invariance in the loadings, and invariance elsewhere.

``` r
library(MIRES)
data(sim_loadings)
head(sim_loadings)
```

```
         x_1        x_2        x_3          x_4        x_5        x_6
1  1.7400023  1.3214588  0.6255326  0.901031733  0.9495002  1.2286485
2 -0.5640201  0.2203121 -1.4399299 -0.571124405 -0.9580103 -0.2434576
3 -1.1030768 -0.4088069  0.7182129 -0.009960373 -1.3954430 -0.1966013
4 -1.6139970 -1.3411869 -2.7745400 -0.818880762 -2.1656676 -1.1999327
5 -0.8192840 -0.9715812 -0.3293007 -0.538967762 -0.7168703 -0.4190734
6  0.3266671  0.8955085  2.4646536  0.596746521  1.5468213  0.4509530
         x_7        x_8 group
1  0.3128417  1.5761853     1
2 -1.1585096 -1.1884185     1
3 -0.6639108 -0.6313458     1
4 -1.4482666 -2.3015851     1
5 -2.7528062 -1.0865613     1
6  3.1724797  2.1475378     1
```

Fit and summarize the model.

``` r
fit <- mires(my_factor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8, group, sim_loadings, iter = 1000)
summary(fit)
```

```
-----
MIRES model object
-----
Specification:
my_factor ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8

-----
Configuration:
	 Factors: Univariate
	 Latent Mean and Variance Identification: Sum-to-zero, product-to-one.
	 Hierarchical inclusion model: Yes
	 Latent Scores Saved: No
	 Inclusion Model Prior Params: 0 ,  0.25
-----
Fixed Effects
-----
Loadings
 Item  Mean Median    SD   L95   U95  Rhat
  x_1 0.757  0.750 0.139 0.502 1.045 1.001
  x_2 0.668  0.656 0.221 0.276 1.154 1.001
  x_3 0.750  0.739 0.143 0.487 1.070 1.000
  x_4 0.579  0.573 0.204 0.193 1.001 1.007
  x_5 0.707  0.693 0.232 0.292 1.215 1.003
  x_6 0.791  0.786 0.229 0.288 1.218 1.005
  x_7 0.594  0.580 0.229 0.159 1.056 1.003
  x_8 1.009  1.003 0.211 0.587 1.429 1.002

Residual Standard Deviations (Unlogged Scale)
 Item  Mean Median    SD   L95   U95  Rhat
  x_1 0.754  0.750 0.061 0.652 0.891 1.001
  x_2 0.743  0.740 0.052 0.633 0.839 1.001
  x_3 0.750  0.747 0.051 0.653 0.854 1.004
  x_4 0.695  0.692 0.051 0.596 0.791 1.001
  x_5 0.687  0.685 0.051 0.590 0.788 1.001
  x_6 0.710  0.709 0.057 0.594 0.819 1.002
  x_7 0.720  0.715 0.051 0.624 0.819 1.005
  x_8 0.788  0.785 0.064 0.679 0.932 1.001

Intercepts
 Item   Mean Median    SD    L95   U95  Rhat
  x_1 -0.040 -0.040 0.083 -0.201 0.121 1.002
  x_2 -0.029 -0.028 0.078 -0.190 0.114 1.003
  x_3 -0.033 -0.035 0.086 -0.198 0.134 1.002
  x_4  0.008  0.008 0.067 -0.112 0.152 1.001
  x_5 -0.107 -0.107 0.079 -0.252 0.057 1.002
  x_6 -0.014 -0.008 0.098 -0.222 0.160 1.005
  x_7 -0.038 -0.037 0.091 -0.219 0.136 1.003
  x_8 -0.025 -0.027 0.092 -0.202 0.156 1.002

-----
Measurement Invariance Assessment
-----
Random Effect Standard Deviations (Unlogged Scale)
   Parameter Item    Factor  Mean Median    SD   L95   U95  Rhat    BF01    BF10 Pr(SD <= 0.1| D) BF(SD <= 0.1)
     Loading  x_1 my_factor 0.219  0.192 0.156 0.000 0.504 1.008   2.107   0.475            0.232         3.046
     Loading  x_2 my_factor 0.437  0.401 0.181 0.182 0.768 1.002   0.001 729.233            0.000         0.000
     Loading  x_3 my_factor 0.225  0.192 0.162 0.000 0.523 1.003   1.933   0.517            0.216         2.777
     Loading  x_4 my_factor 0.388  0.352 0.180 0.104 0.767 1.002   0.017  58.707            0.007         0.071
     Loading  x_5 my_factor 0.460  0.425 0.191 0.161 0.836 1.002   0.004 244.913            0.001         0.010
     Loading  x_6 my_factor 0.473  0.429 0.207 0.146 0.903 1.003   0.023  42.880            0.002         0.025
     Loading  x_7 my_factor 0.467  0.430 0.185 0.173 0.827 1.001   0.002 480.630            0.000         0.000
     Loading  x_8 my_factor 0.407  0.375 0.193 0.076 0.798 1.003   0.131   7.628            0.021         0.222
 Residual SD  x_1           0.114  0.093 0.095 0.000 0.289 1.003   6.617   0.151            0.529        11.355
 Residual SD  x_2           0.092  0.073 0.076 0.000 0.236 1.001   7.097   0.141            0.641        18.052
 Residual SD  x_3           0.079  0.063 0.068 0.000 0.209 1.005   8.589   0.116            0.712        24.934
 Residual SD  x_4           0.083  0.067 0.074 0.000 0.218 0.999   9.658   0.104            0.697        23.202
 Residual SD  x_5           0.090  0.071 0.077 0.000 0.236 1.003   8.955   0.112            0.650        18.776
 Residual SD  x_6           0.100  0.077 0.089 0.000 0.264 1.000   6.815   0.147            0.617        16.287
 Residual SD  x_7           0.087  0.068 0.078 0.000 0.235 1.005 357.644   0.003            0.676        21.094
 Residual SD  x_8           0.084  0.066 0.076 0.000 0.225 1.002   9.266   0.108            0.688        22.294
   Intercept  x_1           0.094  0.072 0.089 0.000 0.258 1.001   8.213   0.122            0.647        18.571
   Intercept  x_2           0.091  0.074 0.077 0.000 0.238 1.003   8.168   0.122            0.636        17.703
   Intercept  x_3           0.096  0.077 0.079 0.000 0.245 1.003   6.632   0.151            0.622        16.601
   Intercept  x_4           0.065  0.050 0.061 0.000 0.173 0.999  13.225   0.076            0.803        41.211
   Intercept  x_5           0.092  0.074 0.078 0.000 0.233 1.004   7.707   0.130            0.644        18.329
   Intercept  x_6           0.126  0.104 0.100 0.000 0.316 1.001   5.705   0.175            0.482         9.426
   Intercept  x_7           0.139  0.121 0.100 0.000 0.328 1.001   2.978   0.336            0.396         6.629
   Intercept  x_8           0.093  0.072 0.085 0.000 0.251 1.001 173.931   0.006            0.642        18.131
```

# Limitations

- Currently only supports one-factor models. Multidimensional models will soon be supported.
- Assumes normality in the indicators. Alternative likelihoods are possible, but not yet implemented.
- Does not allow residual covariance.
- Does not yet support unregularized mixed models. Unregularized (i.e., constant, and independent) priors will soon be supported.

# Details
With three kinds of parameters (loadings, residual SDs, and intercepts), and J items, there are 3J RE-SDs: $\sigma_p, p \in [0,\ldots,3J]$.
The regularizing prior for each RE-SD, $\sigma_p$ is defined as:
$$
\sigma_p \sim \mathcal{N}^+(0, \exp(\tau_p))
$$
If $\exp(\tau_p)$ is very small, then the prior probability that $\sigma_p = 0$ increases greatly, and the parameter corresponding to it is reduced to a fixed effect.

The probability that one parameter is invariant should be informed by other parameters with similar characteristics (e.g., they model the same item, or they are of the same type --- Loading, intercept, or residual SD).
Therefore, the regularization term for each $\sigma_p$ can be dependent on others.
We model the regularization terms $\tau_p$ as:
$$
\tau_p = \tau_0 + \tau_{\text{item}_p} + \tau_{\text{param}_p} + \lambda_p
$$
Therefore, there is a global scaling factor $\tau_0$, an "item" effect ($\tau_{\text{item}, p}), a "parameter" effect ($\tau_{\text{param}, p}), and a unique effect $\lambda_p$.

Therefore, if all parameters are invariant, $\tau_0$ can decrease, and all RE-SDs can jointly approach zero.
Conversely, if some loadings are invariant, then $\tau_{\lambda}$ decreases, and the probability that other loadings are zero increases.
If item $j$'s intercept and loading seems to vary, then $\tau_j$ increases, and the probability that the item's residual SD is invariant decreases.
Finally, the $\lambda_p$ terms provide an "escape hatch" of sorts, if a particular parameter uniquely needs to be invariant or non-invariant.

