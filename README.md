
# biLocPol

<!-- badges: start -->
<!-- badges: end -->

This package features the implementation of the bivariate local
polynomial estimator for covariance kernels up to the second order as
proposed in “Optimal rates for estimating the covariance kernel from
synchronously sampled functional data” by Hajo Holzmann and Max Berger
[arxiv](https://arxiv.org/abs/2407.13641). The package can be installed
from from [GitHub](https://github.com/) via:

``` r
# install.packages("devtools")
devtools::install_github("mbrgr/biLocPol")
# library(biLocPol)
```

The examples in the previous mentiones paper can be found in the
seperate GitHub project
[Optimal-Rates-Covariance-Estimation-in-FDA](https://github.com/mbrgr/Optimal-Rates-Covariance-Kernel-Estimation-in-FDA).

## Structure

The folder “R” contains the .r files. The script
“computational_functions.r” features basic functions for data
transformations, calculations and functions special to the method. The
file “local_polynomial_weights.r” contains the heart of the package. The
function “weight_point” is only a helper function for the main function
“local_polynomial_weights.r” which calculates the weights of the local
polynomial estimator (and its derivatives if demanded). The evaluation
of the estimator is the done by plugging the weights and the transformed
data into the “weights_eval” function. The “cross_validation.r” file
features functions for K-fold and leave one plane out cross validation.

## Minimal example

``` r
set.seed(22)
n = 200
p = 50
x.design = (1:p - 0.5)/p
Y = FDA_observation(n, x.design) # simulation with brownian motion 
p.eval = 50
h = 0.5
m = 2
del = 1

weights = local_polynomial_weights(p, h, p.eval, parallel = F, m, del)
Z = observation_transformation(Y)

estimation = eval_weights(weights, Z)
estimation[,,1] # covariance kernel
plotly::plot_ly(x = ~weights$x.eval, y = ~weights$x.eval, z = ~estimation[,,1]) |> 
                  plotly::add_surface()
estimation[,,2] # derivative in first direction
plotly::plot_ly(x = ~weights$x.eval, y = ~weights$x.eval, z = ~estimation[,,2]) |> 
  plotly::add_surface()
estimation[,,3] # derivative in second direction
plotly::plot_ly(x = ~weights$x.eval, y = ~weights$x.eval, z = ~estimation[,,3]) |> 
  plotly::add_surface()
# Note that the covariance kernel of the brownian motion is not differentiable at the diagonal
```

## Parallel implementation

Note that the calculation of the weights and the simulations can be
paralallized with the “future.apply” and “future” package. Take care
where the parallelization takes place.

## Further Notes

Till now only equidistant grids are allowed. Further grids shall be
implemented soon. The represantation of the calculated weights by
“local_polynomial_weights” is ugly. This shall be fixed soon.
