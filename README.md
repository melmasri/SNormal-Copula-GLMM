# SNormal-Copula-GLMM
An R code implementation of a skew-normal copula-driven generalized linear mixed model (GLMM)

Start with the `main.R` script and read the inline documentations.

## This simplified code is these step:
  * There is a universal `xi` for all units, thus one a universal `Sigma` and also one `Psi`.
  * All units have the same number of observations, that is only `obs`, which is not indexed as `obs[i]`.
  * There is a single skewness vector `lambda` for all units, which is the same size as the number of observations. That also implies a one $\delta$ vector.
  * In the `random.data.R` file, rather then generating multiple skew-normals `Z`'s, one for each unit, we generate one $Z$, with `E[Z]= 0 + delta*sqrt(pi/2)`, where `Omega = Sigma` and `alpha = lambda` in  the `sn` `R` library.
  * Given that there is one `Z`, we have `Z[i] = Z + b_o[i]` for each unit. 

