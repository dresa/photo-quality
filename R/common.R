##
## Common functions for general use
##
## Esa Junttila, 2016-03-23


limit <- function(x, low, high) pmax(low, pmin(high, x))

radianToDegree <- function(rad) rad*(360/(2*pi))

degreeToRadian <- function(degree) degree*(2*pi/360)

binarySearch <- function(func, target, low, high, tol=1e-12) {
  if ( high < low ) { return(NULL) }
  else {
    mid <- (low + high) / 2
    func.val <- func(mid)
    if (abs(func.val - target) < tol) { return(mid) }
    else if (func.val > target) { binarySearch(func, target, low, mid) }
    else if (func.val < target) { binarySearch(func, target, mid, high) }
  }
}

# Misunderstood an article. Obsolete.
# Bessel function of the first kind.
# Approximately correct. It's hard to deduce a good number of terms, now 41.
# Highly restricted search area. The function is monotonic at least in range (-1.84, 1.84).
#besselJ <- function(alpha, x) sum(sapply(0:40, function(k) (-x^2/4)^k/factorial(k)/gamma(alpha+k+1)*(x/2)^alpha))
#invertedBesselJ <- function(alpha, y, low=-1.84, high=1.84) binarySearch(function(x) besselJ(alpha, x), y, low, high, tol=1e-6)

#Evaluates the first and zeroth order Bessel functions of the first kind at a specified non-negative real
#number, and returns the ratio. [Related to circularity]
#See: https://cran.r-project.org/web/packages/circular/circular.pdf
#https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1.R?view=markup&root=circular
A1 <- function(kappa) {
  result <- besselI(kappa, nu=1, expon.scaled = TRUE)/besselI(kappa, nu=0, expon.scaled = TRUE)
  return(result)
}
#Inverse function of the ratio of the first and zeroth order Bessel functions of the first kind.  This
#function is used to compute the maximum likelihood estimate of the concentration parameter of a
#von Mises distribution.
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1inv.R?view=markup&root=circular&pathrev=6
A1inv <- function(x) {
  ifelse (0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5)/6,
          ifelse (x < 0.85, -0.4 + 1.39 * x + 0.43/(1 - x), 1/(x^3 - 4 * x^2 + 3 * x)))
}

