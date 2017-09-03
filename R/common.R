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


# K-means++ clustering algorithm
# Inspiration from mahito-sugiyama/k-meansp2.R
# https://gist.github.com/mahito-sugiyama/ef54a3b17fff4629f106
kmeanspp <- function(x, k, iter.max=10, restarts=1, ...) {
  n <- nrow(x)  # number of data points
  centers <- integer(k) * NA  # row IDs of centers
  # Allocate distances: [i,j] --> distance between point x[i, ] and center x[centers[j], ]
  distances <- matrix(numeric(n*(k-1)), ncol=k-1) * NA
  res.best <- list(tot.withinss=Inf)  # the best result among restarts

  for (rerun in 1:restarts) {
    # K-means++ initialization:
    pr <- rep(1/n, times=n)  # uniform probabilities for sampling first center
    for (i in 1:(k - 1)) {
      centers[i] <- sample.int(n, size=1, prob = pr) # pick the ith center
      # Compute distances from points to center i; use squared Euclidean (no effect):
      distances[, i] <- rowSums((x - rep(x[centers[i], ], each=n))^2)

      # Compute probability for the next sampling
      freq <- apply(distances[, 1:i, drop=FALSE], 1, min)
      pr <- freq / sum(freq)  # probabilities for the next sampling
    }
    centers[k] <- sample.int(n, size=1, prob=pr)  # pick the last (kth) center
    ## Perform k-means clusterin with the obtained centers:
    res <- kmeans(x, x[centers, ], iter.max=iter.max, nstart=1, ...)
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      res$initial.centers <- x[centers, ]
      res.best <- res
    }
  }
  return(res.best)
}

## Test clustering:
#x <- matrix(c(
#  1,2,3,
#  3,6,-2,
#  2,2,2,
#  7,6,5,
#  20,22,27,
#  15,16,17,
#  30,12,20,
#  50,40,50,
#  29,40,45,
#  45,50,38
#), byrow=TRUE, ncol=3)
#res <- kmeanspp(x, 3, iter.max=10, restarts=1)
#y <- res$cluster
#stopifnot(length(unique(y[1:4])) == 1 & length(unique(y[5:7])) == 1 & length(unique(y[8:10])) == 1)

## Speed test
#R <- (4000*3000) / 20 #50000  # 200*3000=600000
#C <- 3 #10
#m <- matrix(rnorm(R*C), ncol=C)
#s <- Sys.time()
#kmeanspp(m, 7, iter.max=20, restarts=10)
#e <- Sys.time()
#print(e-s)

