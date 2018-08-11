##
## Common functions for general use
##
## Esa Junttila, 2016-03-23 (originally)


#' Bound the values of a vector.
#' 
#' Return a new vector where the values of \code{x} are bounded by given
#' limits and otherwise kept intact.
#' @param x Numeric vector whose values we wish to bound.
#' @param low Lower numeric limit.
#' @param high Higher numeric limit.
#' @return bounded vector where for each \code{x_i} in \code{i=1,2,...} we have
#'   \itemize{
#'     \item \code{if low <= x_i <= high}, then \code{x_i}
#'     \item \code{if x_i < low}, then \code{low}
#'     \item \code{if high < x_i}, then \code{high}.
#'   }
#' @export
limit <- function(x, low, high) pmax(low, pmin(high, x))


#' Convert radian angles to degrees.
#' 
#' @param rad angles as numeric radians within \code{[0;2*pi[}
#' @return angles as numeric degrees within \code{[0;360[}
#' @seealso \code{\link{degreeToRadian}}
#' @export
radianToDegree <- function(rad) rad*(360/(2*pi))


#' Convert degree angles to radians.
#' 
#' @param degree angles as numeric degrees within \code{[0;360[}
#' @return angles as numeric radians within \code{[0;2*pi[}
#' @seealso \code{\link{radianToDegree}}
#' @export
degreeToRadian <- function(degree) degree*(2*pi/360)

#' Binary search
#' 
#' Perform a binary search on an increasing function until
#' target value is within an absolute tolerance. The goal is to
#' find such a value \code{x} within search window \code{[low,high]}
#' that \code{y:=func(x)} is closer to \code{target} than
#' given absolute tolerance. The method has logarithmic time complexity.
#' 
#' @param func Function that accepts a single numeric argument within
#'             \code{[low,high]} and is \strong{increasing} in that range.
#' @param target Target numeric value for \code{y:=func(x)}.
#' @param low Minimum numeric value of search window for \code{x}.
#' @param high Maximum numeric value of search window for \code{x}.
#' @param tol [optional] Absolute numeric tolerance; when abolute error is
#'   smaller than tolerance, stop search and return a solution.
#'   Default is \code{1e-12} which may need fine-tuning if function has
#'   very small or large derivatives.
#' @param max.iter [optional] Maximum number of iterations until we give up.
#'   Default is 200 (splits).
#' @return
#'   Return value \code{x} within \code{[low,high]} such that
#'   \code{|func(x) - target| < tol}.
#'   
#'   Return \code{NULL} if target cannot be found (because target value
#'   is not within the search window or function is not increasing).
#' @examples
#' binarySearch(function(x) 2*x - 4, 0, -1, +3)
#' binarySearch(function(x) x^2 - 2, 3, 0, +3)
#' binarySearch(sqrt, 0, 0, +1, tol=1e-6)
#' binarySearch(sqrt, 0, 0, +1)
#' f <- function(x) x^3 + 2*x^2 - 8*x
#' binarySearch(f, 5, -4, +3)  # non-increasing, ambiguous
#' binarySearch(f, 5, -2, +1)  # err: decreasing function, no solution
#' binarySearch(f, -3, 2, 3)  # err: solution not found in window
#' binarySearch(sqrt, 0, 0, +1, tol=1e-18)  # err: floating point precision
#' @export
#'
binarySearch <- function(func, target, low, high, tol=1e-12, max.iter=100) {
  a <- low
  b <- high
  iter <- 1
  while (a <= b && iter <= max.iter) {
      mid <- (a + b) / 2
      func.val <- func(mid)
      if (abs(func.val - target) < tol) { return(mid) }
      else if (a == b) { break }
      else if (func.val > target) { b <- mid }
      else if (func.val < target) { a <- mid }
      iter <- iter + 1
  }
  return(NULL)
}


# Misunderstood an article. Obsolete.
# Bessel function of the first kind.
# Approximately correct. It's hard to deduce a good number of terms, now 41.
# Highly restricted search area. The function is monotonic at least in range (-1.84, 1.84).
#besselJ <- function(alpha, x) sum(sapply(0:40, function(k) (-x^2/4)^k/factorial(k)/gamma(alpha+k+1)*(x/2)^alpha))
#invertedBesselJ <- function(alpha, y, low=-1.84, high=1.84) binarySearch(function(x) besselJ(alpha, x), y, low, high, tol=1e-6)


#' @rdname A1
#' @name A1
#' @title Ratio of Bessel functions of the first kind.
#' @details These functions are related to circularity (estimating
#'   maximum-likelihood on Von Mises distribution),
#'   but I'm not sure about details.
NULL


#' @rdname A1
#' @title Ratio of Bessel functions
#' @description \code{A1} evaluates the first and zeroth order Bessel
#'   functions of the first kind at a specified non-negative real number,
#'   and returns the ratio.
#' @param kappa non-negative concentration arg in Von Mises distribution (?)
#' @return \code{A1} returns ratio of evaluated Bessel functions.
#' @seealso
#'   \url{https://cran.r-project.org/web/packages/circular/circular.pdf}
#'   
#'   \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1.R?view=markup&root=circular}
#'   
A1 <- function(kappa) {
  result <- besselI(kappa, nu=1, expon.scaled = TRUE)/besselI(kappa, nu=0, expon.scaled = TRUE)
  return(result)
}


#' @rdname A1
#' @title Inverse of ratio of Bessel functions.
#' @description \code{A1inv} evaluates the inverse function of \code{A1}.
#'   This function is used to compute the maximum likelihood estimate
#'   of the concentration parameter of a von Mises distribution.
#' @param x numeric distribution argument(?)
#' @return \code{A1inv} returns the estimated kappa (Von Mises concentration).
#' @seealso
#'   \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1inv.R?view=markup&root=circular&pathrev=6}
A1inv <- function(x) {
  ifelse (0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5)/6,
          ifelse (x < 0.85, -0.4 + 1.39 * x + 0.43/(1 - x), 1/(x^3 - 4 * x^2 + 3 * x)))
}


#' Sample a single value from a vector.
#' 
#' Sample a single value from a vector. Sampling takes into account
#' element-wise probabilities or frequencies; if not given, a uniform
#' distribution is used as default.
#'
#' This is a simplistic way of sampling a single number from vector x.
#' For some reason this is faster than using \code{sample.int}
#' -- perhaps it's designed for vectorized use.
#' @param x Vector of elements
#' @param prob Probabilities for sampling each item. Vector \code{prob}
#'   should have same length as \code{x}. Probabilities are normalized,
#'   so frequencies can be used as well.
#' @return sampled element from a vector, according to probabilities
#' @examples
#' sample.one(2) == 2
#' all(abs(table(replicate(3000, sample.one(2:4))) - 1000) < 100)
#' set.seed(1)
#' sample.one(seq(1.23, 2.34, 0.001)) == 1.524
#' is.null(sample.one(c()))
#' tryCatch({sample.one(1:3, c(0.6,0.4)); write("error not catched: x and prob mismatch", stderr())}, error=function(x) TRUE)
#' set.seed(1)
#' all(replicate(20, sample.one(2:5, prob=c(0.1,0.4,0.2,0.3))) == c(3,3,4,5,3,5,5,4,4,2,3,3,4,3,5,3,5,5,3,5))
#' all(abs(table(replicate(10000, sample.one(2:5, prob=c(0.1,0.4,0.2,0.3)))) - 1000*c(1,4,2,3)) < 200)
#' sample.one(0:4, prob=c(0, 1e-12, 0, 1 - 1e-12, 0)) == 3
#' @export
sample.one <- function(x, prob=NULL) {
  if (is.null(prob)) {
    return(x[runif(1) %/% (1 / length(x)) + 1])  # uniform distribution
  } else {
    stopifnot(length(prob) == length(x))
    return(x[findInterval(runif(1), cumsum(prob/sum(prob))) + 1])
  }
}


# Kmeans++ method for center initialization.
# The result may contain NA's if all points are in cluster centers.
initializeCenters <- function(points, k) {
  stopifnot(as.integer(k) == k && length(k) == 1 && k >= 1 && k <= nrow(points))
  n <- nrow(points)
  range <- 1:n                                # point-indices that are available as centers
  pr <- as.numeric(rep(1/n, n))               # select first center-point by uniform sampling
  center.point.idx <- as.integer(rep(NA, k))  # list of center-indices initially empty
  dist.min <- rep(Inf, n)                     # infinite distances to non-existent centers
  total.dist <- Inf
  if (k >= 2) {
    for (i in 1:(k-1)) {                      # add new centers one by one
      c.idx <- sample.one(range, prob=pr)  # sample next center-index from probabilities
      center.point.idx[i] <- c.idx            # append chosen center-index
      c.point <- points[c.idx, ]
      dist.to.center <- rowSums(sweep(points, 2, c.point)^2)  # squared Euclidean distance to new center
      dist.min <- pmin(dist.min, dist.to.center)  # update minimum distances from points to any center
      total.dist <- sum(dist.min)
      if (total.dist == 0) break           # all points are in cluster centers
      pr <- dist.min / total.dist          # update probabilities for next selection
    }
  }
  if (total.dist != 0) center.point.idx[k] <- sample.one(range, prob=pr)  # add last center
  return(center.point.idx)  # return a vector of point indices that refer to chosen centers
}



#' K-means++ clustering algorithm.
#' 
#' K-means++ clustering algorithm, using custom initial seeding,
#' random restarts, and base \code{\link{kmeans}} algorithm.
#' 
#' The K-means++ method was first introduced by David Arthur and
#' Sergei Vassilvitskii in 2007.
#' @param x Data with columns as coordinates and rows as points
#' @param k Number of clusters
#' @param iter.max [optional] Maximum number of iterations; default is 10.
#' @param restarts [optional] Number of separate restarts. The best result
#'   will be used. Default is one.
#' @param ... [optional] Other args passed to standard \code{kmeans} algorithm.
#' @examples
#' x <- matrix(c(
#'    1,2,3,
#'    3,6,-2,
#'    2,2,2,
#'    7,6,5,
#'    20,22,27,
#'    15,16,17,
#'    30,12,20,
#'    50,40,50,
#'    29,40,45,
#'    45,50,38), byrow=TRUE, ncol=3)
#' res <- kmeanspp(x, 3, iter.max=10, restarts=1)
#' y <- res$cluster
#' length(unique(y[1:4])) == 1 & length(unique(y[5:7])) == 1 & length(unique(y[8:10])) == 1
#' @seealso Inspiration from \code{mahito-sugiyama/k-meansp2.R}
#'
#'   K-means++ from Wikipedia.
#'   
#'   \code{\link{kmeans}} base function.
#'   
#'   \url{https://gist.github.com/mahito-sugiyama/ef54a3b17fff4629f106}
#' @export
kmeanspp <- function(x, k, iter.max=10, restarts=1, ...) {
  n <- nrow(x)  # number of data points
  res.best <- list(tot.withinss=Inf)  # the best result among restarts
  for (rerun in 1:restarts) {
    centers.idx <- initializeCenters(x, k)
    ## Perform k-means clustering with the obtained centers:
    C <- as.integer(na.omit(centers.idx))
    init.centers <- x[C, , drop=FALSE]
    # Note: kmeans bug in old R versions when k=1: "Error: number of cluster centres must lie between 1 and nrow(x)"
    res <- kmeans(x, init.centers, iter.max=iter.max, nstart=1, ...)
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      # add extra fields to clustering results
      res$initial.centers <- init.centers
      res$k <- nrow(init.centers)  # may be different than argument 'k'
      res.best <- res
    }
  }
  return(res.best)
}
kmeansppOLD <- function(x, k, iter.max=10, restarts=1, ...) {
  n <- nrow(x)  # number of data points
  # Allocate distances: [i,j] --> distance between point x[i, ] and center x[centers[j], ]
  distances <- matrix(numeric(n*(k-1)), ncol=k-1) * NA
  res.best <- list(tot.withinss=Inf)  # the best result among restarts

  for (rerun in 1:restarts) {
    # K-means++ initialization:
    centers <- integer(k) * NA  # row IDs of centers
    pr <- rep(1/n, times=n)  # uniform probabilities for sampling first center
    if (k >= 2) {
      for (i in 1:(k - 1)) {
        centers[i] <- sample.int(n, size=1, prob = pr) # pick the ith center
        # Compute distances from points to center i; use squared Euclidean (no effect):
        distances[, i] <- rowSums((x - rep(x[centers[i], ], each=n))^2)
        # Compute probability for the next sampling
        freq <- apply(distances[, 1:i, drop=FALSE], 1, min)  # FIXME, make somehow faster!
        if (sum(freq) == 0) break
        pr <- freq / sum(freq)  # probabilities for the next sampling
      }
    }
    if (k == 1 || sum(freq) != 0) centers[k] <- sample.int(n, size=1, prob=pr)  # pick the last (kth) center
    ## Perform k-means clustering with the obtained centers:
    C <- as.integer(na.omit(centers))
    init.centers <- x[C, , drop=FALSE]
    res <- kmeans(x, init.centers, iter.max=iter.max, nstart=1, ...)
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      res$initial.centers <- init.centers
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


## Evidence about R's kmeans bug which still existed in R 3.2.1, bur was fixed until 3.4.3.
## When only one cluster (k=1), passing the initial cluster center (matrix with one row) fails.
## Documented also here: https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=16623
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
#res <- kmeans(x, 3, iter.max=10, nstart=1)
#res <- kmeans(x, x[sample(1:nrow(x), 2), ], iter.max=10, nstart=1)
#res <- kmeans(x, 1, iter.max=10, nstart=1)
#tryCatch(
#  {kmeans(x, x[sample(1:nrow(x), 1), , drop=FALSE], iter.max=10, nstart=1); print('Fixed already?!')},
#  error = function(x) {}
#)
## Error: number of cluster centres must lie between 1 and nrow(x)


#m <- matrix(c(27,12,9,13,8,3,2,25,16,30,11,18,14,4,23,5,10,19,1,24,29,7,17,6,21,26,20,22,15,28), ncol=6)
#clust <- kmeanspp(m, 4)

