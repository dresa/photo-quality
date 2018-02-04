# Test kmeanspp speed, in particular the pp part: initialization of centers
# Esa Junttila
# 2018-02-03

n.points <- as.integer(10000)
n.dims <- as.integer(3)
n.centers <- as.integer(30)
stopifnot(n.centers <= n.points)

points <- matrix(sample(n.points * n.dims), nrow=n.points)

# Method 1: Slow and full
initCenters1 <- function(points, k) {
  n <- nrow(points)
  distances <- matrix(NA, nrow=n, ncol=k)
  pr <- as.numeric(rep(1/n, n))
  centers <- as.integer(rep(NA, k))
  for (i in 1:k) {
    center.idx <- sample(1:n, 1, prob=pr)
    center.point <- points[center.idx, ]  # dropping extra dimension by default
    centers[i] <- center.idx
    distances[ , i] <- rowSums((points - rep(center.point, each=n))^2)
    freq <- apply(distances[ , 1:i, drop=FALSE], 1, min)
    pr <- freq / sum(freq)
  }
  return(centers)
}

# Method 2: Faster last assignment, still very slow
initCenters2 <- function(points, k) {
  n <- nrow(points)
  distances <- matrix(NA, nrow=n, ncol=k)
  pr <- as.numeric(rep(1/n, n))
  centers <- as.integer(rep(NA, k))
  for (i in 1:(k-1)) {
    center.idx <- sample(1:n, 1, prob=pr)
    center.point <- points[center.idx, ]  # dropping extra dimension by default
    centers[i] <- center.idx
    distances[ , i] <- rowSums((points - rep(center.point, each=n))^2)
    freq <- apply(distances[ , 1:i, drop=FALSE], 1, min)
    pr <- freq / sum(freq)
  }
  centers[k] <- sample(1:n, 1, prob=pr)
  return(centers)
}

# Method 3: Faster minimum distances
initCenters3 <- function(points, k) {
  n <- nrow(points)
  pr <- as.numeric(rep(1/n, n))
  centers <- as.integer(rep(NA, k))
  dist.min <- rep(Inf, n)
  for (i in 1:(k-1)) {
    center.idx <- sample(1:n, 1, prob=pr)
    centers[i] <- center.idx
    center.point <- points[center.idx, ]  # dropping extra dimension by default
    dist.to.center <- rowSums((points - rep(center.point, each=n))^2)
    dist.min <- pmin(dist.min, dist.to.center)
    pr <- dist.min / sum(dist.min)
  }
  centers[k] <- sample(1:n, 1, prob=pr)
  return(centers)
}


# Method 4: Faster minimum distances
# A simplistic way of sampling a single number from vector x.
# Faster than using sample.int -- perhaps it's designed for vectorized use.
sample.one <- function(x, prob=NULL) {
  if (is.null(prob)) {
    return(x[runif(1) %/% (1 / length(x)) + 1])  # uniform distribution
  } else {
    stopifnot(length(prob) == length(x))
    return(x[findInterval(runif(1), cumsum(prob/sum(prob))) + 1])
  }
}
## Tests:
#stopifnot(sample.one(2) == 2)
#stopifnot(all(abs(table(replicate(3000, sample.one(2:4))) - 1000) < 100))
#set.seed(1)
#stopifnot(sample.one(seq(1.23, 2.34, 0.001)) == 1.524)
#stopifnot(is.null(sample.one(c())))
#tryCatch({sample.one(1:3, c(0.6,0.4)); write("error not catched: x and prob mismatch", stderr())}, error=function(x) {})
#set.seed(1)
#stopifnot(replicate(20, sample.one(2:5, prob=c(0.1,0.4,0.2,0.3))) == c(3,3,4,5,3,5,5,4,4,2,3,3,4,3,5,3,5,5,3,5))
#stopifnot(abs(table(replicate(10000, sample.one(2:5, prob=c(0.1,0.4,0.2,0.3)))) - 1000*c(1,4,2,3)) < 200)
#stopifnot(sample.one(0:4, prob=c(0, 1e-12, 0, 1 - 1e-12, 0)) == 3)

initCenters4 <- function(points, k) {
  stopifnot(as.integer(k) == k && length(k) == 1 && k >= 1 && k <= nrow(points))
  n <- nrow(points)
  range <- 1:n                                # point-indices that are available as centers
  pr <- as.numeric(rep(1/n, n))               # select first center-point by uniform sampling
  center.point.idx <- as.integer(rep(NA, k))  # list of center-indices initially empty
  dist.min <- rep(Inf, n)                     # infinite distances to non-existent centers
  if (k >= 2) {
    for (i in 1:(k-1)) {                      # add new centers one by one
      c.idx <- sample.one(range, prob=pr)  # sample next center-index from probabilities
      center.point.idx[i] <- c.idx            # append chosen center-index
      c.point <- points[c.idx, ]
      dist.to.center <- rowSums(sweep(points, 2, c.point)^2)  # squared Euclidean distance to new center
      dist.min <- pmin(dist.min, dist.to.center)  # update minimum distances from points to any center
      pr <- dist.min / sum(dist.min)          # update probabilities for next selection
    }
  }
  center.point.idx[k] <- sample.one(range, prob=pr)  # choose the last center
  return(center.point.idx)  # return a vector of point indices that refer to chosen centers
}


## Main test
t0 <- proc.time()
set.seed(1)
centers1 <- initCenters1(points, n.centers)
t1 <- proc.time()
set.seed(1)
centers2 <- initCenters2(points, n.centers)
t2 <- proc.time()
set.seed(1)
centers3 <- initCenters3(points, n.centers)
t3 <- proc.time()
set.seed(1)
centers4 <- initCenters4(points, n.centers)
t4 <- proc.time()
stopifnot(all(centers1 == centers2) && all(centers1 == centers3))
# The centers chosen in Method 4 differ from other Methods,
# because "sample.one" function deviates from "sample.int" function.

Rprof('Rprof.out', line.profiling=TRUE, interval=0.01)
replicate(1000, initCenters4(points, n.centers))
Rprof(NULL)

print(t1-t0)
print(t2-t1)
print(t3-t2)
print(t4-t3)

print(summaryRprof('Rprof.out', lines='show'))


# Results, with 100000 rows, 3 dimensions, 100 centers
#            user  system elapsed
# Method 1: 40.79    0.82   41.64 
# Method 2: 39.45    0.59   40.21 
# Method 3:  3.14    0.09    3.23 
# Method 4:  0.83    0.02    0.84 

