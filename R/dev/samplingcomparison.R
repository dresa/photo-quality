# Speed test for different sampling algorithms

n.points <- 100000
probs <- rep(1/n.points, n.points)


# A simplistic way of sampling a single number from vector x.
# Faster than using sample.int -- perhaps it's designed for vectorized use.
sample.single <- function(x, prob=NULL) {
  if (is.null(prob)) {
    return(x[runif(1) %/% (1 / length(x)) + 1])  # uniform distribution
  } else {
    stopifnot(length(prob) == length(x))
    return(x[findInterval(runif(1), cumsum(prob/sum(prob))) + 1])
  }
}
## Tests:
#stopifnot(sample.single(2) == 2)
#stopifnot(all(abs(table(replicate(3000, sample.single(2:4))) - 1000) < 100))
#set.seed(1)
#stopifnot(sample.single(seq(1.23, 2.34, 0.001)) == 1.524)
#stopifnot(is.null(sample.single(c())))
#tryCatch({sample.single(1:3, c(0.6,0.4)); write("error not catched: x and prob mismatch", stderr())}, error=function(x) {})
#set.seed(1)
#stopifnot(replicate(20, sample.single(2:5, prob=c(0.1,0.4,0.2,0.3))) == c(3,3,4,5,3,5,5,4,4,2,3,3,4,3,5,3,5,5,3,5))
#stopifnot(abs(table(replicate(10000, sample.single(2:5, prob=c(0.1,0.4,0.2,0.3)))) - 1000*c(1,4,2,3)) < 200)
#stopifnot(sample.single(0:4, prob=c(0, 1e-12, 0, 1 - 1e-12, 0)) == 3)



# NOT TESTED!
sample.esa <- function(x, size, prob=NULL, replace=FALSE) {
  if (replace) {
    if (is.null(prob)) {
      return(x[runif(size) %/% (1 / length(x)) + 1])  # uniform distribution
    } else {
      stopifnot(length(prob) == length(x))
      cumu.prob <- cumsum(prob / sum(prob))
      return(x[findInterval(runif(size), cumu.prob) + 1])
    }
  } else {  # without replacement
    if (is.null(prob)) {
      return(x[runif(size) %/% (1 / length(x)) + 1])  # uniform distribution
    } else {
      stopifnot(length(prob) == length(x))
      cumu.prob <- cumsum(prob / sum(prob))
      return(x[findInterval(runif(size), cumu.prob) + 1])
    }
  }
}
#print(table(sample.esa(11:14, 10000)))  
#s <- sample.esa(2:6, 100000, prob=c(0.4, 0.30, 0.15, 0.1, 0.05))
#print(table(s))


repeats <- 10000
set.seed(1)
t0 <- proc.time()
#for (i in 1:repeats) { as.numeric(cut(runif(1), breaks=c(0,cumsum(probs)),include.lowest=FALSE)) }
t1 <- proc.time()
set.seed(1)
for (i in 1:repeats) { sample(1:n.points, 1, prob=probs) }
t2 <- proc.time()
set.seed(1)
range <- 1:n.points
for (i in 1:repeats) { sample(range, 1, prob=probs) }
t3 <- proc.time()
set.seed(1)
for (i in 1:repeats) { sample.int(n.points, 1, prob=probs) }
t4 <- proc.time()
set.seed(1)
for (i in 1:repeats) { sample.single(range, prob=probs) }
t5 <- proc.time()

print(t1-t0)
print(t2-t1)
print(t3-t2)
print(t4-t3)
print(t5-t4)

# Results
#user  system elapsed 
#0       0       0  <-- this is very slow, really 
#user  system elapsed 
#15.28    0.02   15.68 
#user  system elapsed 
#14.29    0.01   14.88 
#user  system elapsed 
#17.09    0.00   17.37
