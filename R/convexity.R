# Convexity
#
# Three primary functions:
# * findConvexHull2D:   Given a set of 2D points, find the convex hull.
# * insideConvexHull2D: Given a point, determine whether it's inside a convex hull.
# * convexHullArea2D:   Compute the area of a 2D convex hull
#
# Esa Junttila
# 2018-02-11

library(grDevices)  # "chull" function for computing convex hull


# Cross product for 2D vectors: two vectors p and q, with x and y coordinates
crossprod2D <- function(px, py, qx, qy) px*qy - py*qx
## Test:
#stopifnot(crossprod2D(-2, 3, 4, -5) == -2)


# Is the turn from 2D vector p to q clockwise? Inclusive.
isClockwiseTurn2D <- function(px, py, qx, qy) {
  return(crossprod2D(px, py, qx, qy) <= 0)
}
# Is the turn from 2D vector p to q counter-clockwise? Inclusive.
isCounterClockwiseTurn2D <- function(px, py, qx, qy) crossprod2D(px, py, qx, qy) >= 0
## Test:
#stopifnot(isClockwiseTurn2D(0, 1, 0, 1))  # extreme: both turns apply
#stopifnot(isClockwiseTurn2D(0, 1, 0, -3))  # extreme: both turns apply
#stopifnot(isCounterClockwiseTurn2D(0, 1, 0, 1))  # extreme: both turns apply
#stopifnot(isCounterClockwiseTurn2D(0, 1, 0, -3))  # extreme: both turns apply
#stopifnot(isClockwiseTurn2D(-2, 3, 4, -6))  # extreme: both turns apply
#stopifnot(isCounterClockwiseTurn2D(-2, 3, 4, -6))  # extreme: both turns apply
#stopifnot(isClockwiseTurn2D(0, 1, 0.1, 1000))  # slightly clockwise
#stopifnot( !isClockwiseTurn2D(0, 1, -0.1, 1000) )  # slightly counter-clockwise
#stopifnot(isClockwiseTurn2D(-2, 3, 4, -5.9))
#stopifnot( !isClockwiseTurn2D(-2, 3, 4, -6.1) )
#stopifnot(isClockwiseTurn2D(3, 4, -2, 1) != isCounterClockwiseTurn2D(3, 4, -2, 1))

# Ordinary Euclidean distance measure
euclideanDistance2D <- function(x1, y1, x2, y2) sqrt((x2-x1)^2 + (y2-y1)^2)

# Distance from a point (x,y) to a line that goes though two
# other points (x1,y1) and (x2,y2).
distanceToLine2D <- function(x, y, x1, y1, x2, y2) {
  stopifnot(length(unique(unlist(lapply(list(x,y,x1,y1,x2,y2), length)))) == 1)
  if (length(x) == 0) return(c())
  dist <- euclideanDistance2D
  norm <- dist(x1, y1, x2, y2)
  d <- ifelse(norm == 0,
    dist(x, y, x1, y1),  # point-to-point distance
    abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1) / norm
  )
  #b <- ((x -x1)*(x2-x1) + (y -y1)*(y2-y1)) / ((x2-x1)^2 + (y2-y1)^2)
  #closest <- c(x1+b*(x2-x1), y1+b*(y2-y1))  # closest (x3,y3) point on the line
  return(d)
}
## Test: 
#stopifnot(distanceToLine2D(0, 1,  0, 0,  0, 2) == 0)
#stopifnot(distanceToLine2D(3, 5,  1, 2,  3, 5) == 0)
#stopifnot(distanceToLine2D(-4, 3,  -2, -1,  -2, 5) == 2)
#stopifnot(distanceToLine2D(-4, 999,  -2, -1,  -2, 5) == 2)
#stopifnot(distanceToLine2D(1, 2,  4, 3,  6, 2) == sqrt(5))
#stopifnot(distanceToLine2D(4.5,2, 6,2, 4,3) == 1.5/sqrt(5))  # closest line-point (4.8, 2.6)


# Distance from a point (x,y) to a line segment that has two
# endpoints as (x1,y1) and (x2,y2). Vectorized; equal-length args expected.
# Inspiratiom from: http://geomalgorithms.com/a02-_lines.html#Distance-to-Ray-or-Segment
distanceToLineSegment2D <- function(x, y, x1, y1, x2, y2) {
  stopifnot(length(unique(unlist(lapply(list(x,y,x1,y1,x2,y2), length)))) == 1)
  # Using two vectors:
  #   w: from line segment start (x1,y1) to reference point (x,y) = (x,y)-(x1,y1)
  #   v: from line segment start (x1,y1) to segment end (x2,y2) = (x2,y2)-(x1,y1)
  c1 <- (x -x1)*(x2-x1) + (y -y1)*(y2-y1)  # dot product between w and v, determine angle
  c2 <- (x2-x1)^2 + (y2-y1)^2  # dot product between v and v, determine angle
  # filters:
  b <- c1 <= 0    # "before" line segment start (x1,y1)
  a <- c1 > 0 & c2 <= c1  # "after" line segment start (x2,y2)
  p <- !b & !a  # projection on line segment (not before nor after endpoints)
  # compute distances:
  dist <- euclideanDistance2D
  d <- numeric(length(x))
  d[b] <- dist(x[b], y[b], x1[b], y1[b])  # distance to segment start point
  d[a] <- dist(x[a], y[a], x2[a], y2[a])  # distance to segment end point
  d[p] <- distanceToLine2D(x[p], y[p], x1[p], y1[p], x2[p], y2[p])  # minimum distance to line segment
  if (any(c1 <= 0 & c2 <= c1)) print(paste('Found special case:', c1, c2))
  return(d)
}
## Test:
#stopifnot(distanceToLineSegment2D(1,2, 4,3, 6,2) == sqrt(10))  # "before"
#stopifnot(distanceToLineSegment2D(1,2, 6,2, 4,3) == sqrt(10))  # "after"
#stopifnot(distanceToLineSegment2D(4.5,2, 6,2, 4,3) == 1.5/sqrt(5))  # "projection"




# Given a set of 2D points in terms of their px and py coordinates,
# determine whether they are inside a given convex hull,
# defined by points ch in clockwise order.
# Inspirtion from: https://salzis.wordpress.com/2014/05/01/convex-hull-how-to-tell-whether-a-point-is-inside-or-outside/
insideConvexHull2D <- function(px, py, ch.x, ch.y, eps=1e-14) {
  n <- length(px)
  k <- length(ch.x)
  stopifnot(length(py) == n && length(ch.y) == k && k >= 1)
  is.degen.point <- k == 1 || sum(abs(c(diff(ch.x), diff(ch.y)))) < 2*k*eps
  if (is.degen.point) return(abs(px - ch.x[1]) < eps & abs(py - ch.y[1]) < eps)
  else if (k == 2) {
    sx <- rep(ch.x[1],n);  sy <- rep(ch.y[1],n)  # facet start replicated
    ex <- rep(ch.x[2],n);  ey <- rep(ch.y[2],n)  # facet end replicated
    return(distanceToLineSegment2D(px, py, sx, sy, ex, ey) < eps)
  }
  else if (k >= 3) {
    # Looking at start--end vector, does start--point vector turn clockwise from it?
    # Check each (point,facet) pair; if TRUE for all facets, then point in inside convex hull.
    facets.x <- diff(c(ch.x, ch.x[1]))  # relative shift w.r.t one facet
    facets.y <- diff(c(ch.y, ch.y[1]))
    inside <- !logical(n)  # initialized with TRUE
    # Conserve memory: one facet at a time. We know k should be small.
    for (i in 1:k) { inside <- inside & isClockwiseTurn2D(facets.x[i], facets.y[i], px - ch.x[i], py - ch.y[i]) }
    return(inside)
  }
}
### Test:
#library(microbenchmark)
#stopifnot(all(insideConvexHull2D(c(1.89, 1.9, 2, 2+1e-15, 2.1), c(3.89, 3.9, 4, 4+1e-15, 4.1), c(1.90, 2.05), c(3.90,4.05)) == c(F,T,T,T,F)))
#test.x <- c(5,4,8,9,4,10,6,4,8,8,2,6,3,4)
#test.y <- c(3,3,1,10,3,8,9,9,5,1,10,1,3,1)
#test.ch <- c(10,14,13,11,4,6)
#test.ch.x <- test.x[test.ch]
#test.ch.y <- test.y[test.ch]
###plot(test.x, test.y)
###lines(c(test.ch.x, test.ch.x[1]), c(test.ch.y, test.ch.y[1]), col='blue')
###points(test.ch.x[1], test.ch.y[1], col="black", pch=19)
#p.x <- c(8,4,3, 2, 9,10,  6, 6, 9,   0, 2,   2,    9.34, 9)
#p.y <- c(1,1,3,10,10, 8,  6, 1, 4.5, 0, 9.9, 10.1, 9.34, 4.4)
#p.expected <- c(T, T, T, T, T, T, T, T, T, F, F, F, F, F)
#stopifnot(all(insideConvexHull2D(p.x, p.y, test.ch.x, test.ch.y) == p.expected))
#print(microbenchmark(insideConvexHull2D(p.x, p.y, test.ch.x, test.ch.y), times=10000))
#x <- runif(10); y <- runif(10); plot(x, y, pch=15); ch <- chull(x,y); lines(x[c(ch, ch[1])], y[c(ch,ch[1])]); points(x[ch[1]], y[ch[1]], col='blue', pch=20)
#test.x <- runif(1000); test.y <- runif(1000); points(test.x, test.y, col=ifelse(insideConvexHull2D(test.x, test.y, x[ch], y[ch]), 'green', 'red'), pch=19)



# Compute the area of a convex hull. We are assuming the
# points x and y are adjacent in the convex hull.
convexHullArea2D <- function(x, y) {
  stopifnot(length(x) == length(y) && length(x) >= 3)
  # Approximate area by finding a mass center and splitting into triangles.
  c.x <- scale(x, scale=FALSE)  # mean-centered
  c.y <- scale(y, scale=FALSE)
  c.next.x <- c(c.x[2:length(x)], c.x[1]) # first as last
  c.next.y <- c(c.y[2:length(y)], c.y[1])
  return(sum(mapply(function(a,b,c,d) abs(a*d-b*c)/2, c.x, c.y, c.next.x, c.next.y)))
}


# Given a set of 2D points (their x and y coordinates), find the convex hull.
# Returns a set of point indices included in the convex hull. The points
# are listed in clock-wise order, and the initial point is not repeated.
findConvexHull2D <- function(x, y) {
  ch <- chull(x,y)  # compute convex hull (indices in clockwise order, initial not repeated)
  return(ch)
}
## Test:
#test.x <- c(5,4,8,9,4,10,6,4,8,8,2,6,3,4)
#test.y <- c(3,3,1,10,3,8,9,9,5,1,10,1,3,1)
#test.ch <- c(10,14,13,11,4,6)
#test.ch.comp <- findConvexHull2D(test.x, test.y)
#stopifnot(all(test.ch.comp == test.ch))
#stopifnot(all(test.x[test.ch.comp] == test.x[test.ch]))
#stopifnot(all(test.y[test.ch.comp] == test.y[test.ch]))


## Stress test
#many.x <- runif(100)
#many.y <- runif(100)
#check.x <- runif(10000)
#check.y <- runif(10000)
#t0 <- proc.time()
#ch <- findConvexHull2D(many.x, many.y)
#t1 <- proc.time()
#stopifnot(all(insideConvexHull2D(many.x, many.y, many.x[ch], many.y[ch])))
#t2 <- proc.time()
#ratio <- sum(insideConvexHull2D(check.x, check.y, many.x[ch], many.y[ch])) / length(check.x)
#t3 <- proc.time()
#print(t1 - t0)
#print(t2 - t1)
#print(t3 - t2)
#print(paste('Randomized ratio', ratio))
#print(paste('Analytical ratio', convexHullArea2D(many.x[ch], many.y[ch])))
#plot(many.x, many.y)
#ch.x <- many.x[ch]
#ch.y <- many.y[ch]
#lines(c(ch.x, ch.x[1]), c(ch.y, ch.y[1]), col='blue')
#points(ch.x[1], ch.y[1], col="black", pch=19)


