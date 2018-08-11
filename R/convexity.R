# Convexity
#
# Three primary functions:
# * findConvexHull2D:   Given a set of 2D points, find the convex hull.
# * insideConvexHull2D: Given a point, determine if it's inside a convex hull.
# * convexHullArea2D:   Compute the area of a 2D convex hull
#
# Esa Junttila
# 2018-02-11

library(grDevices)  # "chull" function for computing convex hull


#' Cross product for 2D vectors.
#' 
#' Computes the cross-product for given 2D vectors. Arguments
#' \code{px}, \code{py}, \code{qx}, and \code{qy} contain the
#' \emph{x} and \emph{y} coordinates for 2D vectors \emph{p} and \emph{q}.
#' The 2D cross product is a special case of 3D Euclidean cross-product,
#' reducing into \emph{p_x * q_y - p_y * q_x}.
#' 
#' @param px \emph{x} coordinates of \emph{p} vectors
#' @param py \emph{y} coordinates of \emph{p} vectors
#' @param qx \emph{x} coordinates of \emph{q} vectors
#' @param qy \emph{y} coordinates of \emph{q} vectors
#' @return 2D cross-product for each pair of \emph{p} and \emph{q}.
#' @examples
#' crossprod2D(-2, 3, 4, -5) == -2
#' all(crossprod2D(c(6,-1), c(2,5), c(-3,-4), c(-7,8)) == c(-36,12))
#' @seealso \url{https://en.wikipedia.org/wiki/Cross_product}
#' @export
crossprod2D <- function(px, py, qx, qy) px*qy - py*qx




#' @rdname vectorTurn2D
#' @name vectorTurn2D
#' @aliases isClockwiseTurn2D
#' @aliases isCounterClockwiseTurn2D
#' @title Direction of vector turns and rotations in 2D space.
#' @description Determine if the turn from a vector \emph{p} to another
#'   vector \emph{q} in Euclidean 2D space is clock-wise or
#'   counter-clockwise.
#' @param px \emph{x} coordinates of vectors \emph{p}
#' @param py \emph{y} coordinates of vectors \emph{p}
#' @param qx \emph{x} coordinates of vectors \emph{q}
#' @param qy \emph{y} coordinates of vectors \emph{q}
#' @details Each argument \code{px}, \code{py}, \code{qx},
#'   and \code{qy} must have equal length.
#' @examples
#' isClockwiseTurn2D(3, 4, -2, 1) != isCounterClockwiseTurn2D(3, 4, -2, 1)
NULL


#' \code{isClockwiseTurn2D} determines if the turn from 2D
#' vector \emph{p} to \emph{q} is clockwise (inclusive).
#' 
#' @return \code{isClockwiseTurn2D} returns \code{TRUE}, if the
#'   turn is clockwise (inclusive); otherwise \code{FALSE};
#'   for all pairs of \each{p} and \each{q}.
#' @examples
#' isClockwiseTurn2D(0, 1, 0, 1)  # extreme: both turns apply
#' isClockwiseTurn2D(0, 1, 0, -3)  # extreme: both turns apply
#' isClockwiseTurn2D(-2, 3, 4, -6)  # extreme: both turns apply
#' isClockwiseTurn2D(0, 1, 0.1, 1000)  # slightly clockwise
#' !isClockwiseTurn2D(0, 1, -0.1, 1000)  # slightly counter-clockwise
#' isClockwiseTurn2D(-2, 3, 4, -5.9)
#' !isClockwiseTurn2D(-2, 3, 4, -6.1)
#' @rdname vectorTurn2D
#' @export
isClockwiseTurn2D <- function(px, py, qx, qy) {
  return(crossprod2D(px, py, qx, qy) <= 0)
}


# Is the turn from 2D vector p to q counter-clockwise? Inclusive.
#' \code{isCounterClockwiseTurn2D} determines if the turn from 2D
#' vector \emph{p} to \emph{q} is counter-clockwise (inclusive).
#' @examples
#' isCounterClockwiseTurn2D(0, 1, 0, 1)  # extreme: both turns apply
#' isCounterClockwiseTurn2D(0, 1, 0, -3)  # extreme: both turns apply
#' isCounterClockwiseTurn2D(-2, 3, 4, -6)  # extreme: both turns apply
#'   
#' @return \code{isCounterClockwiseTurn2D} returns \code{TRUE}, if the
#'   turn is counter-clockwise (inclusive); otherwise \code{FALSE};
#'   for all pairs of \each{p} and \each{q}.
#' @rdname vectorTurn2D
#' @export
isCounterClockwiseTurn2D <- function(px, py, qx, qy) crossprod2D(px, py, qx, qy) >= 0


# Ordinary Euclidean distance measure.
euclideanDistance2D <- function(x1, y1, x2, y2) sqrt((x2-x1)^2 + (y2-y1)^2)



#' @rdname distanceLine2D
#' @name distanceLine2D
#' @aliases distanceToLine2D
#' @aliases distanceToLineSegment2D
#' @title Distance of a point to a line or line segment.
#' @description Compute the minimum distance from a 2D point
#'  (\emph{x},\emph{y}) to a line that goes through two other
#'  points (\emph{x_1},\emph{y_1}) and (\emph{x_2},\emph{y_2}).
#' @param x \emph{x} coordinate of the reference point
#' @param y \emph{y} coordinate of the reference point
#' @param x1 \emph{x} coordinate of the first point that defines the line
#' @param y1 \emph{y} coordinate of the first point that defines the line
#' @param x2 \emph{x} coordinate of the second point that defines the line
#' @param y2 \emph{y} coordinate of the second point that defines the line
#' @return Minimum distance from a point to a line or a line segment.
#' @details Vectorized; each argument \code{x}, \code{y}, \code{x1},
#'   \code{y1}, \code{x2}, and \code{y2} must have equal length.
NULL


#' \code{distanceToLine2D} computes the minimum distance from
#'   (\emph{x},\emph{y}) to a line that goes through the other two points,
#'   and goes on indefinitely.
#' @examples
#' distanceToLine2D(0, 1,  0, 0,  0, 2) == 0  # on the line
#' distanceToLine2D(3, 5,  1, 2,  3, 5) == 0  # on the line
#' distanceToLine2D(-4, 3,  -2, -1,  -2, 5) == 2
#' distanceToLine2D(-4, 999,  -2, -1,  -2, 5) == 2
#' distanceToLine2D(1, 2,  4, 3,  6, 2) == sqrt(5)
#' distanceToLine2D(4.5,2, 6,2, 4,3) == 1.5/sqrt(5)  # closest line-point (4.8, 2.6)
#' @rdname distanceLine2D
#' @export
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


#' \code{distanceToLineSegment2D} computes the minimum distance from 
#' (\emph{x},\emph{y}) to a line segment that spans between the other
#'two points, and stops there.
#' @examples
#' distanceToLineSegment2D(1,2, 4,3, 6,2) == sqrt(10)  # "before"
#' distanceToLineSegment2D(1,2, 6,2, 4,3) == sqrt(10)  # "after"
#' distanceToLineSegment2D(4.5,2, 6,2, 4,3) == 1.5/sqrt(5)  # "projection"
#'   
#' @seealso Inspiration from:
#'   \url{http://geomalgorithms.com/a02-_lines.html#Distance-to-Ray-or-Segment}
#' @rdname distanceLine2D
#' @export
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


#' @rdname convexHull2D
#' @name convexHull2D
#' @aliases insideConvexHull2D
#' @aliases convexHullArea2D
#' @aliases findConvexHull2D
#' @title Convex Hull functions in 2D space.
#' @description Functions that deal with convex hulls on a set of 2D points.
#' @param x \emph{x} coordinates of 2D point set
#' @param y \emph{y} coordinate of 2D point set
#' @param px \emph{x} coordinates of points we test for being inside a convex hull.
#' @param py \emph{y} coordinates of points we test for being inside a convex hull.
#' @param ch.x \emph{x} coordinates of points that define a convex
#' hull (clockwise order), without repeating the initial point.
#' @param ch.y \emph{y} coordinates of points that define a convex
#' hull (clockwise order), without repeating the initial point.
#' @param eps tolerance used for testing whether point is inside a convex hull.
#' @seealso Convex Hull definition:
#' \itemize{
#'   \item \url{https://en.wikipedia.org/wiki/Convex_hull}
#'   \item \url{http://geomalgorithms.com/a10-_hull-1.html}
#' }
NULL


#' \code{findConvexHull2D} computes the convex hull for a set of 2D points,
#' given as \emph{x} and \emph{y} coordinates.
#' @return 
#' \code{findConvexHull2D} returns the point indices
#' (of \code{x} and \code{y}) that define the convex hull.
#' The points are listed in clock-wise order, and the initial point
#' is not repeated.
#' @examples
#' # Compute a convex hull for a set of 2D points.
#' test.x <- c(5,4,8,9,4,10,6,4,8,8,2,6,3,4)
#' test.y <- c(3,3,1,10,3,8,9,9,5,1,10,1,3,1)
#' test.ch <- c(10,14,13,11,4,6)
#' test.ch.comp <- findConvexHull2D(test.x, test.y)
#' all(test.ch.comp == test.ch)
#' all(test.x[test.ch.comp] == test.x[test.ch])
#' all(test.y[test.ch.comp] == test.y[test.ch])
#' 
#' @rdname convexHull2D
#' @export
findConvexHull2D <- function(x, y) {
  ch <- chull(x,y)  # compute convex hull (indices in clockwise order, initial not repeated)
  return(ch)
}


#' \code{insideConvexHull2D} determines whether 2D points \emph{p},
#' given as \emph{x} and \emph{y} coordinates, are inside a convex hull,
#' defined by points \emph{ch} (again as \emph{x} and \emph{y} coordinates)
#' in clockwise order.
#' @return \code{insideConvexHull2D} returns a vector of boolean values
#' that determine whether (\code{px},\code{py}) points are inside
#' given convex hull (inclusive).
#' @examples
#' # Determine which points are inside a convex hull:
#' all(insideConvexHull2D(c(1.89, 1.9, 2, 2+1e-15, 2.1), c(3.89, 3.9, 4, 4+1e-15, 4.1), c(1.90, 2.05), c(3.90,4.05)) == c(F,T,T,T,F))
#' test.x <- c(5,4,8,9,4,10,6,4,8,8,2,6,3,4)
#' test.y <- c(3,3,1,10,3,8,9,9,5,1,10,1,3,1)
#' test.ch <- c(10,14,13,11,4,6)  # incides of points that define a convex hull
#' test.ch.x <- test.x[test.ch]
#' test.ch.y <- test.y[test.ch]
#' plot(test.x, test.y)
#' lines(c(test.ch.x, test.ch.x[1]), c(test.ch.y, test.ch.y[1]), col='blue')
#' points(test.ch.x[1], test.ch.y[1], col="black", pch=19)
#' p.x <- c(8,4,3, 2, 9,10,  6, 6, 9,   0, 2,   2,    9.34, 9)
#' p.y <- c(1,1,3,10,10, 8,  6, 1, 4.5, 0, 9.9, 10.1, 9.34, 4.4)
#' p.expected <- c(T, T, T, T, T, T, T, T, T, F, F, F, F, F)
#' all(insideConvexHull2D(p.x, p.y, test.ch.x, test.ch.y) == p.expected)
#' library(microbenchmark)
#' print(microbenchmark(insideConvexHull2D(p.x, p.y, test.ch.x, test.ch.y), times=10000))
#' x <- runif(10); y <- runif(10); plot(x, y, pch=15); ch <- chull(x,y);
#' lines(x[c(ch, ch[1])], y[c(ch,ch[1])]);
#' points(x[ch[1]], y[ch[1]], col='blue', pch=20)
#' test.x <- runif(1000); test.y <- runif(1000);
#' colors <- ifelse(insideConvexHull2D(test.x, test.y, x[ch], y[ch]), 'green', 'red')
#' points(test.x, test.y, col=colors, pch=19)
#' 
#' @seealso
#' Testing whether a point is inside a convex hull is inspired by
#' \url{https://salzis.wordpress.com/2014/05/01/convex-hull-how-to-tell-whether-a-point-is-inside-or-outside/}.
#' @rdname convexHull2D
#' @export
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


#' \code{convexHullArea2D} computes the area of a given 2D convex hull.
#' We assume that the points \emph{ch} (represented as \emph{x} and
#' \emph{y} coordinates) define a proper convex hull, and the points are
#' adjacent in a convex hull point sequence (either clockwise or
#' counter-clockwise), without repeating the initial point.
#' @return \code{convexHullArea2D} returns the numerical area of given
#' convex hull (geometric 2D object).
#' @examples
#' # Measure the area within a convex hull:
#' convexHullArea2D(c(2,2,4,4), c(6,9,9,6)) == 6  # square
#' convexHullArea2D(c(1,2,3,4,5,5,4,3), c(1,3,3,3,2,1,0,-1)) == 10.5
#' 
#' @rdname convexHull2D
#' @export
convexHullArea2D <- function(ch.x, ch.y) {
  stopifnot(length(ch.x) == length(ch.y) && length(ch.x) >= 3)
  # Approximate area by finding a mass center and splitting the convex hull
  # into triangles, each of which has a point (x0, y0) at origin.
  # Mass center is always inside the convex hull by definition.
  c.x <- scale(ch.x, scale=FALSE)  # mean-centered, origin for triangles
  c.y <- scale(ch.y, scale=FALSE)
  c.next.x <- c(c.x[2:length(ch.x)], c.x[1]) # first as last
  c.next.y <- c(c.y[2:length(ch.y)], c.y[1])
  triangleArea <- function(x1,y1,x2,y2) abs(x1*y2 - y1*x2)/2
  return(sum(mapply(triangleArea, c.x, c.y, c.next.x, c.next.y)))
}



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


