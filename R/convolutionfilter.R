##
## Convolution filters for processing 2D images
##
## Esa Junttila, 2016-03-23


meanFilter <- function(size) {
  if (!(size%%2)) stop('even filter size not possible at the moment')
  f <- matrix(1/(size*size), nrow=size, ncol=size)
  return(f)
}

gaussianFilter1D <- function(n, sigma, as.row=FALSE) {
  variance <- sigma^2
  # Due to normalization, we can omit constants from Gaussian:
  #   1/(sqrt(2*pi)*sigma)*exp(-x^2/(2*sigma^2))
  # This method simply uses the midpoints of the "pixels", which neglects non-linear shapes.
  # Points: -n, -(n-1), ..., 0, 1, 2, ..., n
  pos.indices <- 1:n  # symmetric (and g(0)=1), so only positive points need to be computed
  vals <- exp(-pos.indices^2/(2*variance))  # note that function value at 0 is 1
  m <- matrix(c(rev(vals), 1, vals) / (2*sum(vals)+1))  # reversed, zero-point, values: normalize such that sum=1.0
  return(if (as.row) t(m) else m)
}

# For testing purposes. A more efficient method of applying two 1D filters exists.
gaussianFilter2D <- function(n, sigma) {
  variance <- sigma^2
  # Due to normalization, we can omit constants from Gaussian:
  #   1/(2*pi*sigma^2)*exp(-(x^2+y^2)/(2*sigma^2))
  # This method simply uses the midpoints of the "pixels", which neglects non-linear shapes.
  indices <- -n:n
  idx.grid <- expand.grid(indices, indices)
  x <- idx.grid$Var1
  y <- idx.grid$Var2
  vals <- exp(-(x^2+y^2)/(2*variance))
  vals <- vals / sum(vals)  # normalize such that sum=1.0
  m <- matrix(vals, nrow=length(indices))
  return(m)
}

sobelVertical <- function() {
  return(matrix(c(
    c(-1, 0, +1),
    c(-2, 0, +2),
    c(-1, 0, +1)
  ), nrow=3, byrow=TRUE))
}

