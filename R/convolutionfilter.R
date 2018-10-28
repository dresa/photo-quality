##
## Convolution filters for processing 2D images
##
## Esa Junttila, 2016-03-23


#' @rdname convolutionfilter
#' @name convolutionfilter
#' @aliases meanFilter
#' @aliases gaussianFilter1D
#' @aliases gaussianFilter2D
#' @aliases sobelVertical
#' @title Convolution filters.
#' @description Collection of functions that generate \emph{filter} matrices,
#'   used in 1D and 2D convolution for manipulating matrices and images.
#' @param size [odd integer] size of a square matrix with
#'   dimensions \code{size} x \code{size}.
#' @param n integer length of one tail of Gaussian distribution,
#'   resulting in a vector-like matrix with length \emph{2*n+1}.
#' @param sigma numeric standard deviation of a Gaussian distribution
#'   from which filter weights are generated.
#' @param as.row [boolean] should the result matrix be formatted as
#'   a row vector or as column vector (technically still a matrix).
#'   Default is \code{FALSE} (as column vector).
#' @details The filter matrices are meant to be used together with
#'   convolution methods such as \code{imfilter}. The sum of the values
#'   in a filter matrix is usually 1 (to preserve total mass), but other
#'   sums are also commonly seen for other use cases.
#' @seealso
#'   
#'   Introduction on filters (kernels):
#'   \url{https://en.wikipedia.org/wiki/Kernel_(image_processing)}
#'   
#'   Convolution methods defined in
#'   \code{\link{imageconvolution}}
#' 
NULL


#' \code{meanFilter} generates a filter matrix that converts each point
#'   into a mean value of its immediate neighborhood.
#'   The filter's weights have a uniform value and their sum is 1.
#' @return \code{meanFilter} returns a square matrix of \emph{size} x
#'   \emph{size} where each element has value \emph{1/size^2}.
#' @details Mean filter is explained here:
#'   \url{http://www.librow.com/articles/article-5}
#' @examples
#' meanFilter(1) == 1
#' all(meanFilter(3) == 1/3^2)
#' @rdname convolutionfilter
#' @export
meanFilter <- function(size) {
  if (!(size%%2)) stop('even filter size is not possible at the moment')
  f <- matrix(1/(size*size), nrow=size, ncol=size)
  return(f)
}


#' \code{gaussianFilter1D} generates a filter matrix (more like a vector)
#' that mixes the values of the point's immediate horizontal/vertical
#' neighborhood, with more weight on closer points. The weights are
#' determined by a centered Gaussian distribution and normalized
#' to have sum of 1.
#' @return \code{gaussianFilter1D} returns a matrix of size
#'   \itemize{
#'   \item \emph{1} x (\emph{2*n+1}), if \code{as.row=TRUE}
#'   \item (\emph{2*n+1}) x \emph{1}, if \code{as.row=FALSE} (default).
#'   }
#'   The values are sampled from a centered Gaussian distribution
#'   with variance \code{sigma}^2, and normalized to have sum 1.
#' @details
#'   In the Gaussian filters, with parameter \code{n}, we consider
#'   the closest \code{n} pixels on both sides of center pixel 0, namely
#'   \emph{-n, ..., -1, 0, 1, 2, ..., n}.
#'
#'   The standard Gaussian density is
#'   \code{1/(sqrt(2*pi)*sigma)*exp(-x^2/(2*sigma^2))},
#'   but thanks to normalization, we can omit the common factor
#'   \code{1/(sqrt(2*pi)*sigma)} to use the following simple form instead:
#'   \code{exp(-x^2/(2*sigma^2))}.
#'   
#'   Current implementation samples the Gaussian densities
#'   at the midpoints of the pixels, which neglects non-linear shapes.
#' @examples
#' 
#' gaussianFilter1D(0, 99) == 1
#' all(gaussianFilter1D(1, 1) == dnorm(-1:+1)/sum(dnorm(-1:+1)))
#' all(gaussianFilter1D(2, 3, as.row=TRUE) == dnorm(-2:2,sd=3)/sum(dnorm(-2:2,sd=3)))
#' @rdname convolutionfilter
#' @export
gaussianFilter1D <- function(n, sigma, as.row=FALSE) {
  variance <- sigma^2
  # Compute only positive points because of symmetry and density(0)=1.
  pos.indices <- if (n >= 1) 1:n else c()
  vals <- exp(-pos.indices^2/(2*variance))  # note that value at 0 is 1
  # Build from parts: reversed, zero-point, values. Normalize so that sum is 1.
  m <- matrix(c(rev(vals), 1, vals) / (2*sum(vals)+1))
  return(if (as.row) t(m) else m)
}


#' \code{gaussianFilter2D} generates a filter matrix that mixes
#' the values of the point's immediate neighborhood,
#' with more weight on closer points. The weights are
#' determined by a multivariate 2D Gaussian distribution that is
#' centered, uncorrelated, and normalized to have sum of 1.
#' This method has been implemented only for
#' testing purposes: it is more efficient to apply a 1D filter
#' and then another 1D filter to reach the same result (this is
#' possible because the dimensions are uncorrelated).
#' @return \code{gaussianFilter2D} returns a matrix of size
#'   (\emph{2*n+1}) x (\emph{2*n+1}).
#'   The values are sampled from a centered 2D Gaussian distribution
#'   with variances \code{sigma}^2, and normalized to have sum 1.
#' @details
#'   In a 2D Gaussian filter, with parameter \code{n}, we consider
#'   the closest \code{n} pixels on all sides of center pixel (0,0),
#'   namely from (-n, -n) to (n,n), including (-n, n) and (n, -n).
#'   
#'   The standard uncorrelated 2D Gaussian density (with shared \code{sigma})
#'   is \code{1/(2*pi*sigma^2)*exp(-(x^2+y^2)/(2*sigma^2))},
#'   but thanks to normalization, we can omit the common factor
#'   \code{1/(2*pi*sigma^2)} to use the following simple form instead:
#'   \code{exp(-(x^2+y^2)/(2*sigma^2))}.
#'   
#'   Current implementation samples the 2D Gaussian densities
#'   at the midpoints of the pixels, which neglects non-linear shapes.
#'   
#'   Multivariate Gaussian distibution:
#'   \itemize{
#'   \item \url{https://math.stackexchange.com/questions/803329/multivariate-normal-distribution-independet-iff-uncorrelated}
#'   \item \url{https://www.statisticshowto.datasciencecentral.com/bivariate-normal-distribution/}
#'   }
#'   
#' @examples
#' 
#' library(mvtnorm)  # multivariate normal distribution for testing
#' gaussianFilter2D(0, 99) == 1
#' tol <- 1e-14
#' # Bivariate uncorrelated Gaussian with identical sigmas.
#' g2d <- function(x,s) dmvnorm(x=x, sigma=diag(s^2, 2))
#' testmvg <- function(k, s) {
#'   points <- as.matrix(expand.grid(-k:k, -k:k))
#'   target <- matrix(g2d(points, s)/sum(g2d(points, s)), nrow=2*k+1)
#'   all(abs(gaussianFilter2D(k, s) - target) < tol)
#' }
#' testmvg(1, 1)
#' testmvg(2, 3)
#' @rdname convolutionfilter
#' @export
gaussianFilter2D <- function(n, sigma) {
  variance <- sigma^2
  indices <- -n:n
  idx.grid <- expand.grid(indices, indices)
  x <- idx.grid$Var1
  y <- idx.grid$Var2
  vals <- exp(-(x^2+y^2)/(2*variance))
  vals <- vals / sum(vals)  # normalize such that sum is 1.
  m <- matrix(vals, nrow=length(indices))
  return(m)
}


#' \code{sobelVertical} generates a predefined filter matrix suitable for
#' detecting vertically aligned edges based on neighboring values.
#' The filter's weights have sum 0, and a large deviation from 0 indicates
#' a strong vertically aligned edge. To detect horizontally aligned edges,
#' use the transposed filter matrix instead.
#' @return \code{sobelVertical} returns a predefined square matrix
#'   of size \emph{3} x \emph{3} whose sum is 0.
#' @details Sobel operator is explained here:
#'   \url{https://en.wikipedia.org/wiki/Sobel_operator}
#' @examples
#' 
#' sum(sobelVertical()) == 0
#' 
#' @rdname convolutionfilter
#' @export
sobelVertical <- function() {
  return(matrix(c(
    c(-1, 0, +1),
    c(-2, 0, +2),
    c(-1, 0, +1)
  ), nrow=3, byrow=TRUE))
}

