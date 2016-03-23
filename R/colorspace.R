##
## Color space transformation between color spaces, such as RGB and HSV.
##
## Esa Junttila, 2016-03-23


# Convert RGB image to an HSV image.
# Assumptions:
#   RGB values are between 0 and 1.
#   Hues are degrees between 0 and 360 (or marked as radians);
#   Saturations and Values are between 0 and 1.
toHSV <- function(img.rgb, radians=FALSE, max.value=1) {
  # Extract R,G, and B color channels
  red <- extractRGBChannel(img.rgb, RED)
  green <- extractRGBChannel(img.rgb, GREEN)
  blue <- extractRGBChannel(img.rgb, BLUE)
  # Temporary variables
  dims <- if (is.matrix(red)) dim(red) else c(length(red), 1)
  max.channels <- pmax(red, green, blue)
  diffs <- max.channels - pmin(red, green, blue)
  # Compute H(ue), S(saturation), and V(alue):
  # S(aturation):
  saturations <- diffs / max.channels
  saturations[max.channels == 0] <- 0  # replace Inf by zero (originates from division by zero)
  # H(ue):
  R<-0; G<-2; B<-4  # additive color-shift values in the hue formula
  max.layer <- array(R, dim=dims)  # which color layer has maximum value? R, G, or B?
  max.layer[green > red] <- G
  max.layer[blue > pmax(red, green)] <- B
  red.mask <- max.layer == R
  green.mask <- max.layer == G
  blue.mask <- max.layer == B
  hues <- array(dim=dims)
  angle <- if (radians) pi/3 else 60
  hues[red.mask]   <- angle * ((green[red.mask]  - blue[red.mask])   / diffs[red.mask]   + R)%%6
  hues[green.mask] <- angle * ((blue[green.mask] - red[green.mask])  / diffs[green.mask] + G)%%6
  hues[blue.mask]  <- angle * ((red[blue.mask]   - green[blue.mask]) / diffs[blue.mask]  + B)%%6
  hues[max.channels==0] <- NA  # undefined hue is not needed; NA as default
  hues[diffs==0] <- NA  # undefined hue is not needed; NA as default
  # V(alue):
  values <- max.channels / max.value
  # Return an HSV image as an (M x N x 3) array
  dim.names <- list(NULL, NULL, c('Hue', 'Saturation', 'Value'))
  return(array(c(hues, saturations, values), dim=c(dims[1], dims[2], 3), dimnames=dim.names))
}


# Convert HSV image to a RGB image.
# Assumptions:
#   RGB values are between 0 and 1.
#   Hues are degrees between 0 and 360 (or marked as radians);
#   Saturations and Values are between 0 and 1.
toRGB <- function(img.hsv, radians=FALSE, na.hue=0) {
  # Extract H, S, and V color channels
  hue <- extractHSVChannel(img.hsv, HUE)
  saturation <- extractHSVChannel(img.hsv, SATURATION)
  value <- extractHSVChannel(img.hsv, VALUE)
  # Preprocess hues: replace undefined hues with default red (zero in radians and degrees).
  # Actual hue should not matter since value or saturation should be zero anyway.
  hue[is.na(hue)] <- na.hue
  # Convenience variables
  dims <- if (is.matrix(hue)) dim(hue) else c(length(hue), 1)
  chroma <- value * saturation
  angle <- if (radians) pi/3 else 60  # 60 degrees
  x <- chroma * (1 - abs(((hue / angle) %% 2)- 1))
  m <- value - chroma
  
  # Conversion pre-processing
  red.add <- array(dim=dims)
  green.add <- array(dim=dims)
  blue.add <- array(dim=dims)
  hue.group <- hue %/% angle
  ## proceed by groups
  for (group.id in 0:5) {
    group.px <- hue.group == group.id
    C <- chroma[group.px]
    X <- x[group.px]
    red.add[group.px] <-   switch(group.id + 1, C, X, 0, 0, X, C)
    green.add[group.px] <- switch(group.id + 1, X, C, C, X, 0, 0)
    blue.add[group.px] <-  switch(group.id + 1, 0, 0, X, C, C, X)
  }
  # Wrap it up
  return(createImageRGB(m + red.add, m + green.add, m + blue.add))
}


test <- function() {
  set.seed(1)
  d <- c(4,12,3)
  rgb.vec <- sample(0:255, prod(d), replace=TRUE)
  rgb <- array(rgb.vec, dim=d) / 255
  
  # reference data: use R's own RGB->HSV->RGB conversion
  hsv.ref <- rgb2hsv(255*matrix(c(rgb[,,1], rgb[,,2], rgb[,,3]), nrow=3, byrow=TRUE))
  rgb.ref <- col2rgb(hsv(hsv.ref[1, ], hsv.ref[2, ], hsv.ref[3, ]))
  roundtrip.ok <- all(rgb == array(c(rgb.ref[1, ], rgb.ref[2, ], rgb.ref[3, ]), dim=d) / 255)
  if (!roundtrip.ok) stop('internal test failed: RGB->HSV->RGB conversion')
  
  # custom implementations
  hsv.val <- toHSV(rgb)
  err.hsv <- max(abs(hsv.ref - matrix(c(hsv.val[,,1]/360, hsv.val[,,2], hsv.val[,,3]), nrow=3, byrow=TRUE)))
  print(paste('Matrix RGB->HSV error', err.hsv))
  rgb.roundtrip <- toRGB(hsv.val)
  err.rgb <- max(abs(rgb.roundtrip - rgb))
  print(paste('Matrix HSV->RGB error', err.rgb))
  
  # just one
  rgb <- array(c(200, 100, 150), dim=c(1,1,3)) / 255
  err.hsv <- max(abs(as.vector(toHSV(rgb)) - as.vector(rgb2hsv(200, 100, 150))*c(360,1,1)))
  print(paste('Single HSV error', err.hsv))
  err.rgb <- max(abs(as.vector(toRGB(toHSV(rgb))) - as.vector(rgb)))
  print(paste('Single RGB error', err.rgb))
}

speedTest <- function() {
  set.seed(1)
  d <- c(1000,1500,3)
  rgb.vec <- sample(0:255, prod(d), replace=TRUE)
  rgb <- array(rgb.vec, dim=d) / 255
  
  refFunc <- function(rgb) {
    hsv.ref <- rgb2hsv(matrix(c(rgb[,,1], rgb[,,2], rgb[,,3]), nrow=3, byrow=TRUE))
    rgb.ref <- col2rgb(hsv(hsv.ref[1, ], hsv.ref[2, ], hsv.ref[3, ]))
    return(array(c(rgb.ref[1, ], rgb.ref[2, ], rgb.ref[3, ]), dim=d))
  }
  
  st <- Sys.time(); cmp <- max(abs(toRGB(toHSV(rgb)) - rgb)); et <- Sys.time();
  print(paste('Speed test: maxerror', cmp, ', time', et - st))
  
  st.ref <- Sys.time(); cmp.ref <- max(abs(refFunc(255*rgb) - 255*rgb)); et.ref <- Sys.time();
  print(paste('Reference speed test: maxerror', cmp.ref, ', time', et.ref - st.ref))
}

