#############################################################
#                                                           #
#   circular function                                       #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 19, 2009                                 #
#   Version: 0.8                                            #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################

circular <- function(x, type=c("angles", "directions"), units=c("radians", "degrees", "hours"), template=c("none", "geographics", "clock12", "clock24"), modulo=c("asis", "2pi", "pi"), zero=0, rotation=c("counter", "clock"), names=NULL) {

  type <- match.arg(type)
  units <- match.arg(units)
  template <- match.arg(template) 
  modulo <- match.arg(modulo)
  rotation <- match.arg(rotation)
  
  if (template=="geographics") {
    zero <- pi/2
    rotation <- "clock"
  } else if (template=="clock24") {
    zero <- pi/2
    rotation <- "clock"
  } else if (template=="clock12") {
    zero <- pi/2
    rotation <- "clock"
  }

  if (is.data.frame(x))
    x <- as.matrix(x)
  cl <- class(x)
    
  if (is.matrix(x)) {
    nseries <- ncol(x)
    if (is.null(names)) {
      if (is.null(colnames(x)))
        colnames <- paste("Circular", seq(nseries), sep="")
    } else 
      colnames(x) <- names
  } else {
    nseries <- 1
  }
  if (modulo!="asis") {
    if (modulo=="2pi") {
        ang <- 2
    } else {
        ang <- 1
    }
    if (units=="radians") {
      x <- x %% (ang*pi)
    } else if (units=="degrees") {
      x <- x %% (ang*180)
    } else {
      x <- x %% (ang*12) ## hours
    }
  }

  attr(x, "circularp") <- list(type=type, units=units, template=template, modulo=modulo, zero=zero, rotation=rotation) #-- order is fixed
  if (!inherits(x, "circular"))
    class(x) <- c("circular", cl)
  return(x)
}

#############################################################
#                                                           #
#   is.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: March, 7, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

is.circular <- function (x) inherits(x, "circular")

#############################################################
#                                                           #
#   circularp function                                      #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: March, 7, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

circularp <- function(x) attr(x, "circularp")

#############################################################
#                                                           #
#   circularp<- function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 19, 2009                                 #
#   Version: 0.3                                            #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################

"circularp<-" <- function(x, value) {
  cl <- class(x)
  if (length(value)!=6) stop("value must have six elements")
  
  if (is.list(value)) {
    type <- value$type
    units <- value$units
    template <- value$template
    modulo <- value$modulo
    zero <- value$zero
    rotation <- value$rotation
  } else {
    type <- value[1]
    units <- value[2]
    template <- value[3]
    modulo <- value[4]
    zero <- as.numeric(value[5])
    rotation <- value[6]
    value <- list(type=type, units=units, template=template, modulo=modulo, zero=zero, rotation=rotation)
  }
  
  if (type!="angles" & type!="directions") stop("type (value[1]) must be 'angles', 'directions'")
  
  if (units!="radians" & units!="degrees" & units!="hours") stop("units (value[2]) must be 'radians' or 'degrees' or 'hours'")
  
  if (template!="none" & template!="geographics" & template!="clock24" & template!="clock12") stop("template (value[3]) must be 'none' or 'geographics' or 'clock24' or 'clock12'")
  
  if (modulo!="asis" & modulo!="2pi" &  modulo!="pi") stop("modulo (value[4]) must be 'asis' or 'pi' or '2pi'")
  
  if (!is.numeric(zero)) stop("zero (value[5]) must be numeric")
  
  if (rotation!="clock" & rotation!="counter") stop("rotation (value[6]) must be 'clock' or 'counter'")
  
  attr(x, "circularp") <- value
  if (inherits(x, "circular") && is.null(value))
    class(x) <- if (!identical(cl, "circular")) cl["circular" != cl]
  return(x)
}

#############################################################
#                                                           #
#   as.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 31, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#############################################################

as.circular <- function (x, control.circular=list(), ...) {
  if (is.circular(x))
    return(x)
  else if (!is.null(xcircularp <- circularp(x)))
    circular(x, type=xcircularp$type, units=xcircularp$units, template=xcircularp$template, modulo=xcircularp$modulo, zero=xcircularp$zero, rotation=xcircularp$rotation)
  else {
    warntype <- warnunits <- warntemplate <- warnrotation <- warnmodulo <- warnzero <- ""
    printwarn <- FALSE
    dotc <- list(..., expand=TRUE)
    dc <- control.circular
    if (is.null(dc$type)) {
      if (!is.null(dotc$type)) 
        dc$type <- dotc$type
      else {
        dc$type <- "angles"
        warntype <- "  type: 'angles'\n"
        printwarn <- TRUE
      }
    }
    if (is.null(dc$units)) {
      if (!is.null(dotc$units))
        dc$units <- dotc$units
      else {
        dc$units <- "radians"
        warnunits <- "  units: 'radians'\n"
        printwarn <- TRUE
      }
    }
    if (is.null(dc$template)) {
      if (!is.null(dotc$template))
        dc$template <- dotc$template
      else {
        dc$template <- "none"
        warntemplate <- "  template: 'none'\n"
        printwarn <- TRUE
      }
    }
    if (is.null(dc$modulo)) {
      if (!is.null(dotc$modulo))
        dc$modulo <- dotc$modulo
      else {
        dc$modulo <- "asis"
        warnmodulo <- "  modulo: 'asis'\n"
        printwarn <- TRUE
      }
    }
    if (is.null(dc$zero)) {
      if (!is.null(dotc$zero))
        dc$zero <- dotc$zero
      else {
        dc$zero <- 0
        warnzero <- "  zero: 0\n"
        printwarn <- TRUE
      }
    }
    if (is.null(dc$rotation)) {
      if (!is.null(dotc$rotation))
        dc$rotation <- dotc$rotation
      else {
        dc$rotation <- "counter"
        warnrotation <- "  rotation: 'counter'\n"
        printwarn <- TRUE
      }
    }
    if (printwarn) {
      warn <- paste("an object is coerced to the class 'circular' using default value for the following components:\n", warntype, warnunits, warntemplate, warnmodulo, warnzero, warnrotation, sep="")
      warning(warn, sys.call(-1))
    }
    circular(x, type=dc$type, units=dc$units, template=dc$template, modulo=dc$modulo, zero=dc$zero, rotation=dc$rotation)
  }
}

#############################################################
#                                                           #
#   conversion.circular function                            #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: December, 16, 2009                                #
#   Version: 0.3-1                                          #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################

conversion.circular <- function(x, units=c("radians", "degrees", "hours"), type=NULL, template=NULL, modulo=NULL, zero=NULL, rotation=NULL) {
    units <- match.arg(units)
    if (!is.null(type) && type!="angles" && type!="directions")
       stop("'type' must be 'angles' or 'directions' or NULL")
    if (!is.null(template) && template!="none" && template!="geographics" && template!="clock24" && template!="clock12")
       stop("'template' must be 'none' or 'geographics' or 'clock24' or 'clock12' or NULL")
    if (!is.null(modulo) && modulo!="asis" && modulo!="2pi" && modulo!="pi")
       stop("'modulo' must be 'asis' or 'pi' or '2pi' or NULL")
    if (!is.null(zero) && !is.numeric(zero))
       stop("'zero' must be numeric or NULL")
    if (!is.null(rotation) && rotation!="clock" && rotation!="counter")
       stop("'rotation' must be 'clock' or 'counter' or NULL")
    x <- as.circular(x)
    value <- attr(x, "circularp")
    typep <- value$type
    unitsp <- value$units
    rotationp <- value$rotation
    zerop <- value$zero
    if (!is.null(template)) {
       if (template=="geographics") {
          zero <- pi/2
          rotation <- "clock"
       } else if (template=="clock24") {
          zero <- pi/2
          rotation <- "clock"
       } else if (template=="clock12") {
          zero <- pi/2
          rotation <- "clock"
       }
       value$template <- template 
    }
    if (!is.null(type) && type=="directions" && typep!=type) {
       x <- 2*x
       value$type <- type
    }    
    if (!is.null(units)) {
       if (unitsp=="degrees" & units=="radians") {
          x <- x/180*pi
       } else if (unitsp=="radians" & units=="degrees") {
          x <- x/pi*180
       } else if (unitsp=="degrees" & units=="hours") {
          x <- x/180*12
       } else if (unitsp=="radians" & units=="hours") {
          x <- x/pi*12
       } else if (unitsp=="hours" & units=="degrees") {
          x <- x/12*180
       } else if (unitsp=="hours" & units=="radians") {
          x <- x/12*pi
       }
       value$units <- units
    }

    if (!is.null(zero) && zerop!=zero) {
       if (units=="degrees") {
          zerod <- zero*180/pi
          zeropd <- zerop*180/pi
       } else if (units=="hours") {
          zerod <- zero*12/pi
          zeropd <- zerop*12/pi
       } else {
          zerod <- zero
          zeropd <- zerop
       }
       
       if (rotationp=="counter") {
          x <- x + zeropd - zerod
       } else {
          x <- x - zeropd + zerod
       }
       value$zero <- zero
    }
    
    if (!is.null(rotation) && rotationp!=rotation) {
       x <- -x
       value$rotation <- rotation
    }

    if (!is.null(modulo) && modulo!="asis") {
       if (modulo=="2pi") {
          ang <- 2
       } else {
          ang <- 1
       }
       if (units=="radians") {
         x <- x %% (ang*pi)
       } else if (units=="degrees") {
         x <- x %% (ang*180)
       } else {
         x <- x %% (ang*12) ## time
       }
    }
    if (!is.null(modulo))
       value$modulo <- modulo
    if (!is.null(zero) && zero%%(2*pi)!=pi/2)
       value$template <- "none"
    if (!is.null(rotation) && rotation=="counter")
       value$template <- "none"
    circularp(x) <- value
    return(x)
}

#############################################################
#                                                           #
#   c.circular function                                     #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: March, 30, 2011                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#############################################################

c.circular <- function (..., recursive = FALSE) {
  x <- list(...)
  value <- attr(x[[1]], "circularp")
  n <- length(x)
  if (n>1) {
    for (i in 2:length(x)) {
      x[[i]] <- conversion.circular(x[[i]], type=value$type, units=value$units, template=value$template, modulo=value$modulo, zero=value$zero, rotation=value$rotation)
    }
  }
  x <- structure(c(unlist(lapply(x, unclass))), class = c("circular", "numeric"))
  attr(x, "circularp") <- value
  return(x)
}

#############################################################
#                                                           #
#   [.circular function                                     #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: March, 7, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
"[.circular" <- function(x, i, ...) {
    y <- NextMethod("[", ...)
    class(y) <- class(x)
    attr(y, "circularp") <- attr(x, "circularp")
    return(y)
}

#############################################################
#                                                           #
#   print.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: June, 21, 2003                                    #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
print.circular <- function(x, info=TRUE, ...) {
    x.orig <- x
    x <- as.circular(x)
    if (info) {
        xcircularp <- attr(x, "circularp")
        type <- xcircularp$type
        units <- xcircularp$units
        template <- xcircularp$template
        modulo <- xcircularp$modulo
        zero <- xcircularp$zero
        rotation <- xcircularp$rotation

        cat("Circular Data: \nType =", type,
               "\nUnits =", units,
               "\nTemplate =", template, 
               "\nModulo =", modulo,
               "\nZero =", zero,
               "\nRotation =", rotation, "\n")
    }       
    attr(x, "class") <- attr(x, "circularp") <- attr(x, "na.action") <- NULL
    NextMethod("print", x, quote = FALSE, right = TRUE, ...)
    invisible(x.orig)
}
