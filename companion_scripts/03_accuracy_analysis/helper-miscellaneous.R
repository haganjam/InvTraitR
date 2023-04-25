#'@title: PI() 
#'
#'@description: function to calculate percentile confidence/credible interval
#'from a vector of data. The function is taken from the rethinking package
#'McElreath (2020, https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
#'
#'@param samples vector of data to calculate the percentile interval on
#'@param prob size of the percentile interval (default = 0.89)
#' 

PCI <- function( samples , prob=0.89 ) {
  x <- sapply( prob , function(p) {
    a <- (1-p)/2
    quantile( samples , probs=c(a,1-a) )
  } )
  # now order inside-out in pairs
  n <- length(prob)
  result <- rep(0,n*2)
  for ( i in 1:n ) {
    low_idx <- n+1-i
    up_idx <- n+i
    # lower
    result[low_idx] <- x[1,i]
    # upper
    result[up_idx] <- x[2,i]
    # add names
    a <- (1-prob[i])/2
    names(result)[low_idx] <- paste(round(a*100,0),"%")
    names(result)[up_idx] <- paste(round((1-a)*100,0),"%")
  }
  return(result)
}
PI <- PCI

#'@title: var2() 
#'
#'@description: function to calculate variance that does not use (n-1) as
#'the denominator. The function is taken from the rethinking package
#'McElreath (2020, https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
#'
#'@param x vector of data to calculate the percentile interval on
#' 

var2 <- function( x , na.rm = TRUE ) {
  # use E(x^2) - E(x)^2 form
  mean(x^2) - mean(x)^2
}

#'@title: lm_compare() 
#'
#'@description: function to fit multiple linear models to different sets of
#'predictor variables from a list. The model outputs model coefficients and
#'model comparison metrics like r2 and AIC.
#'
#'@param data data.frame containing the response and predictor variables
#'@param resp name of the response variable
#'@param pred list made of vectors of different predictor variables
#' 

lm_compare <- function(data, resp, pred) {
  
  # set an output list for the model coefficients
  est_lm <- vector("list", length(pred))
  names(est_lm) <- seq(1:length(pred))
  
  # set an output list for the model fit statistics
  fit_lm <- vector("list", length(pred))
  names(fit_lm) <- seq_along(1:length(pred))
  
  # loop over list with different combinations of 
  for (i in 1:length(pred) ) {
    
    # fit model using chosen predictors
    lm_pred <- lm(reformulate(pred[[i]], resp), data = data)
    
    # write coefficients to the est.lm list
    est_lm[[i]] <- broom::tidy(lm_pred)
    
    # write fit statistics to the fit.lm list
    fit_lm[[i]] <- broom::glance(lm_pred)
    
  }
  
  # convert lists to data.frames and join
  mod <- full_join(bind_rows(est_lm, .id = "model"), 
                   bind_rows(fit_lm, .id = "model"),
                   by = "model")
  
  # return the data.frame with relevant information
  return(mod)
  
}

#'@title: seaweed_pal()
#'
#'@description: colour palette used for the four macroalgae species derived
#'from water colour images of these species.
#'

seaweed_pal <- function() {
  
  cols <- c("#D97E46", "#7A5414", "#AF994D", "#EAB20A")
  return(cols)
  
}

### END
