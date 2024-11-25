# Calculate the Freeman-Tukey discrepancy measure 
# Note: the arguments x,y may be matrices

get_ft <- function(x,y) {
  if (any(x<0) | any(y<0)) {
    stop("Error: Freeman-Tukey discrepancy requires all values be non-negative.")
  }
  sum((sqrt(x) - sqrt(y))^2)
}

