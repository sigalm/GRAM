######################################## GRAM HELPER FUNCTIONS: GENERIC ########################################

# Defines all generic functions being used throughout the model.

######################################## generate synthetic sample from microdata
# weights_var is the name of the variable in pop_data containing population weights, not in quotes

generate_synthetic_sample <- function(pop_data, target_size, weights_var = NULL, seed) {
  set.seed(seed)  # Ensure reproducibility for each dataset
  pop_data %>%
    slice_sample(n = target_size, weight_by = weights_var, replace = TRUE)
}

######################################## convert probability to different time

# formula from https://doi.org/10.1007/s40273-020-00937-z

# p = probability
# t_new = targeted time period
# t_old = current time period
# RR = risk ratio (e.g. relative risk or hazard ratio)

f.adjustprobability <- function(p, t_new, t_old, RR) {
  rate1 <- -log(1-p)                  # convert probability to rate
  rate2 <- rate1 * RR * (t_new/t_old) # adjust rate (time and risk ratio)
  out <- 1-exp(-rate2)                # convert rate to probability
  return(out)
}

######################################## sample from categorical variable using pre-defined uniform random values

# p_rand = random values from uniform distribution
# p_cat = probabilities of categorical variable
# values = values corresponding to the probabilities in p_cat (must be in same order)

f.qcat <- function(p_rand, p_cat, values = NULL) {
  
  # see supplemental 'Function f.qcat' for details
  
  if (length(p_rand)==0 | length(p_cat)==0 | any(is.na(p_rand)) | any(is.na(p_cat))) { return(NULL) # if input is NULL
  } else if (sum(p_cat)==0) { return(0) # if input is zero
  } else if (is.vector(p_cat)) { # if p_cat is of type vector
    # checks
    if( (sum(p_cat)<0.99999999 | sum(p_cat)>1.00000001) & sum(p_cat)!=0 & !anyNA(p_cat) ) { warning("vector: p_cat do not sum up to 1"); print(p_cat) }
    if(anyNA(p_cat)) warning("vector: p_cat contains NA")
    if(!is.null(values) & length(p_cat)!=length(values)) stop("vector: length of p_cat does not match length of values")
    if(length(p_cat)<=1) stop("vector: 2 or more p_cat values are required")
    if(is.null(values)) values <- 1:length(p_cat)
    if (any(p_cat == 0) | any(p_cat == 1)) {        # if p_cat contains absolute probabilities, offset them to avoid errors
      p_cat[p_cat == 0] <- 1e-10
      p_cat[p_cat == 1] <- 1 - 1e-10
    }
    
    # function
    breaks <- c(0,cumsum(p_cat)) 
    breaks[duplicated(breaks)] <- breaks[duplicated(breaks)] + 1e-10 
    out <- cut(x = p_rand, breaks = breaks, labels = FALSE, include.lowest = TRUE) # categorize random values based on user-defined breaks
    if(!is.null(values)) out <- values[out] # apply user-defined labels
    return(out)
    
  } else if (is.matrix(p_cat)) { # if p_cat is of type matrix
    # checks
    if(!is.matrix(p_cat)) stop("p_cat is not defined as matrix")
    if(any(is.na(colSums(p_cat))) ) warning("some colsums contain NA")
    if(!all(colSums(p_cat)-1<10E-10) ) warning("matrix: p_cat do not sum up to 1, might be due to NA")
    if(!is.null(values) & nrow(p_cat)!=length(values) ) stop("number of rows of p_cat matrix does not match length of values")
    if(is.null(values) ) values <- 1:nrow(p_cat)
    if(nrow(p_cat)<=1 ) stop("2 or more p_cat values are required (i.e. matrix with 2 or more columns")
    if(length(p_rand)!=ncol(p_cat) ) stop("length of p_rand differs from number of columns of p_cat")
    
    # function
    out <- rep(NA, length(p_rand))
    p_cat.cumsum <- apply(X=p_cat, MARGIN=2, FUN=cumsum)
    n <- nrow(p_cat.cumsum)-1
    out[p_rand>=0 & p_rand<=p_cat.cumsum[1,]] <- values[1]
    for (i in 1:n) {
      out[p_rand>p_cat.cumsum[i,] & p_rand<=p_cat.cumsum[i+1,]] <- values[i+1]
    }
    return(out)
    
  } else stop("p_cat is not vector or matrix")
  
}

######################################## discounting

f.discount <- function(x, discount_rate, n.cycle) {
  # x can be scalar, vector or matrix
  as.matrix(x) / (1 + discount_rate)^(0:(n.cycle - 1))
}

