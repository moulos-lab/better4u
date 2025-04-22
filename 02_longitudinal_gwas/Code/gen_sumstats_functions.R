generate_sumstats <- function(input_data,
                              pheno     = c("BMI", "Weight"),
                              age_group = c("Childhood", "Adulthood", "Elderly"))
  {
  
  if (age_group == "Childhood") {
    max_age <- 18
  } else if (age_group == "Adulthood") {
    min_age <- 18
    max_age <- 65
  } else if (age_group == "Elderly") {
    min_age <- 65
  } # End if 
  
  d <- input_data
  n <- dim(d)[1]
  
  if (pheno == "BMI"){
    
    d$N         <- 0*(1:n)
    d$M_Age     <- 0*(1:n)
    d$S_Age     <- 0*(1:n)
    d$M_BMI     <- 0*(1:n)
    d$S_BMI     <- 0*(1:n)
    d$BMI_Alpha <- 0*(1:n)
    d$BMI_Beta  <- 0*(1:n)
    nCol        <- dim(d)[2]
    
    for (k in 1:n) {
      
      # Split age and BMI strings into vectors
      t <- as.numeric(strsplit(d$Age[k],";")[[1]])
      
      if (age_group == "Childhood") {
        
        bmi <- suppressWarnings(as.numeric(strsplit(d$zBMI[k],";")[[1]]))
        use <- !is.na(t) & !is.na(bmi) & t <= max_age
        t   <- t[use]- max_age # offset so that Alpha is estimate of BMI at age 18
        
      } else if (age_group == "Adulthood") {
        
        bmi <- suppressWarnings(as.numeric(strsplit(d$BMI[k],";")[[1]]))
        use <- !is.na(t) & !is.na(bmi) & t > min_age & t <= max_age
        t   <- t[use] - min_age 
        
      } else if (age_group == "Elderly") {
        bmi <- suppressWarnings(as.numeric(strsplit(d$BMI[k],";")[[1]]))
        use <- !is.na(t) & !is.na(bmi) & t >= min_age
        t   <- t[use] - min_age 
        
      } # End if 
      
      # sub-set data to eligible time points
      bmi    <- bmi[use]
      d$N[k] <- length(t)
      
      if (d$N[k] >= 2) {
        d$M_Age[k] <- mean(t,   na.rm = T)
        d$M_BMI[k] <- mean(bmi, na.rm = T)
        d$S_Age[k] <- sqrt(sum(( t - mean(t, na.rm = T))^2)/( d$N[k] - 1 ))          
        d$S_BMI[k] <- sqrt(sum((bmi - mean(bmi, na.rm = T))^2)/( d$N[k] - 1 ))
        
        f              <- summary( lm( bmi ~ t ) )
        d$BMI_Alpha[k] <- f$coefficients[1,1]  
        
        if (nrow(f$coefficients) > 1) d$BMI_Beta[k]  <- f$coefficients[2,1] else d$BMI_Beta[k] <- 0 
      }# End if
      
    }# End for
    
  } else {
    
    d$N            <- 0*(1:n)
    d$M_Age        <- 0*(1:n)
    d$S_Age        <- 0*(1:n)
    d$M_Weight     <- 0*(1:n)
    d$S_Weight     <- 0*(1:n)
    d$Weight_Alpha <- 0*(1:n)
    d$Weight_Beta  <- 0*(1:n)
    d$Height_med   <- 0*(1:n)
    nCol           <- dim(d)[2]
    
    for (k in 1:n) {
      
      h               <- suppressWarnings(as.numeric(strsplit(d$Height[k],";")[[1]]))
      d$Height_med[k] <- median(h, na.rm = T)
      
      t <- as.numeric(strsplit(d$Age[k],";")[[1]])
      
      if (age_group == "Childhood") {
        
        weight <- suppressWarnings(as.numeric(strsplit(d$Weight[k],";")[[1]]))
        use <- !is.na(t) & !is.na(weight) & t <= max_age
        t   <- t[use]- max_age 
        
      } else if (age_group == "Adulthood") {
        
        weight <- suppressWarnings(as.numeric(strsplit(d$Weight[k],";")[[1]]))
        use <- !is.na(t) & !is.na(weight) & t > min_age & t <= max_age
        t   <- t[use] - min_age 
        
      } else if (age_group == "Elderly") {
        weight <- suppressWarnings(as.numeric(strsplit(d$Weight[k],";")[[1]]))
        use <- !is.na(t) & !is.na(weight) & t >= min_age
        t   <- t[use] - min_age 
        
      } # End if 
      
      # sub-set data to eligible time points
      weight <- weight[use]
      d$N[k] <- length(t)
      
      if (d$N[k] >= 2) {
        
        d$M_Age[k]    <- mean(t, na.rm = T) 
        d$M_Weight[k] <- mean(weight,na.rm = T) 
        d$S_Age[k]    <- sqrt(sum((t - mean(t, na.rm = T))^2)/(d$N[k]-1))              
        d$S_Weight[k] <- sqrt(sum((weight - mean(weight, na.rm = T))^2)/(d$N[k] - 1))
        
        f             <- summary(lm(weight ~ t))
        
        d$Weight_Alpha[k] <- f$coefficients[1,1] 
        if (nrow(f$coefficients) > 1) d$Weight_Beta[k]<- f$coefficients[2,1] else d$Weight_Beta[k] <- 0 # estimator of slope
      
        }# End if
      
    }# End for
  }# End if
  
  # Output:
  return(d)

}# End function