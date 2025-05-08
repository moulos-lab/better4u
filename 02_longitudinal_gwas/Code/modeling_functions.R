modeling_fun <- function(pheno_data,
                         geno_data,
                         SNP_IDs,
                         pheno    = c("BMI", "Weight")#,
                         #numcores = 4
                         )
{
  
  BM2FBM <- function(bm) {
    FBM(nrow        = nrow(bm), 
        ncol        = ncol(bm), 
        type        = typeof(bm),
        backingfile = file.path(bigmemory::dir.name(bm), bigstatsr::sub_bk(file.name(bm))),
        create_bk   = FALSE)
  }
  
  output        <- vector("list", 2)
  names(output) <- c("Effect_sizes", "Variance_ratio")
  
  if (pheno == "BMI"){
    
    #-------------------------------------------------------------------------------------------------------
    # Estimate relative contributions of variance in individual betas from sampling error
    # and inter-individual variation in beta, using a two-stage approech.
    # First fit null model and obtain residuals.
    # Second, calculate index of inverse variance of sampling error and fit linear model
    # square(total residual)*variance(individual beta) ~ variance(individual beta)
    # Calculate sample weights using the ratio of fitted coefficients
    #-------------------------------------------------------------------------------------------------------
    cat("> Estimating variance ratio...")
    
    model_phenotype     <- "BMI_Beta"
    model_covar_nogenos <- paste(c("M_Age*Sex", paste("PC",as.character(1:10), sep = "")), collapse = " + ")
    formula_model0      <- paste(model_phenotype, model_covar_nogenos, sep = " ~ ")
    
    model0              <- lm(formula_model0, data = pheno_data)
    summary_m0          <- summary(model0)
    
    pheno_data$Sample_Weight <- (pheno_data$N-1)*(pheno_data$S_Age^2) # Interaction of gender with age^2
    dd              <- as.data.frame(cbind(pheno_data$Sample_Weight, 
                                           pheno_data$Sample_Weight * (model0$residuals^2)))
    colnames(dd)    <- c("w","y")
    model1          <- MASS::rlm(cbind(0*dd$w + 1,dd$w), dd$y)
    summary_m1      <- summary(model1)
    
    pheno_data$Sample_Weight <- 1 / ( (1/pheno_data$Sample_Weight) + summary_m1$coefficients[2] / summary_m1$coefficients[1] )
    
    output$Variance_ratio <- as.numeric(model1$coefficients[2]/model1$coefficients[1])
    
    cat("> Estimating SNP effects...")
    
    # Extract phenotype and covariates:
    pheno_data$M_Age_Sex <- pheno_data$M_Age * pheno_data$Sex # Interaction term for Mid-age and gender
    
    keep_cols_cols <- c("M_Age_Sex",
                        paste0("PC", 1:10))
    
    y          <- pheno_data$BMI_Beta
    covariates <- as.matrix(pheno_data[, keep_cols_cols])
    
    # Run per-SNP regression adjusting for age and gender
    assoc_results <- bigstatsr::big_univLinReg(X      = geno_data, 
                                               y, 
                                               covar  = covariates#,
                                               #ncores = numcores
                                               )
    
    assoc_results$p.value   <- predict(assoc_results, log10 = FALSE)
    assoc_results$score     <- NULL
    colnames(assoc_results) <- c("Beta","SE","P")
    

   #---- Transform the linear regression coefs and SEs and then calculate the p-value using the X^2 distribution:
    y_w          <- pheno_data$BMI_Beta * sqrt(pheno_data$Sample_Weight)
    covariates_w <- sweep(pheno_data[, keep_cols_cols], 1, sqrt(pheno_data$Sample_Weight), "*")
    covariates_w <- as.matrix(covariates_w)
    
    geno_data_w <- geno_data$bm()
    geno_data_w <- BM2FBM(geno_data_w)

    out <- bigstatsr::big_univLinReg(X      = geno_data_w, 
                                     y_w, 
                                     covar  = covariates_w#,
                                     #ncores = numcores
                                     )
    
    out$p.value <- predict(out, log10 = FALSE)
    out$score <- NULL
    colnames(out) <- c("Beta","SE","P")
    
    rm(assoc_results)
    # head(out)

    output$Effect_sizes <- cbind(SNP_IDs, out)
    
  } else {
    
    #-------------------------------------------------------------------------------------------------------
    # Estimate relative contributions of variance in individual betas from sampling error
    # and inter-individual variation in beta, using a two-stage approech.
    # First fit null model and obtain residuals.
    # Second, calculate index of inverse variance of sampling error and fit linear model
    # square(total residual)*variance(individual beta) ~ variance(individual beta)
    # Calculate sample weights using the ratio of fitted coefficients
    #-------------------------------------------------------------------------------------------------------
    
    cat("> Estimating variance ratio...")
    
    model_phenotype     <- "Weight_Beta"
    model_covar_nogenos <- paste(c("M_Age*Sex", "Height_med", paste("PC",as.character(1:10), sep = "")), collapse = " + ")
    formula_model0      <- paste(model_phenotype, model_covar_nogenos, sep = " ~ ")
    
    model0              <- lm(formula_model0, data = pheno_data)
    summary_m0          <- summary(model0)
    
    pheno_data$Sample_Weight <- (pheno_data$N - 1)*(pheno_data$S_Age^2) # Interaction of gender with age^2
    dd              <- as.data.frame(cbind(pheno_data$Sample_Weight, 
                                           pheno_data$Sample_Weight * (model0$residuals^2)))
    colnames(dd)    <- c("w","y")
    model1          <- MASS::rlm(cbind(0*dd$w + 1,dd$w), dd$y)
    summary_m1      <- summary(model1)
    
    pheno_data$Sample_Weight <- 1 / ( (1/pheno_data$Sample_Weight) + summary_m1$coefficients[2] / summary_m1$coefficients[1] )
    
    output$Variance_ratio <- as.numeric(model1$coefficients[2]/model1$coefficients[1])
    
    cat("> Estimating SNP effects...")
    
    # Extract phenotype and covariates:
    pheno_data$M_Age_Sex <- pheno_data$M_Age * pheno_data$Sex # Interaction term for Mid-age and gender
    
    keep_cols_cols <- c("M_Age_Sex",
                        "Height_med",
                        paste0("PC", 1:10))
    
    y          <- pheno_data$Weight_Beta
    covariates <- as.matrix(pheno_data[, keep_cols_cols])
    
    ###########################
    # Run per-SNP regression adjusting for age and gender
    assoc_results <- bigstatsr::big_univLinReg(X      = geno_data, 
                                               y, 
                                               covar  = covariates#,
                                               #ncores = numcores
    )
    
    assoc_results$p.value   <- predict(assoc_results, log10 = FALSE)
    assoc_results$score     <- NULL
    colnames(assoc_results) <- c("Beta","SE","P")
    
    #---- Transform the linear regression coefs and SEs and then calculate the p-value using the X^2 distribution:
    y_w          <- pheno_data$Weight_Beta * sqrt(pheno_data$Sample_Weight)
    covariates_w <- sweep(pheno_data[, keep_cols_cols], 1, sqrt(pheno_data$Sample_Weight), "*")
    covariates_w <- as.matrix(covariates_w)
    
    geno_data_w <- geno_data$bm()
    geno_data_w <- BM2FBM(geno_data_w)

    out <- bigstatsr::big_univLinReg(X      = geno_data_w, 
                                     y_w, 
                                     covar  = covariates_w#,
                                     #ncores = numcores
    )
    
    out$p.value <- predict(out, log10 = FALSE)
    out$score <- NULL
    colnames(out) <- c("Beta","SE","P")
    
    rm(assoc_results)
    # head(out)
    
    output$Effect_sizes <- cbind(SNP_IDs, out)
    
  }# End if
  
  # Output:
  
  return(output)
  
}# End function