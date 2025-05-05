
# input_data = sumstats_weight_adulthood
# colnames(input_data)[1:15]

modeling_fun <- function(input_data,
                         pheno    = c("BMI", "Weight"),
                         numcores = 6)
{
  
  `%ni%` = Negate(`%in%`)
  
  output        <- vector("list", 2)
  names(output) <- c("Effect_sizes", "Variance_ratio")
  
  # For a linux environment:
  # https://stackoverflow.com/questions/62644063/how-to-run-an-lapply-in-parallel-in-r
  cl            <- makeCluster(numcores, type = "FORK")
  
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
    
    model0              <- lm(formula_model0, data = input_data)
    summary_m0          <- summary(model0)
    
    input_data$Sample_Weight <- (input_data$N-1)*(input_data$S_Age^2) # Interaction of gender with age^2
    dd              <- as.data.frame(cbind(input_data$Sample_Weight, 
                                           input_data$Sample_Weight * (model0$residuals^2)))
    colnames(dd)    <- c("w","y")
    model1          <- MASS::rlm(cbind(0*dd$w + 1,dd$w), dd$y)
    summary_m1      <- summary(model1)
    
    input_data$Sample_Weight <- 1 / ( (1/input_data$Sample_Weight) + summary_m1$coefficients[2] / summary_m1$coefficients[1] )
    
    output$Variance_ratio <- as.numeric(model1$coefficients[2]/model1$coefficients[1])
    
    cat("> Estimating SNP effects...")
    remove_cols <- c("IID",
                     "BMI_Beta", 
                     "M_Age",
                     "S_Age",
                     "N",
                     "Sample_Weight",
                     "Sex", 
                     paste0("PC", 1:10))
    model_covar_genos <- colnames(input_data)[colnames(input_data) %ni% c(remove_cols)]
    
    # Takes a few secs:
    SNP_mapping    <- data.frame(Original_name = model_covar_genos,
                                 Mapping       = paste0("SNP", 1:length(model_covar_genos)))
    
    nonSNP_mapping <- data.frame(Original_name = remove_cols,
                                 Mapping       = remove_cols)
    
    all_mapping    <- rbind(nonSNP_mapping, SNP_mapping)
    
    colnames(input_data) <- all_mapping[match(colnames(input_data), all_mapping[,1]), 2]
    
    #colnames(input_data)[1:20]
    #colnames(input_data)[(length(colnames(input_data))-5):length(colnames(input_data))]
    
    fitAllSNPs <-
      pbapply::pblapply(SNP_mapping$Mapping[1:length(SNP_mapping$Mapping)], 
                        function(SNP, 
                                 input_data,
                                 model_covar_nogenos,
                                 model_phenotype){
                          model1_covars     <- paste(model_covar_nogenos, SNP,       sep = " + ") 
                          formula_model1    <- paste(model_phenotype, model1_covars, sep = " ~ ")
                          
                          fit <- lm(formula_model1,
                                    data    = input_data,
                                    weights = input_data$Sample_Weight
                          )
                          
                          summary_m2    <- summary(fit)$coefficients
                          extract_stats <- summary_m2[rownames(summary_m2)[rownames(summary_m2) %in% SNP],]
                          
                          # If the SNP for all samples is homozygous (2), or heterozygous (1), 
                          # or there is no variant detected (0), there is no variability to be
                          # modelled and the effect estimates are NA. Set them to -9999.
                          tmp_names <- colnames(extract_stats)
                          if( length(extract_stats) == 0)  extract_stats <- rbind(extract_stats,
                                                                                  rep(-9999, 4))
                          
                          return(extract_stats)
                        }, 
                        input_data          = input_data,
                        model_covar_nogenos = model_covar_nogenos,
                        model_phenotype     = model_phenotype,
                        cl = cl)
    
    out           <- data.frame(do.call(rbind, fitAllSNPs))
    out$"t.value" <- NULL
    colnames(out) <- c("Beta","SE","P")
    
    # effect_size_output$AF <- mean(input_data$geno)/2 # Anders should explain the rationale for calculating the Allelic Frequency
    
    out                 <- cbind(SNP_mapping$Original_name, out)
    colnames(out)[1]    <- "SNP"
    output$Effect_sizes <- out
    
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
    
    model0              <- lm(formula_model0, data = input_data)
    summary_m0          <- summary(model0)
    
    input_data$Sample_Weight <- (input_data$N-1)*(input_data$S_Age^2) # Interaction of gender with age^2
    dd              <- as.data.frame(cbind(input_data$Sample_Weight, 
                                           input_data$Sample_Weight * (model0$residuals^2)))
    colnames(dd)    <- c("w","y")
    model1          <- MASS::rlm(cbind(0*dd$w + 1,dd$w), dd$y)
    summary_m1      <- summary(model1)
    
    input_data$Sample_Weight <- 1 / ( (1/input_data$Sample_Weight) + summary_m1$coefficients[2] / summary_m1$coefficients[1] )
    
    output$Variance_ratio <- as.numeric(model1$coefficients[2]/model1$coefficients[1])
    
    cat("> Estimating SNP effects...")
    
    remove_cols <- c("IID", 
                     "Weight_Beta", 
                     "M_Age",
                     "S_Age",
                     "N",
                     "Sample_Weight",
                     "Sex", 
                     "Height_med",
                     paste0("PC", 1:10))
    
    model_covar_genos <- colnames(input_data)[colnames(input_data) %ni% c(remove_cols)]
    
    # Takes a few secs:
    SNP_mapping    <- data.frame(Original_name = model_covar_genos,
                                 Mapping       = paste0("SNP", 1:length(model_covar_genos)))
    
    nonSNP_mapping <- data.frame(Original_name = remove_cols,
                                 Mapping       = remove_cols)
    
    all_mapping    <- rbind(nonSNP_mapping, SNP_mapping)
    
    colnames(input_data) <- all_mapping[match(colnames(input_data), all_mapping[,1]), 2]
    
    #colSums(input_data[c("SNP1", "SNP2", "SNP3", "SNP4")])
    
    fitAllSNPs <-
      pbapply::pblapply(SNP_mapping$Mapping[1:length(SNP_mapping$Mapping)], 
                        function(SNP, 
                                 input_data,
                                 model_covar_nogenos,
                                 model_phenotype){
                          model1_covars     <- paste(model_covar_nogenos, SNP,       sep = " + ") 
                          formula_model1    <- paste(model_phenotype, model1_covars, sep = " ~ ")
                          
                          fit <- lm(formula_model1,
                                    data    = input_data,
                                    weights = input_data$Sample_Weight
                          )
                          
                          summary_m2    <- summary(fit)$coefficients
                          extract_stats <- summary_m2[rownames(summary_m2)[rownames(summary_m2) %in% SNP],]
                          
                          # If the SNP for all samples is homozygous (2), or heterozygous (1), 
                          # or there is no variant detected (0), there is no variability to be
                          # modelled and the effect estimates are NA. Set them to -9999.
                          tmp_names <- colnames(extract_stats)
                          if( length(extract_stats) == 0)  extract_stats <- rbind(extract_stats,
                                                                                  rep(-9999, 4))
                          
                          return(extract_stats)
                        }, 
                        input_data          = input_data,
                        model_covar_nogenos = model_covar_nogenos,
                        model_phenotype     = model_phenotype,
                        cl = cl)
    
    out           <- data.frame(do.call(rbind, fitAllSNPs))
    out$"t.value" <- NULL
    colnames(out) <- c("Beta","SE","P")
    
    # effect_size_output$AF <- mean(input_data$geno)/2 # Anders should explain the rationale for calculating the Allelic Frequency
    
    out                 <- cbind(SNP_mapping$Original_name, out)
    colnames(out)[1]    <- "SNP"
    output$Effect_sizes <- out
    
  }# End if
  
  # Output:
  
  return(output)
  
  stopCluster(cl)
  
}# End function