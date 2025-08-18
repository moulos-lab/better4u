library(parallel)
library(ggplot2)
library(reshape2)

gridSearch <- function(pip=c(0,0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4),
    beta=c(0,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3),metric="prs_r2",prsFile,
    covFile,trait,genoBase,iidCol=2,sum=TRUE,center=FALSE,
    plink=Sys.which("plink"),rc=NULL) {
    # Grid parameters: BETA, PIP, 
    # Write a prsFile for each grid element (griFile)
    # Multiple evalPrs
    # Goal: create arrays of objective function metrics
    
    if (!(metric %in% .metricNames()))
        stop("metric must be one of the known recorded metrics: \n",
            paste(.metricNames(),collapse=", "))
    
    message("Reading PRS file...")
    prs <- read.delim(prsFile)
    #message("Calculating full PRS metrics...")
    #M <- evalPrs(prsFile,covFile,trait,genoBase,iidCol,sum,center,plink)
    
    # Initialize the objective function result arrays
    #obj <- array(NA,dim=c(length(beta)+1,length(pip)+1,length(M$metrics)))
    #rownames(obj) <- c("abs(BETA)>0",paste("abs(BETA)>",beta,sep=""))
    #colnames(obj) <- c("PIP>0",paste("PIP>",pip,sep=""))
    obj <- array(NA,dim=c(length(beta),length(pip),length(.metricNames())))
    rownames(obj) <- paste("abs(BETA)>",beta,sep="")
    colnames(obj) <- paste("PIP>",pip,sep="")
    dimnames(obj)[[3]] <- .metricNames()
    
    #obj["abs(BETA)>0","PIP>0",] <- M$metrics
    
    message("Doing grid search...")
    for (p in pip) {
        message("  PIP > ",p)
        inner <- cmclapply(beta,function(b) {
            message("    abs(BETA) > ",b)
            gprs <- prs[prs$PIP>p & abs(prs$BETA) > b,,drop=FALSE]
            #message("    abs(BETA) > ",b," : ",nrow(gprs)," SNPs")
            griFile <- tempfile()
            write.table(gprs,file=griFile,sep="\t",row.names=FALSE,quote=FALSE)
            return(evalPrs(griFile,covFile,trait,genoBase,iidCol,sum,center,
                plink))
        },rc=rc)
        names(inner) <- paste("abs(BETA)>",beta,sep="")
        for (b in beta)
            obj[paste0("abs(BETA)>",b),paste0("PIP>",p),] <- 
                inner[[paste0("abs(BETA)>",b)]]$metrics
            
        #for (b in beta) {
        #    message("  abs(BETA) > ",b)
        #    gprs <- prs[prs$PIP>p & abs(prs$BETA) > b,drop=FALSE,]
        #    griFile <- tempfile()
        #    write.table(gprs,file=griFile,sep="\t",row.names=FALSE,quote=FALSE)
        #    z <- evalPrs(griFile,covFile,trait,genoBase,iidCol,sum,center,plink)
        #    obj[paste0("abs(BETA)>",b),paste0("PIP>",p),] <- z$metrics
        #}
    }
    
    bestIndex <- .whichMaxLast(obj[,,metric,drop=FALSE])
    message("Best PRS according to ",metric," found at:")
    message("  abs(BETA) > ",beta[bestIndex$ind[1]])
    message("  PIP > ",pip[bestIndex$ind[2]])
    
    # Output some kind of matrix with metrics? R2? PRS R2? PRS p-value
    # Suitable for a heatmap?
    return(list(
        metrics=obj,
        prs=prs[prs$PIP>pip[bestIndex$ind[2]] 
            & abs(prs$BETA) > beta[bestIndex$ind[1]],,drop=FALSE]
    ))
}

gridSearchPlot <- function(obj,what,i=NULL,j=NULL,dec=3) {
    #mat <- as.matrix(obj[,,what,drop=FALSE])
    mat <- as.matrix(obj[,,what])
    rownames(mat) <- rownames(obj[,,what,drop=FALSE])
    colnames(mat) <- colnames(obj[,,what,drop=FALSE])
    
    if (what %in% c("prs_pvalue","prs_pheno_cor_p"))
        mat <- -log10(mat)
    
    n <- dim(mat)[1]
    labs <- matrix(sprintf(paste0("%.",dec,"f"),mat),nrow=n)
    
    # Convert matrix to long data frame
    df <- melt(mat)
    colnames(df) <- c("Row","Col","Value")

    # Optional: melt labs (must be same dim as mat)
    labs_df <- melt(labs)
    colnames(labs_df) <- c("Row","Col","Label")
    df$Label <- labs_df$Label

    # Add text color based on Value (lighter background → black text, 
    # darker → white text)
    fill_palette <- colorRampPalette(c("yellow","grey","blue"))(100)
    fill_mid <- mean(range(df$Value,na.rm=TRUE))
    df$TextColor <- ifelse(df$Value<fill_mid,"black","white")

    # Plot heatmap
    p <- ggplot(df,aes(x=Col,y=Row,fill=Value)) +
        geom_tile(color="white") +
        geom_text(aes(label=Label,color=TextColor),size=4,show.legend=FALSE) +
        scale_fill_gradientn(colors=fill_palette) +
        scale_color_identity() +  # tells ggplot to use color values from data
        theme_minimal(base_size=14) +
        theme(
            axis.text.x=element_text(angle=45,hjust=1),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            legend.position="right",
            plot.margin=margin(5,5,5,5,"pt")
        ) +
        coord_fixed()
  
    if (!is.null(i) && !is.null(j)) {
        highlight <- df[df$Row==rownames(mat)[i] 
            & df$Col==colnames(mat)[j],]
        p <- p + geom_tile(data=highlight,aes(x=Col,y=Row),fill=NA,
            color="red",linewidth=1)
    }
  
    return(p)
}

sanitizePrs <- function(prsFile,genoBase,from=c("sbayesrc","prscs"),pip=0.001) {
    from <- from[1]
    
    # Read SNP data (bim) and initial PRS
    bimFile <- paste0(genoBase,".bim")
    message("Reading BIM ",bimFile)
    bim <- read.delim(bimFile,header=FALSE)
    message("Reading PRS ",prsFile)
    prs <- read.delim(prsFile,header=from=="sbayesrc")

    # Rename BIM columns for clarity
    colnames(bim) <- c("CHR","SNP","CM","BP","A1_bim","A2_bim")
    
    # Assign header names if prscs
    if (from=="prscs")
        colnames(prs) <- c("CHR","SNP","BP","A1","A2","BETA")

    # Merge by SNP
    message("Merging by SNP id")
    merged <- merge(prs,bim[,c("SNP","A1_bim","A2_bim")],by="SNP")
    # Consider reporting coverage somewhere...
    coverage <- nrow(merged)/nrow(prs)
    message("PRS coverage is ",paste0(round(100*coverage,2),"%")," (",
        nrow(merged)," out of ",nrow(prs),") SNPs")

    # Find mismatches
    mismatch <- merged$A1 != merged$A1_bim

    # Flip beta and replace A1 where mismatch
    if (any(mismatch)) {
        message("Found ",length(which(mismatch)),
            " allele mismatches - correcting...")
        merged$BETA[mismatch] <- -merged$BETA[mismatch]
        merged$A1[mismatch] <- merged$A1_bim[mismatch]
    }

    # Output final PRS file (with updated A1 and BETA) - we fake the SE, PIP
    # column in the case of PRS-CS
    if (from=="prscs") {
        merged$SE <- 0
        merged$PIP <- pip
    }
    output <- merged[,c("SNP","A1","BETA","SE","PIP")]
    outFile <- paste0(prsFile,".san")
    write.table(output,file=outFile,quote=FALSE,sep="\t",row.names=FALSE)
    
    message("Sanitized PRS is in ",outFile)
    return(outFile)
}

evalPrs <- function(prsFile,covFile,trait,genoBase,iidCol=2,sum=TRUE,
    center=FALSE,plink=Sys.which("plink")) {
    # Base name for plink score output
    prsName <- sub("\\.[^.]*$","",prsFile)
    
    # Find plink - MUST be in PATH or provided
    if (is.null(plink) || plink == "" || !file.exists(plink))
        stop("plink not found, please provide path to plink 1.90.")
    
    # Read the covariates
    covars <- read.delim(covFile)
    rownames(covars) <- covars[,iidCol]
    covars <- covars[complete.cases(covars),,drop=FALSE]
    
    # Align with genotypes
    fam <- read.table(paste0(genoBase,".fam"))
    rownames(fam) <- fam[,2]
    common <- intersect(rownames(fam),rownames(covars))
    # length(common) can be only smaller or equal to nrow(covars). If smaller, 
    # align with common which is based on fam and possibly re-align later with 
    # remaining fam. If all fam included in common, again, no problem this way
    covars <- covars[common,,drop=FALSE]
    # Continue with fam check
    remFile <- NULL
    if (length(common) < nrow(fam)) { # Can be only smaller or equal
        # File for removal of samples in plink --score
        remove <- setdiff(rownames(fam),common)
        remFile <- tempfile()
        write.table(fam[remove,c(1,2),drop=FALSE],file=remFile,col.names=FALSE,
            row.names=FALSE,quote=FALSE)
        # Nothing to do, covars already aligned with common
    }
    else # All found, final align the covars based on fam
        covars <- covars[rownames(fam),]
    # Remove IID, FID from covars now, not needed and will affect later GLM
    if ("FID" %in% colnames(covars))
        covars <- covars[,-which(colnames(covars)=="FID")]
    if ("IID" %in% colnames(covars))
        covars <- covars[,-which(colnames(covars)=="IID")]
    
    # Ready to run plink --score. If plink found or properly provided, run 
    # --score and exclude samples if required.
    #command <- paste(
    #    paste0(plink," \\"),
    #    paste0("  --bfile ",genoBase," \\"),
    #    paste0("  --score ",prsFile," 1 2 3 header",ifelse(sum," sum",""),
    #        ifelse(center," center","")," \\"),
    #    paste0("  --out ",prsName," \\"),
    #    ifelse(is.null(remFile),"",paste0("  --remove ",remFile," \\")),
    #    "  --silent",
    #    sep="\n"
    #)
    #message("\nExecuting:\n",command)
    args <- c("--bfile",genoBase,"--score",prsFile,"1 2 3 header",
        ifelse(sum,"sum",""),ifelse(center,"center",""),"--out",prsName)
    if (!is.null(remFile))
        args <- c(args,"--remove",remFile)
    args <- c(args,"--silent")
    out <- tryCatch({
        suppressWarnings(system2(plink,args=args))
        TRUE
    },error=function(e) {
        message("Caught error: ",e$message)
        return(FALSE)
    },finally="")

    if (!out)
        stop("Failed to generate score file! Exiting...")
    
    # If all ok, read the score file
    scoreFile <- paste0(prsName,".profile")
    theScore <- read.table(scoreFile,row.names=2,header=TRUE)
    
    # ...and prepare metrics, regressions
    ii <- which(colnames(covars)==trait)
    colnames(covars) <- make.names(colnames(covars))
    covs <- colnames(covars)[-ii]
    ctra <- colnames(covars)[ii]
    if (length(covs) > 0)
        fr <- as.formula(paste(ctra,paste0(covs,collapse="+"),sep="~"))
    else
        fr <- as.formula(paste(ctra,1,sep="~"))
    #message("Null model formula is: ")
    #show(fr)
    nullFit <- glm(fr,data=covars)
    nullModel <- summary(nullFit)
        
    # ...attach score
    pcovars <- cbind(covars,theScore$SCORE)
    colnames(pcovars)[ncol(pcovars)] <- "PRS"
    ii <- which(colnames(pcovars)==trait)
    colnames(pcovars) <- make.names(colnames(pcovars))
    covs <- colnames(pcovars)[-ii]
    ctra <- colnames(pcovars)[ii]
    fm <- as.formula(paste(ctra,paste0(covs,collapse="+"),sep="~"))
    #message("Full model formula is: ")
    #show(fm)
    fullFit <- glm(fm,data=pcovars)
    fullModel <- summary(fullFit)
    
    # - R^2 of the reduced model
    nullR2 <- 1 - nullModel$deviance/nullModel$null.deviance
    
    # - R^2 of the full model
    fullR2 <- 1 - fullModel$deviance/fullModel$null.deviance
    
    # - p-value of the reduced against the full (F test)
    tmp <- tryCatch(anova(nullFit,fullFit,test="F"),error=function(e) {
        return(NULL)
    },finally="")
    if (is.null(tmp))
        P <- 1
    else {
        P <- tmp[["Pr(>F)"]][2]
        # NA p-value? Something bad happens, non-siginificant anyways...
        if (any(is.na(P)))
            P[is.na(P)] <- 1
    }    
    
    # Predictions
    nullPred <- predict(nullFit,pcovars)
    fullPred <- predict(fullFit,pcovars)
    
    # RMSE and MAE
    #nullRmse <- sqrt(sum((nullPred - pcovars[,trait])^2)/nrow(pcovars))
    #fullRmse <- sqrt(sum((fullPred - pcovars[,trait])^2,
    #   na.rm=TRUE)/nrow(pcovars))
    #nullMae <- mean(abs((nullPred - pcovars[,trait])))
    #fullMae <- mean(abs((fullPred - pcovars[,trait])),na.rm=TRUE)
    nullCor <- cor(nullPred,pcovars[,trait])
    if (any(is.na(fullPred)))
        fullCor <- cor(fullPred,pcovars[,trait],use="complete.obs")
    else
        fullCor <- cor(fullPred,pcovars[,trait])
    
    # Now return an object...
    return(list(
        metrics=c(
            null_r2=nullR2,
            full_r2=fullR2,
            prs_r2=fullR2-nullR2,
            prs_pvalue=coef(fullModel)["PRS",4],
            #ftest_pvalue=P,
            #null_rmse=nullRmse,
            #full_rmse=fullRmse,
            #null_mae=nullMae,
            #full_mae=fullMae,
            null_pred_cor=nullCor,
            full_pred_cor=fullCor,
            #null_pred_r2=nullCor^2,
            #full_pred_r2=fullCor^2,
            #prs_pred_r2=fullCor^2-nullCor^2,
            prs_pheno_cor=cor(pcovars[,trait],pcovars$PRS,use="complete.obs"),
            prs_pheno_cor_p=cor.test(pcovars[,trait],pcovars$PRS)$p.value
        ),
        prs=pcovars$PRS#,
        #values=data.frame(
        #    raw_pheno=pcovars[,trait],
        #    null_pheno=nullPred,
        #    full_pheno=fullPred,
        #    prs=pcovars$PRS,
        #    row.names=rownames(pcovars)
        #),
        #reduced_details=coef(nullModel),
        #full_details=coef(fullModel),
        #reduced_ci=confint.default(nullFit),
        #full_ci=confint.default(fullFit),
        #plots=ggs
    ))
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            else 
                m <- FALSE
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

.metricNames <- function() {
    return(c("null_r2","full_r2","prs_r2","prs_pvalue","null_pred_cor",
        "full_pred_cor","prs_pheno_cor","prs_pheno_cor_p"))
}

.whichMaxLast <- function(M) {
  idx <- which(M==max(M),arr.ind=TRUE)
  max_idx <- idx[which.max(idx[,1] * ncol(M) + idx[,2]),]
  return(list(
    ind=c(row=max_idx[1],col=max_idx[2]),
    nam=c(rownames(M)[max_idx[1]],colnames(M)[max_idx[2]])
  ))
}

################################################################################

## ------------------------------------------------------------
## Lambda-by-CI + Gaussian-perturbed betas (summary-statistics)
## ------------------------------------------------------------

# Solve for lambda so that P(|epsilon| < z * SE) = targetProb where 
# epsilon ~ N(0, lambda * SE^2).
# By rescaling, Z = epsilon/SE ~ N(0, lambda), so:
#   P(|Z| < z) = targetProb
# i.e. Φ(z / sqrt(lambda)) - Φ(-z / sqrt(lambda)) = targetProb.
findLambda <- function(targetProb=0.95,zLimit=1.96,lower=1e-4,upper=100) {
    stopifnot(targetProb > 0 && targetProb < 1,zLimit > 0)
    
    f <- function(lambda) {
        pnorm(zLimit/sqrt(lambda)) - pnorm(-zLimit/sqrt(lambda)) - targetProb
    }
    
    return(uniroot(f,lower=lower,upper=upper)$root)
}

# Generate perturbed betas:
# betas_tilde = betas + Normal(0, sqrt(lambda) * SE)
# Optionally draw multiple replicates (columns) for Monte Carlo / bootstraps.
perturbBetas <- function(beta,SE,lambda,N=1,seed=NULL) {
    stopifnot(length(beta)==length(SE),lambda>0,N>=1)
    
    if (!is.null(seed)) # Seed preferably to be set outside this function
        set.seed(seed)
    
    n <- length(beta)
    if (N == 1) {
        return(beta + rnorm(n,mean=0,sd=sqrt(lambda)*SE))
    } 
    else {
        # return an n x N matrix, each column is one perturbation
        eps <- matrix(rnorm(n*N,0,1),nrow=n,ncol=N)
        return(sweep(eps,1,sqrt(lambda)*SE,`*`) + beta)
    }
}

## -------------------------
## Example (drop-in usage)
## -------------------------

# Suppose you have vectors 'beta' and 'SE' from your meta-analysis:
# beta <- gwas$beta
# SE   <- gwas$se

# 1) Get lambda such that ~95% of perturbations fall inside the 95% CI
#lambda_ci95 <- findLambda(targetProb=0.95,zLimit=1.96)
#cat("Lambda (target_prob=0.95, z=1.96):", lambda_ci95, "\n")
# Note: for 95% within 95% CI, lambda will be ~1.0.

# 2) Create one perturbed beta vector
# perturbed <- perturbBetas(beta,SE,lambda=lambda_ci95,seed=123)

# 3) Or create multiple perturbed draws (e.g., 100)
# perturbed_mat <- perturb_betas(beta, SE, lambda = lambda_ci95,
#                                n_draws = 100, seed = 123)
# Each column of 'perturbed_mat' is one perturbed beta set.

## -------------------------
## Optional sanity check
## -------------------------
# # Proportion of a single draw that stays within the original 95% CI:
# ci_low  <- beta - 1.96 * SE
# ci_high <- beta + 1.96 * SE
# in_CI   <- perturbed >= ci_low & perturbed <= ci_high
# mean(in_CI)  # should be close to ~0.95 on average across many runs

#~ # Test
#~ prsFile <- "/media/storage3/playground/b4uprs/work/PRS/baseline/b4u_tgp_sbrc_prs.prs"
#~ prsFile <- "/media/storage3/playground/b4uprs/work/PRS/baseline/b4u_ukb_sbrc_prs.prs"
#~ prsFile <- "/media/storage3/playground/b4uprs/work/PRS/baseline/b4u_tgp_gctb.snpRes.prs"
#~ prsFile <- "/media/storage3/playground/b4uprs/work/PRS/baseline/b4u_ukb_gctb.snpRes.prs"

#~ prsFile = "/media/storage3/playground/b4uprs/work/test/b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto.txt"

#~ covFile <- "/media/storage3/playground/b4uprs/work/HUABB/huabb/HUA_covariates.txt"
#~ trait <- "bmi"
#~ genoBase <- "/media/storage3/playground/b4uprs/work/HUABB/HUA_unrelated_dbsnp"
#~ iidCol <- 2
#~ sum <- TRUE
#~ center <- FALSE
#~ plink <- Sys.which("plink")

#~ # Grid search flow
#~ # 1. Sanitize PRS
#~ #sanFile <- sanitizePrs(prsFile,genoBase)
#~ # 2. Grid search (absolute values of beta)
#~ #obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase)
#~ #pip <- c(0,0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4)
#~ #beta <- c(0,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3)

#~ # PRS coverage is 62.05% (786430 out of 1267425) SNPs
#~ # PRS coverage is 53.81% (3597528 out of 6685102) SNPs
#~ # PRS coverage is 61.48% (563258 out of 916199) SNPs
#~ # PRS coverage is 58.87% (2157515 out of 3664764) SNPs
#~ # GCTB runs not very good... why?
