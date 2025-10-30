library(parallel)
library(ggplot2)
library(reshape2)
library(R.utils)

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
            return(evalPrs(griFile,covFile,trait,genoBase,iidCol=iidCol,sum=sum,
                center=center,plink=plink))
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
    message("  abs(BETA) > ",beta[bestIndex$ind[1]]," (i=",bestIndex$ind[1],")")
    message("  PIP > ",pip[bestIndex$ind[2]]," (j=",bestIndex$ind[2],")")
    
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

formatPrs <- function(prsFile,outFile,from=c("sbayesrc","prscs"),pip=0.001) {
    from <- from[1]
    
    message("Reading PRS ",prsFile)
    prs <- read.delim(prsFile,header=from=="sbayesrc")
    
    # Assign header names if prscs
    if (from=="prscs")
        colnames(prs) <- c("CHR","SNP","BP","A1","A2","BETA")

    # Output final PRS file ready-made for further downstream analysis
    if (from=="prscs") {
        prs$SE <- 0
        prs$PIP <- pip
    }
    
    # Order output by chromosome
    message("Ordering output by chromosome")
    prs <- prs[order(prs$CHR,prs$SNP),]
    
    # Export
    write.table(prs,file=outFile,quote=FALSE,sep="\t",row.names=FALSE)
    
    message("Formatted PRS is in ",outFile)
    return(outFile)
}

# from="ready" to only perform sanitization and exclude SE and PIP from output
# if perChr, BIM files are assumed one per chromosome for chrs variable
# We accept only ending in *{chrSep}{chr}.bim, so if genoBase="COHORT" and
# chrSep="_" then the bim file is COHORT_chr1.bim
sanitizePrs <- function(prsFile,genoBase,perChr=FALSE,chrs=seq(22),chrSep="_",
    from=c("sbayesrc","prscs","ready"),pip=0.001,rc=NULL) {
    from <- from[1]
    
    # Basic check, make sure that PLINK files exist
    if (perChr) {
        if (!file.exists(paste0(genoBase,chrSep,"chr1.bim")))
            stop("PLINK files not found where expected!")
    }
    else {
        if (!file.exists(paste0(genoBase,".bim")))
            stop("PLINK files not found where expected! genoBase is ",genoBase)
    }
    
    # Read SNP data (bim) and initial PRS
    if (perChr) {
        bims <- cmclapply(chrs,function(chr) {
            bimFile <- paste0(genoBase,chrSep,"chr",chr,".bim")
            message("Reading BIM ",bimFile)
            bim <- read.delim(bimFile,header=FALSE)
        },rc=rc)
        bim <- do.call("rbind",bims)
    }
    else {
        bimFile <- paste0(genoBase,".bim")
        message("Reading BIM ",bimFile)
        bim <- read.delim(bimFile,header=FALSE)
    }
    # Rename BIM columns for clarity
    colnames(bim) <- c("CHR","SNP","CM","BP","A1_bim","A2_bim")
    
    message("Reading PRS ",prsFile)
    prs <- read.delim(prsFile,header=from %in% c("sbayesrc","ready"))
    
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
    
    #outCov <- paste0(prsFile,".coverage")
    #message("Writing coverage to ",outCov)
    #writeLines(as.character(round(100*coverage,2)),outCov)
    
    # Output final PRS file (with updated A1 and BETA) - we fake the SE, PIP
    # column in the case of PRS-CS
    if (from=="prscs") {
        merged$SE <- 0
        merged$PIP <- pip
    }
    if (from == "sbayesrc" && !("CHR" %in% names(merged))) # Until I fix in fork
        merged$CHR <- 0L
    if (from %in% c("sbayesrc","prscs"))
        output <- merged[,c("SNP","A1","BETA","CHR","SE","PIP")]
    else
        output <- merged[,c("SNP","A1","BETA","CHR")]
    # Order output by chromosome
    message("Ordering output by chromosome")
    output <- output[order(output$CHR,output$SNP),]
    
    # Export per chromosome if required
    #if (perChr) {
    #    cmclapply(chrs,function(chr) {
    #        outFile <- paste0(prsFile,chrSep,chr,".san")
    #        message("Reading BIM ",bimFile)
    #        bim <- read.delim(bimFile,header=FALSE)
    #    },rc=rc)
    #}
    #else {
    #    outFile <- paste0(prsFile,".san")
    #    write.table(output,file=outFile,quote=FALSE,sep="\t",row.names=FALSE)
    #}
    
    # Export
    outFile <- paste0(prsFile,".san")
    write.table(output,file=outFile,quote=FALSE,sep="\t",row.names=FALSE)
    
    message("Sanitized PRS is in ",outFile)
    return(outFile)
}

# When evaluating by chr, then PRS must be calculated per chr and the individual
# files must be concatenated and then added to the covariates. The same things
# as sanitizePrs apply if calculations are done per chromosome. If the latter,
# PRS file is firstly split and then calculations are made.
# NOTE well that there are tiny differences (range -10e-6, 10e-6 for HUA) 
# between the two approaches which are attributed to floating point rounds and
# calculations.
evalPrs <- function(prsFile,covFile,trait,genoBase,perChr=FALSE,chrs=seq(22),
    chrSep="_",iidCol=2,sum=TRUE,center=FALSE,plink=Sys.which("plink"),
    rc=NULL) {
    # Base name for plink score output
    prsName <- sub("\\.[^.]*$","",prsFile)
    
    # Number of SNPs in PRS
    nsnps <- as.numeric(countLines(prsFile)) - 1
    
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
    if (perChr) {
        # Firstly split prsFile into files per chromosome. It should be 
        # sanitized so it has 4 columns, the 4th is chromosome.
        tmpPrs <- read.delim(prsFile)
        tmpSplit <- split(tmpPrs,tmpPrs$CHR)
        prsSplit <- unlist(cmclapply(names(tmpSplit),function(chr) {
            o <- tempfile()
            o <- paste0(o,"_",chr)
            write.table(tmpSplit[[chr]][,1:3],file=o,sep="\t",row.names=FALSE,
                quote=FALSE)
            return(o)
        },rc=rc))
        names(prsSplit) <- names(tmpSplit)
        
        # Now calculate scores per chromosome
        scoreFiles <- unlist(cmclapply(chrs,function(chr) {
            message("Calculating score with PLINK --score for ",chr)
            bFile <- paste0(genoBase,chrSep,"chr",chr)
            pFile <- prsSplit[chr]
            o <- tempfile()
            o <- paste0(o,"_prs_",chr)
            args <- c("--bfile",bFile,"--score",pFile,"1 2 3 header",
                ifelse(sum,"sum",""),ifelse(center,"center",""),"--out",o)
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
            return(paste0(o,".profile"))
        },rc=rc))
        
        # Then somehow read and combine... essentially for the individuals
        # (row) add #chrs columns and create a data frame with one col SCORE
        chrScores <- cmclapply(scoreFiles,function(f) {
            tmpScore <- read.table(f,row.names=2,header=TRUE)
            return(tmpScore[,"SCORESUM",drop=FALSE])
        })
        chrScores <- do.call("cbind",chrScores)
        theScore <- data.frame(SCORE=rowSums(chrScores))
    }
    else {
        #message("Calculating score with PLINK --score")
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
    }
        
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
    
    # Direct R2 as suggested by Anders
    y <- pcovars[,ii]
    prs <- pcovars[,ncol(pcovars)]
    covs <- pcovars[,-c(ii,ncol(pcovars))]
    g_covs <- nong_covs <- NULL
    hasPC <- grepl("^PC",names(covs),perl=TRUE)
    if (any(hasPC)) {
        g_covs <- covs[,hasPC,drop=FALSE]
        nong_covs <- covs[,!hasPC,drop=FALSE]
    }
    else
        nong_covs <- covs
    craw <- cor(pcovars[,trait],pcovars$PRS,use="complete.obs")
    dr2_raw <- craw^2
    dr2_residy <- directR2(y,prs,nong_covs=nong_covs,g_covs=g_covs)
    dr2_residb <- directR2(y,prs,nong_covs=nong_covs,g_covs=g_covs,
        resid_both=TRUE)
    
    # Now return an object...
    return(list(
        metrics=c(
            null_r2=nullR2,
            full_r2=fullR2,
            prs_r2=fullR2-nullR2,
            prs_pvalue=coef(fullModel)["PRS",4],
            null_pred_cor=nullCor,
            full_pred_cor=fullCor,
            prs_pheno_cor=craw,
            prs_pheno_cor_p=cor.test(pcovars[,trait],pcovars$PRS)$p.value,
            prs_pheno_r2_raw=craw^2,
            prs_pheno_r2_resy=dr2_residy,
            prs_pheno_r2_resb=dr2_residb,
            snps_covered=nsnps#,
            #penal_prs_r2=-log10((fullR2-nullR2)/nsnps),
            #penal_prs_pheno_cor=-log10(craw/nsnps),
            #penal_prs_pheno_r2_resy=-log10(dr2_residy/nsnps),
            #penal_prs_pheno_r2_resb=-log10(dr2_residb/nsnps)
        ),
        prs=pcovars$PRS
    ))
}

# nong_covs: non-genetic covariates, e.g. age, sex etc NO PCs
# g_covs: PCs or other genetic covariates
directR2 <- function(y,prs,nong_covs=NULL,g_covs=NULL,resid_both=FALSE) {
    stopifnot(length(y) == length(prs))
  
    if (is.null(nong_covs) && is.null(g_covs))
        return(cor(y,prs,use="complete.obs")^2)
  
    if (!is.null(nong_covs))
        nong_covs <- as.data.frame(nong_covs)
    if (!is.null(g_covs))
        g_covs <- as.data.frame(g_covs)
        
    if (!is.null(nong_covs) &&  !is.null(g_covs))
        covs <- cbind(nong_covs,g_covs)
    else if (!is.null(nong_covs) && is.null(g_covs))
        covs <- nong_covs
    else if (is.null(nong_covs) && !is.null(g_covs))
        covs <- g_covs
    
    # Residualize y with glm
    fit_y <- glm(y ~ .,data=covs,family="gaussian")
    resid_y <- resid(fit_y,type="response")
  
    if (!resid_both)
        return(cor(resid_y,prs[!is.na(resid_y)],use="complete.obs")^2)
    else {
        fit_prs <- glm(prs~.,data=g_covs,family="gaussian")
        resid_prs <- resid(fit_prs,type="response")
        return(cor(resid_y,resid_prs,use="complete.obs")^2)
    }
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
        "full_pred_cor","prs_pheno_cor","prs_pheno_cor_p","prs_pheno_r2_raw",
        "prs_pheno_r2_resy","prs_pheno_r2_resb","nsnps"))#,"penal_prs_r2",
        #"penal_prs_pheno_cor","penal_prs_pheno_r2_resy",
        #"penal_prs_pheno_r2_resb"))
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

#~ # Direct R2 (GLM, Option B: residual space only)
#~ direct_r2_glm_resid <- function(y, prs, covariates = NULL) {
#~  stopifnot(is.numeric(y), is.numeric(prs))
#~  if (!is.null(covariates)) {
#~      covariates <- as.data.frame(covariates)
#~      if (nrow(covariates) != length(y))
#~          stop("covariates must have same number of rows as y/prs")
#~  }
  
#~  # complete cases
#~  if (is.null(covariates))
#~      ok <- complete.cases(y, prs)
#~  else {
#~      ok <- complete.cases(y, prs, covariates)
#~      covariates <- covariates[ok, , drop = FALSE]
#~  }
#~  y <- y[ok]; prs <- prs[ok]
#~  n <- length(y)
#~  if (n < 3)
#~      stop("Too few complete observations")
  
#~  if (var(y) <= .Machine$double.eps) {
#~      warning("Zero or near-zero variance in y; returning NA")
#~      return(NA_real_)
#~  }
  
#~  if (is.null(covariates)) {
#~      resid_y <- y
#~      resid_prs <- prs
#~  } 
#~  else {
#~      # regress out covariates
#~      fit_y <- glm(y ~ ., data = covariates, family = gaussian())
#~      resid_y <- resid(fit_y)

#~      fit_prs <- glm(prs ~ ., data = covariates, family = gaussian())
#~      resid_prs <- resid(fit_prs)
#~  }
  
#~  # fit resid_y ~ resid_prs
#~  fit_resid <- glm(resid_y ~ resid_prs, family = gaussian())
#~  y_hat_resid <- predict(fit_resid, type = "response")
  
#~  mse <- mean((resid_y - y_hat_resid)^2)
#~  r2  <- 1 - mse / var(y)
    
#~  return(r2)
#~ }

#dr2 <- directR2(pcovars[,ii],pcovars[,ncol(pcovars)],
#    pcovars[,-c(ii,ncol(pcovars))])
#~ directR2_recal <- function(y,prs,covariates=NULL) {
#~   stopifnot(length(y) == length(prs))
  
#~   if (is.null(covariates)) {
#~     # Model: phenotype ~ fixed PRS effect (no refitting weights externally)
#~     fit <- glm(y ~ prs,family="gaussian")
#~     y_hat <- predict(fit)
#~   } else {
#~     # Step 1: regress out covariates from phenotype
#~     fit_cov <- glm(y ~ ., data=data.frame(covariates),family="gaussian")
#~     #resid_y <- resid(fit_cov)
#~     fitted_cov <- predict(fit_cov, type = "response")
    
#~     # Step 2: regress out covariates from PRS
#~     fit_prs <- glm(prs ~ ., data=data.frame(covariates),family="gaussian")
#~     #resid_prs <- resid(fit_prs)
#~     fitted_prs_cov <- predict(fit_prs, type = "response")
    
#~     ## Step 3: predicted values from residualized PRS
#~     #beta_prs <- coef(glm(resid_y ~ resid_prs,family="gaussian"))["resid_prs"]
#~     #y_hat <- beta_prs * resid_prs
    
#~     # residuals
#~     resid_y <- y - fitted_cov
#~     resid_prs <- prs - fitted_prs_cov

#~     denom <- sum(resid_prs^2)
#~     if (denom <= .Machine$double.eps) {
#~       warning("Residualized PRS has (near) zero variance; cannot compute slope")
#~       return(NA_real_)
#~     }

#~     # slope from resid_y ~ resid_prs (OLS)
#~     beta_prs <- sum(resid_y * resid_prs) / denom

#~     # predicted Y on original scale = fitted_cov + beta * resid_prs
#~     y_hat <- fitted_cov + beta_prs * resid_prs
#~   }
  
#~   mse <- mean((y - y_hat)^2)
#~   var_y <- var(y)
  
#~   r2 <- 1 - mse / var_y
#~   return(r2)
#~ }

