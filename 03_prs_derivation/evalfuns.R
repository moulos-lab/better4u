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
# ukb: sanitization against UKB
# permaBim: if bgen, store the converted BIM permanently instead of temp file
sanitizePrs <- function(prsFile,genoBase,perChr=FALSE,chrs=seq(22),chrSep="_",
    bgen=FALSE,ukb=FALSE,permaBim=FALSE,from=c("sbayesrc","prscs","ready"),
    pip=0.001,plink2=Sys.which("plink2"),rc=NULL) {
    if (bgen && !(is.character(plink2) || file.exists(plink2))) 
        stop("PLINK 2.0 is required by sanitizePrs() if you have bgen files!")
        
    from <- from[1]
    
    # Basic check, make sure that PLINK or BGEN files exist
    if (bgen) {
        if (perChr) {
            if (!file.exists(paste0(genoBase,chrSep,"chr1.bgen")))
                stop("BGEN files not found where expected!")
        }
        else {
            if (!file.exists(paste0(genoBase,".bgen")))
                stop("BGEN files not found where expected! genoBase is ",
                    genoBase)
        }
    }
    else {
        if (perChr) {
            if (!file.exists(paste0(genoBase,chrSep,"chr1.bim")))
                stop("PLINK files not found where expected!")
        }
        else {
            if (!file.exists(paste0(genoBase,".bim")))
                stop("PLINK files not found where expected! genoBase is ",
                    genoBase)
        }
    }
    
    # Read SNP data (bim) and initial PRS
    if (bgen) {
        if (perChr) {
            bims <- cmclapply(chrs,function(chr) {
                message("Converting to BIM with PLINK 2.0 for ",chr)
                bFile <- paste0(genoBase,chrSep,"chr",chr,".bgen")
                o <- ifelse(permaBim,paste0(genoBase,chrSep,"chr",chr),
                    tempfile())
                # If UKB, ref-first and a .sample file are required
                bgenArgs <- "ref-unknown"
                if (ukb) {
                    sFile <- paste0(genoBase,chrSep,"chr",chr,".sample")
                    bgenArgs <- paste0("ref-first --sample ",sFile)
                }
                # ref-unknown, sanitization will take care if problem
                args <- c("--bgen",bFile,bgenArgs,"--make-just-bim",
                    "--out",o,"--silent") 
                out <- tryCatch({
                    suppressWarnings(system2(plink2,args=args))
                    TRUE
                },error=function(e) {
                    message("Caught error: ",e$message)
                    return(FALSE)
                },finally="")
                
                if (!out)
                    stop("Failed to generate BIM file for ",chr,
                        "! Exiting...")
                
                message("Reading converted BIM for ",chr)
                bim <- read.delim(paste0(o,".bim"),header=FALSE)
                return(bim)
            },rc=rc)
            bim <- do.call("rbind",bims)
        }
        else {
            message("Converting to BIM with PLINK 2.0")
            bFile <- paste0(genoBase,".bgen")
            o <- ifelse(permaBim,genoBase,tempfile())
            bgenArgs <- "ref-unknown"
            if (ukb) {
                sFile <- paste0(genoBase,".sample")
                bgenArgs <- paste0("ref-first --sample ",sFile)
            }
            args <- c("--bgen",bFile,bgenArgs,"--make-just-bim --out",o,
                "--silent")
            out <- tryCatch({
                suppressWarnings(system2(plink2,args=args))
                TRUE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(FALSE)
            },finally="")

            if (!out)
                stop("Failed to generate BIM file! Exiting...")
            
            # If all ok, read the bim file
            bim <- read.delim(paste0(o,".bim"),header=FALSE)
        }
    }
    else {
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
# calculations. plink2 now required for bgen input.
# As there has been mixups reported regarding plink/plink2 usage, we try to add
# explicit support for both.
# ukb=TRUE controls ref-first or ref-unknown for PLINK2.
# If you wish to keep the file with non-mathcing samples to be removed, give a
# filename to remFile
evalPrs <- function(prsFile,covFile,trait,genoBase,perChr=FALSE,chrs=seq(22),
    chrSep="_",iidCol=2,sum=TRUE,center=FALSE,bgen=FALSE,ukb=FALSE,remFile=NULL,
    plink=Sys.which("plink"),plink2=Sys.which("plink2"),rc=NULL) {
    # Determine plink version(s) to use and for what purpose
    proceed <- .informAboutPlinkVersions(plink,plink2,bgen)
    if (!proceed) return()
    
    if (!is.null(remFile) && !is.character(remFile)) {
        warning("remFile is not a proper file name! A temp file will be used.",
            immediate.=TRUE)
        remFile <- NULL
    }
    if (is.null(remFile))
        remFile <- tempfile()
    
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
    if (bgen) {
        if (perChr)
            fam <- read.table(paste0(genoBase,chrSep,"chr",chrs[1],".sample"),
                header=TRUE)
        else
            fam <- read.table(paste0(genoBase,".sample"),header=TRUE)
        fam <- fam[-1,,drop=FALSE]
    }
    else
        fam <- read.table(paste0(genoBase,".fam"))
    rownames(fam) <- fam[,2]
    common <- intersect(rownames(fam),rownames(covars))
    # length(common) can be only smaller or equal to nrow(covars). If smaller, 
    # align with common which is based on fam and possibly re-align later with 
    # remaining fam. If all fam included in common, again, no problem this way
    covars <- covars[common,,drop=FALSE]
    # Continue with fam check
    if (length(common) < nrow(fam)) { # Can be only smaller or equal
        # File for removal of samples in plink --score
        remove <- setdiff(rownames(fam),common)
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
    # UKB case
    if ("eid" %in% colnames(covars) || "EID" %in% colnames(covars))
        covars <- covars[,-which(colnames(covars) %in% c("EID","eid"))]
    
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
        prsSplit <- .splitPrsPerChr(prsFile,rc=rc)
        
        # Now calculate scores per chromosome
        if (bgen) {
            # For some internal PLINK2 reason, parallelization does not work
            # at all as expected...
            #scoreFiles <- unlist(cmclapply(chrs,function(chr) {
            scoreFiles <- unlist(lapply(chrs,function(chr) {
                message("\nCalculating score with PLINK 2.0 --score for ",chr)
                bFile <- paste0(genoBase,chrSep,"chr",chr,".bgen")
                pFile <- prsSplit[chr]
                o <- tempfile()
                o <- paste0(o,"_prs_",chr)
                # If UKB, ref-first and a .sample file are required
                bgenArgs <- "ref-unknown"
                if (ukb) {
                    sFile <- paste0(genoBase,chrSep,"chr",chr,".sample")
                    bgenArgs <- paste0("ref-first --sample ",sFile)
                }
                args <- c("--bgen",bFile,bgenArgs,"--score",
                    pFile,"1 2 3 header ignore-dup-ids",
                    ifelse(sum,"cols=fid,pheno1,nallele,denom,scoresums",
                        "cols=fid,pheno1,nallele,denom,scoreavgs"),
                ifelse(center,"center",""),"--out",o)
                if (!is.null(remFile))
                    args <- c(args,"--remove",remFile)
                #args <- c(args,"--silent") # Not needed for linear execution
                humanCommand <- .formatPlinkCommand(plink2,args)
                message("Executing:\n",humanCommand)
                out <- tryCatch({
                    suppressWarnings(system2(plink2,args=args))
                    TRUE
                },error=function(e) {
                    message("Caught error: ",e$message)
                    return(FALSE)
                },finally="")
                return(paste0(o,".sscore"))
            }))
        }
        else {
            scoreFiles <- unlist(cmclapply(chrs,function(chr) {
                message("\nCalculating score with PLINK 1.9 --score for ",chr)
                bFile <- paste0(genoBase,chrSep,"chr",chr)
                pFile <- prsSplit[chr]
                o <- tempfile()
                o <- paste0(o,"_prs_",chr)
                args <- c("--bfile",bFile,"--score",pFile,"1 2 3 header",
                    ifelse(sum,"sum",""),ifelse(center,"center",""),"--out",o)
                if (!is.null(remFile))
                    args <- c(args,"--remove",remFile)
                args <- c(args,"--silent")
                humanCommand <- .formatPlinkCommand(plink,args)
                message("Executing:\n",humanCommand)
                out <- tryCatch({
                    suppressWarnings(system2(plink,args=args))
                    TRUE
                },error=function(e) {
                    message("Caught error: ",e$message)
                    return(FALSE)
                },finally="")
                return(paste0(o,".profile"))
            },rc=rc))
        }
        
        # Then somehow read and combine... essentially for the individuals
        # (row) add #chrs columns and create a data frame with one col SCORE
        scoreCol <- ifelse(bgen,ifelse(sum,"SCORE1_SUM","SCORE1_AVG"),
            ifelse(sum,"SCORESUM","SCORE"))
        cc <- ifelse(bgen,"","#")
        chrScores <- cmclapply(scoreFiles,function(f) {
            tmpScore <- read.table(f,row.names=2,header=TRUE,comment.char=cc)
            return(tmpScore[,scoreCol,drop=FALSE])
        })
        chrScores <- do.call("cbind",chrScores)
        theScore <- data.frame(SCORE=rowSums(chrScores))
    }
    else {
        if (bgen) {
            message("\nCalculating score with PLINK 2.0 --score")
            bgenFile <- paste0(genoBase,".bgen")
            bgenArgs <- "ref-unknown"
            if (ukb) {
                sFile <- paste0(genoBase,".sample")
                bgenArgs <- paste0("ref-first --sample ",sFile)
            }
            args <- c("--bgen",bgenFile,bgenArgs,"--score",
                prsFile,"1 2 3 header ignore-dup-ids",
                ifelse(sum,"cols=fid,pheno1,nallele,denom,scoresums",
                    "cols=fid,pheno1,nallele,denom,scoreavgs"),
                ifelse(center,"center",""),"--out",prsName)
            if (!is.null(remFile))
                args <- c(args,"--remove",remFile)
            #args <- c(args,"--silent")
            humanCommand <- .formatPlinkCommand(plink2,args)
            message("Executing:\n",humanCommand)
            out <- tryCatch({
                suppressWarnings(system2(plink2,args=args))
                TRUE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(FALSE)
            },finally="")

            if (!out)
                stop("Failed to generate score file! Exiting...")
            
            # If all ok, read the score file
            scoreFile <- paste0(prsName,".sscore")
            theScore <- read.table(scoreFile,row.names=2,header=TRUE,
                comment.char="")
        }
        else {
            message("\nCalculating score with PLINK 1.9 --score")
            args <- c("--bfile",genoBase,"--score",prsFile,"1 2 3 header",
                ifelse(sum,"sum",""),ifelse(center,"center",""),"--out",prsName)
            if (!is.null(remFile))
                args <- c(args,"--remove",remFile)
            #args <- c(args,"--silent")
            humanCommand <- .formatPlinkCommand(plink,args)
            message("Executing:\n",humanCommand)
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
    hasPC <- grepl("^(PC|genetic_principal)",names(covs),perl=TRUE)
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
            snps_covered=nsnps
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

.splitPrsPerChr <- function(prsFile,outBase=NULL,rc=NULL) {
    message("Reading ",prsFile)
    tmpPrs <- read.delim(prsFile)
    tmpSplit <- split(tmpPrs,tmpPrs$CHR)
    message("Splitting ",prsFile)
    prsSplit <- unlist(cmclapply(names(tmpSplit),function(chr) {
        o <- ifelse(is.null(outBase),tempfile(),outBase)
        o <- paste0(o,"_chr",chr,".prs")
        write.table(tmpSplit[[chr]][,1:3],file=o,sep="\t",row.names=FALSE,
            quote=FALSE)
        return(o)
    },rc=rc))
    names(prsSplit) <- names(tmpSplit)
    return(prsSplit)
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

.informAboutPlinkVersions <- function(plink,plink2,bgen) {
    if (is.character(plink) && file.exists(plink) 
        && is.character(plink2) && file.exists(plink2)) {
        message("Both plink 1.9 and plink 2.0 are found on the system.")
        if (bgen) {
            message("BGEN genotype file input requested.")
            message("plink 2.0 will be used for scoring.")
        }
        else
            message("plink 1.9 will be used for scoring.")
        return(TRUE)
    }
    else if (is.character(plink) && file.exists(plink) 
        && (!is.character(plink2) || !file.exists(plink2))) {
        message("Only plink 1.9 was found on the system.")
        if (bgen) {
            message("BGEN genotype file input requested.")
            message("Please do not set the plink argument to point to plink2 ",
                "executable, the flow will be compromised.")
            message("plink 2.0 is required for scoring. Exiting...")
            return(FALSE)
        }
        else
            message("plink 1.9 will be used for scoring.")
        return(TRUE)
    }
    else if ((!is.character(plink) || !file.exists(plink)) 
        && is.character(plink2) && file.exists(plink2)) {
        message("Only plink 2.0 was found on the system.")
        if (bgen) {
            message("BGEN genotype file input requested.")
            message("plink 2.0 will be used for scoring.")
        }
        else
            message("plink 2.0 will be used for scoring.")
        return(TRUE)
    }
    else {
        message("FATAL! Neither plink 1.9 nor plink 2.0 were found on the ",
            "system! Exiting...")
        return(FALSE)
    }
}

.formatPlinkCommand <- function(plink,args) {
    init <- paste0(plink," ",paste(args,collapse=" "))
    tmp <- strsplit(init," --")[[1]]
    tmp <- trimws(tmp)
    fin <- c(tmp[1],paste("  --",tmp[2:length(tmp)],sep=""))
    fin[1:(length(fin)-1)] <- paste(fin[1:(length(fin)-1)]," \\\n",sep="")
    return(fin)
}

################################################################################

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
