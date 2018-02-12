library(qtl2)

scan12 <- function(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
           intcovar=NULL, weights=NULL, reml=TRUE,
           model=c("normal", "binary"), cores=1, ...){
    # grab dot args
    dotargs <- list(...)
    if("n_perm" %in% names(dotargs))
      stop("You included n_perm as an argument; you probably want to run scan1perm not scan1.")
    
    if(!is.null(kinship)) { # fit linear mixed model
      return(scan1_pg(genoprobs, pheno, kinship, addcovar, Xcovar, intcovar,
                      reml, cores, ...))
    }
    # deal with the dot args
    tol <- grab_dots(dotargs, "tol", 1e-12)

    stopifnot(tol > 0)
    bintol <- grab_dots(dotargs, "bintol", sqrt(tol)) # for model="binary"
    stopifnot(bintol > 0)
    intcovar_method <- grab_dots(dotargs, "intcovar_method", "lowmem",
                                 c("highmem", "lowmem"))
    quiet <- grab_dots(dotargs, "quiet", TRUE)
    max_batch <- grab_dots(dotargs, "max_batch", NULL)
    maxit <- grab_dots(dotargs, "maxit", 100) # for model="binary"
    check_extra_dots(dotargs, c("tol", "intcovar_method", "quiet", "max_batch", "maxit"))
    
    # check that the objects have rownames
    check4names(pheno, addcovar, Xcovar, intcovar)
    
    # force things to be matrices
    if(!is.matrix(pheno))
      pheno <- as.matrix(pheno)
    if(is.null(colnames(pheno))) # force column names
      colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
    if(!is.null(addcovar) && !is.matrix(addcovar))
      addcovar <- as.matrix(addcovar)
    if(!is.null(Xcovar) && !is.matrix(Xcovar))
      Xcovar <- as.matrix(Xcovar)
    if(!is.null(intcovar) && !is.matrix(intcovar))
      intcovar <- as.matrix(intcovar)
    
    # for binary model
    model <- match.arg(model)
    if(model=="binary") {
      if(!is.null(kinship))
        stop("Can't yet account for kinship with model = \"binary\"")
      pheno <- check_binary_pheno(pheno)
    }
    else {
      # square-root of weights (only if model="normal")
      weights <- qtl2::sqrt_weights(weights) # also check >0 (and if all 1's, turn to NULL)
    }
    
    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *all* phenotypes
    ind2keep <- get_common_ids2(genoprobs, addcovar, Xcovar, intcovar,
                               weights, complete.cases=TRUE)
    ind2keep <- get_common_ids2(ind2keep, rownames(pheno)[rowSums(is.finite(pheno)) > 0])
    if(length(ind2keep)<=2) {
      if(length(ind2keep)==0)
        stop("No individuals in common.")
      else
        stop("Only ", length(ind2keep), " individuals in common: ",
             paste(ind2keep, collapse=":"))
    }
    
    # make sure addcovar is full rank when we add an intercept
    addcovar <- drop_depcols(addcovar, TRUE, tol)
    
    # make sure columns in intcovar are also in addcovar
    addcovar <- force_intcovar(addcovar, intcovar, tol)
    
    # drop things from Xcovar that are already in addcovar
    Xcovar <- drop_xcovar(addcovar, Xcovar, tol)
    
    # batch phenotypes by missing values
    phe_batches <- batch_cols(pheno[ind2keep,,drop=FALSE], max_batch)
    
    # drop cols in genotype probs that are all 0 (just looking at the X chromosome)
    genoprob_Xcol2drop <- genoprobs_col2drop(genoprobs)
    is_x_chr <- attr(genoprobs, "is_x_chr")
    if(is.null(is_x_chr)) is_x_chr <- rep(FALSE, length(genoprobs))
    
    # set up parallel analysis
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores)>1) {
      message(" - Using ", n_cores(cores), " cores")
      quiet <- TRUE # make the rest quiet
    }
    # batches for analysis, to allow parallel analysis
    run_batches <- data.frame(chr=rep(seq_len(length(genoprobs)), length(phe_batches)),
                              phe_batch=rep(seq_along(phe_batches), each=length(genoprobs)))
    run_indexes <- seq_len(length(genoprobs)*length(phe_batches))
    
    # the function that does the work
    by_group_func <- function(i) {
      # deal with batch information, including individuals to drop due to missing phenotypes
      chr <- run_batches$chr[i]
      chrnam <- names(genoprobs)[chr]
      phebatch <- phe_batches[[run_batches$phe_batch[i]]]
      phecol <- phebatch$cols
      omit <- phebatch$omit
      these2keep <- ind2keep # individuals 2 keep for this batch
      if(length(omit) > 0) these2keep <- ind2keep[-omit]
      if(length(these2keep)<=2) return(NULL) # not enough individuals
      
      # subset the genotype probabilities: drop cols with all 0s, plus the first column
      Xcol2drop <- genoprob_Xcol2drop[[chrnam]]
      if(length(Xcol2drop) > 0) {
        pr <- genoprobs[[chr]][these2keep,-Xcol2drop,,drop=FALSE]
        pr <- pr[,-1,,drop=FALSE]
      }
      else
        pr <- genoprobs[[chr]][these2keep,-1,,drop=FALSE]
      
      # subset the rest
      ac <- addcovar; if(!is.null(ac)) { ac <- ac[these2keep,,drop=FALSE]; ac <- drop_depcols(ac, TRUE, tol) }
      Xc <- Xcovar;   if(!is.null(Xc)) Xc <- Xc[these2keep,,drop=FALSE]
      ic <- intcovar; if(!is.null(ic)) { ic <- ic[these2keep,,drop=FALSE]; ic <- drop_depcols(ic, TRUE, tol) }
      ph <- pheno[these2keep,phecol,drop=FALSE]
      wts <- weights[these2keep]
      
      # if X chr, paste X covariates onto additive covariates
      # (only for the null)
      if(is_x_chr[chr]) ac0 <- drop_depcols(cbind(ac, Xc), add_intercept=FALSE, tol)
      else ac0 <- ac
      
      if(model=="normal") {
        # FIX_ME: calculating null RSS multiple times :(
        nullrss <- nullrss_clean(ph, ac0, wts, add_intercept=TRUE, tol)
        
        # scan1 function taking clean data (with no missing values)
        rss <- scan1_clean(pr, ph, ac, ic, wts, add_intercept=TRUE, tol, intcovar_method)
        
        # calculate LOD score
        lod <- nrow(ph)/2 * (log10(nullrss) - log10(rss))
      }
      else { # binary traits
        # FIX_ME: calculating null LOD multiple times :(
        nulllod <- null_binary_clean(ph, ac0, wts, add_intercept=TRUE, maxit, bintol, tol)
        
        # scan1 function taking clean data (with no missing values)
        lod <- scan1_binary_clean(pr, ph, ac, ic, wts, add_intercept=TRUE,
                                  maxit, bintol, tol, intcovar_method)
        
        # calculate LOD score
        lod <- lod - nulllod
      }
      
      list(lod=lod, n=nrow(ph)) # return LOD & number of individuals used
    }
 browser()   
    # number of markers/pseudomarkers by chromosome, and their indexes to result matrix
    npos_by_chr <- dim(genoprobs)[3]
    totpos <- sum(npos_by_chr)
    pos_index <- split(seq_len(totpos), rep(seq_len(length(genoprobs)), npos_by_chr))
    
    # object to contain the LOD scores; also attr to contain sample size
    result <- matrix(nrow=totpos, ncol=ncol(pheno))
    n <- rep(NA, ncol(pheno)); names(n) <- colnames(pheno)
    if(totpos==0) { # edge case of no genoprobs
      colnames(result) <- colnames(pheno)
      attr(result, "sample_size") <- n
      class(result) <- c("scan1", "matrix")
      return(result)
    }
    
    if(n_cores(cores)==1) { # no parallel processing
      for(i in run_indexes) {
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        phecol <- phebatch$cols
        
        this_result <- by_group_func(i)
        if(!is.null(this_result)) {
          result[pos_index[[chr]], phecol] <- t(this_result$lod)
          if(chr==1) n[phecol] <- this_result$n
        }
      }
    }
    else {
      # calculations in parallel
      list_result <- cluster_lapply(cores, run_indexes, by_group_func)
      
      # check for problems (if clusters run out of memory, they'll return NULL)
      result_is_null <- vapply(list_result, is.null, TRUE)
      if(any(result_is_null))
        stop("cluster problem: returned ", sum(result_is_null), " NULLs.")
      
      # reorganize results
      for(i in run_indexes) {
        chr <- run_batches$chr[i]
        chrnam <- names(genoprobs)[chr]
        phebatch <- phe_batches[[run_batches$phe_batch[i]]]
        phecol <- phebatch$cols
        
        if(!is.null(list_result[[i]])) {
          result[pos_index[[chr]], phecol] <- t(list_result[[i]]$lod)
          if(chr==1) n[phecol] <- list_result[[i]]$n
        }
      }
    }
    
    pos_names <- unlist(dimnames(genoprobs)[[3]])
    names(pos_names) <- NULL # this is just annoying
    dimnames(result) <- list(pos_names, colnames(pheno))
    
    # add some attributes with details on analysis
    attr(result, "sample_size") <- n
    
    class(result) <- c("scan1", "matrix")
    result
  }
