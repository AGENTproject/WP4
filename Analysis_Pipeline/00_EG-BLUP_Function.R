#' pheno.data: 'Genotype_Matrix', 'Trait'
#' geno.data:  'Marker', 'CHR', 'LOC', 'REF', 'ALT', Genotypes..
#' groups:     'AGENT_ID', 'group'
#' sig_markers: 
#' climate:

agend_egblup <- function(pheno.egblup, geno.egblup, MAF = 0, pred_mode = FALSE, 
                         groups = NULL, sig_markers = NULL, climate = NULL,
                         num.crossvalid = 2, num.folds = 5, nIter = 1000, burnIn = 300) {
  
  if (!is.null(sig_markers) & length(sig_markers) == 0) { sig_markers <- NULL }
  
  geno.subset  <- list()
  pheno.subset <- list()
  
  split_prediction <- list()
  split_PredEB     <- list()
  
  if (!is.null(groups)) {
    total_splits <- max(as.numeric(groups$group))
    
    for(i in 1:total_splits) {
      # TODO: check if the groups has a critical mass of entries
      geno.subset[[i]] <- geno.egblup[, c("Marker", "CHR", "LOC", "REF", "ALT", groups[groups$group == i, 'AGENT_ID'])] 
      pheno.subset[[i]] <- pheno.egblup[pheno.egblup$Genotype_Matrix %in% groups[groups$group == i, 'AGENT_ID'],]
    }
  } else {
    total_splits      <- 1
    geno.subset[[1]]  <- geno.egblup
    pheno.subset[[1]] <- pheno.egblup
  }
  
  for (split_num in 1:total_splits) {
    
    ##############################################################################
    #
    #         Sample R code for EG-BLUP
    #
    #         Code is a derivative of the publication: 
    #         Choosing the right tool: Leveraging of plant genetic resources in wheat (Triticum aestivum L.) 
    #         benefits from selection of a suitable genomic prediction model
    #----------------------------------------------------------------------------------------------------------
    #         Marcel O. Berkner1, Albert W. Schulthess1, Yusheng Zhao1, Yong Jiang1, Markus Oppermann1 and Jochen C. Reif1
    #----------------------------------------------------------------------------------------------------------
    #         1 Leibniz Institute of Plant Genetics and Crop Plant Research (IPK), 06466 Stadt Seeland, Germany
    #         
    #         Corresponding author: Jochen C. Reif, Email: reif@ipk-gatersleben.de
    #----------------------------------------------------------------------------------------------------------
    #         This code was developed with R version 4.0.2
    #
    ##############################################################################
    
    #' Load genotypic data and pure based on minor-allele-frequency
    Geno <- geno.subset[[split_num]]
    
    #' Omit additional information on SNP markers
    Freq <- data.frame(c(1:nrow(Geno)), Geno[,'Marker'], (Geno[,c(6:ncol(Geno))]))
    
    #' Transform SNP matrix into numeric
    Freq[, c(3:ncol(Freq))] <- sapply(Freq[, c(3:ncol(Freq))], as.numeric)
    
    #' Calculate frequencies of both alleles
    Freq <- data.frame(Freq[, c(1,2)],
                       (rowMeans(Freq[, c(3:ncol(Freq))]))/2,
                       1 - ((rowMeans(Freq[, c(3:ncol(Freq))]))/2))
    
    colnames(Freq) <- c('RowNo', 'Marker', 'freq.allele1', 'freq.allele2')
    f <- data.frame(Freq[,c(1,2)],c(1:nrow(Freq)))
    
    for(i in 1:nrow(Freq)) {
      #' Select the minor allele
      if (Freq[i,3] <= 0.5) {
        f[i,3] <- Freq[i,3]
      } else {
        f[i,3] <- Freq[i,4]
      }
    }
    
    colnames(f) <- c('RowNo', 'Marker', 'MAF')
    f <- f[f[,3] >= MAF,]
    
    #' Selected SNP according to minor allele frequency threshold
    #' significant markers should exists regardless MAF threshold
    selmarker <- unique(c(sig_markers, as.character(f[,2])))
    
    #' Create G matrix, marker annotation of -1,0,1
    G <- data.frame(Geno[,c('Marker', 'CHR', 'LOC')], Geno[, c(6:ncol(Geno))]-1)
    
    #' Reduction of matrix to selected markers
    G <- G[c(which(G$Marker %in% selmarker)),]
    
    #' Load phenotypic data
    #' Load phenotype  matrix; phenotype information is given as BLUE
    Pheno <- pheno.subset[[split_num]]
    accession.order <- colnames(G)
    accession.order <- accession.order[4:length(colnames(G))]
    
    #' Bring phenotypic and genotypic information in same order
    Pheno <- Pheno[(match(accession.order, Pheno$Genotype_Matrix)),] 
    
    #' Select trait if multiple traits exist
    dta.P <- Pheno
    
    #' exclude entries with missing value in the phenotypic data BLUEs
    if(pred_mode == FALSE) { dta.P <- dta.P[!(is.na(dta.P$Trait)),] }
    
    p <- as.matrix(dta.P[,'Trait'])
    
    #' EG-BLUP
    #' Addititive Relationship Matrix
    #' First method based on VanRaden, 2008
    
    #' Create marker matrix M with dimension n x m for n genotypes and m markers
    M <- Geno[c(which(Geno$Marker %in% selmarker)), 6:ncol(Geno)]
    rownames(M) <- Geno[c(which(Geno$Marker %in% selmarker)), 'Marker']
    
    M <- M[, c(dta.P$Genotype_Matrix)]
    M <- M - 1
    M <- t(M)
    
    if (!is.null(sig_markers)) {
      Zsig <- M[, sel_markers]
      
      #' remove monomorphic markers 
      M <- M[, !apply(M, 2, sd) == 0]
      
      Mnosig <- M[, -c(which(colnames(M) %in% sel_markers))]
      freq.p <- (colSums(Mnosig == 1) + 1/2 * colSums(Mnosig == 0))/nrow(M)
      
      P <- 2*(t(matrix(rep(freq.p,nrow(Mnosig)), ncol = nrow(Mnosig))) - 0.5)
      Z <- Mnosig - P      
    } else {
      #' remove monomorphic markers 
      M <- M[, !apply(M, 2, sd) == 0]
      
      freq.p <- (colSums(M == 1) + 1/2 * colSums(M == 0))/nrow(M)
      
      P <- 2*(t(matrix(rep(freq.p,nrow(M)), ncol = nrow(M))) - 0.5)
      Z <- M - P
    }
    
    ZZ <- Z %*% t(Z)
    vp <- 2*sum(freq.p*(1 - freq.p))
    A <- ZZ/vp
    
    rm(Z, ZZ, P, vp, freq.p)
    
    #' Covariance Matrix H
    #' Hadamard Product of A as described by Henderson 1985
    H <- A * A
    rownames(H) <- dta.P$Genotype_Matrix
    colnames(H) <- dta.P$Genotype_Matrix
    
    if (is.null(climate)) {
      if (is.null(sig_markers)) {
        ETA.lst <- list(list(K1 = A, model = 'RKHS'), 
                        list(K2 = H, model = 'RKHS'))
      } else {
        ETA.lst <- list(list(K1 = A, model = 'RKHS'),
                        list(K2 = H, model = 'RKHS'),
                        list(X  = Zsig, model = 'FIXED'))  
      }
    } else {
      #' get the climate matrix that restricted to the selected AGENT ids
      E <- climate[rownames(A), colnames(A)]
      
      if (is.null(sig_markers)) {
        ETA.lst <- list(list(K1 = A, model = 'RKHS'), 
                        list(K2 = H, model = 'RKHS'),
                        list(K3 = E, model = 'RKHS'))
      } else {
        ETA.lst <- list(list(K1 = A, model = 'RKHS'),
                        list(K2 = H, model = 'RKHS'),
                        list(K3 = E, model = 'RKHS'),
                        list(X  = Zsig, model = 'FIXED'))  
      }
    }
    
    if(pred_mode == TRUE) {
      EGBLUP <- BGLR::BGLR(y = p,
                           ETA = ETA.lst,
                           nIter = nIter,  # set >= 10000 for large data
                           burnIn = burnIn,  # set >= 3000 for large data
                           verbose = FALSE)  
      
      result <- list(list(data.frame(dta.P$Genotype_Matrix,EGBLUP[1],EGBLUP[14],EGBLUP[15])),
                     list(data.frame(EGBLUP[6],EGBLUP[7],EGBLUP[16],EGBLUP[17],EGBLUP[18],EGBLUP[19])))
      
      split_prediction[[split_num]]  <- result[[1]][[1]]
    } else {
      #' Index-Matrix for runs of five-fold-cross-validation 
      #' Randomly assigned index-matrix
      
      #' select number of five-fold cross-validations
      # num.crossvalid <- 2
      # num.folds      <- 5
      
      #' random split of each round of cross-validation into five parts
      idx.M <- data.frame(as.numeric(row.names.data.frame(dta.P)),
                          c(createFolds(as.numeric(row.names.data.frame(dta.P)), k = num.folds, list = FALSE)))
      
      for (i in 1:(num.crossvalid - 1)) {
        idx.M <- data.frame(idx.M, c(createFolds(as.numeric(row.names.data.frame(dta.P)), k = num.folds, list = FALSE)))
      }
      
      #' column name according to the round of cross-validation
      colnames(idx.M) <- c('RowNo', 1:num.crossvalid)
      
      #' Counterpart to index-matrix for running the models
      x    <- c(1:(num.crossvalid * num.folds))
      run  <- c(rep(1:num.crossvalid, each = num.folds, times = 1))
      fold <- c(rep(1:num.folds, each = 1, times = length(x)/num.folds))
      x.M  <- data.frame(x, run, fold)
      
      rm(run,fold)
      
      #' Create common vector for the evaluation of runs
      run <- c(1:num.crossvalid)
      run <- as.numeric(run)
      cor <- c(1:num.crossvalid)
      cor <- as.numeric(cor)
      
      result <- list()
      
      #' Select number of cores based on the available resources
      no_cores <- detectCores() - 2
      
      cl <- makeCluster(no_cores)
      
      clusterExport(cl = cl, envir = environment(),
                    varlist = list('dta.P', 'ETA.lst', 'idx.M', 'x.M', 'result', 'sig_markers', 'nIter', 'burnIn', 'i', 'para.EGBLUP'))
      
      tryLog(result <- parLapply(cl, x, para.EGBLUP), write.error.dump.file = TRUE)
      stopCluster(cl)
      
      #' Evaluation of model
      Pred.abil.EGBLUP <- data.frame(run, cor)
      PredEB <- data.frame(ncol(6))
      
      #' Combine the results for all folds of the five-fold cross-validations with 
      #' the information of the respective fold of the cross-validation
      for(j in 1:length(x)){
        result.loop <- as.data.frame(result[[j]][1])
        
        test.set <- c(which(is.na(result.loop$y)))
        
        p.test <- data.frame(dta.P,
                             result.loop$yHat,
                             c(rep(j, each = 1, times = nrow(dta.P))),
                             c(rep(x.M[j, 'run'], each = 1, times = nrow(dta.P))),
                             c(rep(x.M[j, 'fold'], each = 1, times = nrow(dta.P))))
        
        p.test <- p.test[test.set,]
        
        PredEB <- rbind(PredEB, p.test)
      }
      
      colnames(PredEB) <- c('Genotype_Matrix', 'Trait', 'pred.Trait', 'x', 'run', 'fold')
      
      crossvalid_result <- aggregate(cbind(Trait, pred.Trait) ~ Genotype_Matrix, data = PredEB, FUN = mean)
      
      split_prediction[[split_num]]  <- crossvalid_result
      split_PredEB[[split_num]] <- PredEB
    }
    
    unlink(c('mu.dat', 'varE.dat', list.files(pattern = 'ETA_.*\\.dat')), force = TRUE)
  }
  
  if (pred_mode == TRUE) {
    predictions <- Reduce(rbind, split_prediction)
  } else {
    crossvalid_result <- Reduce(rbind, split_prediction)
    PredEB            <- Reduce(rbind, split_PredEB)
    
    #' Calculate correlation between all folds of each five-fold cross-validations
    for(i in 1:length(unique(x.M$run))) {
      cor.set <- PredEB[c(which(PredEB$run == i)),]
      Pred.abil.EB <- cor(cor.set[, c('Trait', 'pred.Trait')])
      Pred.abil.EGBLUP[Pred.abil.EGBLUP[, 'run'] == i, 2] <- Pred.abil.EB[1, 2]
    }
  }
  
  model_name <- 'EGBLUP'
  if (!is.null(groups) & is.null(sig_markers) & is.null(climate)) { model_name <- 'Kinship' }
  if (is.null(groups) & !is.null(sig_markers) & is.null(climate)) { model_name <- 'Fixing' }
  if (is.null(groups) & is.null(sig_markers) & !is.null(climate)) { model_name <- 'Climatic' }
  
  if (!is.null(groups) & !is.null(sig_markers) & is.null(climate)) { model_name <- 'Kinship+Fixing' }
  if (!is.null(groups) & is.null(sig_markers) & !is.null(climate)) { model_name <- 'Kinship+Climatic' }
  if (is.null(groups) & !is.null(sig_markers) & !is.null(climate)) { model_name <- 'Fixing+Climatic' }
  
  if (!is.null(groups) & !is.null(sig_markers) & !is.null(climate)) { model_name <- 'Kinship+Fixing+Climatic' }
  
  if (pred_mode == TRUE) {
    return(predictions)
  } else {
    model <- rep(model_name, num.crossvalid)
    model_accuracy <- cbind(model, Pred.abil.EGBLUP)
    
    return(model_accuracy)
  }
}

#' Function of the model
para.EGBLUP <- function(i) {
  #' Link the current run and fold with information which is assigned in the index-matrix
  idx <- idx.M[, c(1, (x.M[i, 'run'] + 1))]
  
  p.train <- dta.P$Trait
  
  #' Replace the phenotypic values of the test set with NA
  p.train[idx[, 2] == x.M[i, 'fold']] <- NA
  
  # set.seed(1002)
  EGBLUP <- BGLR::BGLR(y = p.train,
                       ETA = ETA.lst,
                       nIter = nIter,  # set >= 10000 for large data
                       burnIn = burnIn,  # set >= 3000 for large data
                       verbose = FALSE)
  
  result[[i]] <- list(list(data.frame(dta.P$Genotype_Matrix,EGBLUP[1],EGBLUP[14],EGBLUP[15])),
                      list(data.frame(EGBLUP[6],EGBLUP[7],EGBLUP[16],EGBLUP[17],EGBLUP[18],EGBLUP[19])))
}
