hapMapChar2Numeric <- function(hapMap) {
  # http://adv-r.had.co.nz/C-interface.html
  convertChar2Numeric <- inline::cfunction(signature(SNPsMatrix="character", 
                                                     refAllele="character", 
                                                     rowsNum="integer", 
                                                     colsNum="integer"),
  "int i,j;
  
  // matrix dimentions
  int r = asInteger(rowsNum);
  int c = asInteger(colsNum);
  
  // length of the matrix
  int length = r*c;

  // create matrix of integers with the same size as SNPsMatrix
  SEXP SNPsNum;
  PROTECT(SNPsNum = allocMatrix(INTSXP, r, c));
  
  // convert SNPs codes from the standard IUPAC code (single char) to 
  // numeric 0, 1, or 2 (use 1 for heterozygous)
  for(i = 0; i < r; i++){
    char* x;
    char alleleA;

    // we need to get the reference allele in each SNP (row)
    x = (char*)CHAR(STRING_ELT(refAllele, i));
    alleleA = x[0];
    
    // convert SNPsMatrix to numeric 0,1,2
    // now with alleleA we can convert the genotypes to numeric 0, 1, 2
    for(j = 0; j < c; j++){
      x = (char*)CHAR(STRING_ELT(SNPsMatrix, i*c+j));
      
      // if current SNP is the same of reference allele (alleleA)
      if(x[0] == alleleA){
        // then assign 0 in the SNPsNum matrix
        // take care of the order of the indexes in matrix is by columns
        INTEGER(SNPsNum)[j*r + i] = 0;
      }else if(x[0] == 'A' || x[0] == 'T' || x[0] == 'C' || x[0] == 'G'){
        // if it is homozygous allele [A,T,C,G] 
        // but not alleleA (i.e., minor allele)
        INTEGER(SNPsNum)[j*r + i] = 2;
      }else if(x[0] == 'N'){
        // if it is missing allele [N] 
        INTEGER(SNPsNum)[j*r + i] = -9;
      }else{
        // if it is not (i.e., heterozygous)
        INTEGER(SNPsNum)[j*r + i] = 1;
      }
    }
  }
  UNPROTECT(1);
  return(SNPsNum);")

  hapMap <- as.data.frame(hapMap)
  
  # extract SNP infomation , which is the first 11 columns
  SNPInfo <- hapMap[,1:11]
  
  # remove the first 11 columns
  hapMap <- hapMap[,-c(1:11)]
  
  # convert the hapMap to numeric
  hapMapNumeric <- convertChar2Numeric(unlist(as.matrix(t(hapMap))),
                                       unlist(as.matrix(substr(SNPInfo$alleles,1,1))),
                                       as.integer(nrow(hapMap)),
                                       as.integer(ncol(hapMap)))
  
  # convert to data frame
  hapMapNumeric <- as.data.frame(hapMapNumeric)
  
  # convert -9 values to NA
  hapMapNumeric[hapMapNumeric == -9] <- NA
  
  # get back the column names (accessions)
  colnames(hapMapNumeric) <- colnames(hapMap)
  
  return(cbind(SNPInfo, hapMapNumeric))
}
