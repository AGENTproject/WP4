snps_numeric_recoding <- function(M, ref.allele) {
  # http://adv-r.had.co.nz/C-interface.html
  convertChar2Numeric <- inline::cfunction(signature(SNPsMatrix="character", refAllele="character", rowsNum="integer", colsNum="integer"),
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
  // numeric 0, 1, or 2 (use 1 for heterozygous and 2 for reference allele)
  for (i = 0; i < r; i++) {
    char* x;
    char alleleA;

    // we need to get the reference allele in each SNP (row)
    x = (char*)CHAR(STRING_ELT(refAllele, i));
    alleleA = x[0];
    
    // now with alleleA we can convert SNPsMatrix to numeric values 0, 1, 2, -9 
    // for alternative, heterozygous, reference, and missing allele respectively 
    for (j = 0; j < c; j++) {
      x = (char*)CHAR(STRING_ELT(SNPsMatrix, i*c+j));
      
      if (x[0] == alleleA) {
      
        // if current SNP is the same of reference allele (alleleA)
        // then assign 2 in the SNPsNum matrix
        // take care of the order of the indexes in matrix is by columns
        
        INTEGER(SNPsNum)[j*r + i] = 2;
        
      } else if(x[0] == 'A' || x[0] == 'T' || x[0] == 'C' || x[0] == 'G') {

        // if it is homozygous allele [A,T,C,G] 
        // but not alleleA (i.e., minor allele)
        
        INTEGER(SNPsNum)[j*r + i] = 0;

      } else if(x[0] == 'N') {
      
        // if it is missing allele [N] 
        INTEGER(SNPsNum)[j*r + i] = -9;
        
      } else {
      
        // if it is not (i.e., heterozygous)
        INTEGER(SNPsNum)[j*r + i] = 1;
        
      }
    }
  }
  
  UNPROTECT(1);
  return(SNPsNum);")
  
  M <- as.data.frame(M)

  # convert the hapMap to numeric
  numericM <- convertChar2Numeric(unlist(as.matrix(t(M))),
                                       unlist(as.matrix(ref.allele)),
                                       as.integer(nrow(M)),
                                       as.integer(ncol(M)))
  
  # convert to data frame
  numericM <- as.data.frame(numericM)
  
  # convert -9 values to NA
  numericM[numericM == -9] <- NA
  
  # get back the column names (accessions)
  colnames(numericM) <- colnames(M)
  
  return(numericM)
}
