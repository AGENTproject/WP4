#' Get the genotypic SNPs file 
#' (WP2 output, WP5 data validation, and WP6 data repository/source)

fluidRow(
  style = 'padding: 30px;',
  
  #' Source: Upload (web interface to temp local directory) or URL (optional username/password to access)
  #' Accept *.gz format (7-Zip how-to reference), average genomic file size after compression is 5%
  selectInput(
    inputId = 'geno_input', 
    label   = 'Genotypic SNPs Source*:', 
    choices = list('Upload File' = 'file', 'Copy URL' = 'url'), 
    width   = '200px'
  ),
  fileInput(
    inputId = 'geno_file', 
    label   = NULL, 
    width   = '400px', 
    accept  = c('application/gzip', '.gz', '.txt', '.hmp')
  ),
  textInput(
    inputId = 'geno_url', 
    label   = NULL, 
    value   = '', 
    width   = '400px', 
    placeholder = 'https://urgi.versailles.inrae.fr/fairdom/'
  ),
  checkboxInput(
    inputId = 'geno_example', 
    label = span('Load example ',
                 a('genotypic data', target = '_blank',
                   href = 'example/barley_winter_genotypes_imputed.hmp.txt.gz')), 
    value = FALSE
  ),
  
  fluidRow(
    style = 'padding-right: 0px; padding-left: 0px;',
    box(
      title = 'Notes:',
      width = 12,
      solidHeader = TRUE,
      status = 'success',
      tags$ul(
        tags$li('Accept HapMap format (tab-delimited text file with a header row). 
                The file list SNPs in rows and Accessions (individual samples) 
                in columns. The first 11 columns describe attributes of the SNP, 
                but only the first 4 columns data are required for processing: 
                rs# (SNP id), alleles (e.g., C/G), chrom (chromosome), and pos (position).'),

        tags$li(
          tags$span(
            'Accept numeric coding ([0,1, and 2] or [-1,0, and 1] for reference/major, 
            heterozygous, and alternative/minor alleles respectively), or the ', 
            tags$a('IUPAC single-letter', target = '_blank',
            href = 'https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation'),
            'code (ref. ', tags$a('https://doi.org/10.1093/nar/13.9.3021', target = '_blank',
            href = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC341218/'), ').'
          )
        ),

        tags$li('Position should be in bp (base pairs) not cM (CentiMorgan).'),

        tags$li(
          tags$span(
            'We recommend compressing your HapMap genotypic data using the gzip 
            format (*.gz extension) to significantly reduce file size. On average, 
            the compressed file size is only 5% of the original size. You can use 
            free software such as', tags$a('7-Zip', href = 'https://www.7-zip.org',
            target = '_blank'), 'to perform the compression.'
          )
        ),
      )
    ),
  ),
  
  #' Verify file format: Hapmap file format (with reference link)
  #' highlight that pos unit should be bp not cM
  #' Report summary statistics (#acc, #snps, etc.) 
  #' to let the user verify before proceeding to the next step
  tableOutput('chrom_summary'),
  
  #' Accessions exist in both phenotypic and genotypic files (will be used to train the model)
  #' Accessions have genotypic data but no phenotypic (will predict, add to pheno data file with NA value)
  #' Accessions have phenotypic data but no genotypic (filter them out from the pheno data file)
  verbatimTextOutput('geno_summary'),
  
  #' Select the unassigned chromosome id/name 
  uiOutput('unassigned_chrom'),
)