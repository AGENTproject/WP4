#' Get the settings of the genotypic data filtering 
#' and required pre-processing steps

fluidRow(
  style = 'padding: 30px;',
  
  #' Checkbox for SNPs call rate filtering (proportion with no missing)
  #' if checked ask for ratio (default 10%)
  materialSwitch(
    inputId = 'snp_filtering', 
    label   = 'SNPs call rate filtering (proportion with no missing)', 
    status  = 'primary', 
    right   = TRUE
  ),
  sliderInput(
    inputId = 'snp_rate', 
    label   = NULL, 
    min     = 0, 
    max     = 100, 
    value   = 10, 
    width   = '300px', 
    post    = '%'
  ),
  hr(),
  
  #' Checkbox for Acc. call rate filtering (proportion with no missing)
  #' if checked ask for ratio (default 20%)
  materialSwitch(
    inputId = 'acc_filtering', 
    label   = 'Accessions call rate filtering (proportion with no missing)', 
    status  = 'primary', 
    right   = TRUE
  ),
  sliderInput(
    inputId = 'acc_rate', 
    label   = NULL, 
    min     = 0, 
    max     = 100, 
    value   = 20, 
    width   = '300px', 
    post    = '%'
  ),
  hr(),
  
  #' Checkbox for MAF filtering, if checked ask for ratio (default is 4/#Acc.)
  materialSwitch(
    inputId = 'maf_filtering', 
    label   = 'MAF (Minor Allele Frequency) filtering', 
    status  = 'primary', 
    right   = TRUE
  ),
  sliderInput(
    inputId = 'maf_rate', 
    label   = NULL, 
    min     = 0, 
    max     = 50, 
    value   = 5,
    step    = 0.1,
    width   = '300px', 
    post    = '%'
  ),
  verbatimTextOutput('maf_default'),
  hr(),
  
  actionButton('filter_genodata', 'Calculate active number of SNPs and Accessions'),
  tags$br(), tags$br(),
  verbatimTextOutput('active_genodata'),
  hr(),
  
  #' Conduct imputation option: none or mean 
  #' (possible other options: fixed, random, or beagle)
  selectInput(
    inputId = 'impute', 
    label   = 'Imputation: ', 
    choices = list('None' = FALSE, 'Mean' = TRUE), 
    width   = '200px'
  ),
  
  #' Report the ratio of total missing data (interactive)
  verbatimTextOutput('missing_summary'),
  hr(),
  
  #' Sub-population structure settings: 
  #' none, set (#sub-populations), auto (get best value based on AMOVA)
  selectInput(
    inputId = 'subpop', 
    label   = 'Sub-Population Structure: ', 
    width   = '200px',
    choices = list('None' = FALSE, 
                   'Set (#sub-populations)' = 'set', 
                   'Auto (based on AMOVA)'  = 'auto')
  ),
  numericInput(
    inputId = 'subpop_val', 
    label = '#Sub-Populations', 
    value = 1, 
    width = '200px'
  ),
  hr(),
  
  #' Haplotype SNPrune: yes or no
  materialSwitch(
    inputId = 'snprune', 
    status  = 'primary', 
    right = TRUE,
    label = span('Haplotype SNPrune (ref. ',
                 a('https://doi.org/10.1186/s12711-018-0404-z', target = '_blank',
                   href = 'https://doi.org/10.1186/s12711-018-0404-z'), ')')
  ),
  hr(),
  
  #' LD value: set (get the value in bp), auto (calculate it)
  selectInput(
    inputId = 'ld', 
    choices = list('Auto' = 'auto', 'Set' = 'set'), 
    label = 'Linkage Disequilibrium (LD): ', 
    width = '200px'
  ),
  numericInput(
    inputId = 'ld_val', 
    label = 'LD Value', 
    value = 0, 
    width = '200px'
  ),
)