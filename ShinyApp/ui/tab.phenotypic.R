#' Get the phenotypic BLUEs file 
#' (WP3 output, WP5 data validation, and WP6 data repository/source)

fluidRow(
  style = 'padding: 30px;',
  
  #' Source: upload (web interface to temp local directory) or URL 
  #' Verify file format (CSV comma delimited)
  #' TODO: integrate with AGENT FAIRDOM API (user, pass, data file id & version)
  selectInput(
    inputId = 'pheno_input', 
    label   = 'Phenotypic BLUEs Source*: ', 
    choices = list('Upload File' = 'file', 'Copy URL' = 'url'), 
    width   = '200px'
  ),
  fileInput(
    inputId = 'pheno_file', 
    label   = NULL, 
    width   = '400px', 
    accept  = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
  ),
  textInput(
    inputId = 'pheno_url', 
    label   = NULL, 
    value   = '', 
    width   = '400px', 
    placeholder = 'https://urgi.versailles.inrae.fr/fairdom/'
  ),
  checkboxInput(
    inputId = 'pheno_example', 
    label = span('Load example ',
                 a('phenotypic data', target = '_blank',
                   href = 'example/barley_winter_all_BLUEs_phenotype.csv')), 
    value = FALSE
  ),
  hr(),
  
  #' Set pheno ID column: default is the first column, or user select the 
  #' associated column from a list
  uiOutput('pheno_id'),
  verbatimTextOutput('pheno_id_summary'),
  hr(),
  
  #' TODO: Set coordinates columns (optional): 
  #' Required for FIGS+, latitude & longitude in decimal degree format
  
  #' Set the trait column: List of the CSV file column headers excluding already
  #' selected columns & report summary statistics to let the user verify before 
  #' proceeding to the next step
  uiOutput('pheno_trait'),
  verbatimTextOutput('pheno_trait_summary'),
)