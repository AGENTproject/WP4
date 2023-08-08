fluidRow(
  style = 'padding: 30px;',
  
  textInput(
    inputId = 'email',
    label = 'Your Email Address:',
    placeholder = 'test@example.com',
    width = '400px'
  ),
  hr(),

  materialSwitch(
    inputId = 'rename_files', 
    label   = 'Do you need to change the standard output file names?', 
    status  = 'primary', 
    right   = TRUE
  ),
  
  fluidRow(
    column(
      6, 
      textInput(
        inputId = 'sub_matrix',
        label = 'File to save the sub-genotyping matrix:',
        value = 'sub-genotyping_matrix.csv',
        width = '400px'
      ),
    ),
    
    column(
      6, 
      textInput(
        inputId = 'stats_file',
        label = 'File to save the SNPs stats for each landrace:',
        value = 'snps_stats_for_each_landrace.csv',
        width = '400px'
      ),
    ),
  ),
  
  fluidRow(
    column(
      6, 
      textInput(
        inputId = 'min_grp_acc_file',
        label = 'File to save the minimum group accessions:',
        value = 'minimum_group_accessions.csv',
        width = '400px'
      )
    ),
    
    column(
      6, 
      textInput(
        inputId = 'min_grp_snp_file',
        label = 'File to save the minimum group SNPs:',
        value = 'minimum_group_snps.csv',
        width = '400px'
      ),
    ),
  ),
  
  textInput(
    inputId = 'egblup_results_file',
    label = 'File to save the EGBLUP predictions:',
    value = 'EGBLUP_predictions.csv',
    width = '400px'
  ),

  div(align = 'center', actionButton('submit', label = 'Submit', width = '300px')),
)