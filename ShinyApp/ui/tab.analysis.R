fluidRow(
  style = 'padding: 30px;',
  
  #' Perform GWAS Analysis?
  materialSwitch(
    inputId = 'perform_gwas', 
    label   = 'Perform GWAS Analysis', 
    status  = 'primary', 
    right   = TRUE
  ),

  tags$blockquote(
    id = 'gwas_settings',
    
    #' Significant threshold: set (LOD value)
    numericInput(
      inputId = 'gwas_lod', 
      label = 'Significant threashold (LOD value):', 
      value = 3.5, 
      width = '300px'
    ),
    
    #' Bonferroni correction: 1 - (1 - 0.05)^(1/#snp)
    verbatimTextOutput('alpha_correction'),
    hr(),

    #' Perform stepwise regression: yes or no
    materialSwitch(
      inputId = 'stepwise', 
      label   = 'Perform Stepwise Regression', 
      status  = 'primary', 
      right   = TRUE
    ),
    hr(),
    
    #' Propose to QTLome (get downloadable file to submit)
    materialSwitch(
      inputId = 'QTLome_submit', 
      label   = 'Propose to QTLome', 
      status  = 'primary', 
      right   = TRUE
    ),
    hr(),
    
    #' Get the minimum number of accessions to include all loci: yes or no
    materialSwitch(
      inputId = 'min_subset', 
      label   = 'Get the minimum number of accessions to include all loci', 
      status  = 'primary', 
      right   = TRUE
    ),
  ),
  hr(),
  
  #' Perform Genomic Prediction? Models for testing:
  materialSwitch(
    inputId = 'perform_gs', 
    label   = 'Perform Genomic Prediction', 
    status  = 'primary', 
    right   = TRUE
  ),
  
  tags$blockquote(
    id = 'gs_models',
    
    #' EG-BLUP (default, mandatory with reference to Marcel paper)
    checkboxInput(
      inputId = 'eg_blup', 
      label = span('EG-BLUP (ref. ',
                   a('https://doi.org/10.1007/s00122-022-04227-4', target = '_blank',
                     href = 'https://doi.org/10.1007/s00122-022-04227-4'), ')'), 
      value = TRUE
    ),
    
    #' Kinship split into sub-population (only active if sub-population option was not none)
    disabled(
      checkboxInput(
        inputId = 'kinship_split', 
        label = 'Kinship split into sub-population', 
        value = FALSE
      )
    ),
    
    #' Markers as fixed (check with Marcel for BGLR implementation or use ASReml-R v4.2)
    disabled(
      checkboxInput(
        inputId = 'fixed_markers', 
        label = 'Markers as fixed effect', 
        value = FALSE
      )
    ),
    
    #' Inclusion of climatic matrix (cords required, use QBMS & TerraClimate, FIGS+ codes)
    disabled(
      checkboxInput(
        inputId = 'climate_model', 
        label = 'Inclusion of climatic matrix', 
        value = FALSE
      )
    ),
  ),
)