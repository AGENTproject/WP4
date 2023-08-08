observeEvent(
  input$geno_input, 
  if (input$geno_input == 'file') {
    shinyjs::show('geno_file')
    shinyjs::hide('geno_url')
  } else {
    shinyjs::hide('geno_file')
    shinyjs::show('geno_url')
  }
)

geno_data <- reactive({
  if (input$geno_input == 'file') {
    if (is.null(input$geno_file)) return(NULL)
    snps_file <- input$geno_file$datapath
  } else {
    if (input$geno_url == '') return(NULL)
    snps_file <- input$geno_url
  }

  # library(vcfR)
  # vcf <- read.vcfR(file.choose())
  # 
  # SNPs_info <- vcfR2tidy(vcf, info_only = TRUE)$fix
  # View(SNPs_info)
  # 
  # gt <- extract.gt(vcf, as.numeric = TRUE)
  # View(gt[1:100,1:100])
  
  show_modal_spinner('fading-circle', text = 'Loading...')
  df <- as.data.frame(data.table::fread(snps_file, sep = '\t', header = TRUE))
  remove_modal_spinner()
  
  hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 
                       'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode')

  if (!all(colnames(df)[1:11] == hapmap_snp_attr)) {
    show_alert(title = 'Error !!', text = 'Not a valid HapMap file format :-(', type = 'error')
    return(NULL)
  }
  
  first_row   <- df[1, -c(1:11)]
  valid_IUPAC <- c('A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N')

  #' IUPAC single-letter code 
  if (all(first_row %in% valid_IUPAC)) {
    
    show_modal_spinner('fading-circle', text = 'Converting...')
    df <- hapMapChar2Numeric(df)
    remove_modal_spinner()
    
  #' -1, 0, 1 numeric coding
  } else if (min(as.numeric(first_row), na.rm = TRUE) == -1 & 
             max(as.numeric(first_row), na.rm = TRUE) == 1) {

    df <- cbind(df[, 1:11], 
                data.frame(apply(df[, -c(1:11)], 2, function(x) 1 + as.numeric(as.character(x)))))
    
  #' 0, 1, 2 numeric coding
  } else if (min(as.numeric(first_row), na.rm = TRUE) == 0 & 
             max(as.numeric(first_row), na.rm = TRUE) == 2) {
    
    df <- cbind(df[, 1:11], 
                data.frame(apply(df[, -c(1:11)], 2, function(x) as.numeric(as.character(x)))))
  
  #' something else!
  } else {
    show_alert(title = 'Error !!', text = 'Not a valid HapMap file format :-(', type = 'error')
    return(NULL)
  }
  
  return(df)
})

output$geno_summary <- renderText({
  if (!is.null(geno_data()) & input$pheno_id != '') {
    paste(
      " Data Integrity Checks:\n",
      
      sum(colnames(geno_data()) %in% pheno_data()[, input$pheno_id]),
      "Accessions exist in both phenotypic and genotypic files (will be used to train the model)\n",
      
      sum(!colnames(geno_data()) %in% pheno_data()[, input$pheno_id]),
      "Accessions have genotypic data but no phenotypic (will be predicted, add to pheno data file with NA value)\n",
      
      sum(!pheno_data()[, input$pheno_id] %in% colnames(geno_data())),
      'Accessions have phenotypic data but no genotypic (will filtered them out from the pheno dataset)'
    )
  }
})

output$chrom_summary <- renderTable({
  if (!is.null(geno_data())) {
    data.frame(
      chrom = unique(geno_data()[,'chrom']),
      min_pos = aggregate(pos ~ chrom, data = geno_data(), FUN = min)[,2],
      max_pos = aggregate(pos ~ chrom, data = geno_data(), FUN = max)[,2],
      snps_count = aggregate(pos ~ chrom, data = geno_data(), FUN = length)[,2]
    )
  }
})

output$geno_snps <- renderInfoBox({
  infoBox(
    title = 'Geno SNPs', 
    value = nrow(geno_data()), 
    icon  = icon('dna'), 
    color = 'light-blue', 
  )
})

output$geno_indv <- renderInfoBox({
  infoBox(
    title = 'Geno Acc.', 
    value = ncol(geno_data()) - 11, 
    icon  = icon('vial'), 
    color = 'light-blue', 
  )
})

output$unassigned_chrom <- renderUI({
  selectInput(
    inputId = 'unassigned_chrom', 
    label   = 'Unassigned Chromosome:', 
    width   = '200px', 
    choices = as.list(c('', unique(geno_data()[,'chrom'])))
  )
})

observeEvent(
  input$geno_example, 
  if (input$geno_example) {
    updateSelectInput(session, 'geno_input', selected = 'url')
    
    geno_example_url <-  paste0(session$clientData$url_protocol, '//',
                                session$clientData$url_hostname, ':',
                                session$clientData$url_port,
                                session$clientData$url_pathname,
                                'example/barley_winter_genotypes_imputed.hmp.txt.gz')
    
    updateTextInput(session, 'geno_url', value = geno_example_url)
    
    shinyjs::hide('geno_file')
    shinyjs::show('geno_url')
  } else {
    updateSelectInput(session, 'geno_input', selected = 'file')
    updateTextInput(session, 'geno_url', value = '')
    
    shinyjs::show('geno_file')
    shinyjs::hide('geno_url')
  }
)
