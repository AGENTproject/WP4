job_id   <- paste0('AGENT_T4_JOB_', format(Sys.time(), '%Y%m%d%H%M%S'))
job_file <- paste0(job_id, '.RData')

observeEvent(
  input$perform_gwas | input$rename_files, 
  if (input$perform_gwas & input$rename_files) {
    shinyjs::show('sub_matrix') 
    shinyjs::show('stats_file') 
  } else {
    shinyjs::hide('sub_matrix') 
    shinyjs::hide('stats_file') 
  }
)

observeEvent(
  input$min_subset | input$rename_files, 
  if (input$min_subset & input$rename_files) {
    shinyjs::show('min_grp_acc_file') 
    shinyjs::show('min_grp_snp_file') 
  } else {
    shinyjs::hide('min_grp_acc_file') 
    shinyjs::hide('min_grp_snp_file') 
  }
)

observeEvent(
  input$perform_gs | input$rename_files, 
  if (input$perform_gs & input$rename_files) {
    shinyjs::show('egblup_results_file') 
  } else {
    shinyjs::hide('egblup_results_file') 
  }
)

observeEvent(
  input$submit, 
  {
    BLUE <- pheno_data()
    if (is.null(BLUE)) {
      show_alert(title = 'Error !!', text = 'No Phenotypic Data...', type = 'error')
      updateTabsetPanel(session, 'tabset', selected = 'Tab1')
      return(NULL)
    }
    
    config <- list()

    config$acc_id <- input$pheno_id
    if (input$pheno_id == '') {
      show_alert(title = 'Error !!', text = 'No Phenotypic ID Selected :-(', type = 'error')
      updateTabsetPanel(session, 'tabset', selected = 'Tab1')
      return(NULL)
    }

    config$trait <- input$pheno_trait
    if (input$pheno_trait == '') {
      show_alert(title = 'Error !!', text = 'No Phenotypic Trait Selected :-(', type = 'error')
      updateTabsetPanel(session, 'tabset', selected = 'Tab1')
      return(NULL)
    }
    
    SNPs <- geno_data()
    if (is.null(SNPs)) {
      # show_alert(title = 'Error !!', text = 'No Genotypic Data :-(', type = 'error')
      # updateTabsetPanel(session, 'tabset', selected = 'Tab2')
      # return(NULL)
      geno_alert <- 'We will use the default genotypic dataset AGENT_barley_gbs_datashare_230511'
    } else {
      geno_alert <- ''
    }
    
    config$email <- input$email
    if (!grepl('\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>', 
               as.character(input$email), ignore.case = TRUE)) {
      show_alert(title = 'Error !!', text = 'No Valid Email Address :-(', type = 'error')
      updateTabsetPanel(session, 'tabset', selected = 'Tab5')
      return(NULL)
    }
    
    config$snp_filtering <- input$snp_filtering
    config$snp_rate      <- input$snp_rate
    
    config$acc_filtering <- input$acc_filtering
    config$acc_rate      <- input$acc_rate
    
    config$maf_filtering <- input$maf_filtering
    config$maf_rate      <- input$maf_rate
    
    config$impute        <- input$impute
    config$snprune       <- input$snprune
    
    config$subpop        <- input$subpop
    config$subpop_val    <- input$subpop_val
    
    config$ld            <- input$ld
    config$ld_val        <- input$ld_val

    config$perform_gwas  <- input$perform_gwas
    config$gwas_lod      <- input$gwas_lod

    config$perform_gs    <- input$perform_gs
    config$eg_blup       <- input$eg_blup
    config$kinship_split <- input$kinship_split
    config$fixed_markers <- input$fixed_markers
    config$climate_model <- input$climate_model

    config$stepwise      <- input$stepwise
    config$min_subset    <- input$min_subset
    config$QTLome_submit <- input$QTLome_submit
    
    config$sub_matrix <- input$sub_matrix
    config$stats_file <- input$stats_file
    
    config$min_grp_acc_file <- input$min_grp_acc_file
    config$min_grp_snp_file <- input$min_grp_snp_file
    
    config$egblup_results_file <- input$egblup_results_file

    config$unassigned_chrom <- input$unassigned_chrom

    show_modal_spinner('fading-circle', text = 'Saving...')
    save(config, SNPs, BLUE, file = job_file, compress = TRUE, compression_level = 6)
    remove_modal_spinner()
    
    show_alert(
      title = 'Success !!',
      text = span(geno_alert, br(), br(), 'Job ID: ', downloadLink('download_job', label = job_id)), 
      type = 'success',
    )
  }
)

output$download_job <- downloadHandler(
  filename <- function() {
    job_file
  },
  
  content <- function(file) {
    file.copy(paste0(getwd(), '/', job_file), file)
  },
  
  contentType = 'application/zip'
)
