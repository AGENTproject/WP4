observeEvent(
  input$snp_filtering, 
  if (input$snp_filtering) {
    shinyjs::show('snp_rate') 
  } else {
    shinyjs::hide('snp_rate')
  }
)

observeEvent(
  input$acc_filtering, 
  if (input$acc_filtering) {
    shinyjs::show('acc_rate')
  } else {
    shinyjs::hide('acc_rate')
  }
)

observeEvent(
  input$maf_filtering, 
  if (input$maf_filtering) {
    shinyjs::show('maf_rate') 
    shinyjs::show('maf_default') 
  } else {
    shinyjs::hide('maf_rate')
    shinyjs::hide('maf_default')
  }
)

output$maf_default <- renderText({
  if (!is.null(geno_data())) {
    paste(
      'Suggested MAF ratio (4/#Acc.): ',
      round(100 * 4 / (ncol(geno_data())-11), digits = 1),
      '%'
    )
  }
})


observeEvent(
  input$filter_genodata, {
    # snp.attr <- df[,1:11]
    # geno.M   <- df[,-c(1:11)]
    # rownames(geno.M) <- snp.attr[,1]
    # geno.M <- t(geno.M)
    # 
    # qc <- qc.filtering(M = geno.M, marker.callrate = 0.1, ind.callrate = 0.2, maf = 0.05, impute = FALSE)
    # 
    # qc.snp.attr <- merge(colnames(qc$M.clean), snp.attr, by.x = 1, by.y = 1)
    # colnames(qc.snp.attr) <- colnames(snp.attr)
    # 
    # df2 <- cbind(qc.snp.attr, t(qc$M.clean))
    # rownames(df2) <- NULL

    show_modal_spinner('fading-circle', text = 'Calculating...')
    
    geno.M <- geno_data()[, -c(1:11)]
    rownames(geno.M) <- geno_data()[, 1]
    geno.M <- t(geno.M)
    
    qc <- qc.filtering(
      M = geno.M, 
      marker.callrate = ifelse(input$snp_filtering, input$snp_rate/100, 1), 
      ind.callrate = ifelse(input$acc_filtering, input$acc_rate/100, 1), 
      maf = ifelse(input$maf_filtering, input$maf_rate/100, 0), 
      message = FALSE
    )
    
    output$active_genodata <- renderText({ 
      paste(
        'SNPs call rate:',
        ifelse(input$snp_filtering, input$snp_rate/100, 1),
        
        '| Acc. call rate:',
        ifelse(input$acc_filtering, input$acc_rate/100, 1),
        
        '| MAF threshold:',
        ifelse(input$maf_filtering, input$maf_rate/100, 0),
        
        "\nActive number of Acc.:", 
        dim(qc$M.clean)[1],
        
        "\nActive number of SNPs:",
        dim(qc$M.clean)[2]
      )
    })
    
    remove_modal_spinner()
  }
)

output$missing_summary <- renderText({
  if (!is.null(geno_data())) {
    paste(
      'Missing ratio:',
      round((100 * sum(is.na(geno_data()))) / (nrow(geno_data()) * (ncol(geno_data())-11)), digits = 2),
      '%'
    )
  }
})

observeEvent(
  input$subpop, 
  if (input$subpop == 'set') {
    shinyjs::show('subpop_val')
  } else {
    shinyjs::hide('subpop_val')
  }
)

observeEvent(
  input$ld, 
  if (input$ld == 'set') {
    shinyjs::show('ld_val') 
  } else {
    shinyjs::hide('ld_val')
  }
)
