output$alpha_correction <- renderText({
  if (!is.null(geno_data())) {
    paste(
      'Suggested LOD value (Bonferroni correction): ',
      round(-log10(1 - (1 - 0.05)^(1/nrow(geno_data()))), digits = 1)
    )
  }
})

observeEvent(
  input$perform_gwas, 
  if (input$perform_gwas) {
    shinyjs::show('gwas_settings') 
  } else {
    shinyjs::hide('gwas_settings')
  }
)

observeEvent(
  input$perform_gs, 
  if (input$perform_gs) {
    shinyjs::show('gs_models') 
  } else {
    shinyjs::hide('gs_models') 
  }
)
