observeEvent(
  input$pheno_input, 
  if (input$pheno_input == 'file') {
    shinyjs::show('pheno_file')
    shinyjs::hide('pheno_url')
  } else {
    shinyjs::hide('pheno_file')
    shinyjs::show('pheno_url')
  }
)

pheno_data <- reactive({
  if (input$pheno_input == 'file') {
    if (is.null(input$pheno_file)) return(NULL)
    return(read.csv(input$pheno_file$datapath))
  } else {
    if (input$pheno_url == '') return(NULL)
    return(read.csv(input$pheno_url))
  }
})

output$pheno_id <- renderUI({
  selectInput(
    inputId = 'pheno_id', 
    label   = 'Entry ID*:', 
    width   = '200px', 
    choices = as.list(c('', colnames(pheno_data())))
  )
})

output$pheno_trait <- renderUI({
  header <- colnames(pheno_data())
  selectInput(
    inputId = 'pheno_trait', 
    label   = 'Phenotypic Trait*:', 
    width   = '200px', 
    choices = as.list(c('', header[header != input$pheno_id]))
  )
})

output$pheno_id_summary <- renderText({
  if (!is.null(input$pheno_id) & input$pheno_id != '') {
    paste(
      'First values: ',
      paste(head(pheno_data()[, input$pheno_id]), collapse = ', '),
      ", etc.\nDistinct values: ",
      length(unique(pheno_data()[, input$pheno_id]))
    )
  }
})

output$pheno_trait_summary <- renderPrint({
  if (!is.null(input$pheno_trait) & input$pheno_trait != '') {
    summary(pheno_data()[, input$pheno_trait])
  }
})

output$pheno_rows <- renderInfoBox({
  infoBox(
    title = 'Pheno Records', 
    value = nrow(pheno_data()), 
    icon  = icon('list'), 
    color = 'light-blue', 
  )
})

output$pheno_cols <- renderInfoBox({
  infoBox(
    title = 'Pheno Columns', 
    value = ncol(pheno_data()), 
    icon  = icon('seedling'), 
    color = 'light-blue', 
  )
})

observeEvent(
  input$pheno_example, 
  if (input$pheno_example) {
    updateSelectInput(session, 'pheno_input', selected = 'url')
    
    pheno_example_url <-  paste0(session$clientData$url_protocol, '//',
                                 session$clientData$url_hostname, ':',
                                 session$clientData$url_port,
                                 session$clientData$url_pathname,
                                 'example/barley_winter_all_BLUEs_phenotype.csv')
    
    updateTextInput(session, 'pheno_url', value = pheno_example_url)
    
    shinyjs::hide('pheno_file')
    shinyjs::show('pheno_url')
  } else {
    updateSelectInput(session, 'pheno_input', selected = 'file')
    updateTextInput(session, 'pheno_url', value = '')
    
    shinyjs::show('pheno_file')
    shinyjs::hide('pheno_url')
  }
)
