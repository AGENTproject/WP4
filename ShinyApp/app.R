library(shiny)
library(shinyjs)
library(shinybusy)
library(shinyWidgets)
library(shinydashboard)

library(utils)
library(inline)
library(data.table)

#' https://vsni.co.uk/free-software/asrgenomics
library(ASRgenomics)

#' https://vsni.co.uk/free-software/asrgwas
#' library(ASRgwas)

# Setup the Shiny app header ---------------------------------------------------

header <- dashboardHeader(disable = FALSE, title = 'T4 Job Requester')

# Setup the Shiny app sidebar --------------------------------------------------

sidebar <- dashboardSidebar(disable = FALSE, collapsed = FALSE)

# Setup the Shiny app body -----------------------------------------------------

body <- dashboardBody(
  useShinyjs(),
  tags$head(tags$link(rel = 'icon', type = 'image/x-icon', href = 'favicon.ico')),
  tags$img(src = 'logo.png', height = 168, style = 'padding: 20px;'),
  box(
    width = 12,
    infoBoxOutput('pheno_rows', width = 3),
    infoBoxOutput('pheno_cols', width = 3),
    infoBoxOutput('geno_snps', width = 3),
    infoBoxOutput('geno_indv', width = 3),
    
    tabBox(
      width = 12, 
      id    = 'tabset',
           
      tabPanel(
        title = 'Phenotypic BLUEs',
        value = 'Tab1',
        source('./ui/tab.phenotypic.R', local = TRUE)$value,
      ),
      tabPanel(
        title = 'Genotypic SNPs',
        value = 'Tab2',
        source('./ui/tab.genotypic.R', local = TRUE)$value,
      ),
      tabPanel(
        title = 'Pre-Processing',
        value = 'Tab3', 
        source('./ui/tab.pre-processing.R', local = TRUE)$value,
      ),
      tabPanel(
        title = 'Analysis Settings',
        value = 'Tab4',
        source('./ui/tab.analysis.R', local = TRUE)$value,
      ),
      tabPanel(
        title = 'Job Submit',
        value = 'Tab5',
        source('./ui/tab.submit.R', local = TRUE)$value,
      ),
      tabPanel(
        title = 'Download Results',
        value = 'Tab6',
        source('./ui/tab.results.R', local = TRUE)$value,
      ),
    ),
      
    uiOutput('navigate')
  )
)

# Setup the Shiny app UI components --------------------------------------------

ui <- dashboardPage(header, sidebar, body)

# Setup the Shiny app back-end components --------------------------------------

server <- function(input, output, session) {
  #' Automatically stop a Shiny app when closing the browser tab
  session$onSessionEnded(stopApp)
  
  #' Increase the Shiny limits file uploads to 64MB per file.
  options(shiny.maxRequestSize = 264 * 1024 ^ 2)
  
  source('./server/hapMapChar2Numeric.v2.R', local = TRUE)$value
  source('./server/run.phenotypic.R', local = TRUE)$value
  source('./server/run.genotypic.R', local = TRUE)$value
  source('./server/run.pre-processing.R', local = TRUE)$value
  source('./server/run.analysis.R', local = TRUE)$value
  source('./server/run.submit.R', local = TRUE)$value
  source('./server/run.results.R', local = TRUE)$value
  
  back_bn  <- actionButton('prev_tab', 'Back')
  next_bn  <- actionButton('next_tab', 'Next')
  tab_list <- c('Tab1', 'Tab2', 'Tab3', 'Tab4', 'Tab5', 'Tab6') 
  
  output$navigate <- renderUI({
    div(align = 'center',
      column(1, if (which(tab_list == input$tabset) != 1) back_bn),
      column(1, offset = 10, if (which(tab_list == input$tabset) != length(tab_list)) next_bn),
    )
  })
  
  observeEvent(input$prev_tab,
    {
      n <- which(tab_list == input$tabset)
      updateTabsetPanel(session, 'tabset', selected = tab_list[n - 1])
    }
  )

  observeEvent(input$next_tab,
    {
      n <- which(tab_list == input$tabset)
      updateTabsetPanel(session, 'tabset', selected = tab_list[n + 1])
    }
  )
}

# Render the Shiny app ---------------------------------------------------------

shinyApp(ui, server)
