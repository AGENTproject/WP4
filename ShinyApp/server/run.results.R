observeEvent(
  input$result_example, 
  if (input$result_example) {
    updateTextInput(session, 'job_id', value = 'AGENT_T4_JOB_20230508060723')
  } else {
    updateTextInput(session, 'job_id', value = '')
  }
)

observeEvent(
  input$download, 
  {
    if (input$job_id == '') {
      show_alert(title = 'Error !!', 
                 text = 'Insert Your Job ID First...', 
                 type = 'error')
      return(NULL)
    } else if (!file.exists(paste0('./www/output/', input$job_id, '.zip'))) {
      show_alert(title = 'Warning !!', 
                 text = 'Your Analysis Results Not Ready Yet...', 
                 type = 'warning')
      return(NULL)
    } else {
      show_alert(title = 'Success !!', 
                 text = span('Find Your Analysis Results ',
                             a('Here', target = '_blank',
                               href = paste0('output/', input$job_id, '.zip'))), 
                 type = 'success')
    }
  }
)