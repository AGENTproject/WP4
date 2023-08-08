fluidRow(
  style = 'padding: 30px;',
  
  textInput(
    inputId = 'job_id',
    label = 'Your Job ID:',
    placeholder = 'e.g., AGENT_T4_JOB_20230508060723',
    width = '400px',
  ),

  checkboxInput(
    inputId = 'result_example', 
    label = 'Get example job ID', 
    value = FALSE,
  ),
  
  div(align = 'center', 
      actionButton('download', 
                   label = 'Check if your results are ready to download', 
                   width = '300px')),
)
