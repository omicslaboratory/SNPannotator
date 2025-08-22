## logger setup function
setupLogOptions <- function(log.file.path) {

  if(is.null(log.file.path))
    return(NULL)

  if(file.exists(log.file.path))
    file.remove(log.file.path)

  flog.appender(appender.file(log.file.path), name='SNPannotator_logger')
  flog.threshold(INFO)
}

## the main logger function
print_and_log <- function(a_message,
                          level = "info",
                          display = TRUE,
                          LF=TRUE) {

  if(is.null(a_message))
    return()

  #tryCatch({
    # 1- display in console
    if (.SNPannotator$verbose==TRUE && display == TRUE)
      message(a_message, appendLF = LF)

    # 2- save in file

    # remove new line from the message
    if(a_message != "")
      a_message <- gsub(pattern = '\n',x = a_message,replacement = '')

    if (level == "info") {
      flog.info(msg = a_message, name = "SNPannotator_logger")
    } else if (level == "warning") {
      flog.warn(msg = a_message, name = "SNPannotator_logger")
    } else if (level == "fatal") {
      flog.fatal(msg = a_message, name = "SNPannotator_logger")
      stop( "Pipeline stopped.", call. = FALSE)
    } else {
      message('unknown log level is defined.')
    }
  # },
  # warning = function(x){
  #   message(paste("Warning in logger: " ,x$message))
  # },
  # error = function(x){
  #   message(paste("Error in logger: " ,x$message))
  # })
}
