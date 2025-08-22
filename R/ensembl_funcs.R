
#' Checks if the service is alive
#'
#' This function test whether the Ensembl server is accessible or not
#'
#' @param server name of the server. "https://rest.ensembl.org" can be used for GRCh38
#' and "https://grch37.rest.ensembl.org" for GRCh37.
#' @return a message is displayed to the user
#'
pingEnsembl <- function(server)
{
  accessible <- TRUE

  tryCatch(
    {
      if(exists('print_and_log'))
        print_and_log("\nPinging Ensembl server ... ",LF = F)
      else
        message("\nPinging Ensembl server ... ",appendLF = F)

      ext <- "/info/ping?"
      #r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
      r <- get_response_from_url(server, ext)
      #stop_for_status(r)
      checkResponseStatusCode_stop(r)

      response <- fromJSON(toJSON(content(r)))

      if(exists('print_and_log'))
      {
        if(response$ping == 1)
        {
          print_and_log("accessible.")
        } else {
          print_and_log("not accessible.", level='fatal')
          accessible <- FALSE
        }
      }
      else
      {
        if(response$ping == 1)
        {
          message("accessible.")
        } else{
          stop("not accessible.")
          accessible <- FALSE
        }

      }

    },
    error = function(err) stop(err$message,call. = FALSE)
  )

  return(accessible)
}

#' List population from human database (1000 Genomes project)
#'
#' This function list the name, description and size of the available populations
#' in 1000 Genomes project database. This database will be used for returning variables in high LD
#' with the target SNP.
#'
#' @param build Genome build. Either 37 or 38. default: 38
#' @return A data table is returned which includes the name, description and size of the available populations
#' in 1000 Genomes project database.
#' @export
#'
EnsemblDatabases <- function(build = 38)
{
  server <- NULL

  if(build == 37)
    server <- "https://grch37.rest.ensembl.org"
    else if (build == 38)
      server <- "https://rest.ensembl.org"

  if(is.null(server))
    message('Wrong build selected, available options are 37, 38.')

  #========================================================

  ext <- "/info/variation/populations/human?"
  #r <- GET(paste0(server, ext), content_type("application/json"))
  #stop_for_status(r)
  r <- get_response_from_url(server, ext)
  checkResponseStatusCode_stop(r)


  response <- fromJSON(toJSON(content(r)))

  setDT(response)

  return(response[grepl(x=name,pattern = '1000'),list(name,size,description)])

}

#' data release available on this REST server.
#'
#' Shows the data releases available on this REST server.
#' May return more than one release (infrequent non-standard Ensembl configuration).
#'
#' @param build Genome build. Either 37 or 38. default: 38
#' @return a message is displayed to the user
#' @export
#'
EnsemblReleases <- function(build = 38)
{

  server <- NULL

  if(build == 37)
    server <- "https://grch37.rest.ensembl.org"
  else if (build == 38)
    server <- "https://rest.ensembl.org"

  if(is.null(server))
    message('Wrong build selected, available options are 37, 38.')

  #====================================================================

  ext <- "/info/data?"
  #r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #stop_for_status(r)
  r <- get_response_from_url(server, ext)
  checkResponseStatusCode_stop(r)


  return( fromJSON(toJSON(content(r))))

}


# this function replaces GET to check if server is available or not
get_response_from_url <- function(server,ext)
{

  tryCatch({
    r <- GET(paste0(server, ext), content_type("application/json"))
    return(r)
  },
  error = function(e) {stop('An error has occured trying to access the server: ',
                            server,
                            '\nmessage: ', e$message,
                            '\nCheck your internet connection and try agian later.',call.=FALSE)})

}

fetch <- function(ext,server)
{
  counter <- 1
  notFound <- TRUE

  while(counter < 6 && notFound)
  {

    r <- GET(trimws(paste0(server, ext)), content_type("application/json"))

    if(checkStatusCode(r))
      notFound <- FALSE
    else
      counter <- counter + 1

  }

  if(!notFound)
    return(fromJSON(toJSON(content(r))))
  else
  {
    return(NULL)
  }
}


checkStatusCode <- function(response)
{
  code <- response$status_code
  #https://github.com/Ensembl/ensembl-rest/wiki/HTTP-Response-Codes

  if(code == 200L)
    return(TRUE)
  else if(code == 400L)
  {
    #cat('Bad Request.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 403L) {
    message('You are submitting far too many requests.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 404L)
  {
    message('Indicates a badly formatted request. Check your URL.', fill=TRUE)
    return(FALSE)
  }
  else if (code== 408L)
  {
    message('Timeout.', fill=TRUE)
    return(FALSE)
  }
  else if(code == 429L)
  {
    message('You have been rate-limited.', fill=TRUE)
    return(FALSE)
  }
  else if (grepl(x=code, pattern='^5.+')){
    message('Server-side error.', fill=TRUE)
    return(FALSE)
  }
  else
  {
    message('unknown error. code = ',code)
    return(FALSE)
  }


}

# this is similar to the above function
# but is used for responses rather than variant information
# will stop the package if not successful
checkResponseStatusCode_stop <- function(response)
{
  code <- response$status_code
  response.url <- response$url

  if(code ==200L)
  {}
  else if(code == 400L)
  {
    stop('Bad Request is sent to server (400).\nCheck this message ',  response.url)
  }
  else if (code== 403L) {
    stop('You are submitting far too many requests to ',response.url)
  }
  else if (code== 404L)
  {
    stop('Badly formatted request was used (404). Check your URL.\n',  response.url )
  }
  else if (code== 408L)
  {
    stop('Timeout occured whe trying to access ',  response.url ,'\nCheck you internet connection or try later.')
  }
  else if(code == 429L)
  {
    stop('You have been rate-limited because of too many requests (429).\nTry again in an hour.')

  }
  else if (code== 503L){
    stop('The service is temporarily unavailable (503) at ', response.url,'\nCheck your internet connection or try later.')
  }
  else
  {
    stop('An unknown error has occured trying to access ', response.url,'\nError code: ', code)
  }
}

checkRemainingLimit <- function(response)
{
  if(r$headers$`x-ratelimit-remaining` < 1)
    stop(sprintf('your requests are limited. Try again in %s seconds.'),r$headers$`x-ratelimit-reset`)
}


#' Query Ensembl for variant information based on genomic position
#'
#' This function retrieves variant information from Ensembl based on the specified genomic position.
#' It takes the chromosome number, start position, and end position as input parameters and searches
#' for variants within this window, using the specified genomic build.
#' If only the start position is provided, the function automatically sets the end position equal
#' to the start position. This is particularly relevant for SNP variants, where the start and
#' end positions are the same. The function returns all variants found within the defined window.
#'
#' @param chromosome Numeric, specifying the chromosome number.
#' @param start_position Numeric, specifying the starting base pair position.
#' @param end_position Numeric, specifying the ending base pair position.
#' @param build Numeric, specifying the genomic build, default value is 38.
#' @param file_path character, path to a file for saving results as Excell spreadsheet.
#'
#' @return A `data.table` containing variant information including:
#' - `id`: variant id in rsID format
#' - `alleles`: variant alleles
#' - `seq_region_name`: chromosome number
#' - `start`: starting base pair
#' - `end`: ending base pair
#'
#' @export
findRSID <- function(chromosome, start_position,end_position=NULL, build = "38",file_path=NULL) {

  if(is.null(chromosome) ||
     is.null(start_position) ||
     !is.numeric(chromosome) ||
     !is.numeric(start_position))
    stop("Chromosome, position not in correct format.")

  if(build =='37')
    server <- .SNPannotator$ENSEMBL_API_37
  else if(build =='38')
    server <- .SNPannotator$ENSEMBL_API_38
  else
    stop("Assmebly version is wrong (should be either 37 or 38).")

  # use start position for end position of not specified OR if lower than start_position
  if(is.null(end_position) || end_position < start_position)
    end_position <- start_position


  output.tbl <- NULL

  region <- paste0(chromosome, ":", start_position, "-", end_position)
  query <- paste0("/overlap/region/homo_sapiens/", region, "?feature=variation")

  query <- paste0(server, query)

  r <- tryCatch({
    GET(query, content_type("application/json"))
  }, error = function(cond) {
    stop(paste('Error occured:',cond$message))
    return(NULL)
  }, warning = function(cond) {
    message(paste('Warning occured:',cond$message))
  })


  if(!is.null(r) && r$status_code == 200)
  {
    results <- fromJSON(content(r, "text", encoding = "UTF-8"))

    if (!is.null(results) &&
        length(results) > 0 &&
        "id" %in% names(results))
    {
      setDT(results)
      output.tbl <- results[,list(id,alleles,seq_region_name, start,end)]
    }
    else
    {
      message("Variant not found.")
      return(NULL)
    }
  }

  # return if no file is defined
  if(is.null(file_path))
    return(output.tbl)


  # save as excel
  dir_path <- dirname(file_path)

  # Check if directory exists; if not, create it
  if (!dir.exists(dir_path)) {
    message("Directory does not exist. File will not be saved.")
    file_path <- NULL
  }

  # Try-catch block for exporting Excel file
  if(!is.null(file_path))
  {
    tryCatch({
      wb <- createWorkbook()
      addWorksheet(wb, "Results")
      writeData(wb, "Results", output.tbl)

      saveWorkbook(wb, file_path, overwrite = TRUE)

      message("File saved successfully: ", file_path)
    }, error = function(e) {
      message("Error saving file: ", e$message)
    })
  }

  return(output.tbl)
}
