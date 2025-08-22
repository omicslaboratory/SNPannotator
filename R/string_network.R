
convert_logical <- function(x)
{
  x_logical <- suppressWarnings(as.logical(x))

  if (is.na(x_logical) || length(x_logical) != 1) {
    stop(sprintf("%s must be coercible to a single logical value TRUE or FALSE.",deparse(substitute(x))))
  } else{
    return(x_logical)
  }
}

convert_numeric <- function(x)
{
  x_numeric <- suppressWarnings(as.numeric(x))

  if (is.na(x_numeric) || length(x_numeric) != 1) {
    stop(sprintf("%s must be coercible to a single number.",deparse(substitute(x))))
  } else{
    return(x_numeric)
  }
}


#' Analyze STRING DB Interactions and perform functional enrichment
#'
#' This function takes a vector of gene symbols, retrieves their interaction partners
#' from STRING DB, and performs functional enrichment analysis.
#'
#' @param name A character string specifying a unique identifier for this analysis run.
#' @param gene_list A character vector of gene symbols (e.g., HGNC symbols or Ensembl gene IDs).
#' @param required_score Threshold of significance to include an interaction, a number between 0 and 1000.
#' @param limit Limits the number of interaction partners retrieved per protein, a number between 0 and 100.
#' @param ... Additional arguments passed to downstream functions for extended customization.
#'
#' @return set of report files, including images, text and excel files containing functional enrichment analysis results.
#' @export
#'
run_stringdb_annotation <- function(name, gene_list, required_score = 700, limit = 0, ...)
{

  ## check function parameters
  args <- list(...)

  if(is.null(name))
    stop('A name for the job is required.')

  if(is.null(gene_list) || length(gene_list) == 0)
    stop('No genes selected.')

  required_score <- convert_numeric(required_score)
  if ( required_score < 0 || required_score > 1000) {
    stop("required_score must be a single number between 0 and 1000.")
  }

  limit <- convert_numeric(limit)
  if ( limit < 0 || limit > 100) {
    stop("limit must be a single number between 0 and 100.")
  }


  network_image <- ifelse(is.null(args$network_image), TRUE, convert_logical(args$network_image))
  enrichment_image <- ifelse(is.null(args$enrichment_image), TRUE, convert_logical(args$enrichment_image))
  functional_enrichment <- ifelse(is.null(args$functional_enrichment), TRUE, convert_logical(args$functional_enrichment))
  webLink <- ifelse(is.null(args$webLink), TRUE, convert_logical(args$webLink))
  description <- ifelse(is.null(args$description), TRUE, convert_logical(args$description))

  imagetype <- ifelse(!is.null(args$imagetype) && is.element(args$imagetype, c('image','svg','highres_image')) ,
                      args$imagetype,
                      'image')

  ## generate parameters list
  STRING_config_list <- list()
  STRING_config_list[['network_image']] = network_image
  STRING_config_list[['enrichment_image']] = enrichment_image
  STRING_config_list[['functional_enrichment']] = functional_enrichment
  STRING_config_list[['webLink']] = webLink
  STRING_config_list[['description']] = description
  STRING_config_list[['required_score']] = required_score
  STRING_config_list[['limit']] = limit
  STRING_config_list[['imagetype']] = imagetype

  options(error=NULL)
  options(warn = -1)

  ## empty the env on exit
  on.exit({

    if(exists('start.time') && !is.null(start.time))
    {
      if(exists('print_and_log'))# END LOG
      {
        print_and_log("")
        print_and_log(sprintf("Run time: %s",timetaken(start.time)))
      }
      else
        message(sprintf("\nRun time: %s",timetaken(start.time)))
    }

    # remove the environment variable
    if(exists('.SNPannotator'))
      rm(.SNPannotator)

  })


  #=========================

  string.server = .SNPannotator$STRINGDB_API

  # setup the verbose option for this run
  .SNPannotator$verbose <- TRUE

  # add other config parameters
  STRING_config_list[['network_type']] = 'functional'
  STRING_config_list[['network_flavor']] = 'evidence'
  STRING_config_list[['add_color_nodes']] = 0
  STRING_config_list[['add_white_nodes']] = 0
  STRING_config_list[['hide_node_labels']] = 0
  STRING_config_list[['hide_disconnected_nodes']] = 0
  STRING_config_list[['show_query_node_labels']] = 0
  STRING_config_list[['block_structure_pics_in_bubbles']] = 0
  STRING_config_list[['center_node_labels']] = 0
  STRING_config_list[['flat_node_design']] = 0
  STRING_config_list[['custom_label_font_size']] = 12

  STRING_config_list[['functional_annotation']] = FALSE


  # paths for saving files
  output.xlsx.file = paste(name,'xlsx',sep = '.')
  output.html.file = paste(name,'html',sep = '.')
  output.log.file = paste(name,'log',sep = '.')

  output.string.enrich.file = paste(name,'enrichment_report.csv',sep="_")

  output.string.network.image.file = paste(name,'string_network.png',sep="_")
  output.string.enrichment.image.file.function = paste(name,'string_function.png',sep="_")
  output.string.enrichment.image.file.process = paste(name,'string_process.png',sep="_")
  output.string.enrichment.image.file.component = paste(name,'string_component.png',sep="_")
  output.string.enrichment.image.file.diseases = paste(name,'string_diseases.png',sep="_")
  output.string.enrichment.image.file.tissues = paste(name,'string_tissues.png',sep="_")


  # set start time
  start.time <-  proc.time()

  #===========================================================
  # setup a logger
  setupLogOptions(output.log.file)

  #################################
  ####### STRING DB search ########
  #################################


  string.report.list <- NULL
  pngContent_network <- NULL
  pngContent_enrichment_function <- NULL
  pngContent_enrichment_process <- NULL
  pngContent_enrichment_component <- NULL
  pngContent_enrichment_diseases <- NULL
  pngContent_enrichment_tissues <- NULL
  enrich.report <- NULL
  annot.report <- NULL
  string.description <- NULL
  string.webLink <- NULL


  # start the run
  print_and_log("=================================================",display=FALSE)
  print_and_log(sprintf("======== SNPannotator package v%s =========",.SNPannotator$script.version),display=FALSE)
  print_and_log("=================================================",display=FALSE)
  print_and_log("",display=FALSE)
  print_and_log("============= run time parameters ===============",display=FALSE)
  print_and_log(sprintf("Output foler = %s",getwd()),display=FALSE)
  print_and_log(sprintf("Project name = %s",name),display=FALSE)
  print_and_log(sprintf("required_score = %s",required_score),display=FALSE)
  print_and_log(sprintf("limit = %s",limit),display=FALSE)

  print_and_log('STRING DB search started ...\n')


  string.report.list <- do.string.analysis(gene_list,
                                           string.server,
                                           STRING_config_list,
                                           standalone = TRUE)

  if(!is.null(string.report.list)){
    pngContent_network <- string.report.list[['pngContent_network']]
    pngContent_enrichment_function <- string.report.list[['pngContent_enrichment_function']]
    pngContent_enrichment_process <- string.report.list[['pngContent_enrichment_process']]
    pngContent_enrichment_component <- string.report.list[['pngContent_enrichment_component']]
    pngContent_enrichment_diseases <- string.report.list[['pngContent_enrichment_diseases']]
    pngContent_enrichment_tissues <- string.report.list[['pngContent_enrichment_tissues']]
    enrich.report <- string.report.list[['enrichment']]
    annot.report <- string.report.list[['annotation']]
    string.webLink <- string.report.list[['webLink']]
    string.description <- string.report.list[['description']]
  }

  # save network image
  if(!is.null(pngContent_network))
    png::writePNG(image = pngContent_network,target = output.string.network.image.file,dpi = 300)

  if(!is.null(pngContent_enrichment_function))
    png::writePNG(image = pngContent_enrichment_function,target = output.string.enrichment.image.file.function,dpi = 300)

  if(!is.null(pngContent_enrichment_process))
    png::writePNG(image = pngContent_enrichment_process,target = output.string.enrichment.image.file.process,dpi = 300)

  if(!is.null(pngContent_enrichment_component))
    png::writePNG(image = pngContent_enrichment_component,target = output.string.enrichment.image.file.component,dpi = 300)

  if(!is.null(pngContent_enrichment_diseases))
    png::writePNG(image = pngContent_enrichment_diseases,target = output.string.enrichment.image.file.diseases,dpi = 300)

  if(!is.null(pngContent_enrichment_tissues))
    png::writePNG(image = pngContent_enrichment_tissues,target = output.string.enrichment.image.file.tissues,dpi = 300)



  # save excel file - string description
  if(!is.null(string.description) && nrow(string.description) > 0)
    appendXLSXfile_stringdbstandalone(string.description,thisSheetName = 'Description',fileName = output.xlsx.file, addFirst = FALSE)


  # save excel file - string (func enrich)
  if(!is.null(enrich.report) && nrow(enrich.report) > 0)
    appendXLSXfile_stringdbstandalone(enrich.report,thisSheetName = 'Func. enrichment.',fileName = output.xlsx.file, addFirst = FALSE)


}

