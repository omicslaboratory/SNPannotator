merge_data <- function(file_names) {

  obj_list <- lapply(file_names, readRDS)

  # config list
  merged_configs <- stats::setNames(
    lapply(obj_list, function(x) x$config),
    paste0(tools::file_path_sans_ext(basename(file_names)), "_config")
  )

  # input list
  input_lists <- lapply(obj_list, function(x) x$input)
  merged_inputs <- do.call(c, input_lists)


  # main data table
  tables <- lapply(obj_list, function(x) flatten_list_columns(x$mainData))

  index_offset <- 0
  for (i in seq_along(tables)) {
    idx_col <- tables[[i]][, 1]
    if (i > 1) {
      tables[[i]][, 1] <- idx_col + index_offset
    }
    index_offset <- max(tables[[i]][, 1])
  }
  merged_table <- do.call(rbind, tables)

  var_ids <- unique(merged_table[,c(1,3)])
  names(var_ids) <- c('#gSNP','id')

  # GWAScatalog data table

  GWASCatalog_tables <- extract_and_flatten_tables(obj_list, "GWASCatalog")
  merged_GWASCatalog_tables <- data.table::rbindlist(GWASCatalog_tables, fill=TRUE)
  merged_GWASCatalog_tables <- update_list_ids(var_ids, merged_GWASCatalog_tables)

  # eQTL data table

  eqtl_tables <- extract_and_flatten_tables(obj_list, "eQTL")
  merged_eqtl_table <- data.table::rbindlist(eqtl_tables, fill=TRUE)
  merged_eqtl_table <- update_list_ids(var_ids, merged_eqtl_table)

  # sQTL data table

  sqtl_tables <- extract_and_flatten_tables(obj_list, "sQTL")
  merged_sqtl_table <- data.table::rbindlist(sqtl_tables, fill=TRUE)
  merged_sqtl_table <- update_list_ids(var_ids, merged_sqtl_table)

  # make the output list
  merged_obj <- list()
  merged_obj$input <- merged_inputs
  merged_obj$configs <- merged_configs
  merged_obj$mainData <- merged_table
  merged_obj$GWASCatalog <- merged_GWASCatalog_tables
  merged_obj$eQTL <- merged_eqtl_table
  merged_obj$sQTL <- merged_sqtl_table

  return(merged_obj)
}

flatten_list_columns <- function(dt) {
  for (colname in names(dt)) {
    if (is.list(dt[[colname]])) {
      dt[, (colname) := sapply(dt[[colname]], function(x) if(is.null(x)) NA else x)]
    }
  }
  return(dt)
}

update_list_ids <- function(dt1,dt2){
  if(!is.element('#gSNP',names(dt1)) || !is.element('#gSNP',names(dt2)))
    stop('Var ID column not found for merging datasets.')

  if(!is.element('id', names(dt2)) && is.element('rsid',names(dt2)))
    data.table::setnames(dt2,'rsid','id')

  dt2 <- merge(
    dt2[, -'#gSNP', with=FALSE],
    dt1,
    by = "id",
    all.x = TRUE
  )
  data.table::setcolorder(dt2,'#gSNP')
  data.table::setorder(dt2,'#gSNP')
  return(dt2)
}

extract_and_flatten_tables <- function(obj_list, column_name) {
  tables <- lapply(obj_list, function(x) {
    tbl <- x[[column_name]]
    if (!is.null(tbl)) flatten_list_columns(tbl) else NULL
  })
  tables <- Filter(Negate(is.null), tables)
  return(tables)
}

build_string_list <- function(dots)
{

  if (!is.null(dots$required_score)) {
    required_score <- dots$required_score
  } else {
    required_score <- 700
  }


  if (!is.null(dots$limit)) {
    limit <- dots$limit
  } else {
    limit <- 0
  }


  if (!is.null(dots$description)) {
    description <- as.logical(dots$description)
  } else {
    description <- TRUE
  }

  if (!is.null(dots$webLink)) {
    webLink <- as.logical(dots$webLink)
  } else {
    webLink <- TRUE
  }


  if (!is.null(dots$enrichment_image)) {
    enrichment_image <- as.logical(dots$enrichment_image)
  } else {
    enrichment_image <- TRUE
  }

  if (!is.null(dots$functional_enrichment)) {
    functional_enrichment <- as.logical(dots$functional_enrichment)
  } else {
    functional_enrichment <- TRUE
  }


  if (!is.null(dots$network_image)) {
    network_image <- as.logical(dots$network_image)
  } else {
    network_image <- TRUE
  }

  if (!is.null(dots$imagetype)) {
    imagetype <- dots$imagetype
  } else {
    imagetype <- 'image'
  }

  if (!is.null(dots$hide_node_labels)) {
    hide_node_labels <- dots$hide_node_labels
  } else {
    hide_node_labels <- 0
  }

  if (!is.null(dots$hide_disconnected_nodes)) {
    hide_disconnected_nodes <- dots$hide_disconnected_nodes
  } else {
    hide_disconnected_nodes <- 0
  }

  if (!is.null(dots$show_query_node_labels)) {
    show_query_node_labels <- dots$show_query_node_labels
  } else {
    show_query_node_labels <- 0
  }

  if (!is.null(dots$block_structure_pics_in_bubbles)) {
    block_structure_pics_in_bubbles <- dots$block_structure_pics_in_bubbles
  } else {
    block_structure_pics_in_bubbles <- 0
  }

  if (!is.null(dots$center_node_labels)) {
    center_node_labels <- dots$center_node_labels
  } else {
    center_node_labels <- 0
  }

  if (!is.null(dots$custom_label_font_size)) {
    custom_label_font_size <- dots$custom_label_font_size
  } else {
    custom_label_font_size <- 0
  }


  STRING_config_list <- list()
  STRING_config_list[['network_image']] = network_image
  STRING_config_list[['enrichment_image']] = enrichment_image
  STRING_config_list[['functional_enrichment']] = functional_enrichment
  STRING_config_list[['webLink']] = webLink
  STRING_config_list[['description']] = description
  STRING_config_list[['imagetype']] = imagetype
  STRING_config_list[['required_score']] = required_score
  STRING_config_list[['hide_node_labels']] = hide_node_labels
  STRING_config_list[['limit']] = limit
  STRING_config_list[['hide_disconnected_nodes']] = hide_disconnected_nodes
  STRING_config_list[['show_query_node_labels']] = show_query_node_labels
  STRING_config_list[['block_structure_pics_in_bubbles']] = block_structure_pics_in_bubbles
  STRING_config_list[['center_node_labels']] = center_node_labels
  STRING_config_list[['custom_label_font_size']] = custom_label_font_size

  return(STRING_config_list)
}

#' Merge multiple inSilico annotation output objects and export to Excel
#'
#' This function loads multiple result objects produced by SNPannotator package.
#' It merges the relevant tables (e.g., input variants, configuration, ...) across all objects
#' and exports the combined results to a single Excel file.
#' Optionally, it can trigger annotation and enrichment analyses using the STRING database.
#'
#' @param file_names Character vector. Paths to the annotation output \code{.rds} files to merge.
#'        Each file should be a result object containing tables with variant, input, and configuration data.
#' @param project_name Character. The project name used for naming the output Excel file and the log file.
#' @param STRING Logical. If \code{TRUE}, performs annotation and enrichment analysis using STRING for the merged data.
#' @param ... Optional. Additional parameters passed to internal functions.
#' @return Writes an object file, graphs and an Excel file with merged annotation results to disk as a side effect.
#'         The output file contains separate sheets for different data types.
#'
#' @examples
#' \dontrun{
#' merge_annotations(
#'   file_names = c("output1.rds", "output2.rds", "output3.rds"),
#'   project_name = 'project1' )
#' }
#'
#' @export
merge_annotations <- function(file_names, project_name, STRING=FALSE, ...)
{

  # setup the verbose option for this run
  .SNPannotator$verbose <- TRUE

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


  # read the dot arguments in the function
  dots <- list(...)

  if (!is.null(dots$string_server)) {
    string_server <- dots$string_server
  } else {
    string_server <- "https://string-db.org/api"
  }

  if (!is.null(dots$graph_layout)) {
    graph_layout <- eval(parse(text = dots$graph_layout))
  } else {
    graph_layout <- igraph::layout_nicely
  }


  #================================================================
  # set file names
  output.log.file = file.path( paste0(project_name,'.log'))
  output.object.file = file.path( paste0(project_name,'.rds'))
  output.xlsx.file = file.path( paste0(project_name,'.xlsx'))
  output.eqtl.file = file.path( paste0(project_name,'._eqtl.png'))
  output.sqtl.file = file.path( paste0(project_name,'_sqtl.png'))
  output.traits.file = file.path( paste0(project_name,'_traits.png'))
  output.string.enrich.file = file.path( paste0(project_name,'_enrichment.tsv'))
  output.string.network.image.file = file.path( paste0(project_name,'_network.png'))
  output.string.enrichment.image.file.function = file.path( paste0(project_name,'_enrichment_function.png'))
  output.string.enrichment.image.file.process = file.path( paste0(project_name,'_enrichment_process.png'))
  output.string.enrichment.image.file.component = file.path( paste0(project_name,'_enrichment_component.png'))
  output.string.enrichment.image.file.diseases = file.path( paste0(project_name,'_enrichment_diseases.png'))
  output.string.enrichment.image.file.tissues = file.path( paste0(project_name,'_enrichment_tissues.png'))

  # STOP if output file already exists
  if (file.exists(output.xlsx.file))
    stop(
      paste0(
        "The output file '", output.xlsx.file, "' already exists in the current folder. ",
        "Please choose a different project name and try again."
      ),
      call. = FALSE
    )

  # set start time
  start.time <-  proc.time()
  # start logger function
  setupLogOptions(output.log.file)

  # get a merged list of all data from input files
  merged_list <- merge_data(file_names)

  # Save output object
  saveRDS(object = merged_list, file = output.object.file)
  print_and_log('Output object file is saved.')

  # get the required output data
  output <- NULL
  ens.gwascat.output <- NULL
  eqtl.output <- NULL
  sqtl.output <- NULL
  output.index.table <- NULL
  eqtl.cluster.report <- NULL
  sqtl.cluster.report <- NULL
  trait.cluster.report <- NULL

  output <- merged_list$mainData
  ens.gwascat.output <- merged_list$GWASCatalog
  eqtl.output <- merged_list$eQTL
  sqtl.output <- merged_list$sQTL
  output.index.table <- get.input.index.table(output)
  setDT(eqtl.output)
  setDT(sqtl.output)
  data.table::setnames(eqtl.output,'id','rsid')
  data.table::setnames(sqtl.output,'id','rsid')

  # add graphs
  if(is.element('Phenotype', names(output)) && output[Phenotype != '', .N] > 0)
  {
    graph.analysis.output <- make.graph.for.alldataset.clumped(output,
                                                               graph_layout,
                                                               output.traits.file )

    graph.analysis.output.tbl <- graph.analysis.output[['graph.analysis.output']]
    if(!is.null(graph.analysis.output.tbl) & nrow(graph.analysis.output.tbl) > 0 )
      trait.cluster.report <- graph.analysis.output.tbl
  }


  if(nrow(eqtl.output) > 0)
    eqtl.cluster.report <- make.eQTL.graphs.ebi.for.alldataset.clumped(eqtl.output,
                                                                       output.index.table,
                                                                       graph_layout,
                                                                       output.eqtl.file)


  if(nrow(sqtl.output) > 0)
    sqtl.cluster.report <- make.eQTL.graphs.ebi.for.alldataset.clumped(sqtl.output,
                                                                       output.index.table,
                                                                       graph_layout,
                                                                       output.sqtl.file)

  ##### EXCEL #####


  # generate project info table
  project_info_tbl <- NULL
  project_info_tbl <- data.table("Project name" = project_name,
                                 "Output folder" = getwd(),
                                 "Input files" = paste(file_names,collapse = ';'))

  if(STRING == TRUE)
    project_info_tbl = cbind(project_info_tbl , data.table("STRING server URL" = string_server))



  # start writing to excel
  if(nrow(output) > 0)
    writeXLSXfile(output,project_info_tbl,fileName = output.xlsx.file)


  if(!is.null(ens.gwascat.output) && nrow(ens.gwascat.output) > 0)
    appendXLSXfile(ens.gwascat.output, thisSheetName = 'GWASCatalog',fileName = output.xlsx.file, addFirst = FALSE)


  # save excel file - eQTL
  if(!is.null(eqtl.output) && nrow(eqtl.output) > 0)
    appendXLSXfile(eqtl.output,thisSheetName = 'eQTL',fileName = output.xlsx.file, addFirst = FALSE)



  # save excel file - eQTL report
  if(!is.null(eqtl.output) && nrow(eqtl.output) > 0)
    appendXLSX_qtl_report(eqtl.output,'eQTL report',  eqtl.cluster.report , output.xlsx.file)


  # save excel file - sQTL
  if(!is.null(sqtl.output) && nrow(sqtl.output) > 0)
    appendXLSXfile(sqtl.output,thisSheetName = 'sQTL',fileName = output.xlsx.file, addFirst = FALSE)


  # save excel file - sQTL report
  if(!is.null(sqtl.output) && nrow(sqtl.output) > 0)
    appendXLSX_qtl_report(sqtl.output,'sQTL report',  sqtl.cluster.report , output.xlsx.file)




  ##### STRING #####

  if(STRING == FALSE)
    return(invisible(NULL))

  string_gene_list <- NULL
  string.report.list <- NULL
  pngContent_network <- NULL
  pngContent_enrichment_function <- NULL
  pngContent_enrichment_process <- NULL
  pngContent_enrichment_component <- NULL
  pngContent_enrichment_diseases <- NULL
  pngContent_enrichment_tissues <- NULL
  enrich.report <- NULL
  annot.report <- NULL
  trait.cluster.report <- NULL
  eqtl.cluster.report <- NULL
  sqtl.cluster.report <- NULL
  string.description <- NULL
  string.webLink <- NULL
  output.index.table <- NULL # table including gSNP number and rsID. used for plotting
  PingSTRING <- NULL


  # generate string configuration list from function inputs in dots
  STRING_config_list <- build_string_list(dots)



  # extract genes
  string_gene_list <- filter.variants.for.STRING(output)
  output.index.table <- get.input.index.table(output)


  PingSTRING <- PingSTRING(string_server)

  if( !is.null(string_gene_list) && PingSTRING == TRUE)
  {

    string.report.list <- do.string.analysis(string_gene_list,
                                             string_server,
                                             STRING_config_list)

    if(!is.null(string.report.list)){
      pngContent_network <- string.report.list[['pngContent_network']]
      pngContent_enrichment_function <- string.report.list[['pngContent_enrichment_function']]
      pngContent_enrichment_process <- string.report.list[['pngContent_enrichment_process']]
      pngContent_enrichment_component <- string.report.list[['pngContent_enrichment_component']]
      pngContent_enrichment_diseases <- string.report.list[['pngContent_enrichment_diseases']]
      pngContent_enrichment_tissues <- string.report.list[['pngContent_enrichment_tissues']]
      enrich.report <- string.report.list[['enrichment']]
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


    ## append excel file
    # string description
    if(!is.null(string.description) && nrow(string.description) > 0)
      appendXLSXfile(string.description,thisSheetName = 'Description',fileName = output.xlsx.file, addFirst = FALSE)


    # string (func enrich)
    if(!is.null(enrich.report) && nrow(enrich.report) > 0)
      appendXLSXfile(enrich.report,thisSheetName = 'Func. enrichment.',fileName = output.xlsx.file, addFirst = FALSE)



  }

}
