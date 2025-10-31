#' Copy a sample configuration file
#'
#' This function provides a sample configuration file. The user can modify the parameters as desired
#' @param dir.path The existing folder for copying the file.
#' @export
#'
getConfigFile <- function(dir.path)
{
  if(missing(dir.path))
  {
    stop('No directory specified.',
         call. = FALSE)
  }

  if(!is(dir.path,"character"))
    stop('Argument is not in a correct format. A string parameter is required.', call. = FALSE)
  else
    dir.path <- gsub('/+$', '' , dir.path)

  # check if folder exists
  if (!dir.exists(dir.path))
    stop(sprintf('Directory not found \'%s\'', dir.path), call. = FALSE)

  # check if folder is writable
  if (file.access(dir.path, 2) != 0)
    stop(sprintf('Permission denied to write at \'%s\'', dir.path),
         call. = FALSE)

  # check if config file exists in package
  configFile <-
    system.file("extdata", "inSilicoSeqConfig.ini", package = "SNPannotator")

  if (!file.exists(configFile))
    stop('Config file not found in the package!', call. = FALSE)


  # save file
  # exit if file already exists

  file.name = 'inSilicoSeqConfig.ini'
  file.path = file.path(dir.path , file.name)

  config.file.exists  <- TRUE

  while (config.file.exists) {
    if (file.exists(file.path))
    {
      file.name = sprintf('config_%s.ini', sample(1:100000, 1))
      file.path = file.path(dir.path , file.name)
    } else
      config.file.exists <- FALSE

  }


  saveResult <- file.copy(from = configFile,
                          to = file.path,
                          overwrite = FALSE)

  # double check if file is saved
  if (saveResult &&
      file.exists(file.path))
  {
    message(sprintf('Sample config file saved: \'%s\'!', file.path))
    return(invisible(file.path))
  }
  else
    warning(sprintf('Could not save sample config file at: \'%s\'!', file.path),
            call. = FALSE)


}


loadConfigurationFile = function(path, demo = FALSE)
{
  if(!file.exists(path)) {
    message("Configuration file does not exist \"",path,"\"")
    return(NULL)
  }

  tryCatch({
    configs = read.ini(filepath = path)

    ## GENERAL
    ##
    if(is.null(configs$general$projectName))
      configs$general$projectName='myProject'

    ##
    if(is.null(configs$general$outputFolder))
      configs$general$outputFolder='.'

    configs$general$outputFolder = normalizePath(configs$general$outputFolder, mustWork = FALSE)

    if(!dir.exists(configs$general$outputFolder))
    {
      message('Output directory does not exist.')
      return(NULL)
    }

    ##
    # DEPRECATED
    # if(!is.null(configs$general$geneNames_File))
    # {
    #   configs$general$geneNames_File = normalizePath(configs$general$geneNames_File, mustWork = FALSE)
    #   configs$general$geneNames_File = doesFileExist_or_STOP(configs$general$geneNames_File)
    # }

    if(demo == TRUE)
    {
      configs$general$inputFile = system.file("extdata", "demo_input.txt", package = "SNPannotator")
    }

    if(!is.null(configs$general$inputFile))
    {
      configs$general$inputFile = normalizePath(configs$general$inputFile,mustWork = FALSE)
      configs$general$inputFile = doesFileExist_or_STOP(configs$general$inputFile)
      configs$general$inputFile = load_RSlist(configs$general$inputFile)
    } else {
      message("Input file argument not found in config file.")
      return(NULL)
    }




    configs$general$generateHTML = as.logical(toupper(trimws(configs$general$generateHTML)))
    if(is.na(configs$general$generateHTML))
      configs$general$generateHTML <- TRUE
    #stop('Unknown input found for generateHTML parameter.')

    # build 37 OR 38
    configs$general$build = is_character_or_default(configs$general$build ,
                                                    'build',
                                                    c('grch37','grch38'),
                                                    'grch37',
                                                    ignore_error = FALSE)
    # ##

    configs$general$nearestGene_type = is_character_or_default(configs$general$nearestGene_type,
                                                               'nearestGene_type',
                                                               c("all","pseudogene","lincRNA","protein_coding","antisense",
                                                                 "processed_transcript","snRNA","sense_intronic","miRNA",
                                                                 "misc_RNA","snoRNA","rRNA","3prime_overlapping_ncRNA",
                                                                 "polymorphic_pseudogene","sense_overlapping","IG_V_gene",
                                                                 "IG_C_gene","IG_J_gene","IG_V_pseudogene","TR_C_gene",
                                                                 "TR_J_gene","TR_V_gene","TR_V_pseudogene","IG_C_pseudogene",
                                                                 "TR_D_gene","TR_J_pseudogene","IG_J_pseudogene","IG_D_gene",
                                                                 "transcribed_unprocessed_pseudogene","unprocessed_pseudogene",
                                                                 "processed_pseudogene","transcribed_processed_pseudogene","TEC",
                                                                 "transcribed_unitary_pseudogene","scaRNA","rRNA_pseudogene",
                                                                 "unitary_pseudogene","scRNA","sRNA","ribozyme",
                                                                 "translated_processed_pseudogene","vault_RNA",
                                                                 "IG_pseudogene","artifact"),
                                                               'all',
                                                               ignore_error = FALSE)


    ## ENSEMBL server

    configs$general$server <- switch (configs$general$build,
                                      'grch37' = .SNPannotator$ENSEMBL_API_37,
                                      'grch38' = .SNPannotator$ENSEMBL_API_38)

    if(is.null(configs$general$db))
    {
      message('Database for Ensembl not provided.')
      return(NULL)
    }

    ##
    if(!is.null(configs$general$window_size) &&
       !is.na(as.numeric(configs$general$window_size)) &&
       as.numeric(configs$general$window_size) <= 500 &&
       as.numeric(configs$general$window_size) > 0)
      configs$general$window_size = as.numeric(configs$general$window_size)
    else
      configs$general$window_size = 500

    ##
    if(!is.null(configs$general$r2) &&
       !is.na(as.numeric(configs$general$r2)) &&
       as.numeric(configs$general$r2) <= 1 &&
       as.numeric(configs$general$r2) > 0)
      configs$general$r2 = as.numeric(configs$general$r2)
    else
      configs$general$r2 = 0.8

    # check that r2  is above 0.5
    checkR2(configs$general$r2)
    ##

    configs$general$LDlist = as.logical(toupper(trimws(configs$general$LDlist)))
    NA_or_false(configs$general$LDlist,'Unknown input found for LDlist parameter.')
    ##

    if(!is.null(configs$general$start_index) &&
       !is.na(as.numeric(configs$general$start_index)))
      configs$general$start_index = as.numeric(configs$general$start_index)
    else
      configs$general$start_index = 1

    # DECPRECATED - parallel proccessing
    #
    #     if(!is.null(configs$general$cores) && !is.na(as.numeric(configs$general$cores)))
    #       configs$general$cores = as.numeric(configs$general$cores)
    #     else
    #       configs$general$cores = 1


    ##
    ## GTEX
    ## remove GTEX for now
    # if(!is.null(configs$GTEX$server) && !startsWith(tolower(configs$GTEX$server),'http'))
    #   stop('URL for GTEx server is not correct.')
    #
    # if(!is.null(configs$GTEX$server) && is.null(configs$GTEX$version))
    #   stop('GTEx version not provided.')
    #
    # if(!is.null(configs$GTEX$server) && is.null(configs$GTEX$tissueID))
    #   stop('GTEx tissue_id not provided.')
    #
    #
    # if(!is.null(configs$GTEX$pvalue_threshold) &&
    #    !is.na(as.numeric(configs$GTEX$pvalue_threshold)) &&
    #    as.numeric(configs$GTEX$pvalue_threshold) < 1 &&
    #      as.numeric(configs$GTEX$pvalue_threshold) > 0)
    #   configs$GTEX$pvalue_threshold = as.numeric(configs$GTEX$pvalue_threshold)
    # else
    #   configs$GTEX$pvalue_threshold = 5E-8

    ##
    ## QTL
    if(is.null(configs$QTL$eQTL))
      configs$QTL$eQTL <- FALSE
    else
      configs$QTL$eQTL = as.logical(toupper(trimws(configs$QTL$eQTL)))

    if(is.na(configs$QTL$eQTL))
    {
      message('Unknown input found for eQTL parameter.')
      return(NULL)
    }

    if(is.null(configs$QTL$sQTL))
      configs$QTL$sQTL <- FALSE
    else
      configs$QTL$sQTL = as.logical(toupper(trimws(configs$QTL$sQTL)))

    if(is.na(configs$QTL$sQTL))
    {
      message('Unknown input found for sQTL parameter.')
      return(NULL)
    }

    # gtex dataset id , gtex_v8 or gtex_v10

    if(is.null(configs$QTL$dataset))
      configs$QTL$dataset <- 'gtex_v10'
    else
      configs$QTL$dataset = tolower(trimws(configs$QTL$dataset))

    configs$QTL$dataset <- is_character_or_default(parameter = configs$QTL$dataset,
                                                   parameter_name = 'GTEx portal dataset',
                                                   req_list = c('gtex_v8','gtex_v10'),
                                                   default_value = 'gtex_v10',
                                                   ignore_error = TRUE)
    ##
    configs$QTL$server = NULL

    if(configs$QTL$eQTL || configs$QTL$sQTL)
      configs$QTL$server <- switch (configs$general$build,
                                      'grch37' = .SNPannotator$EBI_eqtl_API,
                                      'grch38' = .SNPannotator$GTEx_eqtl_API)

    if(configs$QTL$sQTL && configs$general$build =='grch37') ## GTEx only accepts variant id in build 38 , based on genomic position
    {
      stop('sQTL information is currently available for build38 only.\n')
    }


    # DEPRECATED - might want to search all variants
    # if(configs$QTL$eQTL && is.null(configs$QTL$eQTL_group))
    # {
    #   message('eQTL_group not provided.')
    #   return(NULL)
    # }


    if(!is.null(configs$QTL$pvalue_threshold) &&
       !is.na(as.numeric(configs$QTL$pvalue_threshold)) &&
       as.numeric(configs$QTL$pvalue_threshold) < 1 &&
       as.numeric(configs$QTL$pvalue_threshold) > 0)
      configs$QTL$pvalue_threshold = as.numeric(configs$QTL$pvalue_threshold)
    else
      configs$QTL$pvalue_threshold = 5E-2

    ## SRTING DB
    if(is.null(configs$STRING_DB$network_image))
      configs$STRING_DB$network_image <- FALSE
    else
      configs$STRING_DB$network_image = as.logical(toupper(trimws(configs$STRING_DB$network_image)))


    if(is.null(configs$STRING_DB$enrichment_image))
      configs$STRING_DB$enrichment_image <- FALSE
    else
      configs$STRING_DB$enrichment_image = as.logical(toupper(trimws(configs$STRING_DB$enrichment_image)))



    if(is.null(configs$STRING_DB$functional_annotation))
      configs$STRING_DB$functional_annotation <- FALSE
    else
      configs$STRING_DB$functional_annotation = as.logical(toupper(trimws(configs$STRING_DB$functional_annotation)))




    if(is.null(configs$STRING_DB$functional_enrichment))
      configs$STRING_DB$functional_enrichment <- FALSE
    else
      configs$STRING_DB$functional_enrichment = as.logical(toupper(trimws(configs$STRING_DB$functional_enrichment)))




    if(is.null(configs$STRING_DB$webLink))
      configs$STRING_DB$webLink <- FALSE
    else
      configs$STRING_DB$webLink = as.logical(toupper(trimws(configs$STRING_DB$webLink)))




    if(is.null(configs$STRING_DB$description))
      configs$STRING_DB$description <- FALSE
    else
      configs$STRING_DB$description = as.logical(toupper(trimws(configs$STRING_DB$description)))

    NA_or_false(configs$STRING_DB$network_image,'Unknown input found for network_image parameter.')
    NA_or_false(configs$STRING_DB$enrichment_image,'Unknown input found for enrichment_image parameter.')
    NA_or_false(configs$STRING_DB$functional_annotation,'Unknown input found for functional_annotation parameter.')
    NA_or_false(configs$STRING_DB$functional_enrichment,'Unknown input found for functional_enrichment parameter.')
    NA_or_false(configs$STRING_DB$webLink,'Unknown input found for webLink parameter.')
    NA_or_false(configs$STRING_DB$description,'Unknown input found for description parameter.')

    # network image type
    configs$STRING_DB$imagetype = is_character_or_default(configs$STRING_DB$imagetype,
                                                          'imagetype',
                                                          c('image','highres_image','svg'),
                                                          'image')

    configs$STRING_DB$network_type = is_character_or_default(configs$STRING_DB$network_type,
                                                             'network_type',
                                                             c('functional','physical'),
                                                             'functional')

    configs$STRING_DB$network_flavor = is_character_or_default(configs$STRING_DB$network_flavor,
                                                               'network_flavor',
                                                               c('evidence', 'confidence', 'actions'),
                                                               'evidence')

    configs$STRING_DB$add_color_nodes = is_numeric_or_default(configs$STRING_DB$add_color_nodes,0)
    configs$STRING_DB$add_white_nodes = is_numeric_or_default(configs$STRING_DB$add_white_nodes,0)
    configs$STRING_DB$required_score = is_numeric_or_default(configs$STRING_DB$required_score,0)
    configs$STRING_DB$hide_node_labels = is_numeric_or_default(configs$STRING_DB$hide_node_labels,0)
    configs$STRING_DB$hide_disconnected_nodes = is_numeric_or_default(configs$STRING_DB$hide_disconnected_nodes,0)
    configs$STRING_DB$show_query_node_labels = is_numeric_or_default(configs$STRING_DB$show_query_node_labels,0)
    configs$STRING_DB$block_structure_pics_in_bubbles = is_numeric_or_default(configs$STRING_DB$block_structure_pics_in_bubbles,0)
    configs$STRING_DB$center_node_labels = is_numeric_or_default(configs$STRING_DB$center_node_labels,0)
    configs$STRING_DB$flat_node_design = is_numeric_or_default(configs$STRING_DB$flat_node_design,0)
    configs$STRING_DB$custom_label_font_size = is_numeric_or_default(configs$STRING_DB$custom_label_font_size,12)
    configs$STRING_DB$limit = is_numeric_or_default(configs$STRING_DB$limit,0)

    if( configs$STRING_DB$custom_label_font_size < 6 )
      configs$STRING_DB$custom_label_font_size = 6
    else if ( configs$STRING_DB$custom_label_font_size > 14)
      configs$STRING_DB$custom_label_font_size = 14


    configs$STRING_DB$server = NULL
    if(configs$STRING_DB$network_image ||
       configs$STRING_DB$enrichment_image ||
       configs$STRING_DB$webLink ||
       configs$STRING_DB$description ||
       configs$STRING_DB$functional_annotation ||
       configs$STRING_DB$functional_enrichment )
      configs$STRING_DB$server = .SNPannotator$STRINGDB_API

    ##
    ## DEPRECATED - GWAS_Catalogue
    #
    # if(!is.null(configs$GWAS_Catalogue$graph_file))
    # {
    #   configs$GWAS_Catalogue$graph_file = normalizePath(configs$GWAS_Catalogue$graph_file, mustWork = FALSE)
    #   configs$GWAS_Catalogue$graph_file = doesFileExist_or_STOP(configs$GWAS_Catalogue$graph_file)
    # }
    #
    #
    # if(!is.null(configs$GWAS_Catalogue$granularity) &&
    #    is.element(configs$GWAS_Catalogue$granularity,c('1','2','3')))
    #   configs$GWAS_Catalogue$granularity = as.numeric(configs$GWAS_Catalogue$granularity)
    # else
    #   configs$GWAS_Catalogue$granularity = 2

    ##
    ## Graphs
    configs$graphs$generate_graph = as.logical(toupper(trimws(configs$graphs$generate_graph)))
    NA_or_false(configs$graphs$generate_graph,'Unknown input found for generate_graph parameter.')

    if(!is.null(configs$graphs$layout) &&
       is.element(configs$graphs$layout,c('layout_in_circle' ,
                                          'layout_nicely',
                                          'layout_as_star',
                                          'layout_on_grid',
                                          'layout_on_sphere',
                                          'layout_randomly',
                                          'layout_with_dh',
                                          'layout_with_dri',
                                          'layout_with_fr',
                                          'layout_with_gem',
                                          'layout_with_graphopt',
                                          'layout_with_kk',
                                          'layout_with_lgl',
                                          'layout_with_mds')))
    configs$graphs$layout = eval(parse(text = configs$graphs$layout))
    else
      configs$graphs$layout = igraph::layout_in_circle



    ## file paths
    # file for saving output object as plain text
    configs$paths$output.tsv.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_output.tsv'))

    configs$paths$output.xlsx.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'.xlsx'))

    # file for saving output object
    configs$paths$output.rds.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'.rds'))

    # file for saving eQTL data from EBI as plain text
    configs$paths$output.eqtl.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_eQTL.tsv'))
    configs$paths$output.eqtl.graph = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_eQTL.png'))
    configs$paths$output.sqtl.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_sQTL.tsv'))
    configs$paths$output.sqtl.graph = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_sQTL.png'))

    configs$paths$output.gwascatalog.graph = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_gwascatalog.png'))

    # file for saving HTML results
    configs$paths$output.html.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'.html'))

    # file for saving not found variants
    configs$paths$output.notFound.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_not_found.txt'))

    # file for GO table
    configs$paths$output.GO.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_GO_terms.tsv'))

    # file for log
    configs$paths$output.log.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'.log'))

    configs$paths$output.string.enrich.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment.tsv'))
    configs$paths$output.string.annot.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_annotation.tsv'))
    configs$paths$output.string.network.image.file = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_network.png'))
    configs$paths$output.string.enrichment.image.file.function = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment_function.png'))
    configs$paths$output.string.enrichment.image.file.process = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment_process.png'))
    configs$paths$output.string.enrichment.image.file.component = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment_component.png'))
    configs$paths$output.string.enrichment.image.file.diseases = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment_diseases.png'))
    configs$paths$output.string.enrichment.image.file.tissues = file.path(configs$general$outputFolder, paste0(configs$general$projectName,'_enrichment_tissues.png'))


    configs$paths$output.log.file = normalizePath(configs$paths$output.log.file, mustWork = FALSE)
    configs$paths$output.notFound.file = normalizePath(configs$paths$output.notFound.file, mustWork = FALSE)

    configs$paths$output.tsv.file = normalizePath(configs$paths$output.tsv.file, mustWork = FALSE)
    configs$paths$output.xlsx.file = normalizePath(configs$paths$output.xlsx.file, mustWork = FALSE)
    configs$paths$output.rds.file = normalizePath(configs$paths$output.rds.file, mustWork = FALSE)

    configs$paths$output.eqtl.file = normalizePath(configs$paths$output.eqtl.file, mustWork = FALSE)
    configs$paths$output.eqtl.graph = normalizePath(configs$paths$output.eqtl.graph, mustWork = FALSE)
    configs$paths$output.sqtl.file = normalizePath(configs$paths$output.sqtl.file, mustWork = FALSE)
    configs$paths$output.sqtl.graph = normalizePath(configs$paths$output.sqtl.graph, mustWork = FALSE)

    configs$paths$output.gwascatalog.graph = normalizePath(configs$paths$output.gwascatalog.graph, mustWork = FALSE)

    configs$paths$output.html.file = normalizePath(configs$paths$output.html.file, mustWork = FALSE)

    configs$paths$output.string.enrich.file = normalizePath(configs$paths$output.string.enrich.file, mustWork = FALSE)
    configs$paths$output.string.annot.file = normalizePath(configs$paths$output.string.annot.file, mustWork = FALSE)

    configs$paths$output.string.network.image.file = normalizePath(configs$paths$output.string.network.image.file, mustWork = FALSE)
    configs$paths$output.string.enrichment.image.file.function = normalizePath(configs$paths$output.string.enrichment.image.file.function, mustWork = FALSE)
    configs$paths$output.string.enrichment.image.file.process = normalizePath(configs$paths$output.string.enrichment.image.file.process, mustWork = FALSE)
    configs$paths$output.string.enrichment.image.file.component = normalizePath(configs$paths$output.string.enrichment.image.file.component, mustWork = FALSE)
    configs$paths$output.string.enrichment.image.file.diseases = normalizePath(configs$paths$output.string.enrichment.image.file.diseases, mustWork = FALSE)
    configs$paths$output.string.enrichment.image.file.tissues = normalizePath(configs$paths$output.string.enrichment.image.file.tissues, mustWork = FALSE)

    # check if file exists
    if(file.exists(configs$paths$output.tsv.file))
    {
      message('Output file already exists!  ', configs$paths$output.tsv.file)
      return(NULL)

      # craete a random file if already exists
      #projectName = sprintf('%s_%s',createRandString(),projectName)
      #outputPath = sprintf("%s/%s.xlsx",outputFolder,projectName)
    }

    return(configs)

  }, error = function(cond) {
    message('Error occured: ')
    message(cond)
    return(NULL)
  },
  warning = function(cond) {
    message('Warning occured: ')
    message(cond)
    return(NULL)
  })
}

load_RSlist <- function(fileName)
{
  tryCatch({
    data = fread(fileName,
                 header=FALSE,
                 strip.white = TRUE,
                 blank.lines.skip = TRUE)
    # data = read.table(file = fileName,
    #            header = FALSE,
    #            strip.white = TRUE,
    #            blank.lines.skip = TRUE,
    #            comment.char = '#',)

    if(!is.data.table(data))
      setDT(data)

    if(nrow(data) == 0)
    {
      message('Input file error. No data found.')
      return(NULL)
    }

    if(ncol(data) > 1)
    {
      message('Input file error. Input file has more than one column.')
      return(NULL)
    }

    dup_count <- data[duplicated(V1),.N]
    if(dup_count > 0)
    {
      message(sprintf('Removed %s duplicated variants from input file.',dup_count))
      data <- data[!duplicated(V1),]
    }


    return(tolower(trimws(data$V1)))

  },
  error = function(x) {
    stop(x$message, call. = FALSE)
  })

}

checkR2 <- function(r2)
{
  if(r2 < 0.5)
    stop('Selected r2 parameter is very low. Consider values >= 0.5.')
}


doesFileExist_or_NULL <- function(filePath)
{

  if(file.exists(filePath))
    return(filePath)
  else
  {
    message(paste0('File not found! (',filePath,')'))
    return(NULL)
  }

}

doesFileExist_or_STOP <- function(filePath)
{

  if(file.exists(filePath))
    return(filePath)
  else
  {
    message(paste0('File not found! (',filePath,')'))
    stop(call. = FALSE)
  }

}

NA_or_false = function(x,x_message)
{
  if(is.na(x))
  {
    message(x_message)
    return(NULL)
  }
}

is_numeric_or_default = function(parameter,default_value)
{

  if(is.null(parameter) || is.na(as.numeric(trimws(parameter))))
    return(default_value)
  else
    return(as.numeric(trimws(parameter)))

}

is_character_or_default = function(parameter,parameter_name,req_list,default_value,ignore_error=TRUE)
{

  if(is.null(parameter) || !is.element(trimws(tolower(parameter)),req_list))
  {
    if(ignore_error == TRUE)
      return(default_value)
    else
      stop(paste('Error in parameter:',parameter_name))
  }
  else
    return(trimws(tolower(parameter)))

}
