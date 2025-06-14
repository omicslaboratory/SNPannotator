#' Run the annotation pipeline on a list of variants
#'
#' This function receives the path to the configuration file.
#' A list of variants is received and their information is checked on various API servers.
#'
#' @param configurationFilePath Character. The path to the configuration file.
#' @param verbose Logical. Whether to display messages in the console.
#' @return A data table containing all variant information is returned based on the user's selected
#' specifications and parameters. Additionally, numerous other report files in various formats,
#' including text, HTML, Excel, and image, are saved in the output folder.
#' @export
#'
run_annotation <- function(configurationFilePath, verbose = TRUE) {

  options(error=NULL)
  options(warn = -1)

  # setup the verbose option for this run
  .SNPannotator$verbose <- verbose

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

  #===========================================================
  # read and validate configuration file
  config.list = loadConfigurationFile(configurationFilePath)
  if(is.null(config.list))
    stop('Could not validate the configuration file.',call. = FALSE)

  #===========================================================


  #  set environment variables from config file
  rslist = config.list$general$inputFile
  skipHTML = !config.list$general$generateHTML # generateHTML is set in the config file which sould be inverted
  projectName = config.list$general$projectName
  outputFolder = config.list$general$outputFolder
  build = config.list$general$build

  ensembl.server = config.list$general$server
  db = config.list$general$db
  window_size = config.list$general$window_size
  r2 = config.list$general$r2
  LDlist = config.list$general$LDlist
  start_index = config.list$general$start_index

  #DEPRECATED
  #geneNames.file = config.list$general$geneNames_File
  #cores = config.list$general$cores
  nearestGene_type = config.list$general$nearestGene_type

  # STRINGDB
  string.server = config.list$STRING_DB$server
  STRING_config_list <- list()
  STRING_config_list[['network_image']] = config.list$STRING_DB$network_image # T/F
  STRING_config_list[['enrichment_image']] = config.list$STRING_DB$enrichment_image # T/F
  STRING_config_list[['functional_annotation']] = config.list$STRING_DB$functional_annotation # T/F
  STRING_config_list[['functional_enrichment']] = config.list$STRING_DB$functional_enrichment # T/F
  STRING_config_list[['webLink']] = config.list$STRING_DB$webLink # T/F
  STRING_config_list[['description']] = config.list$STRING_DB$description # T/F

  STRING_config_list[['imagetype']] = config.list$STRING_DB$imagetype
  STRING_config_list[['network_type']] = config.list$STRING_DB$network_type
  STRING_config_list[['network_flavor']] = config.list$STRING_DB$network_flavor
  STRING_config_list[['add_color_nodes']] = config.list$STRING_DB$add_color_nodes
  STRING_config_list[['add_white_nodes']] = config.list$STRING_DB$add_white_nodes
  STRING_config_list[['required_score']] = config.list$STRING_DB$required_score
  STRING_config_list[['hide_node_labels']] = config.list$STRING_DB$hide_node_labels
  STRING_config_list[['hide_disconnected_nodes']] = config.list$STRING_DB$hide_disconnected_nodes
  STRING_config_list[['show_query_node_labels']] = config.list$STRING_DB$show_query_node_labels
  STRING_config_list[['block_structure_pics_in_bubbles']] = config.list$STRING_DB$block_structure_pics_in_bubbles
  STRING_config_list[['center_node_labels']] = config.list$STRING_DB$center_node_labels
  STRING_config_list[['flat_node_design']] = config.list$STRING_DB$flat_node_design
  STRING_config_list[['custom_label_font_size']] = config.list$STRING_DB$custom_label_font_size
  STRING_config_list[['limit']] = config.list$STRING_DB$limit

  # eQTL catalogue
  ebi.server = config.list$eQTL_Catalogue$server
  ebi.eQTL_group = config.list$eQTL_Catalogue$eQTL_group
  ebi.p.value.threshold = config.list$eQTL_Catalogue$pvalue_threshold

  # DEPRECATE - GWAS catalogue
  #gwas.catalogue.graph.file = config.list$GWAS_Catalogue$graph_file
  granularity = 1 #config.list$GWAS_Catalogue$granularity

  # graph features
  skip_graph = !config.list$graphs$generate_graph
  graph_layout = config.list$graphs$layout


  # paths
  output.xlsx.file = config.list$paths$output.xlsx.file
  output.rds.file = config.list$paths$output.rds.file
  output.tsv.file  = config.list$paths$output.tsv.file
  output.ebi.file = config.list$paths$output.ebi.file
  output.ebi.graph = config.list$paths$output.ebi.graph
  output.gwascatalog.graph = config.list$paths$output.gwascatalog.graph

  output.html.file = config.list$paths$output.html.file
  output.notFound.file = config.list$paths$output.notFound.file
  output.log.file = config.list$paths$output.log.file

  output.string.enrich.file = config.list$paths$output.string.enrich.file
  output.string.annot.file = config.list$paths$output.string.annot.file

  output.string.network.image.file = config.list$paths$output.string.network.image.file
  output.string.enrichment.image.file.function = config.list$paths$output.string.enrichment.image.file.function
  output.string.enrichment.image.file.process = config.list$paths$output.string.enrichment.image.file.process
  output.string.enrichment.image.file.component = config.list$paths$output.string.enrichment.image.file.component
  output.string.enrichment.image.file.diseases = config.list$paths$output.string.enrichment.image.file.diseases
  output.string.enrichment.image.file.tissues = config.list$paths$output.string.enrichment.image.file.tissues

  ####
  stopifnot('RS list must be a vector.' = is.vector(rslist),
            'Range for r2 is 0-1' = (r2 >0 & r2<=1),
            'Range for window size is 100-500' = (window_size > 100 & window_size<=500)
  )

  #===========================================================

  # set start time
  start.time <-  proc.time()

  #===========================================================
  # setup a logger
  setupLogOptions(output.log.file)


  #===========================================================
  # write run parameters
  print_and_log("=================================================",display=FALSE)
  print_and_log(sprintf("======== SNPannotator package v%s =========",.SNPannotator$script.version),display=FALSE)
  print_and_log("=================================================",display=FALSE)
  print_and_log("",display=FALSE)
  print_and_log("============= run time parameters ===============",display=FALSE)
  print_and_log(sprintf("Output foler = %s",outputFolder),display=FALSE)
  print_and_log(sprintf("Project name = %s",projectName),display=FALSE)
  print_and_log(sprintf("LDlist = %s",as.character(LDlist)),display=FALSE)
  print_and_log(sprintf("Build = %s",build),display=FALSE)
  print_and_log(sprintf("Population = %s",db),display=FALSE)
  print_and_log(sprintf("R2 threshold = %s",r2),display=FALSE)
  print_and_log(sprintf("Window = %s",window_size),display=FALSE)

  if(!is.null(ebi.server))
  {
    if(!is.null(ebi.eQTL_group))
      print_and_log(sprintf("eQTL group = %s",ebi.eQTL_group),display=FALSE)
    print_and_log(sprintf("eQTL p_value threshold = %s",ebi.p.value.threshold),display=FALSE)
  }


  print_and_log("",display=FALSE)
  print_and_log(paste("started at :", date()), display=FALSE)
  print_and_log("===================================================",display=FALSE)
  print_and_log("",display=FALSE)

  #===========================================================

  # check if system has the number of selected cores
  # check.core.count(cores)


  #===========================================================
  # check FGRepo version package
  #check_package_version('FGRepo')
  log.package.version()

  #===========================================================

  #===========================================================

  # load gene name file - to add gene information for intergenic variations
  # based on grch37 or 38

  g.set <- loadGeneNameData(build)


  #===========================================================

  # load cytoband file
  c.set <- loadCytobandData(build)

  #===========================================================

  ## load the GWAS catalog graph

  traits.graph=loadGWAScatGraph()

  #===========================================================

  # initiate objects
  output <- data.table() # final dataset for traits
  ebi.output <- data.table() # final dataset for eQTL

  notFound.list = c()# not found variants
  output_list <- list('config' = config.list) # this list keeps all environment variables, to be saved as one RDS object



  ######################
  ##### ping APIs ######
  ######################

  if(!is.null(ensembl.server))
    pingEnsembl(ensembl.server)

  if(!is.null(ebi.server))
    PingEBI(ebi.server)

  if(!is.null(string.server))
    PingSTRING(string.server)


  ######################
  ##### fetch data #####
  ######################
  print_and_log("\nFetching variant information ...")
  ### add header only for the first variant that was found to the output file
  add.header.tsv.file = TRUE
  add.header.ebi.file = TRUE

  for(i in 1:length(rslist))
  {
    rs= rslist[i]
    ebi.tab <- data.table()
    gwascatalog.tab <- data.table()


    #===========================================================

    # get variant info + variants in high LD
    varInfo <- tryCatch(
      {
        inSilicoSeqPipeline(rs,
                            ensembl.server,
                            db ,
                            window_size ,
                            r2 ,
                            i,
                            length(rslist),
                            LDlist)
      },
      error = function(err) {
        print_and_log(paste0("Error in fetching variant information: ", err$message), level='warning')
        return(NULL)
      }
    )

    if(is.null(varInfo))
    {
      print_and_log("Variant not found.", level='warning')
      notFound.list = c(notFound.list, rs) # add the variant to NOT_FOUND list
      next # skip to the next variant
    }

    #===========================================================
    ## check if this variants is already present in the output
    ## maybe the provided variants are not independent with each other
    checkifOutputIncludesRS(rs,output)


    # TO BE DONE: rsid is checked vs synonym id
    # if(!is.null(varInfo$LDList))
    #   varInfo$LDList <- checkReturnedVariantData(varInfo)


    # the variants is found
    tab <-  returnVariantDatatable(i, varInfo,db)

    # bypass this variant if errors occured during parsing
    if(is.null(tab))
      next



    #===========================================================

    # fill missing genes
    if(!is.null(c.set))
    {
      print_and_log("Checking cytoband information ...", LF = FALSE)

      tab[Cytoband=='',
          Cytoband := find.band(Chr,Pos,c.set),
          by=list(Chr,Pos)]
      print_and_log('done.')
    }

    # fill missing genes
    if(!is.null(g.set))
    {
      print_and_log("Checking gene information ...", LF = FALSE)

      tab[Gene=='',
          c('GeneId','Gene') := find.nearest.gene(Chr,Pos,g.set,nearestGene_type),
          by=list(Chr,Pos)]


      print_and_log('done.')
    }

    # find nearest gene
    # deprecated - only the nearest is reported now
    # tab[, Nearest_gene := mapply(keep_nearest_gene_in_row,Gene,GeneId)]


    #===========================================================

    ##################################
    ### GWAS catalog associations  ###
    ##################################
    print_and_log('Checking Phenotype associations ...',LF=FALSE)
    if(!is.null(traits.graph))
    {

      gwascatalog.tab <- find.associations.offline(traits.graph,tab)

      print_and_log('done.')

      if(!is.null(gwascatalog.tab) && nrow(gwascatalog.tab) > 0 )
      {
        gwascatalog.tab[,Phenotype:=paste(Phenotype,collapse = ';'),by=SNP]
        gwascatalog.tab = gwascatalog.tab[!duplicated(SNP),]

        old.name.order <- names(tab)
        tab <- merge(x=tab,
                     y=gwascatalog.tab,
                     by.x = 'Linked_SNP',
                     by.y = 'SNP',
                     all.x = TRUE,
                     sort = FALSE)

        setcolorder(tab,old.name.order)

      } else {
        # add an empty column for consistency
        tab[, Phenotype := ""]
        print_and_log('No phenotype association found.',display=FALSE)

      }
    }else
    {
      print_and_log("skipped.")
      tab[, Phenotype := ""]
    }

    #===========================================================

    #############
    ### save  ###
    #############

    # add a sheet for each variant
    # appendXLSXfile(tab,sheetName = rs,fileName = outputPath)

    # save current data
    if(!is.null(tab))
    {

      fwrite(tab,
             file = output.tsv.file ,
             quote = FALSE,
             append = TRUE,
             sep='\t',
             row.names = FALSE,
             col.names = add.header.tsv.file)

      # only add one column name row to the file
      if(add.header.tsv.file == TRUE)
        add.header.tsv.file <- FALSE

    }

    #===========================================================

    ############
    ### EBI  ###
    ############

    print_and_log('Checking eQTL data ...', LF = FALSE)
    if(!is.null(ebi.server))
    {
      ebi.tab <- find.eqtl.ebi(ebi.server,
                               tab[tab$LD >0.8,]$Linked_SNP, # only find eQTL for top SNP and linked snps with LD > 0.8
                               ebi.eQTL_group,
                               ebi.p.value.threshold)

      if(!is.null(ebi.tab) && nrow(ebi.tab) > 0)
      {
        ebi.tab <- add.gene.name(ebi.tab , g.set)
        #ebi.tab <- clean.ebi.output(output,ebi.tab)
        ebi.tab$`#gSNP` <- i

        setcolorder(ebi.tab,'#gSNP')

        # save current data
        fwrite(ebi.tab,
               file = output.ebi.file ,
               quote = FALSE,
               append = TRUE,
               sep='\t',
               row.names = FALSE,
               col.names = add.header.ebi.file)

        # only add one column name row to the file
        if(add.header.ebi.file == TRUE)
          add.header.ebi.file <- FALSE
      } else {
        print_and_log('No eQTL association was found.',display=FALSE)
      }

      print_and_log('Checking eQTL data ... done.')

    }else
    {
      print_and_log('skipped.')
    }

    #===========================================================


    ###################################
    # finished analyzing this variant #
    ###################################

    # add current variant info to the rest
    if(!is.null(tab))
    {
      output <- rbind(output,tab,fill=TRUE)
      output_list[['input']][[rs]] <- varInfo

      # add eQTL terms to the final list
      ebi.output <- rbind(ebi.output ,ebi.tab,fill=TRUE)

    }
  }

  #===========================================================
  #===========================================================
  #===========================================================

  #################################
  ####### STRING DB search ########
  #################################

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
  string.description <- NULL
  string.webLink <- NULL
  output.index.table <- NULL # table including gSNP number and rsID. used for plotting


  if(nrow(output) > 0)
  {
    print_and_log("===========================================", display=FALSE)
    print_and_log("")
    #print_and_log("Checking functional annotations ...")

    # gene names separated by %0d
    string_gene_list <- filter.variants.for.STRING(output)

    output.index.table <- get.input.index.table(output)
  }


  if(!is.null(string.server) && !is.null(string_gene_list))
  {

    string.report.list <- do.string.analysis(string_gene_list,
                                             string.server,
                                             STRING_config_list)

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



    # save enrichment analysis report
    if(!is.null(enrich.report) && nrow(enrich.report) > 0)
    {
      # save current data
      fwrite(enrich.report,
             file = output.string.enrich.file ,
             quote = FALSE,
             append = TRUE,
             sep='\t',
             row.names = FALSE)

    }

    # save annotation analysis report
    if(!is.null(annot.report) && nrow(annot.report) > 0)
    {
      # save current data
      fwrite(annot.report,
             file = output.string.annot.file ,
             quote = FALSE,
             append = TRUE,
             sep='\t',
             row.names = FALSE)
    }
  }else{
    print_and_log('Checking functional annotations ... skipped.')
  }

  #===========================================================
  #===========================================================
  #===========================================================

  ######################
  ####### igraphs ######
  ######################

  print_and_log('Plotting ...',LF=FALSE)
  if(!skip_graph)
  {

    print_and_log('Plotting association graphs ...',display=FALSE)
    if(!is.null(traits.graph) &&
       is.element('Phenotype', names(output)) &&
       output[Phenotype != '', .N] > 0)
    {
      graph.analysis.output <- make.graph.for.alldataset.clumped(traits.graph,
                                                                 output,
                                                                 graph_layout,
                                                                 granularity,
                                                                 output.gwascatalog.graph )

      graph.analysis.output.tbl <- graph.analysis.output[['graph.analysis.output']]
      if(!is.null(graph.analysis.output.tbl) & nrow(graph.analysis.output.tbl) > 0 )
        trait.cluster.report <- graph.analysis.output.tbl
    }


    if(nrow(ebi.output) > 0)
    {
      print_and_log('Plotting eQTL graphs ...',display=FALSE)
      setDT(ebi.output)
      # plot the eqtl graph and return cluster table
      eqtl.cluster.report <- make.eQTL.graphs.ebi.for.alldataset.clumped(ebi.output,
                                                                         output.index.table,
                                                                         graph_layout,
                                                                         output.ebi.graph)
    }

    print_and_log('done.')
  }else
  {
    print_and_log('skipped.')
  }

  #################################
  ####### save final reports ######
  #################################

  if(nrow(output) > 0)
  {
    print_and_log("===========================================", display=FALSE)
    print_and_log("")
    print_and_log("Saving final reports ...", LF=FALSE)
    # output=setorderv(output,c("#gSNP","LD"),c(1,-1))
    output_list <- c(output_list,list('mainData' = output))

    # save excel file - main
    appendXLSXfile(output,thisSheetName = 'All Variants',fileName = output.xlsx.file, addFirst = FALSE)
  }




  #===========================================================

  # save excel file - eQTL EBI
  if(!is.null(ebi.output) && nrow(ebi.output) > 0)
  {
    appendXLSXfile(ebi.output,thisSheetName = 'eQTL',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('eQTL'= ebi.output))

  }


  # save excel file - eQTL EBI report
  if(!is.null(ebi.output) && nrow(ebi.output) > 0)
  {
    appendXLSX_eqtl_report(ebi.output,  eqtl.cluster.report , output.xlsx.file)
  }

  #===========================================================

  # save excel file - string description
  if(!is.null(string.description) && nrow(string.description) > 0)
  {
    appendXLSXfile(string.description,thisSheetName = 'Description',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('String_description'= string.description))
  }

  # save excel file - string (func annot)
  if(!is.null(annot.report) && nrow(annot.report) > 0)
  {
    appendXLSXfile(annot.report,thisSheetName = 'Func. annot.',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('Func_annot'= annot.report))

  }

  # save excel file - string (func enrich)
  if(!is.null(enrich.report) && nrow(enrich.report) > 0)
  {
    appendXLSXfile(enrich.report,thisSheetName = 'Func. enrichment.',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('Func_enrich'= enrich.report))
  }


  # save excel file - string weblibk
  if(!is.null(string.webLink))
  {
    output_list <- c(output_list,list('STRING_webLink'= string.webLink))
  }

  #===========================================================

  if(length(notFound.list) > 0 )
  {
    # save variants that are not found
    write.table(notFound.list,
                file= output.notFound.file,
                col.names = FALSE,
                row.names = FALSE,
                quote=FALSE)

    appendXLSXfile(notFound.list,thisSheetName = 'Not found',fileName = output.xlsx.file, addFirst = FALSE)

    output_list <- c(output_list,list('notFound' = notFound.list))

  }


  #===========================================================

  if(!is.null(trait.cluster.report) )
  {
    # write.table(trait.cluster.report,
    #             file= trait.cluster.file,
    #             col.names = FALSE,
    #             row.names = FALSE,
    #             quote=FALSE)

    appendXLSXfile(trait.cluster.report,
                   thisSheetName = 'Associated trait clusters',
                   fileName = output.xlsx.file,
                   addFirst = FALSE)

    output_list <- c(output_list,list('traitClusters' = trait.cluster.report))

  }
  #===========================================================

  # save RDS
  saveRDS(output_list,file=output.rds.file)

  #===========================================================

  # save html file
  if(skipHTML)
  {
    print_and_log('HTML reports are skipped.', display=FALSE)
  } else if(rmarkdown::pandoc_available())
  {
    if(nrow(output) > 0 && LDlist== TRUE)
    {

      if(generate.report.file(output_list,outputFolder, output.html.file))
        print_and_log(sprintf('File saved: %s',output.html.file), display=FALSE)
    }
  }else {
    print_and_log('Pandoc is not available. Skipping HTML report.',level='warning', display=FALSE)
  }

  #===========================================================

  #################
  ### finished! ###
  #################

  print_and_log("done.")
  invisible(output)
}
