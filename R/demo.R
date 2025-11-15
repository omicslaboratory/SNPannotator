#' Demo run of the annotation pipeline
#'
#' This function is a demo of the annotation algorithm.
#'
#' @return A data table containing the variant information for testing is returned.
#' Report files are also saved in the current working directory.
#'
#' @export
#'
demo_annotation <- function() {

  configFile = system.file("extdata", "demo_config.ini", package = "SNPannotator")

  if(is.null(configFile) || !file.exists(configFile))
    stop('demo configuration file is missing. re-install the package.')

  options(error=NULL)
  options(warn = -1)

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

  #===========================================================
  # read and validate configuration file

  config.list = loadConfigurationFile(configFile, demo = TRUE)

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

  # make a table for output excel file (1st sheet)
  project_info_tbl <- NULL
  if(LDlist == TRUE)
    project_info_tbl <- data.table("Project name" = projectName,
                                   "Output folder" = outputFolder,
                                   "Genomic build" = build,
                                   "DB" = db,
                                   "Window Size" = window_size,
                                   "r2" = r2,
                                   "RS list" = paste(rslist,collapse = ';'))
  else
    project_info_tbl <- data.table("Project name" = projectName,
                                   "Output folder" = outputFolder,
                                   "Genomic build" = build,
                                   "RS list" = paste(rslist,collapse = ';'))

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

  # QTL
  qtl.server = config.list$QTL$server
  qtl.p.value.threshold = config.list$QTL$pvalue_threshold
  eqtl.run <- config.list$QTL$eQTL
  sqtl.run <- config.list$QTL$sQTL

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
  output.eqtl.file = config.list$paths$output.eqtl.file
  output.eqtl.graph = config.list$paths$output.eqtl.graph
  output.sqtl.file = config.list$paths$output.sqtl.file
  output.sqtl.graph = config.list$paths$output.sqtl.graph
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

  if(!is.null(qtl.server))
  {
    print_and_log(sprintf("QTL p_value threshold = %s",qtl.p.value.threshold),display=FALSE)
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

  g.set <- getGeneFile_v2(build)


  #===========================================================

  # load cytoband file
  c.set <- getCytobandFile_v2(build)


  #===========================================================

  # initiate objects
  output <- data.table() # final dataset for traits
  eqtl.output <- data.table() # final dataset for eQTL
  sqtl.output <- data.table() # final dataset for sQTL
  ens.clinvar.output <- data.table() # final dataset for clinvar
  ens.gwascat.output <- data.table() # final dataset for gwascat from ensembl
  ens.phenotype.output <- data.table() # final dataset for ALL phenotype data from ensembl

  notFound.list = c()# not found variants
  output_list <- list('config' = config.list) # this list keeps all environment variables, to be saved as one RDS object



  ######################
  ##### ping APIs ######
  ######################
  .SNPannotator$pingEnsembl <- NULL
  .SNPannotator$PingSTRING <- NULL
  .SNPannotator$PingQTL <- NULL

  if(!is.null(ensembl.server))
    .SNPannotator$pingEnsembl <- pingEnsembl(ensembl.server)

  if(!is.null(qtl.server))
    .SNPannotator$PingQTL <- PingQTL(qtl.server,build)

  if(!is.null(string.server))
    .SNPannotator$PingSTRING <- PingSTRING(string.server)


  ######################
  ##### fetch data #####
  ######################
  print_and_log("\nFetching variant information ...")
  ### add header only for the first variant that was found to the output file
  add.header.tsv.file = TRUE
  add.header.eqtl.file = TRUE
  add.header.sqtl.file = TRUE

  for(i in 1:length(rslist))
  {
    rs= rslist[i]
    eqtl.tab <- data.table()
    sqtl.tab <- data.table()
    gwascatalog.tab <- data.table()

    tab <- NULL
    clinvar_tab <- NULL
    gwascat_tab <- NULL
    phenotype_tab <- NULL


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
    variant_data_list <-  returnVariantDatatable(i, varInfo,db)

    tab <- variant_data_list['variant_table'][[1]]
    clinvar_tab <- variant_data_list['clinvar_table'][[1]]
    gwascat_tab <- variant_data_list['gwascat_table'][[1]]
    phenotype_tab <- variant_data_list['phenotype_table'][[1]]

    ## add gSNP number to the first column
    if(!is.null(clinvar_tab) && nrow(clinvar_tab) > 0)
    {
      clinvar_tab$`#gSNP` <- i
      setcolorder(clinvar_tab,'#gSNP')
    }

    if(!is.null(gwascat_tab) && nrow(gwascat_tab) > 0)
    {
      gwascat_tab$`#gSNP` <- i
      setcolorder(gwascat_tab,'#gSNP')
    }

    if(!is.null(phenotype_tab) && nrow(phenotype_tab) > 0)
    {
      phenotype_tab$`#gSNP` <- i
      setcolorder(phenotype_tab,'#gSNP')
    }

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

    if(!is.null(qtl.server) && .SNPannotator$PingQTL == TRUE)
    {
      if(eqtl.run==TRUE)
      {


        temp_tab <- tab[LD >0.8,] # only find eQTL for top SNP and linked snps with LD > 0.8
        temp_tab[,b38_id := sprintf("chr%s_%s_%s_%s_b38",Chr,Pos,Ref_Allele,Alt_Allele)]
        b38_ids <- unique(temp_tab$b38_id)

        ## eqtl from GTEx
        print_and_log('Checking eQTL data ...',LF=FALSE)
        eqtl.tab <- find.qtl.GTEx(server = qtl.server,
                                  qtl='eQTL',
                                  variant_IDs = b38_ids,
                                  pvalue_threshold = qtl.p.value.threshold)

        if(!is.null(eqtl.tab) && nrow(eqtl.tab) > 0)
        {
          eqtl.tab$`#gSNP` <- i
          setcolorder(eqtl.tab,'#gSNP')

          # save current data
          fwrite(eqtl.tab,
                 file = output.eqtl.file ,
                 quote = FALSE,
                 append = TRUE,
                 sep='\t',
                 row.names = FALSE,
                 col.names = add.header.eqtl.file)

          # only add one column name row to the file
          if(add.header.eqtl.file == TRUE)
            add.header.eqtl.file <- FALSE

        } else {
          print_and_log('No eQTL association was found.',display=FALSE)
        }
      }

      if(sqtl.run == TRUE)
      {


        ## sqtl from GTEx
        print_and_log('Checking sQTL data ...',LF=FALSE)
        sqtl.tab <- find.qtl.GTEx(server = qtl.server,
                                  qtl='sQTL',
                                  variant_IDs = b38_ids,
                                  pvalue_threshold = qtl.p.value.threshold)

        if(!is.null(sqtl.tab) && nrow(sqtl.tab) > 0)
        {
          sqtl.tab$`#gSNP` <- i

          setcolorder(sqtl.tab,'#gSNP')

          # save current data
          fwrite(sqtl.tab,
                 file = output.sqtl.file ,
                 quote = FALSE,
                 append = TRUE,
                 sep='\t',
                 row.names = FALSE,
                 col.names = add.header.sqtl.file)

          # only add one column name row to the file
          if(add.header.sqtl.file == TRUE)
            add.header.sqtl.file <- FALSE

        } else {
          print_and_log('No sQTL association was found.',display=FALSE)
        }
      }
    }

    ###################################
    # finished analyzing this variant #
    ###################################

    # add current variant info to the rest
    if(!is.null(tab))
    {
      output <- rbind(output,tab,fill=TRUE)
      output_list[['input']][[rs]] <- varInfo

      if(!is.null(clinvar_tab))
        ens.clinvar.output <- rbind(ens.clinvar.output ,clinvar_tab,fill=TRUE)

      if(!is.null(gwascat_tab))
        ens.gwascat.output <- rbind(ens.gwascat.output ,gwascat_tab,fill=TRUE)

      if(!is.null(phenotype_tab))
        ens.phenotype.output <- rbind(ens.phenotype.output ,phenotype_tab,fill=TRUE)

      # add eQTL terms to the final list
      if(!is.null(eqtl.tab))
        eqtl.output <- rbind(eqtl.output ,eqtl.tab,fill=TRUE)
      # add sQTL terms to the final list
      if(!is.null(sqtl.tab))
        sqtl.output <- rbind(sqtl.output ,sqtl.tab,fill=TRUE)

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
  sqtl.cluster.report <- NULL
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


  if(!is.null(string.server) && !is.null(string_gene_list) && .SNPannotator$PingSTRING == TRUE)
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

    ## association graph
    print_and_log('Plotting association graphs ...',display=FALSE)
    if(is.element('Phenotype', names(output)) &&
       output[Phenotype != '', .N] > 0)
    {
      graph.analysis.output <- make.graph.for.alldataset.clumped(output,
                                                                 graph_layout,
                                                                 output.gwascatalog.graph )

      graph.analysis.output.tbl <- graph.analysis.output[['graph.analysis.output']]
      if(!is.null(graph.analysis.output.tbl) & nrow(graph.analysis.output.tbl) > 0 )
        trait.cluster.report <- graph.analysis.output.tbl
    }

    ## eqtl graph
    if(nrow(eqtl.output) > 0)
    {
      print_and_log('Plotting eQTL graphs ...',display=FALSE)
      setDT(eqtl.output)
      # plot the eqtl graph and return cluster table
      eqtl.cluster.report <- make.eQTL.graphs.ebi.for.alldataset.clumped(eqtl.output,
                                                                         output.index.table,
                                                                         graph_layout,
                                                                         output.eqtl.graph)
    }

    ## sqtl graph
    if(nrow(sqtl.output) > 0)
    {
      print_and_log('Plotting sQTL graphs ...',display=FALSE)
      setDT(sqtl.output)
      # plot the eqtl graph and return cluster table
      sqtl.cluster.report <- make.eQTL.graphs.ebi.for.alldataset.clumped(sqtl.output,
                                                                         output.index.table,
                                                                         graph_layout,
                                                                         output.sqtl.graph)
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
    writeXLSXfile(output,project_info_tbl,fileName = output.xlsx.file)
    # appendXLSXfile(output,thisSheetName = 'All Variants',fileName = output.xlsx.file, addFirst = FALSE)
  }


  #===========================================================

  # save excel file - ENSEMBL ClinVar
  if(!is.null(ens.clinvar.output) && nrow(ens.clinvar.output) > 0)
  {
    appendXLSXfile(ens.clinvar.output, thisSheetName = 'ClinVar',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('ClinVar'= ens.clinvar.output))

  }

  # save excel file - ENSEMBL GWASCatalog
  if(!is.null(ens.gwascat.output) && nrow(ens.gwascat.output) > 0)
  {
    appendXLSXfile(ens.gwascat.output, thisSheetName = 'GWASCatalog',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('GWASCatalog'= ens.gwascat.output))

  }

  # save excel file - ENSEMBL Phenotype
  if(!is.null(ens.phenotype.output) && nrow(ens.phenotype.output) > 0)
  {
    appendXLSXfile(ens.phenotype.output, thisSheetName = 'Phenotypes',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('Phenotypes'= ens.phenotype.output))

  }
  #===========================================================

  # save excel file - eQTL
  if(!is.null(eqtl.output) && nrow(eqtl.output) > 0)
  {
    setcolorder(eqtl.output,'#gSNP')

    appendXLSXfile(eqtl.output,thisSheetName = 'eQTL',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('eQTL'= eqtl.output))

  }


  # save excel file - eQTL report
  if(!is.null(eqtl.output) && nrow(eqtl.output) > 0)
  {
    appendXLSX_qtl_report(eqtl.output,'eQTL report',  eqtl.cluster.report , output.xlsx.file)
  }

  #===========================================================

  # save excel file - sQTL
  if(!is.null(sqtl.output) && nrow(sqtl.output) > 0)
  {
    setcolorder(sqtl.output,'#gSNP')

    appendXLSXfile(sqtl.output,thisSheetName = 'sQTL',fileName = output.xlsx.file, addFirst = FALSE)
    output_list <- c(output_list,list('sQTL'= sqtl.output))

  }


  # save excel file - sQTL report
  if(!is.null(sqtl.output) && nrow(sqtl.output) > 0)
  {
    appendXLSX_qtl_report(sqtl.output,'sQTL report',  sqtl.cluster.report , output.xlsx.file)
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

    # appendXLSXfile(trait.cluster.report,
    #                thisSheetName = 'Associated trait clusters',
    #                fileName = output.xlsx.file,
    #                addFirst = FALSE)

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
