#' Run the annotation pipeline on a list of variants from shiny app
#'
#' This function should not be used outside shiny app
#' A list of variants and parameters are received and their information is checked on various API servers.
#'
#' @param config.list List. A list of variants andconfiguration parameters.
#' @return a data table with all variant information is returned.
#' @export
#'
annotate_shiny <- function(config.list) {


  #  set environment variables from config file
  rslist = config.list$rslist
  build = config.list$build
  db = config.list$db
  window_size = config.list$window_size
  r2 = config.list$r2
  LDlist = config.list$LDlist
  nearestGene_type = config.list$nearestGene_type
  ensembl.server = switch(build,
                          'grch37' = .SNPannotator$ENSEMBL_API_37,
                          'grch38' = .SNPannotator$ENSEMBL_API_38)

  .SNPannotator$verbose = FALSE
  #setupLogOptions('snpannotator.log')

  #===========================================================
  g.set <- getGeneFile_v2(build)

  c.set <- getCytobandFile_v2(build)

  #===========================================================

  # initiate objects
  output <- data.table() # final dataset for traits



  ######################
  ##### fetch data #####
  ######################

  for(i in 1:length(rslist))
  {
    rs= rslist[i]

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
        return(NULL)
      }
    )

    if(is.null(varInfo))
    {

      next # skip to the next variant
    }

    #===========================================================
    ## check if this variants is already present in the output
    ## maybe the provided variants are not independent with each other
    checkifOutputIncludesRS(rs,output)

    # the variants is found
    variant_data_list <-  returnVariantDatatable(i, varInfo,db)

    tab <- variant_data_list['variant_table'][[1]]

    # bypass this variant if errors occured during parsing
    if(is.null(tab))
      next



    #===========================================================

    if(!is.null(c.set))
    {

      tab[Cytoband=='',
          Cytoband := find.band(Chr,Pos,c.set),
          by=list(Chr,Pos)]
    }

    if(!is.null(g.set))
    {

      tab[Gene=='',
          c('GeneId','Gene') := find.nearest.gene(Chr,Pos,g.set,nearestGene_type),
          by=list(Chr,Pos)]

      }


    #===========================================================



    ###################################
    # finished analyzing this variant #
    ###################################

    # add current variant info to the rest
    if(!is.null(tab))
    {
      output <- rbind(output,tab,fill=TRUE)
    }
  }

  #===========================================================
  #===========================================================
  #===========================================================

  invisible(output)
}
