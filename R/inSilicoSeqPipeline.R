inSilicoSeqPipeline <- function(rsID,
                                server,
                                database,
                                window_size,
                                r2,
                                i,
                                variantCount,
                                addLinkedVars)
{

  rsID <- trimws(tolower(rsID))
  if(!startsWith(rsID,"rs"))
  {
    print_and_log("rsID is not correct.",level='warning')
    return(NULL)
  }
  print_and_log("===========================================", display=FALSE)
  print_and_log(sprintf("\nFetching data for target SNP (#%s/%s): %s ",i,variantCount,rsID))


  # initiate a counter
  pb <- NULL

  varInfo <- getVariantInfo(rsID,server,pb)

  # check if varInfo has data
  if(is.null(varInfo))
    return(NULL)


  #
  if(addLinkedVars)
  {
    VarLDList <- getVariantLDs(rsID,server,database, window_size, r2)
  }
  else
  {
    print_and_log('Proxy variants... skipped.')
    return(varInfo) # return if LDlist is not asked for
  }

  #=========================================
  # if LDlist is wanted and found something
  if(!exists('VarLDList') || is.null(VarLDList))
  {
    varInfo$LDcount <- 0
    varInfo$LDlistFull <- NULL
    print_and_log('Proxy variants: 0')
  }else
  {
    VarLDList[,r2 := as.numeric(r2)]
    VarLDList <- VarLDList[order(-r2),]
    varInfo$LDList <- VarLDList

    varInfo$LDcount <- nrow(varInfo$LDList)
    varInfo$LDlistFull <- list()

    if(.SNPannotator$verbose==TRUE)
    {
      pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = varInfo$LDcount)
      pb$tick(0)
    }


    print_and_log(sprintf('Proxy variants: %s', varInfo$LDcount))
    print_and_log(sprintf('Proxy list: %s', paste(varInfo$LDList$variation2,collapse = ";")), display=FALSE)


    #####
    ##### DEPRECATED ####
    ##### searching variants one by one
    # for(j in seq_len(varInfo$LDcount))
    # {
    #   proxy.rsid <- varInfo$LDList[j,variation2]
    #   proxy.varInfo <- getVariantInfo(proxy.rsid , server, pb)
    #
    #   if(is.null(proxy.varInfo))
    #     print_and_log(paste('Proxy variant not found:',proxy.rsid),level='warning',display=FALSE)
    #   else
    #     varInfo$LDlistFull[[j]] <- proxy.varInfo
    # }

    ##### CURRENT METHOD ####
    ##### searching variants in one step ####
    proxy.IDs <- VarLDList$variation2

    #### method 1 ####
    ## no separation
    # proxylist <- getMultiVariantInfo(proxy.IDs,server)


    #### method 2 ####
    # Split IDs into chunks of 100
    id_chunks <- split_list(proxy.IDs, 150)

    # Process each chunk separately and combine results
    proxylist <- Reduce(c, lapply(seq_along(id_chunks), function(i) {
      print_and_log(sprintf("Processing set %d of %d ...", i, length(id_chunks)))
      getMultiVariantInfo(id_chunks[[i]], server)
    }))


    for(j in seq_len(varInfo$LDcount))
    {
      proxy.rsid <- varInfo$LDList[j,variation2]
      proxy.varInfo <- proxylist[[proxy.rsid]]

      if(is.null(proxy.varInfo))
        print_and_log(paste('Proxy variant not found:',proxy.rsid),level='warning',display=FALSE)
      else
        varInfo$LDlistFull[[j]] <- proxy.varInfo
    }
  }

  #varInfo$LDList <- correctSynonymIDs(varInfo)

  return(varInfo)
}


## separate to chunks of desired size
split_list <- function(x, n) {
  split(x, ceiling(seq_along(x) / n))
}
