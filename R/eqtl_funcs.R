
################
### Ensembl ####
################

PingEBI <- function(server)
{
  print_and_log("Pinging EBI server ... ",LF = FALSE)

  accessible <- TRUE

  tryCatch(
    {
      r <- GET(sprintf('%s/datasets/QTD000600',server), content_type("application/json"))

      if (!is.null(r) && r$status_code == 200)
      {
        print_and_log('accessible.')
      }else
      {
        print_and_log('not accessible.', level='fatal')
        accessible <- FALSE
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in EBI ping.',cond$message), level='fatal')
      accessible <- FALSE
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in EBI ping.',cond$message), level='warning')
      accessible <- FALSE
    }
  )

  return(accessible)
}

get.EBI.genes <- function(server,rsid,eqtl_group,ebi.p.value.threshold)
{
  r <- GET(sprintf('%s/associations/%s?qtl_group=%s&p_upper=%s&size=100',server,rsid,eqtl_group, ebi.p.value.threshold),
           content_type("application/json"))

  out <- data.table()

  if(r$status_code != 200)
    return(out)

  out=data.table(t(sapply(content(r)$`_embedded`$associations,
                          function(x) return(c(x$rsid,x$gene_id,x$pvalue,x$qtl_group)))))

  if(!is.null(out) && nrow(out) > 0 && ncol(out) ==4)
  {
    names(out) <- c('rsId','Gene','pValue','eQTL_group')
    return(out)
  }
  else
    return(data.table())
}

get.EBI.genes.all.tissues <- function(server,rsid)
{
  r <- GET(sprintf('%s/associations?rsid=%s&size=10',server,rsid),
           content_type("application/json"))

  out <- data.table()

  if(r$status_code != 200)
    return(out)

  out=rbindlist(content(r), fill=TRUE)
  setDT(out)

  if(!is.null(out) && nrow(out) > 0 )
  {
    #names(out) <- c('rsId','Ref','Alt','MAF','GeneId','eQTL_group','Beta','SE','Median_tpm','pValue','Study')
    out <- out[,list(rsid,alt,gene_id,beta,pvalue,se)]
    out <- out[order(rsid,beta,pvalue),]
    return(out)
  }
  else
    return(data.table())
}

find.eqtl.ebi <- function(server,
                          snpIDs,
                          p.value.threshold)
{

  out=data.table()

  tryCatch(
    {

      pb <- NULL

      if(.SNPannotator$verbose==TRUE)
      {
        pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(snpIDs))
        pb$tick(0)
      }

      #===================
      for (snpID in snpIDs)
      {
        # if(is.null(eqtl_group))
        out = rbind(out,
                    get.EBI.genes.all.tissues(server,
                                              snpID),
                    fill=TRUE)
        # else
        #   out = rbind(out,
        #               get.EBI.genes(server,
        #                             snpID,
        #                             eqtl_group,
        #                             ebi.p.value.threshold),
        #               fill=TRUE)

        if(!is.null(pb))
          pb$tick(1)
      }

      setDT(out)

      # remove unneccessary columns
      if(is.element('Ref', names(out)))
        out$Ref <- NULL

      if(is.element('MAF', names(out)))
        out$MAF <- NULL
      #

      # log row and column count of table
      print_and_log(sprintf('eQTL data fetched with %s rows and %s columns.',nrow(out),ncol(out)),display=FALSE)

      out[pvalue < p.value.threshold,]
      # log row and column count of table after filtering
      print_and_log(sprintf('eQTL data fetched with %s rows and %s columns after filtering (P < %s).',
                            nrow(out),ncol(out),p.value.threshold),display=FALSE)

      return(out)

    }, error = function(cond) {
      print_and_log(paste('Error occured: ',cond$message), level='fatal')
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured: ',cond$message), level='warning')
      return(NULL)
    }
  )
}

clean.ebi.output <- function(output,ebi.output)
{
  ## add gSNP number to the rows
  setDT(ebi.output)
  setkey(ebi.output,rsId)

  output2 <- output[,list(`#gSNP`,Linked_SNP)]
  setkey(output2,Linked_SNP)

  ebi.output2 <- output2[ebi.output,]


  #

  return(ebi.output2)
}


add.gene.name <- function(data, gdata)
{
  if(!is.data.table(gdata))
    setDT(gdata)

  data <- merge(x=data,y=gdata[,list(id,name,type)],by.x='gene_id',by.y='id',all.x=TRUE,sort=FALSE)
  # setcolorder(data,c('rsId','Ref','Alt','MAF','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'))

  # if(all(is.element(c('rsId','Alt','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'), names(data))))
  #   setcolorder(data,c('rsId','Alt','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'))

  if(is.element('name',names(data)))
    setnames(data,'name','Gene')

  if(is.element('type',names(data)))
    setnames(data,'type','Type')

  data[is.na(Gene) | Gene =='', Gene := gene_id]
  setcolorder(data,c('rsid','alt','gene_id','Gene','Type','beta','se','pvalue'))

  return(data)
}

################
### GTEx ####
################


PingGTEx <- function(server)
{
  print_and_log("Pinging GTEx server ... ",LF = FALSE)

  accessible <- TRUE
  tryCatch(
    {
      r <- GET(server, content_type("application/json"))

      if (!is.null(r) && r$status_code == 200)
      {
        print_and_log('accessible.')
      }else
      {
        print_and_log('not accessible.', level='fatal')
        accessible <- FALSE
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in GTEx ping.',cond$message), level='fatal')
      accessible <- FALSE
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in GTEx ping.',cond$message), level='warning')
      accessible <- FALSE
    }
  )

  return(accessible)
}


get.GTEx.tissues <- function(datasetId)
{
  ext="dataset/tissueSiteDetail"
  server= .SNPannotator$GTEx_eqtl_API
  output.tbl <- NULL

  r <- tryCatch({
    query=sprintf("%s/%s?format=json&datasetId=%s",server,ext,datasetId)

    GET(query, content_type("application/json"))

  }, error = function(cond) {
    print_and_log(paste('Error occured: ',cond$message), level='fatal')
    return(NULL)
  }, warning = function(cond) {
    print_and_log(paste('Warning occured: ',cond$message), level='warning')
    return(NULL)
  })

  if(!is.null(r) && r$status_code == 200)
  {
    output.tbl <- content(r,as = "parsed", simplifyVector = TRUE)
    output.tbl <- as.data.table(output.tbl$data)
  }

  return(output.tbl)
}

find.qtl.GTEx <- function(server,datasetId=NULL,tissueSiteDetailId=NULL,qtl,variant_IDs,pvalue_threshold){

  #print_and_log(sprintf('Getting %s data from GTEx ...',qtl), LF=FALSE)

  ext <- switch(qtl,
                'eQTL' = "association/singleTissueEqtl",
                'sQTL' = "association/singleTissueSqtl",
                'ieQTL' = "association/singleTissueIEqtl",
                'isTL' = "association/singleTissueISqtl",
  )

  output.tbl <- NULL

  query <- sprintf("%s/%s?format=json",server,ext)

  if(!is.null(tissueSiteDetailId))
    query <- sprintf("%s&tissueSiteDetailId=%s",query,tissueSiteDetailId)

  if(!is.null(datasetId))
    query <- sprintf("%s&datasetId=%s",query,datasetId)

  variant_IDs <- paste0("variantId=", variant_IDs, collapse = "&")
  query <- sprintf("%s&%s",query,variant_IDs)


  r <- tryCatch({
    r <- GET(query, content_type("application/json"))
    print_and_log('done.')
    r
  }, error = function(cond) {
    print_and_log(paste('Error occured: ',cond$message), level='fatal')
    return(NULL)
  }, warning = function(cond) {
    print_and_log(paste('Warning occured: ',cond$message), level='warning')
    return(NULL)
  })

  if(!is.null(r) && r$status_code == 200)
  {
    output.tbl <- content(r,as = "parsed", simplifyVector = TRUE)
    output.tbl <- as.data.table(output.tbl$data)

    # log row and column count of table before filtering
    print_and_log(sprintf('%s data fetched with %s rows and %s columns.',qtl,nrow(output.tbl),ncol(output.tbl)),display=FALSE)

    desired_cols <- intersect(c("snpId" ,"geneSymbol","gencodeId","tissueSiteDetailId", "nes","pValue"), names(output.tbl))

    if("pValue" %in% names(output.tbl))
      output.tbl <- output.tbl[pValue < pvalue_threshold, ..desired_cols , drop=F]
    else
      output.tbl <- output.tbl[, ..desired_cols , drop=F]

    # log row and column count of table after filtering
    print_and_log(sprintf('%s data fetched with %s rows and %s columns after filtering (P < %s).',
                          qtl,nrow(output.tbl),ncol(output.tbl),pvalue_threshold),display=FALSE)

    if("snpId" %in% names(output.tbl))
      setnames(output.tbl, "snpId", "rsid")

    if("geneSymbol" %in% names(output.tbl))
      setnames(output.tbl, "geneSymbol", "Gene")

    if("geneSymbol" %in% names(output.tbl))
      setnames(output.tbl, "gencodeId", "gene_id")

    if("tissueSiteDetailId" %in% names(output.tbl))
      setnames(output.tbl, "tissueSiteDetailId", "Tissue")

    if("nes" %in% names(output.tbl))
      setnames(output.tbl, "nes", "beta")

    if("pValue" %in% names(output.tbl))
      setnames(output.tbl, "pValue", "pvalue")

    if(all(is.element(c("rsid","pvalue","beta"),names(output.tbl))))
      output.tbl <- output.tbl[order(rsid,beta,pvalue),]
  }

  return(output.tbl)
}

# find variant info (rsid, varid) from within the package
find.id.GTEx <- function(server,ID){

  ext="dataset/variant"
  output.tbl <- NULL
  output <- NULL
  query <- NULL

  if(startsWith(x = tolower(ID),prefix = 'rs'))
  {
    print_and_log('Converting rsID to variantID from GTEx ...', LF=FALSE)
    query <- sprintf("%s/%s?format=json&snpId=%s", server,ext,tolower(ID))
  } else if(startsWith(x = toupper(ID),prefix = 'CHR')) {
    print_and_log('Converting variantID to rsID from GTEx ...', LF=FALSE)
    query <- sprintf("%s/%s?format=json&variantId=%s", server,ext,snpID)
  } else {
    print_and_log('ID not in correct format.',level = 'fatal')
  }


  r <- tryCatch({
    GET(query, content_type("application/json"))
    print_and_log('done.')
  }, error = function(cond) {
    print_and_log(paste('Error occured: ',cond$message), level='fatal')
    return(NULL)
  }, warning = function(cond) {
    print_and_log(paste('Warning occured: ',cond$message), level='warning')
    return(NULL)
  })

  if(!is.null(r) && r$status_code == 200)
  {
    output.tbl <- content(r,as = "parsed", simplifyVector = TRUE)
    output.tbl <- as.data.table(output.tbl$data)
  }

  ## output columns
  #"snpId" ==> rsID
  #"b37VariantId" ==> GTEx_ID b37
  #"variantId" ==> GTEx_ID b38

  return(output.tbl)
}

#' Query GTEx portal for Variant's genomic position based on rsID
#' Retrieves variant information from the GTEx portal using either
#' an rsID or a variant ID formatted as `CHR_POS_REF_ALT`.
#' If an rsID is provided, the function returns the corresponding
#' genomic positions in both GRCh37 and GRCh38 builds.
#' When searching for an rsID based on genomic position, the position
#' parameter should be specified according to the GRCh38 reference genome.
#'
#' @param id Character string representing the rsID (e.g., `"rs12345"`) or the variant ID
#' in the format `"CHR_POS_REF_ALT"` (e.g., `"1_1234567_A_T"`), depending on `type`.
#' @param type Character string specifying the type of query. Must be either `"rsid"` or `"varid"`.
#' @param file_path character, path to a file for saving results as Excell spreadsheet.
#'
#' @return A `data.table` containing variant information including:
#' - `rsid`: variant id in rsID format
#' - `chromosome`: chromosome number
#' - `position_b37`: genomic position
#' - `position_b38`: genomic position
#' - `ref`: reference allele
#' - `alt`: alternate allele
#'
#' @export
findGenomicPos <- function(id , type="rsid", file_path = NULL){

  if(is.null(id) || is.null(type) || !type %in% c("rsid","varid"))
    stop('Wrong input parameters',call. = FALSE)


  server = .SNPannotator$GTEx_eqtl_API
  ext="dataset/variant"
  output.tbl <- NULL
  query <- NULL

  id <- trimws(id)

  if(type=="rsid" && startsWith(x = tolower(id),prefix = 'rs'))
  {
    query <- sprintf("%s/%s?format=json&snpId=%s", server,ext,tolower(id))
  } else if(grepl(pattern = "^\\d+_\\d+_[ATCG]+_[ATCG]+$", x = id)) {
    query <- sprintf("%s/%s?format=json&variantId=chr%s_b38", server,ext,id)
  } else {
    stop(sprintf('ID not in correct format: %s',id),call. = FALSE)
  }


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
    output.tbl <- content(r,as = "parsed", simplifyVector = TRUE)
    output.tbl <- as.data.table(output.tbl$data)

    if(nrow(output.tbl) == 0 )
    {
      message("Variant not found.")
      return(NULL)
    }

    output.tbl[,position_b37 := strsplit(b37VariantId,"_")[[1]][2]]

    setnames(output.tbl,'pos',"position_b38")
    output.tbl <-  output.tbl[,c("snpId","chromosome","position_b37","position_b38","ref","alt")]
  }

  # return table if no file is specified
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

################################
############ run time functions#
################################


PingQTL <- function(server,build)
{
  accessible <- ifelse(build=='grch37', PingEBI(server), PingGTEx(server))

  return(accessible)
}
