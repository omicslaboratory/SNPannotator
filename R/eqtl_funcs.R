
################
### Ensembl ####
################

PingEBI <- function(server)
{
  print_and_log("Pinging EBI server ... ",LF = FALSE)

  tryCatch(
    {
      r <- GET(server, content_type("application/json"))

      if (!is.null(r) && r$status_code == 200)
      {
        print_and_log('accessible.')
      }else
      {
        print_and_log('not accessible.', level='fatal')
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in EBI ping.',cond$message), level='fatal')
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in EBI ping.',cond$message), level='warning')
    }
  )
}

get.EBI.tissues <- function(server)
{
  r <- GET(sprintf('%s/tissues?size=-1',servre), content_type("application/json"))

  tissues <- data.table(fromJSON(toJSON(content(r)))$`_embedded`$tissues)

  return(tissues)
}

#' Get the list of available tissue data from EBI website.
#'
#' This function can be used to check the available tissue data from EBI API server.
#'
#' @param server the address for EBI API server.
#' @return A data.table containing the results.
#' @export
#'
getEBIgroups <- function()
{
  r <- GET(sprintf('%s/qtl_groups?size=150',.SNPannotator$EBI_eqtl_API ), content_type("application/json"))

  qtl_groups <- data.table(fromJSON(toJSON(content(r)))$`_embedded`$qtl_groups)

  return(qtl_groups)
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

get.EBI.genes.all.tissues <- function(server,rsid,ebi.p.value.threshold)
{
  r <- GET(sprintf('%s/associations/%s?p_upper=%s&size=1000',server,rsid,ebi.p.value.threshold),
           content_type("application/json"))

  out <- data.table()

  if(r$status_code != 200)
    return(out)

  out=data.table(t(sapply(content(r)$`_embedded`$associations,
                          function(x) return(c(x$rsid,x$ref,x$alt,x$maf,x$gene_id,x$qtl_group,x$beta,x$se,x$median_tpm,x$pvalue,x$study_id)))))


  if(!is.null(out) && nrow(out) > 0 && ncol(out) == 11)
  {
    names(out) <- c('rsId','Ref','Alt','MAF','GeneId','eQTL_group','Beta','SE','Median_tpm','pValue','Study')
    return(out)
  }
  else
    return(data.table())
}

find.eqtl.ebi <- function(server,
                          snpIDs,
                          eqtl_group,
                          ebi.p.value.threshold)
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
        if(is.null(eqtl_group))
          out = rbind(out,
                      get.EBI.genes.all.tissues(server,
                                                snpID,
                                                ebi.p.value.threshold),
                      fill=TRUE)
        else
          out = rbind(out,
                      get.EBI.genes(server,
                                    snpID,
                                    eqtl_group,
                                    ebi.p.value.threshold),
                      fill=TRUE)

        if(!is.null(pb))
          pb$tick(1)
      }

      # remove unneccessary columns
      if(is.element('Ref', names(out)))
        out$Ref <- NULL

      if(is.element('MAF', names(out)))
        out$MAF <- NULL
      #

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

  data <- merge(x=data,y=gdata[,list(id,name,type)],by.x='GeneId',by.y='id',all.x=TRUE,sort=FALSE)
 # setcolorder(data,c('rsId','Ref','Alt','MAF','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'))

  if(all(is.element(c('rsId','Alt','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'), names(data))))
    setcolorder(data,c('rsId','Alt','GeneId','name','type','eQTL_group','Beta','SE','Median_tpm','pValue','Study'))

  if(is.element('name',names(data)))
    setnames(data,'name','Gene')

  if(is.element('type',names(data)))
    setnames(data,'type','Type')

  data[is.na(Gene) | Gene =='', Gene := GeneId]

  return(data)
}
#####DEPRECATED
# library(httr)
# library(jsonlite)
#
# server='https://gtexportal.org/api/v2'
# gtex.version = 'gtex_v8'
# snpID = 'RS698'
# snpID = 'chr4_99339632_T_C_b38'
# tissueID='Thyroid'
#
# r <- GET(sprintf("%s/association/singleTissueEqtl?format=json&variantId=%s&tissueSiteDetailId=%s&datasetId=%s",
#                  server,snpID,tissueID,gtex.version), content_type("application/json"))
# content(r)$data
#
# r <- GET(sprintf("%s/association/singleTissueEqtlByLocation?format=json&tissueSiteDetailId=%s&datasetId=%s&chromosome=chr4&start=99339600&end=99339650",
#                  server,tissueID,gtex.version), content_type("application/json"))
#
#
# r <- GET(sprintf("%s/association/singleTissueEqtl?format=json&snpId=%s&tissueSiteDetailId=%s&datasetId=%s",
#                  server,snpID,tissueID,gtex.version), content_type("application/json"))
#
#
# a <- get.Gtex.tissues(server,gtex.version,T)
#
