returnVariantDatatable <- function(varNum , varInfo, database )
{

  gVarId <- unique(unlist(varInfo$id))
  tab <- data.table()

  # create an empty data table where the first row is the target SNP
  varData <-  getVariantData(varNum,gVarId,varInfo,database, LD = 1)
  if(is.null(varData))
    return(NULL)

  tab <- rbind(tab , varData, fill =TRUE)

  # fill the above data table if any variant found in high LD
  if(!is.null(varInfo$LDList) && !is.null(ncol(varInfo$LDList)))
  {
    LDTable <- setDT(varInfo$LDList)

    for(i in seq_len(length(varInfo$LDlistFull)))
    {
      if(!is.null(varInfo$LDlistFull[[i]]$id) )#&& length(varInfo$LDlistFull[[i]]$population) > 0)
      {

        LD = LDTable[variation2 == unique(unlist(varInfo$LDlistFull[[i]]$id)) , r2]

        thisVarInfo <- NULL
        thisVarInfo <-  getVariantData(varNum,gVarId,varInfo$LDlistFull[[i]],database,LD)
        if(!is.null(thisVarInfo))
        {
          varInfo$LDlistFull[[i]] <-  thisVarInfo
          tab <- rbind(tab, thisVarInfo, fill=TRUE)
        }
      }

    }
  }

  # tab <- tab[order(-LD),]
  suppressWarnings(tab[Chr =='X', Chr := '23'])
  suppressWarnings(tab[Chr =='Y', Chr := '24'])
  suppressWarnings(tab[,Chr := as.numeric(Chr)])
  suppressWarnings(tab[,Pos := as.numeric(Pos)])

  tab[is.na(Chr), Chr := ""]
  tab[is.na(Pos), Pos := ""]

  return(tab)
}


getVariantData <- function(varNum, gVarId, varInfoList, database , LD)
{

  getString <- function(x) return(as.character(unique(unlist(x))))

  ## transcript consequence table ----
  name <- ""
  colocated_name <- ""
  population <- ""
  chr <- ""
  pos <- ""
  band <- ""
  strand <- ""
  most_severe_consequence <- ""
  geneNames <- ""
  geneIDs <- ""
  band <- ""
  geneBiotype <- ""
  reg_cadd <- ""
  reg_id <- ""
  reg_biotype <- ""
  minor_all <- ""
  minor_all_frq <- ""
  alleles <- ""
  trans_consq_tbl <- data.table()
  coloc_var_tbl <- data.table()
  reg_feat_tbl <- data.table()
  frq_tbl <- data.table()
  allele_string <- ""
  sift_prediction <- ""
  polyphen_prediction <- ""
  variant_class <- ""
  # go_tbl <- data.table()
  #=============================

  list.index <- 1
  if(length(varInfoList$seq_region_name) > 1)
  {
    regions <- suppressWarnings(sapply(varInfoList$seq_region_name , as.numeric))
    list.index <- which(!is.na(regions))
  }

  name <- getString(varInfoList$id[list.index])
  chr <- getString(varInfoList$seq_region_name[list.index])
  pos <- getString(varInfoList$start[list.index])
  strand <- getString(varInfoList$strand[list.index])
  most_severe_consequence <- getString(varInfoList$most_severe_consequence[list.index])

  if(!is.null(varInfoList$variant_class))
    variant_class <- varInfoList$variant_class[[1]]

  if(is.na(chr) || is.na(pos) || chr == '' || pos == '')
    return(NULL)


  out <- tryCatch(
    {

      alleles_string <- getString(varInfoList$allele_string[list.index])

      alleles <- unlist(strsplit(alleles_string,split = '/'))
      minor_all <- alleles[2]
      # TODO: check other alleles
      # minor_all = if(length(alleles) ==2)
      # {
      #   alleles[2]
      # } else{
      #   alleles[3]
      # }

      population <- if(grepl(x = database, pattern = 'eur',ignore.case = T))
      {'eur'} else if( grepl(x = database, pattern = 'afr',ignore.case = T))
      {'afr'} else if(grepl(x = database, pattern = 'sas',ignore.case = T))
      { 'sas'} else if(grepl(x = database, pattern = 'eas',ignore.case = T))
      {'eas'}else if(grepl(x = database, pattern = 'amr',ignore.case = T))
      { 'amr'} else
      {'af'}



      #=============================
      var.index <- 1 # which of the returned data to be used baase on deq_region_name
      if(!is.null(varInfoList$colocated_variants[[list.index]]))
      {
        colocated_variants_table <- varInfoList$colocated_variants[list.index]
        # some variants have more than one row, select the one with nmeric chromosome
        #regions <- suppressWarnings(sapply(varInfoList$colocated_variants, function(x) as.numeric(x$seq_region_name)))
        same_allele <- sapply(colocated_variants_table, function(x) {
          if(is.element('allele_string',names(x)))
            x$allele_string == alleles_string
        })
        # find the row with the same alleles
        if(any(same_allele))
          var.index <- which(same_allele)

        # check if the colocated variantID is the same as original variantID
        # upadte the RS number if not the same
        colocated_name <- trimws(tolower(colocated_variants_table[[1]][var.index,]$id))
        if(!is.null(colocated_name) &&
           length(colocated_name) == 1 &&  # some chr positions map to more than one variant "rs2003944"    "rs2100210096" "rs2100210107"
           startsWith(colocated_name,'rs') &&
           (colocated_name != name))
        {
          print_and_log(sprintf('Variant name updated from %s to %s', name,colocated_name),
                        display=FALSE,
                        level='warning')

          name <- colocated_name
        }

        # find the row with minor allele
        frq_tbl <- colocated_variants_table[[1]][var.index,]$frequencies

        # which alleles have frq information
        avail_alleles <- names(frq_tbl)
        # find the reported frq in population
        avail_frq <- sapply(avail_alleles, function(x) frq_tbl[[x]][[population]])
        if(!is.null(avail_frq))
        {

          #remove those without frq
          avail_frq <- Filter(Negate(is.null),avail_frq)

          frq_vector <- unlist(avail_frq)
          if(!is.null(frq_vector) && length(frq_vector) > 0)
          {

            frq_vector <- frq_vector[which.min(frq_vector)]

            if(!is.null(names(frq_vector)))
              minor_all <- names(frq_vector)

            minor_all_frq <- frq_vector
          }

        }

        #if(!is.null(frq_tbl) && !is.element(minor_all , names(frq_tbl)))
        #  minor_all <- names(frq_tbl)[1]

        # if(!is.null(frq_tbl))
        # {
        #   #minor_all_frq <- getString(as.data.table(varInfoList$colocated_variants[[1]]$frequencies[[minor_all]])[index,][[population]])
        #   minor_all_frq <- getString(frq_tbl[[minor_all]][[population]])
        # }

      }
      #=============================
      if(!is.null(varInfoList$transcript_consequences[list.index]))
      {
        trans_consq_tbl <-  as.data.table(varInfoList$transcript_consequences[list.index][[1]])

        # we dont; want the gene with most sever consequence !!
        if('consequence_terms' %in% names(trans_consq_tbl) &&
           trans_consq_tbl[consequence_terms == 'intron_variant', .N] > 0)
        {
          trans_consq_tbl <- trans_consq_tbl[consequence_terms == 'intron_variant',]

        } else if('consequence_terms' %in% names(trans_consq_tbl) &&
                  most_severe_consequence != "" &&
                  trans_consq_tbl[consequence_terms == most_severe_consequence, .N] > 0)
        {
          trans_consq_tbl <- trans_consq_tbl[consequence_terms == most_severe_consequence,]
        }

        # check if there are any transcripts for another gene, based on distance
        if('distance' %in% names(trans_consq_tbl) &&
           trans_consq_tbl[distance  == "NULL", .N] > 0 &&
           'variant_allele' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all & distance  == "NULL",]
        else if('variant_allele' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all,]

        # TODO: find a better sort
        if(is.element('go', names(trans_consq_tbl)))
        {
          go_count <- apply(trans_consq_tbl,1, function(x) {
            tbl <- as.data.table(x$go)
            nrow(tbl)
          })
          go_order <- order(go_count,decreasing = TRUE)
          trans_consq_tbl[go_order,]
        }

        ### GO Table ----
        # DEPRECATED - USE STRING DB
        # go_tbl <- as.data.table(trans_consq_tbl[1]$go)
        #
        # if(!is.null(go_tbl) && nrow(go_tbl) >0)
        # {
        #   go_tbl <- go_tbl[, lapply(.SD, function(x) as.character(unlist(x)))]
        #   go_tbl$`#gSNP` <- varNum
        #   go_tbl$rsID <- name
        #   go_tbl$Gene <- unique(unlist(trans_consq_tbl[1]$gene_symbol))
        #   go_tbl$GeneID <- unique(unlist(trans_consq_tbl[1]$gene_id))
        #   setcolorder(go_tbl,c('#gSNP','rsID','Gene','GeneID'))
        #
        # }


        if('gene_symbol' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[!duplicated(gene_symbol ), ]

      }

      if(nrow(trans_consq_tbl) > 0)
      {

        # make sure the nearest gene is reported
        # this table contains the gene with the transcript that might not be the nearest
        # if distance ==NULL ==> intronic
        # if NOT, the data is not reliable
        if('distance' %in% names(trans_consq_tbl) && trans_consq_tbl[distance  == "NULL", .N] > 0)
        {
          geneNames <- paste(unique(trans_consq_tbl$gene_symbol), collapse = ",")

          geneIDs <- paste(unique(trans_consq_tbl$gene_id), collapse = ",")
          geneBiotype <- paste(unique(trans_consq_tbl$biotype), collapse = ",")
        }

        reg_cadd <- unlist(unique(trans_consq_tbl$cadd_phred))

        ## TODO
        if('sift_prediction' %in% names(trans_consq_tbl))
          sift_prediction <- unlist(unique(trans_consq_tbl$sift_prediction))

        if('polyphen_prediction' %in% names(trans_consq_tbl))
          polyphen_prediction  <- unlist(unique(trans_consq_tbl$polyphen_prediction ))
      }

      #=============================

      #TODO perhaps use [list.index]
      if(!is.null(varInfoList$regulatory_feature_consequences[list.index]))
      {
        reg_feat_tbl <-  as.data.table(varInfoList$regulatory_feature_consequences[list.index][[1]])

        if('variant_allele' %in% names(reg_feat_tbl))
          reg_feat_tbl <- reg_feat_tbl[variant_allele == minor_all ,]
      }

      if(nrow(reg_feat_tbl) > 0)
      {
        reg_id <- paste(unique(reg_feat_tbl$regulatory_feature_id), collapse = ",")
        reg_biotype <- paste(unique(reg_feat_tbl$biotype), collapse = ",")
        if(is.null(reg_cadd) || reg_cadd=='')
          reg_cadd <- paste(unique(reg_feat_tbl$cadd_phred), collapse = ",")
      }


      # for intergenic variants
      if((is.null(reg_cadd) || reg_cadd == '') && !is.null(varInfoList$intergenic_consequences[list.index]))
      {
        intergenic_consequences = as.data.table(varInfoList$intergenic_consequences[list.index])


        if(nrow(intergenic_consequences) > 0 &
           'variant_allele' %in% names(intergenic_consequences) &
           'cadd_phred' %in% names(intergenic_consequences))
          reg_cadd = getString(intergenic_consequences[variant_allele == minor_all,]$cadd_phred)
      }

      # make sure CADD has a value
      if(is.null(reg_cadd) || identical(reg_cadd , character(0)))
        reg_cadd <- ""

      if(identical(minor_all_frq , character(0)))
        minor_all_frq <- ""

      if(identical(most_severe_consequence , character(0)))
        most_severe_consequence <- ""

      if(identical(geneNames , character(0)))
        geneNames <- ""

      if(identical(geneIDs , character(0)))
        geneIDs <- ""


      if(identical(geneBiotype , character(0)))
        geneBiotype <- ""

      if(identical(reg_id , character(0)))
        reg_id <- ""

      if(identical(reg_biotype , character(0)))
        reg_biotype <- ""

      if(identical(sift_prediction , character(0)))
        sift_prediction <- ""

      if(identical(polyphen_prediction , character(0)))
        polyphen_prediction <- ""


      ## create table ----
      tab <- data.table(`#gSNP` = varNum,
                        gSNP = gVarId,
                        Linked_SNP = name,
                        Chr = chr,
                        Pos = pos,
                        Class = variant_class,
                        Allele = minor_all,
                        Allele_frequency = minor_all_frq,
                        LD = LD ,
                        Cytoband = band,
                        Most_severe_consequence = most_severe_consequence,
                        Gene = geneNames,
                        GeneId = geneIDs,
                        geneBiotype = geneBiotype,
                        RegId = reg_id,
                        RegType = reg_biotype,
                        CADD = reg_cadd,
                        SIFT_prediction = sift_prediction,
                        Polyphen_prediction = polyphen_prediction
      )

      # varInfoList[['variant_data']]= tab
      # varInfoList[['GO_terms']]= go_tbl
      return(tab)
    },
    error=function(x){
      print_and_log(sprintf('Error in %s: %s', name,x$message),
                    level='warning',
                    LF = TRUE,
                    display=FALSE)
      return(NULL)
    }
  )

  return(out)
}


getVariantInfo <- function(rsID,server,pb=NULL)
{
  on.exit({
    if(!is.null(pb))
      pb$tick(1)
  })

  var_vep <- NULL

  # add variant info
  # ext <- sprintf("/variation/human/%s?phenotypes=1;pops=1",rsID)
  # var <- fetch(ext,server)
  #
  # if(is.null(var) || length(var$mappings) == 0)
  # {
  #   message(sprintf("no data found"))
  #   return(NULL)
  # }

  # add VEP info
  # full list
  ## "/vep/human/id/rs72844527?CADD=1&GO=1&Phenotypes=1&canonical=1&merged=1&minimal=1&per_gene=1&protein=1&variant_class=1&mutfunc=1&AlphaMissense=1&DisGeNET=1&EVE=1&Enformer=1&Geno2MP=1&LoF=1&Mastermind=1"
  ext_vep <- sprintf("/vep/human/id/%s?CADD=1&mutfunc=1&variant_class=1",rsID)
  var_vep <- fetch(ext_vep,server)

  if(!is.null(var_vep))
    var_vep <- as.list(var_vep)

  return(var_vep)

}

getMultiVariantInfo <- function(IDs, server)
{

  ext <- "/vep/human/id"
  r <- POST(paste(server, ext, sep = ""),
            content_type("application/json"),
            accept("application/json"),
            body = sprintf('{"CADD":"1", "mutfunc":"1","variant_class":"1", "ids" : %s }',paste0('["', paste(IDs, collapse = '", "'), '"]')))

  stop_for_status(r)

  tbl <- fromJSON(toJSON(content(r)))
  tbl$seq_region_name = suppressWarnings(as.numeric(tbl$seq_region_name))
  tbl <- tbl[!is.na(tbl$seq_region_name),]


  tbl.list <- list()

  for(i in 1:nrow(tbl))
  {
    rs <- unlist(tbl[i,]$id)
    tbl.list[[rs]] <- as.list(tbl[i,])
  }

  Sys.sleep(2)
  return(tbl.list)
}

