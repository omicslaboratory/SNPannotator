returnVariantDatatable <- function(varNum , varInfo, database )
{

  overall_clinvar_table <- data.table()
  overall_gwascat_table <- data.table()
  overall_variant_table <- data.table()
  overall_phenotype_table <- data.table()
  overall_output_tables_list <- list()



  gVarId <- unique(unlist(varInfo$id))

  # create an empty data table where the first row is the target SNP
  varData <-  getVariantData(varNum,gVarId,varInfo,database, LD = 1)

  this_variant_table <- varData['variant_table'][[1]]
  this_clinvar_table <- varData['clinvar_table'][[1]]
  this_gwascat_table <- varData['gwascat_table'][[1]]
  this_phenotype_table <- varData['phenotype_table'][[1]]

  if(is.null(this_variant_table))
    return(NULL)

  overall_variant_table <- rbind(overall_variant_table , this_variant_table, fill =TRUE)
  overall_clinvar_table <- rbind(overall_clinvar_table , this_clinvar_table, fill =TRUE)
  overall_gwascat_table <- rbind(overall_gwascat_table , this_gwascat_table, fill =TRUE)
  overall_phenotype_table <- rbind(overall_phenotype_table , this_phenotype_table, fill =TRUE)

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

        this_proxy_variant_table <- thisVarInfo['variant_table'][[1]]
        this_proxy_clinvar_table <- thisVarInfo['clinvar_table'][[1]]
        this_proxy_gwascat_table <- thisVarInfo['gwascat_table'][[1]]
        this_proxy_phenotype_table <- thisVarInfo['phenotype_table'][[1]]

        if(!is.null(this_proxy_variant_table))
        {
          varInfo$LDlistFull[[i]] <-  this_proxy_variant_table

          overall_variant_table <- rbind(overall_variant_table , this_proxy_variant_table, fill =TRUE)
          overall_clinvar_table <- rbind(overall_clinvar_table , this_proxy_clinvar_table, fill =TRUE)
          overall_gwascat_table <- rbind(overall_gwascat_table , this_proxy_gwascat_table, fill =TRUE)
          overall_phenotype_table <- rbind(overall_phenotype_table , this_proxy_phenotype_table, fill =TRUE)

        }
      }

    }
  }

  # tab <- tab[order(-LD),]
  suppressWarnings(overall_variant_table[Chr =='X', Chr := '23'])
  suppressWarnings(overall_variant_table[Chr =='Y', Chr := '24'])
  suppressWarnings(overall_variant_table[,Chr := as.numeric(Chr)])
  suppressWarnings(overall_variant_table[,Pos := as.numeric(Pos)])

  overall_variant_table[is.na(Chr), Chr := ""]
  overall_variant_table[is.na(Pos), Pos := ""]

  if(nrow(overall_phenotype_table) > 0)
  {
    overall_phenotype_table$seq_region_name <- NULL
    overall_phenotype_table$start <- NULL
    overall_phenotype_table$variation_names <- NULL
    overall_phenotype_table$strand <- NULL
    overall_phenotype_table$type <- NULL
    overall_phenotype_table$end <- NULL
    overall_phenotype_table$external_id <- NULL
    overall_phenotype_table$attrib_type <- NULL
    overall_phenotype_table$associated_gene <- NULL
    overall_phenotype_table$datelastevaluated <- NULL
    overall_phenotype_table$submitter_name <- NULL
    overall_phenotype_table$mim <- NULL
    overall_phenotype_table$mutation_consequence <- NULL
    overall_phenotype_table$g2p_confidence <- NULL
    overall_phenotype_table$inheritance_type <- NULL
    overall_phenotype_table$clinvar_clin_sig <- NULL
    overall_phenotype_table$datelastevaluated <- NULL
    overall_phenotype_table$review_status <- NULL
    overall_phenotype_table$samples_tested <- NULL
    overall_phenotype_table$so_accession <- NULL
    overall_phenotype_table$mutated_samples <- NULL
    overall_phenotype_table$samples_mutation <- NULL

    setcolorder(overall_phenotype_table,c('id','phenotype','source'))
  }

  overall_output_tables_list <- list('variant_table' = overall_variant_table,
                                     'clinvar_table' = overall_clinvar_table,
                                     'gwascat_table' = overall_gwascat_table,
                                     'phenotype_table' = overall_phenotype_table)

  return(overall_output_tables_list)
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
  ref_all <- ""
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
  clinvar_table <- data.table()
  gwascat_table <- data.table()
  var_table <- data.table()
  output_table_list <- list()
  gwascat_phenotypes <- ""
  phenotype_table <- data.table()
  phenotype_table_filtered <- data.table()
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
      ref_all <- alleles[1]
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

      #===========================
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


        ##=============================
        # clinvar & gwas catalog

        if ("phenotypes" %in% colnames(trans_consq_tbl)) {

          phenotype_table <- as.data.table(trans_consq_tbl[["phenotypes"]][1])

          if(!is.null(phenotype_table))
          {
            # we dont want phenptypes associated with genes, ONLY variants
            if(all(is.element(c('type','phenotype'),names(phenotype_table))))
              phenotype_table <-  phenotype_table[!is.na(phenotype) & phenotype !='' & type=='Variation']

            if(nrow(phenotype_table) > 0 &&
               phenotype_table[grepl('clinvar',x = source,ignore.case = T) &
                               !grepl('phenotype not specified',x = phenotype),.N] > 0)
            {

              cols_to_check_clinvar <- c("id","phenotype","risk_allele","associated_gene","clinvar_clin_sig","pubmed_id")
              existing_cols_clinvar <- intersect(cols_to_check_clinvar, colnames(phenotype_table))

              clinvar_table <- phenotype_table[ grepl('clinvar',x = source,ignore.case = T) &
                                                  !grepl('phenotype not specified',x = phenotype) ,
                                                ..existing_cols_clinvar,drop=FALSE]
            }

            if(nrow(phenotype_table) > 0 &&
               phenotype_table[grepl('gwas_cat',x = source,ignore.case = T) &
                               !grepl('phenotype not specified',x = phenotype),.N] > 0)
            {
              cols_to_check_gwascat <- c("id","phenotype","risk_allele","p_value","beta_coef","odds_ratio")
              existing_cols_gwascat <- intersect(cols_to_check_gwascat, colnames(phenotype_table))

              gwascat_table <- phenotype_table[ grepl('gwas_cat',x = source,ignore.case = T) &
                                                  !grepl('phenotype not specified',x = phenotype),
                                                ..existing_cols_gwascat, drop = FALSE]

              # remove additional string about ukb
              gwascat_table[,phenotype := sub(pattern = " ukb data field.+",
                                              replacement = "",x = phenotype,ignore.case = TRUE)]

              # get the phenotype string
              phenos <- trimws(gwascat_table$phenotype)
              phenos <- tolower(phenos)
              phenos_unique <- unique(phenos)
              gwascat_phenotypes <- paste(phenos_unique, collapse = "; ")

            }

            # keep other phenotypes after extracting ClinVar and GwasCatalog
            if('source' %in% names(phenotype_table))
              phenotype_table_filtered <- phenotype_table[!grepl('gwas_cat',x = source,ignore.case = T) &
                                                            !grepl('clinvar',x = source,ignore.case = T),]
          }

        }

        #======================

        # # check if there are any transcripts for another gene, based on distance
        # if('distance' %in% names(trans_consq_tbl) &&
        #    trans_consq_tbl[distance  == "", .N] > 0 &&
        #    'variant_allele' %in% names(trans_consq_tbl))
        #   trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all & distance  == "",]
        # else if('variant_allele' %in% names(trans_consq_tbl))
        #   trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all,]
        #
        # if(!is.null(trans_consq_tbl) && is.element('variant_allele',names(trans_consq_tbl)))
        #   reg_cadd = getString(trans_consq_tbl[variant_allele == minor_all,]$cadd_phred)
        #

        # if('gene_symbol' %in% names(trans_consq_tbl))
        #   trans_consq_tbl <- trans_consq_tbl[!duplicated(gene_symbol ), ]

      }

      #======================
      # the intergenic_consequences sometime exists instead of transcript_consequences

      if(!is.null(varInfoList$intergenic_consequences[list.index][[1]]) & nrow(phenotype_table) == 0)
      {
        trans_consq_tbl <-  as.data.table(varInfoList$intergenic_consequences[list.index][[1]])
        {

          phenotype_table <- as.data.table(trans_consq_tbl[["phenotypes"]][1])

          if(!is.null(phenotype_table))
          {

            # we dont want phenptypes associated with genes, ONLY variants
            if(all(is.element(c('type','phenotype'),names(phenotype_table))))
              phenotype_table <-  phenotype_table[!is.na(phenotype) & phenotype !='' & type=='Variation']

            if(nrow(phenotype_table) > 0 &&
               phenotype_table[grepl('clinvar',x = source,ignore.case = T) &
                               !grepl('phenotype not specified',x = phenotype),.N] > 0)
            {

              cols_to_check_clinvar <- c("id","phenotype","risk_allele","associated_gene","clinvar_clin_sig","pubmed_id")
              existing_cols_clinvar <- intersect(cols_to_check_clinvar, colnames(phenotype_table))

              clinvar_table <- phenotype_table[ grepl('clinvar',x = source,ignore.case = T) &
                                                  !grepl('phenotype not specified',x = phenotype) ,
                                                ..existing_cols_clinvar,drop=FALSE]
            }

            if(nrow(phenotype_table) > 0 &&
               phenotype_table[grepl('gwas_cat',x = source,ignore.case = T) &
                               !grepl('phenotype not specified',x = phenotype),.N] > 0)
            {
              cols_to_check_gwascat <- c("id","phenotype","risk_allele","p_value","beta_coef","odds_ratio")
              existing_cols_gwascat <- intersect(cols_to_check_gwascat, colnames(phenotype_table))

              gwascat_table <- phenotype_table[ grepl('gwas_cat',x = source,ignore.case = T) &
                                                  !grepl('phenotype not specified',x = phenotype),
                                                ..existing_cols_gwascat, drop = FALSE]

              # remove additional string about ukb
              gwascat_table[,phenotype := sub(pattern = " ukb data field.+",
                                              replacement = "",x = phenotype,ignore.case = TRUE)]

              # get the phenotype string
              phenos <- trimws(gwascat_table$phenotype)
              phenos <- tolower(phenos)
              phenos_unique <- unique(phenos)
              gwascat_phenotypes <- paste(phenos_unique, collapse = "; ")

            }

            # keep other phenotypes after extracting ClinVar and GwasCatalog
            if('source' %in% names(phenotype_table))
              phenotype_table_filtered <- phenotype_table[!grepl('gwas_cat',x = source,ignore.case = T) &
                                                            !grepl('clinvar',x = source,ignore.case = T),]
          }

        }
      }

      #===================================
      #===================================

      if(nrow(trans_consq_tbl) > 0)
      {

        # make sure the nearest gene is reported
        # this table contains the gene with the transcript that might not be the nearest
        # if distance ==NULL ==> intronic
        # if NOT, the data is not reliable


        # the distance column is as a list => convert to char
        unwanted_genes <- NULL

        if('distance' %in% names(trans_consq_tbl))
        {
          trans_consq_tbl <- trans_consq_tbl[, distance := sapply(distance, function(x) {
            if (length(x) == 0) "" else as.character(x)
          })]
          unwanted_genes <- unlist(unique(trans_consq_tbl[distance != "",]$gene_id ))
        }

        if('gene_symbol' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[, gene_symbol := sapply(gene_symbol, function(x) {
            if (length(x) == 0) "" else as.character(x)
          })]

        if('gene_id' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[, gene_id := sapply(gene_id, function(x) {
            if (length(x) == 0) "" else as.character(x)
          })]


        # remove genes with distance !="" that are duplicates
        if(!is.null(unwanted_genes))
          trans_consq_tbl <- trans_consq_tbl[!is.element(gene_id,unwanted_genes), ]

        # put gene_id in place of gene_symbol if missing
        if(all(is.element(c('gene_symbol','gene_id') , names(trans_consq_tbl))))
          trans_consq_tbl[gene_symbol=="", gene_symbol := gene_id]



        if('distance' %in% names(trans_consq_tbl) &&
           trans_consq_tbl[distance  == "", .N] > 0 &&
           'variant_allele' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all & distance  == "",]
        else if('variant_allele' %in% names(trans_consq_tbl))
          trans_consq_tbl <- trans_consq_tbl[variant_allele == minor_all,]


        if(!is.null(trans_consq_tbl) && is.element('variant_allele',names(trans_consq_tbl)))
        {
          geneNames <- paste(unique(trans_consq_tbl[variant_allele == minor_all,]$gene_symbol), collapse = ",")

          geneIDs <- paste(unique(trans_consq_tbl[variant_allele == minor_all,]$gene_id), collapse = ",")
          geneBiotype <- paste(unique(trans_consq_tbl[variant_allele == minor_all,]$biotype), collapse = ",")
          reg_cadd = getString(trans_consq_tbl[variant_allele == minor_all,]$cadd_phred)

        }


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




      # for intergenic variants
      if((is.null(reg_cadd) || all(reg_cadd == '')) && !is.null(varInfoList$intergenic_consequences[list.index]))
      {
        intergenic_consequences = as.data.table(varInfoList$intergenic_consequences[list.index])


        if(nrow(intergenic_consequences) > 0 &
           'variant_allele' %in% names(intergenic_consequences) &
           'cadd_phred' %in% names(intergenic_consequences))
          reg_cadd = getString(intergenic_consequences[variant_allele == minor_all,]$cadd_phred)
      }

      # if CADD is still empty
      if(nrow(reg_feat_tbl) > 0)
      {
        reg_id <- paste(unique(reg_feat_tbl$regulatory_feature_id), collapse = ",")
        reg_biotype <- paste(unique(reg_feat_tbl$biotype), collapse = ",")
        if(is.null(reg_cadd) || reg_cadd=='' || length(reg_cadd) ==0)
          reg_cadd <- paste(unique(reg_feat_tbl$cadd_phred), collapse = ",")
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
      var_table <- data.table(`#gSNP` = varNum,
                              gSNP = gVarId,
                              Linked_SNP = name,
                              Chr = chr,
                              Pos = pos,
                              Class = variant_class,
                              Ref_Allele=ref_all,
                              Alt_Allele = minor_all,
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
                              Polyphen_prediction = polyphen_prediction,
                              Phenotype = gwascat_phenotypes
      )

      # varInfoList[['variant_data']]= tab
      # varInfoList[['GO_terms']]= go_tbl

      output_table_list <- list('variant_table' = var_table,
                                'clinvar_table' = clinvar_table,
                                'gwascat_table' = gwascat_table,
                                'phenotype_table' = phenotype_table_filtered)
      return(output_table_list)
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
  ext_vep <- sprintf("/vep/human/id/%s?CADD=1&mutfunc=1&variant_class=1&Phenotypes=1",rsID)
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
            body = sprintf('{"CADD":"1","Phenotypes":"1", "mutfunc":"1","variant_class":"1", "ids" : %s }',paste0('["', paste(IDs, collapse = '", "'), '"]')))

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

getRSID <- function(chr,pos1,pos2)
{
  id_query <- sprintf("/overlap/region/human/%s:%s-%s?feature=variation",chr,pos1,pos2)
  id_tbl <- fetch(id_query,server)

  if(!is.null(id_tbl))
    var_vep <- as.list(var_vep)

  return(var_vep)

}
