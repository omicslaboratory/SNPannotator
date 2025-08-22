# this function is called from withni the pipeline
find.gene.partners<- function(geneList,
                              server,
                              STRING_config_list)
{
  geneList.count <- length(unlist(strsplit(geneList, "%0d")))

  # already too many
  if(geneList.count > 250)
  {
    print_and_log("No partrners added from STRING.",display=FALSE)
    return(geneList)
  }


  tryCatch(
    {
      query <- sprintf("%s/tsv/interaction_partners?identifiers=%s&species=9606",server,geneList)

      if(STRING_config_list[['required_score']] == 0 )
      {
        query=paste(sprintf("%s&%s",query,'required_score=700'))
        print_and_log("STRING parameter for STRING interaction partners: required_score = 700",display=FALSE)
      }


      if(geneList.count > 100 )
      {
        query=paste(sprintf("%s&%s",query,'limit=2'))
        print_and_log("STRING parameter for STRING interaction partners: limit = 2",display=FALSE)
      } else if(geneList.count > 50 )
      {
        query=paste(sprintf("%s&%s",query,'limit=5'))
        print_and_log("STRING parameter for STRING interaction partners: limit = 5",display=FALSE)
      } else
      {
        query=paste(sprintf("%s&%s",query,'limit=10'))
        print_and_log("STRING parameter for STRING interaction partners: limit = 10",display=FALSE)
      }


      r <- GET(query)

      if (!is.null(r) && r$status_code == 200)
      {
        string.partner.tb <- as.data.table(content(r, show_col_types=FALSE))

        if(nrow(string.partner.tb) > 0)
        {
          geneNames=unique(c(string.partner.tb$preferredName_A,
                             string.partner.tb$preferredName_B))

          print_and_log(paste('Gene list for STRING with partners:',paste(geneNames,collapse = ';')),display=FALSE)
          print_and_log(paste('Gene list for STRING with partners (count):',length(geneNames)),display=FALSE)

          return(paste(geneNames,collapse='%0d'))
        }
        else
        {
          print_and_log(paste('No interaction partner data from STRING DB was found.',content(r,encoding = "UTF-8"),sep = '\n'),display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log(paste('error in interaction partner analysis from STRING.',content(r,encoding = "UTF-8"),sep = '\n'), level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING interaction partner analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING interaction partner analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

#this function is called when only String function is called
find.gene.partners.standalone<- function(geneList,
                                         server,
                                         STRING_config_list)
{

  print_and_log(paste('Gene list for STRING without partners:',paste(geneList,collapse = ';')),display=FALSE)
  print_and_log(paste('Gene list for STRING without partners (count):',length(geneList)),display=FALSE)


  tryCatch(
    {
      query <- sprintf("%s/tsv/interaction_partners?identifiers=%s&species=9606",server, paste(geneList,collapse = "%0d"))


      query=paste(sprintf("%s&required_score=%s",query,STRING_config_list[['required_score']]))
      print_and_log(sprintf("STRING parameter for STRING interaction partners: required_score = %s",STRING_config_list[['required_score']]),display=FALSE)

      query=paste(sprintf("%s&limit=%s",query,STRING_config_list[['limit']]))
      print_and_log(sprintf("STRING parameter for STRING interaction partners: limit = %s",STRING_config_list[['limit']]),display=FALSE)





      r <- GET(query)

      if (!is.null(r) && r$status_code == 200)
      {
        string.partner.tb <- as.data.table(content(r, show_col_types=FALSE))

        if(nrow(string.partner.tb) > 0)
        {
          geneNames=unique(c(string.partner.tb$preferredName_A,
                             string.partner.tb$preferredName_B))

          print_and_log(paste('Gene list for STRING with partners:',paste(geneNames,collapse = ';')),display=FALSE)
          print_and_log(paste('Gene list for STRING with partners (count):',length(geneNames)),display=FALSE)

          return(paste(geneNames,collapse='%0d'))
        }
        else
        {
          print_and_log(paste('No interaction partner data from STRING DB was found.',content(r,encoding = "UTF-8"),sep = '\n'),display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log(paste('error in interaction partner analysis from STRING.',content(r,encoding = "UTF-8"),sep = '\n'), level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING interaction partner analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING interaction partner analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

do.string.analysis <- function(geneList,
                               server,
                               STRING_config_list,
                               standalone = FALSE)
{

  #===========
  print_and_log("Checking STRING interaction partners ...", LF = FALSE)

  string.gene.partner.list <- NULL

  if(standalone == TRUE)
  {
    string.gene.partner.list <- find.gene.partners.standalone(geneList, server,STRING_config_list)
  }
  else
  {
    string.gene.partner.list <- find.gene.partners(geneList, server,STRING_config_list)
  }

  if(!is.null(string.gene.partner.list))
    print_and_log("Checking STRING interaction partners ...done")



  #===========
  string.report.list <- list()
  string.report.list[['pngContent_network']] <- NULL
  string.report.list[['pngContent_enrichment_function']] <- NULL
  string.report.list[['pngContent_enrichment_process']] <- NULL
  string.report.list[['pngContent_enrichment_component']] <- NULL
  string.report.list[['pngContent_enrichment_diseases']] <- NULL
  string.report.list[['pngContent_enrichment_tissues']] <- NULL
  string.report.list[['enrichment']] <- NULL
  string.report.list[['annotation']] <- NULL
  string.report.list[['webLink']] <- NULL
  string.report.list[['description']] <- NULL

  pngContent <- NULL
  enrich.report <- NULL
  annot.report <- NULL
  webLink <- NULL
  description <- NULL


  if(is.null(string.gene.partner.list))
  {
    print_and_log("No genes found for STRING DB analysis  ...")
    return(string.report.list)
  }
  #===========

  print_and_log("Checking functional enrichment ...", LF = FALSE)

  if(STRING_config_list[['functional_enrichment']] == TRUE)
  {
    enrich.report <- get.string.functional.enrichment(server, string.gene.partner.list)
    print_and_log('Checking functional enrichment ...done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(enrich.report))
    string.report.list[['enrichment']] <- enrich.report


  #===========

  # print_and_log("Checking functional annotation ...", LF = FALSE)
  # if(STRING_config_list[['functional_annotation']] == TRUE)
  # {
  #   annot.report <- get.string.annot.enrichment(server, string.gene.partner.list)
  #   print_and_log('Checking functional annotation ...done.')
  #   Sys.sleep(1)
  # }else{
  #   print_and_log('skipped.')
  # }
  #
  # if(!is.null(annot.report))
  #   string.report.list[['annotation']] <- annot.report

  #===========

  print_and_log("Checking descriptions ...", LF = FALSE)

  if(STRING_config_list[['description']] == TRUE)
  {
    description <- get.string.description(server, string.gene.partner.list)
    print_and_log('Checking description ...done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(description))
    string.report.list[['description']] <- description

  #===========

  print_and_log("Checking weblink ...", LF = FALSE)

  if(STRING_config_list[['webLink']] == TRUE)
  {
    webLink <- get.string.webLink(server, string.gene.partner.list)
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(webLink))
    string.report.list[['webLink']] <- webLink

  #===========

  print_and_log("Checking network image ...", LF = FALSE)

  if(STRING_config_list[['network_image']] == TRUE)
  {
    pngContent <- get.string.network.image(server, string.gene.partner.list,STRING_config_list)
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_network']] <- pngContent

  #===== enrichment image function ======

  print_and_log("Checking enrichment image (Function) ...", LF = FALSE)

  pngContent <- NULL

  if(STRING_config_list[['enrichment_image']] == TRUE)
  {
    pngContent <- get.string.enrichment.image(server, string.gene.partner.list,STRING_config_list,"Function")
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_enrichment_function']] <- pngContent

  #===== enrichment image process ======

  print_and_log("Checking enrichment image (Process) ...", LF = FALSE)

  pngContent <- NULL

  if(STRING_config_list[['enrichment_image']] == TRUE)
  {
    pngContent <- get.string.enrichment.image(server, string.gene.partner.list,STRING_config_list,"Process")
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_enrichment_process']] <- pngContent


  #===== enrichment image component ======

  print_and_log("Checking enrichment image (Component) ...", LF = FALSE)

  pngContent <- NULL

  if(STRING_config_list[['enrichment_image']] == TRUE)
  {
    pngContent <- get.string.enrichment.image(server, string.gene.partner.list,STRING_config_list,"Component")
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_enrichment_component']] <- pngContent


  #===== enrichment image DISEASES ======

  print_and_log("Checking enrichment image (Diseases) ...", LF = FALSE)

  pngContent <- NULL

  if(STRING_config_list[['enrichment_image']] == TRUE)
  {
    pngContent <- get.string.enrichment.image(server, string.gene.partner.list,STRING_config_list,"DISEASES")
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_enrichment_diseases']] <- pngContent

  #===== enrichment image TISSUES ======

  print_and_log("Checking enrichment image (Tissues) ...", LF = FALSE)

  pngContent <- NULL

  if(STRING_config_list[['enrichment_image']] == TRUE)
  {
    pngContent <- get.string.enrichment.image(server, string.gene.partner.list,STRING_config_list,"TISSUES")
    print_and_log('done.')
    Sys.sleep(1)
  }else{
    print_and_log('skipped.')
  }

  if(!is.null(pngContent))
    string.report.list[['pngContent_enrichment_tissues']] <- pngContent

  #=================================

  return(string.report.list)
}


PingSTRING <- function(server)
{
  print_and_log("Pinging STRING DB server ... ",LF = FALSE)
  accessible <- TRUE

  tryCatch(
    {
      r <- GET(paste(server,"json/version",sep = '/'), content_type("application/json"))
      #stop_for_status(r, "connect to STRING API.")

      if (!is.null(r) && r$status_code == 200)
      {
        print_and_log('accessible.')

        string_vs <- fromJSON(content(r))
        if(!is.null(string_vs) && !is.null(string_vs$string_version))
          print_and_log(sprintf('STRING db version: %s', string_vs$string_version), display=FALSE)

      }else
      {
        print_and_log('not accessible.', level='fatal')
        accessible <- FALSE

      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING ping.',cond$message), level='fatal')
      accessible <- FALSE

    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING ping.',cond$message), level='warning')
      accessible <- FALSE

    }
  )

  return(accessible)
}


get.string.network.image <- function(server,geneList,STRING_config_list){

  tryCatch(
    {

      geneCount = length(unlist(strsplit(x=geneList,split = '%0d')))
      query <- paste(sprintf("%s/%s/network?identifiers=%s&species=9606",
                             server,
                             STRING_config_list[['imagetype']],
                             geneList))

      # not used any more because genes are increased with partners
      # if(STRING_config_list[['add_color_nodes']] == 0 )
      # {
      #   query=paste(sprintf("%s&%s=%s",query,'add_color_nodes',geneCount))
      #   print_and_log(sprintf("STRING parameter: add_color_nodes = %s",geneCount),display=FALSE)
      # }
      #
      #
      # if(STRING_config_list[['add_white_nodes']] == 0 )
      # {
      #   query=paste(sprintf("%s&%s=%s",query,'add_white_nodes',geneCount))
      #   print_and_log(sprintf("STRING parameter: add_white_nodes = %s",geneCount),display=FALSE)
      # }
      #
      #
      # if(STRING_config_list[['required_score']] == 0 )
      # {
      #   query=paste(sprintf("%s&%s",query,'required_score=700'))
      #   print_and_log("STRING parameter: required_score = 700",display=FALSE)
      # }
      #
      #
      # if(STRING_config_list[['network_type']] != 'functional' )
      #   query=paste(sprintf("%s&%s=%s",query,'network_type',STRING_config_list[['network_type']]))
      #
      #
      # if(STRING_config_list[['network_flavor']] != 'evidence' )
      #   query=paste(sprintf("%s&%s=%s",query,'network_flavor',STRING_config_list[['network_flavor']]))
      #
      #
      if(STRING_config_list[['hide_node_labels']] != 0 )
        query=paste(sprintf("%s&%s=%s",query,'hide_node_labels',STRING_config_list[['hide_node_labels']]))


      if(STRING_config_list[['hide_disconnected_nodes']] != 0 )
        query=paste(sprintf("%s&%s=%s",query,'hide_disconnected_nodes',STRING_config_list[['hide_disconnected_nodes']]))


      if(STRING_config_list[['show_query_node_labels']] != 0 )
        query=paste(sprintf("%s&%s=%s",query,'show_query_node_labels',STRING_config_list[['show_query_node_labels']]))


      if(STRING_config_list[['block_structure_pics_in_bubbles']] != 0 )
        query=paste(sprintf("%s&%s=%s",query,'block_structure_pics_in_bubbles',STRING_config_list[['block_structure_pics_in_bubbles']]))


      if(STRING_config_list[['center_node_labels']] != 0 )
        query=paste(sprintf("%s&%s=%s",query,'center_node_labels',STRING_config_list[['center_node_labels']]))

      if(STRING_config_list[['custom_label_font_size']] != 12 )
        query=paste(sprintf("%s&%s=%s",query,'custom_label_font_size',STRING_config_list[['custom_label_font_size']]))

      ## ======
      r <- GET(query)
      # message_for_status(r, "get netowrk image from string.")

      if (!is.null(r) && r$status_code == 200)
      {
        return(content(r))
      }else
      {
        print_and_log('error in network image from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING network image function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING network image function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

get.string.enrichment.image <- function(server,geneList,STRING_config_list,enrichmentfigure_category){

  # Format   	  Description
  # image 	        network PNG image with alpha-channel
  # highres_image 	high resolution network PNG image with alpha-channel
  # svg 	          vector graphic format (SVG)


  # Parameter 	Description
  # identifiers 	required parameter for multiple items, e.g. DRD1_HUMAN%0dDRD2_HUMAN
  # species 	NCBI/STRING taxon (e.g. 9606 for human, or STRG0AXXXXX see: STRING organisms).
  # category 	term category (e.g., KEGG, WikiPathways, etc. See the table below for all category keys. Default is Process)
  # group_by_similarity 	threshold for visually grouping related terms on the plot, ranging from 0.1 to 1, in steps of 0.1 (e.g. 0.8), with no grouping applied by default.
  # color_palette 	color palette to represent FDR values (mint_blue, lime_emerald, green_blue, peach_purple, straw_navy, yellow_pink, default is mint_blue).
  # number_of_term_shown 	maximum number of terms displayed on the plot (default is 10).
  # x_axis 	specifies the order of the terms and the variable on the X-axis (signal, strength, FDR, gene_count, default is signal).


  # Category ID 	Description
  # Process 	Biological Process (Gene Ontology)
  # Function 	Molecular Function (Gene Ontology)
  # Component 	Cellular Component (Gene Ontology)
  # Keyword 	Annotated Keywords (UniProt)
  # KEGG 	KEGG Pathways
  # RCTM 	Reactome Pathways
  # HPO 	Human Phenotype (Monarch)
  # MPO 	The Mammalian Phenotype Ontology (Monarch)
  # DPO 	Drosophila Phenotype (Monarch)
  # WPO 	C. elegans Phenotype Ontology (Monarch)
  # ZPO 	Zebrafish Phenotype Ontology (Monarch)
  # FYPO 	Fission Yeast Phenotype Ontology (Monarch)
  # Pfam 	Protein Domains (Pfam)
  # SMART 	Protein Domains (SMART)
  # InterPro 	Protein Domains and Features (InterPro)
  # PMID 	Reference Publications (PubMed)
  # NetworkNeighborAL 	Local Network Cluster (STRING)
  # COMPARTMENTS 	Subcellular Localization (COMPARTMENTS)
  # TISSUES 	Tissue Expression (TISSUES)
  # DISEASES 	Disease-gene Associations (DISEASES)
  # WikiPathways 	WikiPathways

  tryCatch(
    {

      geneCount = length(unlist(strsplit(x=geneList,split = '%0d')))
      query <- paste(sprintf("%s/%s/enrichmentfigure?identifiers=%s&species=9606&category=%s&group_by_similarity=0.8",
                             server,
                             STRING_config_list[['imagetype']],
                             geneList,
                             enrichmentfigure_category))


      ## ======
      r <- GET(query)

      if (!is.null(r) && r$status_code == 200)
      {
        return(content(r))
      }else
      {
        print_and_log('error in enrichment image from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING enrichment image function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING enrichment image function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

get.string.functional.enrichment <- function(server,geneList){

  string.enrich.tb <- NULL

  tryCatch(
    {
      query <- sprintf("%s/tsv/enrichment?identifiers=%s&species=9606",server,geneList)

      r <- GET(query)
      #message_for_status(r, "get functional enrichment analysis from string.")

      if (!is.null(r) && r$status_code == 200)
      {
        string.enrich.tb <- as.data.table(content(r, show_col_types=FALSE))

        if(nrow(string.enrich.tb) > 0)
        {
          if("ncbiTaxonId" %in% names(string.enrich.tb))
            string.enrich.tb <- string.enrich.tb[,-c('ncbiTaxonId')] # remove the taxon id , because it is always 9606

          string_names <- c("category","term","number_of_genes","number_of_genes_in_background",
                            "inputGenes","preferredNames", "p_value","fdr",  "description" )

          if(all(is.element(string_names,names(string.enrich.tb))))
            setcolorder(string.enrich.tb,c("category","term",  "description","p_value","fdr"))

          return(string.enrich.tb)
        }
        else
        {
          print_and_log('No functional enrichment data from STRING DB was found.',display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log('error in functional enrichment analysis from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING functional enrichment analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING functional enrichment analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

get.string.annot.enrichment <- function(server,geneList){

  string.annot.tb <- NULL

  tryCatch(
    {
      query <- sprintf("%s/tsv/functional_annotation?identifiers=%s&species=9606",server,geneList)

      r <- GET(query)
      #message_for_status(r, "get annotation from string.")

      if (!is.null(r) && r$status_code == 200)
      {
        string.annot.tb <- as.data.table(content(r, show_col_types=FALSE))

        if(nrow(string.annot.tb) > 0)
        {
          return(string.annot.tb)
        }
        else
        {
          print_and_log('No functional annotation data from STRING DB was found.',display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log('error in annotation analysis from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING annotation analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING annotation analysis function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

get.string.description <- function(server,geneList){

  string.annot.tb <- NULL

  tryCatch(
    {
      query <- sprintf("%s/tsv/get_string_ids?identifiers=%s&species=9606&echo_query=1",server,geneList)

      r <- GET(query)
      #message_for_status(r, "get annotation from string.")

      if (!is.null(r) && r$status_code == 200)
      {
        string.annot.tb <- as.data.table(content(r, show_col_types=FALSE))

        if(nrow(string.annot.tb) > 0)
        {
          if('stringId' %in% names(string.annot.tb) &
             'preferredName' %in% names(string.annot.tb) &
             'annotation' %in% names(string.annot.tb))
            string.annot.tb <- string.annot.tb[,list(stringId,preferredName,annotation)]

          return(string.annot.tb)
        }
        else
        {
          print_and_log('No description data from STRING DB was found.',display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log('error in description call from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING description function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING description function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}

get.string.webLink <- function(server,geneList){

  webLink <- NULL

  tryCatch(
    {
      query <- sprintf("%s/json/get_link?identifiers=%s&species=9606",server,geneList)

      r <- GET(query)
      #message_for_status(r, "get annotation from string.")

      if (!is.null(r) && r$status_code == 200)
      {
        webLink <- fromJSON(content(r))

        if(!is.null(webLink))
        {
          print_and_log(paste('STRING weblink:',webLink),display=FALSE)
          return(webLink)
        }
        else
        {
          print_and_log('No weblink data was returned from STRING DB.',display=FALSE)
          return(NULL)
        }

      } else {
        print_and_log('error in weblink call from STRING.', level='warning',display=FALSE)
        return(NULL)
      }
    }, error = function(cond) {
      print_and_log(paste('Error occured in STRING weblink function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in STRING weblink function.',cond$message), level='warning',display=FALSE)
      return(NULL)
    }
  )
}


filter.variants.for.STRING <- function(data)
{

  if(!is.element('Gene', names(data)))
    return(NULL)

  if(!is.data.table(data))
    setDT(data)

  geneList <- data[gSNP == Linked_SNP |
                     LD > 0.8 |
                     (Most_severe_consequence == 'missense_variant' & LD > 0.5),]

  #geneList <- geneList[, processed_geneId_col := mapply(keep_nearest_gene_in_row,Gene,GeneId)]
  #geneList <-  unique(gsub(trimws(unlist(sapply(unlist(geneList$GeneId), function(x) strsplit(x,',|/')))), pattern= '(.+)\\(.+)',replacement = '\\1'))

  geneList <- unique(unlist(strsplit(geneList$GeneId,split = ',')))

  if(length(geneList) == 0)
  {
    return(NULL)
  }else{
    print_and_log(paste('Gene list for STRING:',paste(geneList,collapse = ";")), display=FALSE)
    print_and_log(paste('Gene list for STRING (count):',length(geneList)), display=FALSE)
    return(paste(geneList,collapse = "%0d"))
  }

}

