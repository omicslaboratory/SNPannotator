#' @import data.table
NULL

#' @import httr
NULL

#' @import jsonlite
NULL

#' @import tools
NULL

#' @import xml2
NULL

#' @import openxlsx
NULL

#' @import progress
NULL

#' @import ggplot2
NULL

#' @import futile.logger
NULL

#' @import methods
NULL

#' @import rmarkdown
NULL

#' @import ini
NULL

#' @import igraph
NULL

#' @import ggraph
NULL

#' @import png
NULL

#' @importFrom kableExtra %>%
NULL

#' @importFrom kableExtra column_spec
NULL

#' @importFrom kableExtra kable_styling
NULL

#' @importFrom kableExtra footnote
NULL


.onAttach <-
  function(libname, pkgname) {

    welcomeString = sprintf("initializing SNPannotator package [v%s]",
                            .SNPannotator$script.version)

    packageStartupMessage(welcomeString)

  }


.onLoad <- function(libname, pkgname) {

}


utils::globalVariables(c('r','variation2','variation3','syn_id','rs_id','checkSynonyms',
                         'seq_region_name', '.',
                         'variation1','chr','start','Chr',
                         'end','on.gene','type','dist1',
                         'dist2','i','population','feature_type',
                         'external_name','gene_id','id','d_prime',
                         'name','size','description','r2','Gene',
                         'Pos','Cytoband','gSNP','Linked_SNP','number','Ref_Allele',
                         'Alt_Allele','Deleteriousness','#gSNP', '..desired_cols', '..existing_cols_clinvar',
                         '..existing_cols_gwascat', 'Associations', 'GeneId', 'Imports', 'LD',
                         'Most_severe_consequence', 'Package', 'Phenotype', 'Tissue', 'V1', 'Version',
                         'alleles', 'alt', 'annotation', 'b37VariantId', 'b38_id', 'cluster',
                         'consequence_terms', 'cytoband_grch37', 'cytoband_grch38', 'distance',
                         'eQTL_group', 'from', 'geneSymbol', 'gene_names_grch37', 'gene_names_grch38',
                         'gene_symbol', 'item', 'pValue', 'phenotype', 'position_b37', 'preferredName',
                         'pvalue', 'rsId', 'rsid', 'se', 'server', 'snpID', 'stringId', 'variant_allele'))


