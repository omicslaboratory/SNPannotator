# load gene name file

# gene data file had to be put outside the package due to size limit from CRAN
getGeneFile <- function(server)
{

  grch37 <- grepl(pattern = '37',x = server)

  if(grch37)
    gene.File <- system.file("extdata", "Gene_Names_Ensembl_104_GRCh37.rds", package = "SNPannotator")
  else
    gene.File <- system.file("extdata", "Gene_Names_Ensembl_104_GRCh38.rds", package = "SNPannotator")


  if (!file.exists(gene.File))
    stop('Gene file not found in the package!', call. = FALSE)
  else
    return(readRDS(gene.File))
}

getGeneFile_v2 <- function(build)
{
  g.File <- ''

  if(build == 'grch37')
    g.File <- system.file("extdata", "gene_names_grch37.rda", package = "SNPannotator")
  else if(build == 'grch38')
    g.File <- system.file("extdata", "gene_names_grch38.rda", package = "SNPannotator")

  if (!file.exists(g.File))
    print_and_log('Cytoband information file not found in the package!', level='fatal')
  else
    load(g.File)

  if(build == 'grch37')
    return(gene_names_grch37)
  else if(build == 'grch38')
    return(gene_names_grch38)
}

getGeneFile.external <- function(geneFilePath,build)
{

  gene.File <- NULL
  build <- toupper(build)

  if(file.exists(geneFilePath))
    gene.File <- tryCatch({
      gene.File <- readRDS(geneFilePath)
      gene.File <- gene.File[[build]]
    },
    error = function(er){
      print_and_log(paste0('Error reading file. (',geneFilePath,')',er$message),level='warning')
    }
    )
  else
    print_and_log(paste0('File not found! (',geneFilePath,')'),level='fatal')


  return(gene.File)

}

# load regulatory file
# regulatory data file had to be put outside the package due to size limit from CRAN
getRegulatoryFile <- function(server)
{

  reg.File <- system.file("extdata", "homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.rds", package = "SNPannotator")


  if (!file.exists(reg.File))
    print_and_log('Regulatory file not found in the package!', level='fatal')
  else
    return(readRDS(reg.File))
}

getRegulatoryFile.external <- function(regFilePath)
{
  reg.File <- NULL

  if(file.exists(regFilePath))
    reg.File <- tryCatch({
      readRDS(regFilePath)
    },
    error = function(er){
      print_and_log(paste0('Error reading file. (',regFilePath,')',er$message),level='warning')

    }
    )
  else
    print_and_log(paste0('File not found! (',regFilePath,')'),level='warning')



  return(reg.File)
}

# load cytoband file
getCytobandFile <- function(build)
{
  c.File <- ''

  if(build == 'grch37')
    c.File <- system.file("extdata", "cytoband_GRCh37.rds", package = "SNPannotator")
  else if(build == 'grch38')
    c.File <- system.file("extdata", "cytoband_GRCh38.rds", package = "SNPannotator")

  if (!file.exists(c.File))
    print_and_log('Cytoband information file not found in the package!', level='fatal')
  else
    return(readRDS(c.File))
}

getCytobandFile_v2 <- function(build)
{
  c.File <- ''

  if(build == 'grch37')
    c.File <- system.file("extdata", "cytoband_GRCh37.rda", package = "SNPannotator")
  else if(build == 'grch38')
    c.File <- system.file("extdata", "cytoband_GRCh38.rda", package = "SNPannotator")

  if (!file.exists(c.File))
    print_and_log('Cytoband information file not found in the package!', level='fatal')
  else
    load(c.File)

  if(build == 'grch37')
    return(cytoband_grch37)
  else if(build == 'grch38')
    return(cytoband_grch38)
}
