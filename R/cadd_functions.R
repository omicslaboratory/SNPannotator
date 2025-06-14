convert.CADD <- function(cadd)
{
  tryCatch({

    cadd =ifelse(grepl(pattern = "NA", x=cadd), # maybe NA is returned as the output
                 NA,
                 gsub(pattern = '(.+=\\s)(.+)',replacement = '\\2',x = cadd))

    cadd = suppressWarnings(as.numeric(cadd))

    if(is.na(cadd))
      return('')
    else
      return(sprintf('top %#.1f%%',(10^(-cadd/10))*100))
  },error=function(cond){
    return(cadd)
  })

}


select.CADD.scores <- function(alleles,cadd.scores)
{
  tryCatch(
    {
      if (alleles == "" | !grepl(x = cadd.scores, pattern = ','))
        return(cadd.scores)

      alleles <- unlist(strsplit(alleles,'/'))
      cadd.scores <- unlist(strsplit(cadd.scores,' = '))
      cadd.scores1 <- unlist(strsplit(cadd.scores[1],','))
      cadd.scores2 <- unlist(strsplit(cadd.scores[2],','))

      i <- which(is.element(cadd.scores1,alleles))[1]

      return(paste(paste(cadd.scores1[i],collapse = ','),'=',paste(cadd.scores2[i],collapse = ',')))
    },
    error = function(cond) {
      return(cadd.scores)
    }
  )

}
