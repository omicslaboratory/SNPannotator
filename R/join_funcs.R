joinResults <- function(filesDir)
{
  
  files <- list.files(path = filesDir,pattern = 'rds$',full.names = TRUE)
  outputFile <- data.table()
  
  for(thisfile in files)
  {
    thisData=readRDS(thisfile)
    
    if(!is.null(thisData[['mainData']]))
    {
      thisData=thisData[['mainData']]
      outputFile <- data.table::rbindlist(list(outputFile,thisData),fill = TRUE)
    }
    
  }
  
  outputFile[,Chr := as.numeric(Chr)]
  outputFile[,Pos := as.numeric(Pos)]
 # outputFile=outputFile[order(outputFile$Chr,outputFile$Pos)]
  varList=unique(outputFile$gSNP)
  varList=data.table(gSNP = varList , number = seq(1,length(varList)))
  
  outputFile = merge.data.table(outputFile,varList,by='gSNP',all.x = TRUE,sort = FALSE)
  outputFile$`#gSNP` <- NULL
  setnames(outputFile,'number','#gSNP')
  setcolorder(outputFile,'#gSNP')
  
  
  openxlsx::write.xlsx(outputFile,file=file.path(filesDir,'mergedData.xlsx'))
  
  invisible(outputFile)
}