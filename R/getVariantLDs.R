#' Finds variants in high LD
#'
#' This function returns a list of variables that are in high LD with a list of selected variants using data from the Ensembl website.
#'
#' @param rslist A vector of rs numbers.
#' @param file Path to the Excel file for saving search results.
#' @param build Genome build. Either 37 or 38. default: 38
#' @param db The population database for calculating LD scores.
#' This can be found using `Ensembl.Databases()` function. default: 1000GENOMES:phase_3:EUR
#' @param window_size Number of base pairs around the variant for checking LD scores (max = 500kb). default: 500
#' @param r2 The minimum LD threshold for selecting variants around the target SNP. default: 0.8.
#' @return A data table with variant information.
#' @export
#'
findProxy <- function(rslist,file=NULL,build='38',db="1000GENOMES:phase_3:EUR", window_size=500, r2=0.8)
{
  start.time <-  proc.time()

  on.exit({
    if(exists('start.time') && !is.null(start.time))
      message(sprintf("Run time: %s",timetaken(start.time)),appendLF = TRUE)# END LOG
  })

  #==============================================
  # some cleaning
  rslist <- trimws(rslist)
  # remove duplicates
  rslist <- rslist[!duplicated(rslist)]
  # remove those not starting with rs
  rslist <- rslist[startsWith(x = rslist,prefix = 'rs')]
  # count them
  rsCount <- length(rslist)

  if(rsCount==0)
    stop('No variants specified.',call. = FALSE)
  else if(rsCount > 50)
    stop('Too many variants are selected. A maximum of 50 are strongly advised due to Ensembl server policy.',call. = FALSE)

  if(build == '37')
    server = .SNPannotator$ENSEMBL_API_37
  else if(build=='38')
    server = .SNPannotator$ENSEMBL_API_38
  else
    stop('Build parameter should be either 37 or 38.',call. = FALSE)


  if(r2<0 || r2 >1)
    stop('r2 parameter should be between 0-1.',call. = FALSE)


  #===============================================
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = rsCount)
  pb$tick(0)

  output <- data.frame()
  message("Searching for proxy variants ...",appendLF = TRUE)
  message(sprintf("Build: %s, min-r2: %s, window-size: %s, db: %s",build,r2,window_size,db),appendLF = TRUE)

  for(i in 1:rsCount)
  {
    rs= rslist[i]

    # get variants in high LD
    ldlist <- getVariantLDs(rs,server,db, window_size, r2,file.log = FALSE)

    #
    if(!is.null(ldlist))
    {
      # add variant number in the first column
      ldlist$number = i

      # check if  RSID has been changed
      returned_rs <- unique(ldlist$variation1)

      if(returned_rs !="" && returned_rs != rs)
        message(sprintf('Variant name updated from %s to %s', rs,returned_rs),appendLF = TRUE)

      # bind the tables
      l <- list(output,ldlist)
      output <- rbindlist(l,use.names=TRUE,fill=TRUE)
    }

    pb$tick(1)
  }

  #==============================================

  if(!is.null(output) && nrow(output)>0)
  {
    setDT(output)
    output <- output[order(-r2),
                     list(variation2,r2),
                     keyby=list(number,variation1)]
    names(output) <- c('number','variant1','variant2','R2')


    # check the number of found proxies
    message(sprintf("Found %s proxies for %s variants.",nrow(output), rsCount),appendLF = TRUE)

    # save to excel
    if(!is.null(file))
    {
      tryCatch({
        dir_path <- dirname(file)
        if (!dir.exists(dir_path)) {
          message("Directory does not exist. File will not be saved.")
          file <- NULL
        }
      }, error = function(e) {
        message("Error in file path: ", e$message)
        file <- NULL
      })
    }

    if(!is.null(file))
      save_xlsx_report(file, rslist,build,db,r2,window_size,output)

    # return the table
    return(output)
  }else
  {
    message("No proxies were found.",appendLF = TRUE)
    return(NULL)
  }
}

save_xlsx_report = function(file_path, rslist,build,db,r2,window_size,output_tbl)
{
  tbl1 <- data.table(
    'input' = paste(rslist,collapse = ";"),
    'build' = build,
    'db' = db,
    'r2' = r2,
    'window_size' = window_size
  )

  tryCatch({
    wb <- createWorkbook()

    addWorksheet(wb, "Parameters")
    addWorksheet(wb, "Data")

    # Write data tables to the sheets
    writeData(wb, sheet = "Parameters",
              x = sprintf("SNPannotator report v%s",.SNPannotator$script.version ),
              startRow = 1, colNames = FALSE)

    writeData(wb, sheet = "Parameters",
              x = data.frame(DateTime = as.character(Sys.time())),
              startRow = 2, colNames = FALSE)

    writeData(wb, sheet = "Parameters", t(tbl1), startRow = 4,rowNames = TRUE,colNames = FALSE)


    writeData(wb, sheet = "Data", output_tbl)

    # Save the workbook to the file
    saveWorkbook(wb, file=file_path, overwrite = TRUE)

    message('Excel file saved.',appendLF = TRUE)
  },
  error = function(x) message(paste("Error in saving the Excel report:", x$message)))
}


getVariantLDs <- function(rsID,server,db, window_size, r2, file.log=TRUE)
{

  ext <- sprintf("/ld/human/%s/%s?window_size=%s;r2=%s",
                 rsID,
                 db,
                 window_size,
                 r2)
  resp <- fetch(ext,server)
  LDlist <- setDT(resp)

  if(!is.null(LDlist) && nrow(LDlist) > 0)

  {
    LDlist <- LDlist[, lapply(.SD, as.character)]
    non.rs <- LDlist[!startsWith(variation2 ,prefix = "rs"), .N]


    # remove variants that don't start with rs , e.g. esv569117
    if(LDlist[!startsWith(variation2 ,prefix = "rs"), .N] > 0)
    {

      noRS_message <- sprintf("Inappropriate IDs: %s",
                              paste(LDlist[!startsWith(variation2 ,prefix = "rs"),variation2],
                                    collapse ="|"))
      if(file.log == TRUE)
        print_and_log(noRS_message, level='warning')
      else
        message(noRS_message,appendLF = TRUE)
      LDlist <- LDlist[startsWith(variation2 ,prefix = "rs"),]


    }

    LDlist <- LDlist[ ,r2:= as.numeric(r2)]
    LDlist <- LDlist[ ,d_prime:= as.numeric(d_prime)]
    LDlist[order(-r2)]
    return(LDlist)
  }
  else
  {
    # if(file.log == TRUE)
    #   print_and_log("No proxies were found.",LF = TRUE)
    # else
    #   message("\nNo proxies were found.",appendLF = TRUE)

    return(NULL)
  }

}



#' Computes and returns LD values between the given variants.
#'
#' This function returns a data frame of LD values between the given variants in a selected population.
#'
#' @param rsList A vector of rs numbers.
#' @param file Path to the Excel file for saving search results.
#' @param pairwise If TRUE, compute pairwise LD between all elements of a list.
#' If FALSE, computes the LD between first and other elements of the list. default: FALSE
#' @param build Genome build. Either 37 or 38. default: 38
#' @param db The population database for calculating LD scores.
#' This can be found using `Ensembl.Databases()` function. default: "1000GENOMES:phase_3:EUR"
#' @param r2 Only return pairs of variants whose r-squared value is equal to or greater than the value provided. default: 0.1.
#' @return A data table with variant information.
#' @export
#'
findPairwiseLD <- function(rsList,file=NULL, pairwise = FALSE,build=38,db="1000GENOMES:phase_3:EUR", r2=0.1)
{

  start.time <-  proc.time()

  on.exit({
    if(exists('start.time') && !is.null(start.time))
      message(sprintf("Run time: %s",timetaken(start.time)),appendLF = TRUE)# END LOG
  })

  # some cleaning
  rsList <- trimws(rsList)
  # remove duplicates
  rsList <- rsList[!duplicated(rsList)]
  # remove those not starting with rs
  rsList <- rsList[startsWith(x = rsList,prefix = 'rs')]

  if(length(rsList) < 2)
    stop("At least two variants should be provided.",call. = FALSE)
  else if(length(rsList) > 5)
    stop('Too many variants are selected. A maximum of 5 are strongly advised due to Ensembl server policy.',call. = FALSE)

  if(r2<0 || r2 >1)
    stop('r2 parameter should be between 0-1.',call. = FALSE)

  if(build == '37')
    server = .SNPannotator$ENSEMBL_API_37
  else if(build=='38')
    server = .SNPannotator$ENSEMBL_API_38
  else
    stop('Build parameter should be either 37 or 38.',call. = FALSE)

  #=====================================================
  output_ld_list <- list()
  output_ld_tbl <- data.table()

  message("Searching for LD data ...",appendLF = TRUE)
  message(sprintf("Build: %s, db: %s",build,db),appendLF = TRUE)

  if(pairwise == TRUE)
  {
    rsList_pairs <- combn(x = rsList, m =2)
    output_ld_list <- apply(rsList_pairs, 2, function(pair) getPiarwiseLD(pair[1], pair[2],server,db,r2))
  }else{
    output_ld_list <- lapply(2:length(rsList), function(i) getPiarwiseLD(rsList[1], rsList[i],server,db,r2))
  }

  #=====================================================

  if(length(output_ld_list) > 0)
  {
    output_ld_list <- Filter(Negate(is.null), output_ld_list)

    output_ld_tbl <- do.call(rbind,output_ld_list)
    setDT(output_ld_tbl)
    if(is.element("variation1", names(output_ld_tbl)))
      setcolorder(output_ld_tbl,c("variation1", "variation2","r2", "d_prime"))
  }

  #=====================================================
  if(nrow(output_ld_tbl) > 0)
  {

    # save to excel
    if(!is.null(file))
    {
      tryCatch({
        dir_path <- dirname(file)
        if (!dir.exists(dir_path)) {
          message("Directory does not exist. File will not be saved.")
          file <- NULL
        }
      }, error = function(e) {
        message("Error in file path: ", e$message)
        file <- NULL
      })
    }

    # save to excel
    if(!is.null(file))
      save_xlsx_report(file, rsList,build,db,r2,'NA',output_ld_tbl)


    return(output_ld_tbl)
  }
  else
    message("No data available.")

}

getPiarwiseLD <- function(rs1,rs2,server,db,r2)
{

  ext <- sprintf("/ld/human/pairwise/%s/%s?content-type=application/json;population_name=%s&r2=%s",rs1,rs2,db,r2)
  resp <- fetch(ext,server)
  tbl <- as.data.table(resp)

  if(nrow(tbl) > 0)
  {
    if(tbl$variation1 != rs1)
      message(sprintf('Variant name updated from %s to %s', rs1,tbl$variation1),appendLF = TRUE)

    if(tbl$variation2 != rs2)
      message(sprintf('Variant name updated from %s to %s', rs2,tbl$variation2),appendLF = TRUE)
    return(tbl)
  }else
  {
    return(data.table())
  }

}
