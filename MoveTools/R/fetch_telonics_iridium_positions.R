#' Batch download and process Telonics Iridium GPS collar data
#'
#' Written by Jerod Merkle. Last updated August 2023.
#'
#' @param username A character specifying the telonics username
#' @param password A character specifying the telonics password
#' @param fldr_out  An empty folder for writing temporary files. Must use \\ slashes to denote. And the folder must be empty!
#' @param TDC_path Character vector of the file pathway to TDC.exe. Must use \\ slashes and must include the TDC.exe at the end of the string. Cand download the TDC here: https://www.telonics.com/software/tdc.php
#' @param keep.reports If FALSE (the default), all temporary files and reports will be removed from your temp folder
#' @param start_date a posixct object denoting the earliest data that should be returned. If NULL, then ALL data are returned. Default is Sys.time()-10*24*3600
#'
#' @return Returns an sf point dataframe of the relocations. Collars are denoted by CollarSerialNumber, and the GPS.Fix.Time column represents the fix date/time in UTC.
#'
#' @examples
#' # To come

#' @export

fetch_telonics_iridium_positions <- function(
    username = "username",
    password = "password",
    fldr_out = "C:/Users/lolson10/Documents/SourceCode/TelonicsTemp",   # must be \\ slashes
    TDC_path = "C:\\Program Files (x86)\\Telonics\\Data Converter\\TDC.exe",
    keep.reports = FALSE,
    start_date = Sys.time()-25*24*3600  # just grab the previous 10 days of data. If NULL, then will return all data
    #used to be 10 instead of 35  
){
  
  if(all(c("processx","sf","readr","dplyr") %in% installed.packages()[, 1]) == FALSE)
    stop("You must install the following packages: processx and sf")
  
  # get packages aligned
  require(processx)
  require(sf)
  require(dplyr)
  require(readr)
  
  # some tests
  if(length(dir(fldr_out))>0)
    stop("Your fldr_out must be empty to proceed!")
  
  # create a reports folder
  fldr_reports <- paste0(fldr_out, "\\reports")
  dir.create(fldr_reports)
  
  # create the xml file
  txt <- paste("<BatchSettingsV2>",
               "\t<Iridium>",
               paste0("\t\t<Username>",username,"</Username>"),
               paste0("\t\t<Password>",password,"</Password>"),
               "\t</Iridium>",
               "\t<DownloadData>true</DownloadData>",
               "\t<ConvertAllData />",
               paste0("\t<BatchLog>",fldr_out,"\\BatchLog.txt</BatchLog>"),
               paste0("\t<OutputFolder>",fldr_reports,"</OutputFolder>"),
               "<ReportFileMode>overwrite</ReportFileMode>",
               "</BatchSettingsV2>",
               sep="\n")
  
  Batch_path <- paste0(fldr_out, "\\TelonicsBatchFile.xml")
  # save the xml file
  cat(txt, file=Batch_path)
  
  print("Downloading data from Telonics")
  Batch_path = paste0("/batch:", Batch_path)  # create new batch path for processx
  processx::run(TDC_path, Batch_path)  # TDC should be closed on your computer
  
  print("Importing CSV files")
  #Import the csv files from batch ####

  fixes<-readr::read_csv(list.files(fldr_reports, "Condensed.csv$", full.names = TRUE), 
                         col_select = c("Acquisition Time", 
                                        "GPS Fix Attempt", 
                                        "GPS Latitude", 
                                        "GPS Longitude"),
                         id = 'CollarSerialNumber', 
                         skip = 22, 
                         show_col_types = FALSE,
                         ) %>% 
    dplyr::mutate(CollarSerialNumber = substr(CollarSerialNumber, nchar(CollarSerialNumber)-20, nchar(CollarSerialNumber)-14)) %>% 
    dplyr::filter(`GPS Fix Attempt` == 'Succeeded') %>% 
    dplyr::mutate(`Acquisition Time` = as.POSIXct(`Acquisition Time`, format="%Y.%m.%d %H:%M:%S", tz = "UTC"))

    if(is.null(start_date)==FALSE){
      # reduce to specified start date
      fixes <- fixes %>% 
        filter(`Acquisition Time` >= start_date)
    }
  
  # order by serial number and then by date
  fixes <- fixes[order(fixes$CollarSerialNumber, fixes$`Acquisition Time`),]
  
  # convert data to an sf object
  fixes <- sf::st_as_sf(fixes,
                        coords = c("GPS Longitude", "GPS Latitude"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # remove all the temporary files if keep.reports = FALSE
  if(keep.reports==FALSE){
    fls = dir(fldr_out, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    unlink(fls, force=TRUE, recursive = TRUE)
  }
  
  return(fixes)  # return the sf dataframe
  
}  # end of function
