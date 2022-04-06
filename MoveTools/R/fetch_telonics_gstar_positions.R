#' Batch download and process Telonics Globalstar GPS collar data
#'
#' Written by Jerod Merkle. Last updated March 2022.
#'
#' @param username A character specifying the telonics globalstar username
#' @param password A character specifying the telonics globalstar password
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

fetch_telonics_gstar_positions <- function(
  username = "username",
  password = "password",
  fldr_out = "C:\\Users\\jmerkle\\Desktop\\teltest",   # must be \\ slashes
  TDC_path = "C:\\Program Files (x86)\\Telonics\\Data Converter\\TDC.exe",
  keep.reports = FALSE,
  start_date = Sys.time()-10*24*3600   # just grab the previous 10 days of data. If NULL, then will return all data
){

  if(all(c("processx","sf") %in% installed.packages()[, 1]) == FALSE)
    stop("You must install the following packages: processx and sf")

  # get packages aligned
  require(processx)
  require(sf)

  # some tests
  if(length(dir(fldr_out))>0)
    stop("Your fldr_out must be empty to proceed!")

  # create a reports folder
  fldr_reports <- paste0(fldr_out, "\\reports")
  dir.create(fldr_reports)

  # create the xml file
  txt <- paste("<BatchSettingsV2>",
               "\t<Globalstar>",
               paste0("\t\t<Username>",username,"</Username>"),
               paste0("\t\t<Password>",password,"</Password>"),
               "\t</Globalstar>",
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
  # Create a list of the data you will import
  fls <- list.files(fldr_reports, ".csv$")

  ## Run a loop that goes over the list, cleans and merges the data
  # Create an empty data frame where all the individuals will be merged in
  fixes <- do.call(rbind, lapply(1:length(fls), function(i){
    # The skip parameter is because there is some meta information above where the recordings begin
    df.i = read.csv(paste0(fldr_reports,"/",fls[i]), skip = 22, header = TRUE)

    # Get the ID and add it as a column (I am using the name the file is saved under and extracting the
    # component that will match with the way it is saved in my meta data column)
    df.i$CollarSerialNumber <- substr(fls[i], 1, 7)

    # Isolate the cases with a successful fix attempt
    df.i <- df.i[which(df.i$GPS.Fix.Attempt=="Succeeded"),]

    print(paste0(nrow(df.i), " rows of data importanted for individual: ", substr(fls[i], 1, 7)))

    # Work on the DateTime stamp. It is a character string so I will first convert it to POSIXct
    # I always try to avoid deleting raw data (you never know when you will need it) so I will create a new DateTime Column
    df.i$GPS.Fix.Time = as.POSIXct(df.i$GPS.Fix.Time, format="%Y.%m.%d %H:%M:%S", tz = "UTC")


    if(is.null(start_date)==FALSE){
      # reduce to specified start date
      df.i <- df.i[df.i$GPS.Fix.Time >= start_date,]
    }

    return(df.i)
  }))

  # order by serial number and then by date
  fixes <- fixes[order(fixes$CollarSerialNumber, fixes$GPS.Fix.Time),]

  # remove any rows with NAs in fix category
  fixes <- fixes[is.na(fixes$GPS.Longitude)==FALSE,]

  # convert data to an sf object
  fixes <- sf::st_as_sf(fixes,
                        coords = c("GPS.Longitude", "GPS.Latitude"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # remove all the temporary files if keep.reports = FALSE
  if(keep.reports==FALSE){
    fls = dir(fldr_out, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    unlink(fls, force=TRUE, recursive = TRUE)
  }

  return(fixes)  # return the sf dataframe

}  # end of function

