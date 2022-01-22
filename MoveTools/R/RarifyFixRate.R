#' Rarify GPS collar point data
#'
#' This function takes a dataset of animal locations and subsets it to the step interval of interest.  It will maximize the number of steps available at the ideal fix rate.  It can deal with complicated collar schedules.Written by Jerod Merkle. Last updated January 2021.
#'
#' @param data An ordered dataframe or sf POINT dataframe with a posix column and a burst or animal id column
#' @param burst_name A character specifying the name of the column representing burst from CalcBurst, or simply animal id.
#' @param date_name  A character specifying the name of the column representing date and time stamps of the locations.
#' @param ideal.fix.interval Your ideal fix rate in hours
#' @param fix.interval.var Variation (in hours) around your ideal fix interval that you are willing to allow.
#' @param maximize.pts Do you want to maxmimize the number of fixes with your idea fix interval? If TRUE, you will also have many fixes that are less than your ideal fix interval.
#' @param max.cpus Maximum number of cpus you want this function to use during the parallelization process. Default is detectCores()-1

#'
#' @return Returns the same data.frame entered with the addition of the following columns: dist (distance between steps in meters), dt (time elapsed bewteen steps in seconds), speed (meters/second), abs.angle (angle moved relative no north, in degrees), rel.angle (angle moved relative to the direction of the previous step, in degrees), and StepFlag (logical, denoting whether a step is connected and has all movement params calculated)
#'
#' @examples
#' # To come

#' @export

RarifyFixRate <- function(data = data,
                          burst_name="burst", # This is AID or a burst column
                          date_name="date",
                          ideal.fix.interval = 3,
                          fix.interval.var = 0.5,
                          maximize.pts=TRUE,
                          max.cpus=detectCores()-1
){

  #manage packages
  if(all(c("sf","parallel") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: sf, parallel")
  require(sf)
  require(parallel)

  #some checks to start
  if(any(colnames(data) == date_name) == FALSE)
    stop(print("Your date_name is not correct."))
  if(any(colnames(data) == burst_name) == FALSE)
    stop(print("Your burst_name is not correct."))
  if("date1234" %in% colnames(data))
    stop("The column named date1234 cannot be in your data!")
  if("id1234" %in% colnames(data))
    stop("The column named id1234 cannot be in your data!")
  orig <- data   # save original database for later
  if(inherits(data, "sf")){
    data <- sf::st_drop_geometry(data)
  }
  if(inherits(data, "tbl_df") == TRUE){  # take out of dumb tibble!
    data <- data.frame(data)
  }
  data <- data[,c(burst_name, date_name)]  # IT WILL ALWASYS BE BURST FIRST COLUMN AND DATE SECOND COLUMN
  colnames(data) <- c("id1234","date1234")  # rename columns
  if(!inherits(data$date1234, "POSIXct"))
    stop(print("date column is not POSIXct"))
  if(any(is.na(data$date1234) == TRUE))
    stop("You have NAs in your date column")
  if(any(is.na(data$id1234) == TRUE))
    stop("You have NAs in your burst column")
  key <- 1:nrow(data)
  key2 <- key[order(data$id1234, data$date1234)]
  if(all(key==key2)==FALSE)
    stop(print("Your data are not ordered correctly"))
  rm(key, key2)
  if(any(duplicated(data[c("id1234", "date1234")])) == TRUE)
    stop("You have duplicate burst_dates in your database!")

  # some check to see how many diffs are actually possibly within that range.
  #perhaps someone specified an ideal fix interval that is too small for the data

  data$date1234 <- as.numeric(data$date1234)/3600  # turn date column into hours

  dt <- diff(data$date1234)  # a quick calculation of the date differences
  flags <- c(diff(data$id1234),1)   # identify when bursts switch
  dt[flags == 1] <- NA  # add NA when burst switches

  if(mean(dt>ideal.fix.interval, na.rm=T) > 0.50)
    print("Warning! More than half of your current fix interval is already larger than your ideal.fix.interval! Thus, this function will not remove many points, if any.")
  rm(dt, flags)

  data$row <- 1:nrow(data)  # add a row column to keep track of rows
  bursts <- unique(data$id1234)  # identify unique bursts

  # identify cores
  no_cores <- ifelse(max.cpus > length(bursts), length(bursts), max.cpus)
  # Setup cluster
  clust <- parallel::makeCluster(no_cores)
  # export the objects you need for your calculations from your environment to each node's environment
  parallel::clusterExport(clust, varlist=c("data","bursts","fix.interval.var",
                                           "ideal.fix.interval","maximize.pts"),envir=environment())

  # now for the actual loop
  rows2keep <- do.call(c, parallel::clusterApplyLB(clust, 1:length(bursts), function(i){

    data.b <- data[data$id1234 == bursts[i],]
    rows.keep.b <- c(1)  # initialize vector of rows to keep for each burst subset

    for(e in 2:nrow(data.b)){ # loop over the rows, starting at 2

      if(max(rows.keep.b) >= e){  # check to see if we skipped previous point
        next
      }

      diffs <- data.b$date1234[e:nrow(data.b)]-(data.b$date1234[max(rows.keep.b)]+ideal.fix.interval)

      goods <- which(diffs <= fix.interval.var & diffs >= -fix.interval.var)

      if(length(goods)==0){  # if there are no points for next step
        if(maximize.pts == TRUE){  # if maximize is true, go to next step
          rows.keep.b <- c(rows.keep.b,e)
        }else{  # if maxmize.pts is FALSE, go to next point largar then max.fix.interval
          rows.keep.b <- c(rows.keep.b,min(which(diffs > fix.interval.var))+e-1)
        }
        next
      }  # end of if length(goods)==0

      if(length(goods)==1){
        rows.keep.b <- c(rows.keep.b,goods+e-1)
        next
      } # end of if length(goods)==1

      # finally, if length of goods is greater than 1
      # figure out the step that is closest to the ideal fix rate.
      best.fix <- which.min(abs(diffs[goods]))
      rows.keep.b <- c(rows.keep.b,goods[best.fix]+e-1)

    }  # end of for loop over the rows

    return(data.b$row[rows.keep.b])  # return the row numbers

  }))
  parallel::stopCluster(clust)   # you must stop the parallelization process


  print(paste0("You removed ", round((1-length(rows2keep)/nrow(orig))*100,1), "% of your data!"))

  dt <- diff(data$date1234[rows2keep])  # a quick calculation of the date differences

  print(paste0("Your new median fix rate is approximately ", round(median(dt, na.rm=T),1), " hours!"))

  return(orig[rows2keep,])  # return only the rows of the origional dataframe

} #end of RarifyFixRate function
