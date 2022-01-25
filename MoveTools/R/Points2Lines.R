#' Transform points into lines
#'
#' Can either create a line string for each id, or a line for each connected point in the database. Written and updated by Jerod Merkle. Last updated January 2022.
#'
#' @param data A ordered sf POINT dataframe with a posix column and an animal id column
#' @param id_name A character specifying the name of the column representing animal ID.
#' @param date_name  A character specifying the name of the column representing date and time stamps of the locations.
#' @param byid If TRUE, returns a sf line objective with nrow equal to length(unique(data[,id])). If FALSE, returns a sf objective with nrow equal to nrow(data)
#' @param no_cores How many processor cores would you like to use? e.g., detectCores()-1
#'
#' @return Returns a linestring sf datafame. If byid=TRUE, nrow will be the number of unique ids. If byid==FALSE, nrow of value will be nrow(data)-length(unique(data[,id_name]))
#'
#' @examples
#' #To come
#'
#' @export

Points2Lines <- function(data=data,
                         date_name="date",
                         id_name="id",
                         byid=TRUE,
                         no_cores=4
){

  #manage packages
  if(all(c("sf","parallel") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: sf and parallel")
  require(sf)
  require(parallel)

  # checks
  if("sf" %in% class(data) == FALSE | "sfc_POINT" %in% class(sf::st_geometry(data)) == FALSE)
    stop("data must be a sfc_POINT")
  if(any(is.na(st_drop_geometry(data)[,date_name])))
    stop("You have NAs in your date column!")
  if(any(is.na(st_drop_geometry(data)[,id_name])))
    stop("You have NAs in your id column!")
  if("id1234" %in% names(data))
    stop("You need to remove the column named id1234!")
  if("date1234" %in% names(data))
    stop("You need to remove the column named date1234!")

  # add new columns
  data$id1234 <- sf::st_drop_geometry(data)[,id_name]
  data$date1234 <- sf::st_drop_geometry(data)[,date_name]

  # make sure the data are ordered
  data <- data[order(data$id1234, data$date1234),]

  if(byid == TRUE){

    u <- unique(data$id1234)
    lns <- do.call(c, lapply(1:length(u), function(e){
      return(sf::st_cast(sf::st_combine(data[data$id1234 == u[e],]), "LINESTRING"))
    }))

    lns <- data.frame(mig=u,
                      firstdate=do.call(c, lapply(u, function(e){min(data$date1234[data$id1234==e], na.rm=TRUE)})),
                      lastdate=do.call(c, lapply(u, function(e){max(data$date1234[data$id1234==e], na.rm=TRUE)})),
                      geometry=lns)
    lns <- sf::st_as_sf(lns, sf_column_name = "geometry")
    names(lns)[1] <- id_name    # rename the id column to id_name

  }else{   # when byid = FALSE

    u <- unique(data$id1234)

    # identify cores (use 1 less than you have)
    no_cores <- ifelse(length(u) > no_cores, no_cores, length(u))
    # Setup cluster
    clust <- parallel::makeCluster(no_cores)
    # export the objects you need for your calculations from your environment to each node's environment
    parallel::clusterExport(clust, varlist=c("data","u"),envir=environment())


    lns <- do.call(rbind, parallel::clusterApplyLB(clust, 1:length(u), function(e){

      library(sf)
      dat.id <- data[data$id1234 == u[e],]

      dat.id.xy1 <- sf::st_coordinates(dat.id)
      dat.id.xy2 <- dat.id.xy1[2:nrow(dat.id.xy1),]
      dat.id.xy1 <- dat.id.xy1[1:(nrow(dat.id.xy1)-1),]

      dataL <- lapply(1:nrow(dat.id.xy1), function(i){
        return(sf::st_linestring(rbind(dat.id.xy1[i,],
                                       dat.id.xy2[i,]),dim="XY"))
      })
      dataL <- sf::st_as_sfc(dataL, crs=sf::st_crs(data))
      dataL <- data.frame(sf::st_drop_geometry(dat.id)[1:(nrow(dat.id)-1),], geometry=dataL)
      dataL <- sf::st_as_sf(dataL, sf_column_name = "geometry")
      return(dataL)
    }))
    stopCluster(clust)   # you must stop the parallelization process

    # nrow(data); nrow(lns)+length(u)  # these of course should be exactly the same!
    lns$id1234 <- NULL
    lns$date1234 <- NULL
  }  # end of section for by_id=FALSE
  return(lns)
} # end of function
