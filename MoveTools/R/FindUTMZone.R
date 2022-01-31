#' Identify UTM zone from longitude and latitude
#'
#' Can be a single lat/long value or a vector of lat/long values. Written by Jerod Merkle with some help from the internet. Last updated January 2022.
#'
#' @param longitude a numeric vector of longitude values
#' @param latitude a numeric vector of latitude values
#'
#' @return Returns numeric value or vector of values representing the UTM zone of your input data.
#'
#' @examples
#' # none

#' @export


FindUTMZone <- function(longitude, latitude) {

  if("numeric" %in% class(longitude) == FALSE)
    stop("longitude is not a numeric value or vector")
  if("numeric" %in% class(latitude) == FALSE)
    stop("latitude is not a numeric value or vector")
  if(length(longitude)!=length(latitude))
    stop("Your longitude and latitude are not the same length")

  if(length(longitude)==1){
    # Special zones for Svalbard and Norway
    if (latitude >= 72.0 && latitude < 84.0 )
      if (longitude >= 0.0  && longitude <  9.0)
        return(31);
    if (longitude >= 9.0  && longitude < 21.0)
      return(33)
    if (longitude >= 21.0 && longitude < 33.0)
      return(35)
    if (longitude >= 33.0 && longitude < 42.0)
      return(37)
    (floor((longitude + 180) / 6) %% 60) + 1
  }else{
    # loop through each value
    return(do.call(c, lapply(1:length(longitude), function(i){
      # Special zones for Svalbard and Norway
      if (latitude[i] >= 72.0 && latitude[i] < 84.0 )
        if (longitude[i] >= 0.0  && longitude[i] <  9.0)
          return(31);
      if (longitude[i] >= 9.0  && longitude[i] < 21.0)
        return(33)
      if (longitude[i] >= 21.0 && longitude[i] < 33.0)
        return(35)
      if (longitude[i] >= 33.0 && longitude[i] < 42.0)
        return(37)
      (floor((longitude[i] + 180) / 6) %% 60) + 1
    })))
  }
}
