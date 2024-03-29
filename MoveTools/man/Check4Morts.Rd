% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Check4Morts.R
\name{Check4Morts}
\alias{Check4Morts}
\title{Function to identify if and when you have mortalies, or periods of time when you animal was not moving, given a dist (in m) and time (hrs) threshold.}
\usage{
Check4Morts(
  data = data,
  dist_thresh = 50,
  time_thresh = 24,
  id_name = "id",
  date_name = "date"
)
}
\arguments{
\item{data}{A dataframe or sf POINT dataframe with a posix column and an animal id column (optional)}

\item{dist_thresh}{Distance threshold (in meters) for how far the animal can move without being considered alive.}

\item{time_thresh}{Time threshold (in hrs) for how long the animal needs to move less than the dist_thresh without being considered alive.}

\item{id_name}{A character specifying the name of the column representing animal ID.}

\item{date_name}{A character specifying the name of the column representing date and time stamps of the locations.}
}
\value{
Returns a datafame denoting all mortlities, including a column for animal ID, start time and end time of the mortality event.
}
\description{
Data must be ordered and have columns dist, id, date, from mov.param(). Written and updated by Jerod Merkle. Last updated January 2021.
}
\examples{
#To come

}
