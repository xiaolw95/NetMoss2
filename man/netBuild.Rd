% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netBuild.R
\name{netBuild}
\alias{netBuild}
\title{Building microbial networks}
\usage{
netBuild(case_dir, control_dir, method = "sparcc")
}
\arguments{
\item{case_dir}{string.The directory of diseased data set or a single diseased file.}

\item{control_dir}{string.The directory of healthy data set or a single healthy file.}

\item{method}{string. The method chosen to build microbial networks which can be "sparcc" or "pearson".}
}
\value{
microbial networks of case and control data set.
}
\description{
This is a function to build microbial networks using abundance matrix.
}
\examples{
data(testData)
case_dir = mydata[[1]]
control_dir = mydata[[2]]
netBuild(case_dir = case_dir,
     control_dir = control_dir,
     method = "sparcc")


}
