% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetMoss.R
\name{NetMoss}
\alias{NetMoss}
\title{Identifying Differential Bacteria Based on a Network Workflow.}
\usage{
NetMoss(
  case_dir,
  control_dir,
  net_case_dir,
  net_control_dir,
  scaled = TRUE,
  deepSplit = 4,
  minModuleSize = 20,
  soft_pow = 5
)
}
\arguments{
\item{case_dir}{string.The directory of diseased data set or a single diseased file.}

\item{control_dir}{string.The directory of healthy data set or a single healthy file.}

\item{net_case_dir}{string.The directory of network correlation of diseased data set or a single diseased network file.}

\item{net_control_dir}{string.The directory of network correlation of healthy data set or a single healthy network file.}

\item{scaled}{logical.If TURE then the NetMoss score are scaled.}

\item{deepSplit}{numerical. For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method, a finer control can be achieved via maxCoreScatter and minGap below.}

\item{minModuleSize}{numerical. Minimum module size. Set to 20 by default.}

\item{soft_pow}{numerical. Soft power for module division. Set to 5 by default.}
}
\value{
The output of the NetMoss function is a list contains NetMoss score and integrated results, which could be used in the downstream analysis.
The first list: nodes_result[\link{1}] is the NetMoss score and p value of each taxon;
The second and the third list: nodes_result[\link{2}] and nodes_result[\link{3}] are integrated networks constructed from case and control data;
The 4th and thd 5th list: nodes_result[\link{4}] and nodes_result[\link{5}] are case and control data from the input.
}
\description{
This is a function to identify differential bacteria between case and control data set.
}
\examples{
data(testData)
nodes_result = NetMoss(case_dir = mydata[[1]],
     control_dir = mydata[[2]],
     net_case_dir = mydata[[3]],
     net_control_dir = mydata[[4]])

}
