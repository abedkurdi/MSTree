\name{makeGraphFromChewBBACA}
\alias{makeGraphFromChewBBACA}
\title{Construct an igraph object for minimum spanning tree plotting from 
the output of ChewBBACA pipeline}
\usage{makeGraphFromChewBBACA(chewbbaca_ExtractCgMLST_output, 
max_allelic_difference=-1)
}
\description{
Given a data frame or the full path of the table from ChewBBACA output, 
this function will construct an igraph object that will be used to generate 
a minimum spanning tree (MST) that visualize the distances between the 
bacterial isolates. 
}

\arguments{
    \item{chewbbaca_ExtractCgMLST_output}{Full path of the file name or 
                a data frame}

    \item{max_allelic_difference}{numeric; This is an optional argument, can 
                be used to connect nodes that have a distance less than or 
                equal to the specified number, useful to show clusters of 
                isolates. The default is set to -1 which means no clustering 
                is applied}
}

\details{
The function \emph{makeGraphFromChewBBACA} will make an igraph object using 
the ChewBBACA ExtractCgMLST output. The generated object will be used in a 
second step to plot the minimum spanning tree.

\emph{chewbbaca_ExtractCgMLST_output} is a full path of the file or a 
data frame if the user has loaded the file.

\emph{max_allelic_difference} allows the user to choose whether to connect 
nodes with a maximum distance. Nodes that have a distance equal or less the 
indicated number, they will form a cluster. Default is -1 which means 
no clustering is applied.

}

\examples{
my_graph <- makeGraphFromChewBBACA("inst/extdata/cgMLST95.csv", 
    max_allelic_difference=-1)
}

\value{
The returned value is a igraph object containing the nodes and edges between nodes.
}

\author{

    Abdullah El-Kurdi

}