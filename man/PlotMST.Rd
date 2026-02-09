\name{PlotMST}
\alias{PlotMST}
\title{Generate a minimum spanning tree from the output of 
makeGraphFromChewBBACA function}
\usage{PlotMST(igraph_object, metadata = NULL, show_clustering = TRUE,
                    MST_edges_color = "black", cluster_edges_color = "grey",
                    MST_edges_width = 0.75, cluster_edges_width = 0.5,
                    edge_label_size = 3, edge_label_dodge = 2,
                    node_color = "red", node_size = 3, node_label_size = 0,
                    show_legend=TRUE, title = NULL)
}
\description{
Using the igraph object generated from makeGraphFromChewBBACA function, 
PlotMST function will plot a minimum spanning tree.
}

\arguments{
    \item{igraph_object}{igraph object, output of the makeGraphFromChewBBACA 
    function}

    \item{metadata}{data frame; This is an optional argument. If given, 
                the function will assume the first column contains the 
                isolates names and will add the remaining columns as 
                metadata to the igraph object, useful for generating 
                custom coloring of the nodes i.e. Year, Source, etc. 
                The default is NULL}

    \item{show_clustering}{logical; If the user has used a value for 
                max_allelic_difference when running makeGraphFromChewBBACA, 
                the user has the option to show the clustering of nodes based 
                on that threshold is set to TRUE, or FALSE to show the 
                minimum spanning tree without any clustering between nodes, as 
                if the user did not apply any clustering. The default is TRUE}

    \item{MST_edges_color}{character; This argument will allow the user to set 
                the color of the edges that connect the nodes. The default is 
                "black"}

    \item{cluster_edges_color}{character; This argument will allow the user to 
                set the color of the edges that connect the clustered nodes. 
                The default is "grey"}

    \item{MST_edges_width}{numeric; This argument will allow the user to set 
                the width of the edges that connect the nodes. 
                The default is 0.75}

    \item{cluster_edges_width}{numeric; This argument will allow the user to set 
                the width of the edges that connect the clustered nodes. 
                The default is 0.5}

    \item{edge_label_size}{numeric; This argument will allow the user to set 
                the label size of the distances on the edges. 
                The default is 3}

    \item{edge_label_dodge}{numeric; This argument will allow the user control 
                the position of the distance on the edges. 
                The default is 2}

    \item{node_color}{character; This argument will allow the user to choose 
                the color of the nodes. The color can be specified by its name 
                or using hexadecimal representation. Also, if the metadata 
                argument is given, the user can specify a column name from the 
                metadata data frame to color the nodes i.e. "Year". 
                The default is "red"}

    \item{node_size}{numeric; This argument will allow the user to control the 
                node size. The default is 3}

    \item{node_label_size}{numeric; This argument will allow the user to 
                control the node label size. If set to zero, no labels will be 
                shown. The default is 0}

    \item{show_legend}{logical; This argument will allow the user to 
                show or hide the legend of the generated plot. 
                The default is TRUE}

    \item{title}{character; This argument will allow the user to 
                a title to the generated plot. 
                The default is NULL}
}

\details{
The function \emph{PlotMST} will make a minimum spannig tree using 
the output of makeGraphFromChewBBACA function. Many arguments are 
available to customize the minimum spanning tree.
}

\examples{
PlotMST(my_graph, show_clustering = TRUE, show_legend=FALSE, 
    MST_edges_color="#b97b29", node_color = "#3b17db", 
    node_label_size = 3, title = "MST")
}

\value{
The returned value is a ggraph object showing the minimum spanning tree.
}

\author{

    Abdullah El-Kurdi

}