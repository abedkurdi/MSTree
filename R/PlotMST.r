#' @name PlotMST
#' @title MSTree R package
#' 
#' @description Using the igraph object generated from makeGraphFromChewBBACA 
#' function, PlotMST function will plot a minimum spanning tree.
#' 
#' @details The function \emph{PlotMST} will make a minimum spannig tree using 
#' the output of makeGraphFromChewBBACA function. Many arguments are 
#' available to customize the minimum spanning tree.
#' 
#' @param igraph_object igraph object, output of the makeGraphFromChewBBACA 
#' function
#' @param metadata data frame; This is an optional argument. If given, 
#'            the function will assume the first column contains the 
#'            isolates names and will add the remaining columns as 
#'            metadata to the igraph object, useful for generating 
#'            custom coloring of the nodes i.e. Year, Source, etc. 
#'            The default is NULL
#' @param show_clustering logical; If the user has used a value for 
#'            max_allelic_difference when running makeGraphFromChewBBACA, 
#'            the user has the option to show the clustering of nodes based 
#'            on that threshold if set to TRUE, or FALSE to show the 
#'            minimum spanning tree without any clustering between nodes, as 
#'            if the user did not apply any clustering. The default is TRUE
#' @param MST_edges_color character; This argument will allow the user to set 
#'            the color of the edges that connect the nodes. The default is 
#'            "black" 
#' @param cluster_edges_color character; This argument will allow the user to 
#'            set the color of the edges that connect the clustered nodes. 
#'            The default is "grey"
#' @param MST_edges_width numeric; This argument will allow the user to set 
#'            the width of the edges that connect the nodes. 
#'            The default is 0.75 
#' @param cluster_edges_width numeric; This argument will allow the user to 
#'            set the width of the edges that connect the clustered nodes. 
#'            The default is 0.5 
#' @param edge_label_size numeric; This argument will allow the user to set 
#'            the label size of the distances on the edges. 
#'            The default is 3
#' @param edge_label_dodge numeric; This argument will allow the user control 
#'            the position of the distance on the edges. 
#'            The default is 2 
#' @param node_color character; This argument will allow the user to choose 
#'            the color of the nodes. The color can be specified by its name 
#'            or using hexadecimal representation. Also, if the metadata 
#'            argument is given, the user can specify a column name from the 
#'            metadata data frame to color the nodes i.e. "Year". 
#'            The default is "red" 
#' @param node_size numeric; This argument will allow the user to control the 
#'            node size. The default is 3 
#' @param node_label_size numeric; This argument will allow the user to 
#'            control the node label size. If set to zero, no labels will be 
#'            shown. The default is 0 
#' @param show_legend logical; This argument will allow the user to 
#'            show or hide the legend of the generated plot. 
#'            The default is TRUE
#' @param title character; This argument will allow the user to 
#'            add a title to the generated plot. 
#'            The default is NULL
#' @param basic logical; This argument will allow the user to 
#'            generate a very basic minimum spanning tree layout using the 
#'            the package NetPathMiner. 
#'            The default is FALSE

#' @return Graph object.
#' @importFrom ggraph scale_edge_color_manual
#' @importFrom ggraph geom_edge_link
#' @importFrom igraph E
#' @importFrom igraph E<-
#' @importFrom igraph subgraph.edges
#' @importFrom igraph set_vertex_attr
#' @importFrom ggraph ggraph
#' @importFrom ggraph scale_edge_width_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom ggraph theme_graph
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom NetPathMiner plotNetwork
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 unit
#' @importFrom methods is
#' 
#' @examples
#' my_graph <- makeGraphFromChewBBACA(system.file("extdata", "cgMLST95.csv", 
#'    package = "MSTree"), max_allelic_difference=-1)
#' mst <- PlotMST(my_graph, show_clustering = TRUE, show_legend=FALSE, 
#'    MST_edges_color="#b97b29", node_color = "#3b17db", 
#'    node_label_size = 3, title = "MST")
#' 
#' @export PlotMST
globalVariables(c(".data"))

PlotMST <- function(igraph_object, metadata = NULL, show_clustering = TRUE,
                    MST_edges_color = "black", cluster_edges_color = "grey",
                    MST_edges_width = 0.75, cluster_edges_width = 0.5,
                    edge_label_size = 3, edge_label_dodge = 2,
                    node_color = "red", node_size = 3, node_label_size = 0,
                    show_legend=TRUE, title = NULL, basic = FALSE){

    ready_graph <- igraph_object

    `%notin%` <- Negate(`%in%`)

    # either show the clustering of nodes or not
    if(!show_clustering) {
        selected_eids <- which(E(ready_graph)$is_mst == TRUE)
        ready_graph <- subgraph.edges(ready_graph, eids = selected_eids, 
            delete.vertices = TRUE)
    }

    # check if metadata data frame is supplied or not - and color accordingly
    if(is.null(metadata)){
        p <- ggraph(ready_graph, layout = "stress") +

        geom_edge_link(aes(label = allele_diff,
        color = is_mst,
        width = is_mst),
        label_size = edge_label_size, 
        angle_calc = 'along', 
        label_dodge = unit(edge_label_dodge, 'mm'), show.legend=show_legend) +

        scale_edge_color_manual(values = c("TRUE" = MST_edges_color,
                                            "FALSE" = cluster_edges_color),
                                labels = c("TRUE" = "MST", 
                                            "FALSE" = "Threshold"),
                                name = "Edge Type") +

        scale_edge_width_manual(values = c("TRUE" = MST_edges_width, 
                                            "FALSE" = cluster_edges_width),
                        labels = c("TRUE" = "MST", "FALSE" = "Threshold"),
                        name = "Edge Type")+

        geom_node_point(color = node_color, size = node_size) +
        
        ggtitle(title) +
        
        theme_graph()

    } else if(!is.null(metadata) && is(metadata,"data.frame")){
        metadata <- metadata[match(V(ready_graph)$name, metadata[,1]),]

        for(i in 2:ncol(metadata)){
            ready_graph <- set_vertex_attr(ready_graph, 
                colnames(metadata)[i], 
                value = as.character(metadata[[colnames(metadata)[i]]]))
        }

        p <- ggraph(ready_graph, layout = "stress") +

        geom_edge_link(aes(label = allele_diff,
            color = is_mst,
            width = is_mst),
            label_size = edge_label_size, 
            angle_calc = 'along', 
            label_dodge = unit(edge_label_dodge,'mm'),show.legend=show_legend)+

        scale_edge_color_manual(values = c("TRUE" = MST_edges_color,
                                            "FALSE" = cluster_edges_color),
                                labels = c("TRUE" = "MST", 
                                            "FALSE" = "Threshold"),
                                name = "Edge Type") +

        scale_edge_width_manual(values = c("TRUE" = MST_edges_width, 
                                            "FALSE" = cluster_edges_width),
                        labels = c("TRUE" = "MST", "FALSE" = "Threshold"),
                        name = "Edge Type")+

        ggtitle(title) +
        
        theme_graph()

        if(node_color %notin% colnames(metadata)[2:ncol(metadata)]){
            message("\nset node_color parameter from the metadata\n")
            stop("the provided colum name does not exist")
        } else {
            p <- p + geom_node_point(aes(color = .data[[node_color]]), 
                size = node_size)
        }       
    }       
        
    if(node_label_size != 0 && node_label_size > 0){
        p <- p + geom_node_text(aes(label = name),
                repel = TRUE,
                size = node_label_size,
                box.padding = unit(1.5, "lines"),
                point.padding = unit(1.5, "lines"),
                segment.size = 0.5,
                segment.color = "grey50")
    }

    if(basic){
        p <- plotNetwork(ready_graph)
    }

    return(p)
}