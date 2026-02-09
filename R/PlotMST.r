##' @title MSTree

PlotMST <- function(igraph_object, metadata = NULL, show_clustering = TRUE,
                    MST_edges_color = "black", cluster_edges_color = "grey",
                    MST_edges_width = 0.75, cluster_edges_width = 0.5,
                    edge_label_size = 3, edge_label_dodge = 2,
                    node_color = "red", node_size = 3, node_label_size = 0,
                    show_legend=TRUE, title = NULL){

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

        scale_edge_width_manual(values = c("TRUE" = MST_edges_with, 
                                            "FALSE" = cluster_edges_width),
                          labels = c("TRUE" = "MST", "FALSE" = "Threshold"),
                          name = "Edge Type")+

        geom_node_point(color = node_color, size = node_size) +
        
        ggtitle(title) +
        
        theme_graph()

    } else if(!is.null(metadata) && class(metadata) == "data.frame"){
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
            label_dodge = unit(edge_label_dodge, 'mm'), show.legend=show_legend) +

        scale_edge_color_manual(values = c("TRUE" = MST_edges_color,
                                            "FALSE" = cluster_edges_color),
                                labels = c("TRUE" = "MST", 
                                            "FALSE" = "Threshold"),
                                name = "Edge Type") +

        scale_edge_width_manual(values = c("TRUE" = MST_edges_with, 
                                            "FALSE" = cluster_edges_width),
                          labels = c("TRUE" = "MST", "FALSE" = "Threshold"),
                          name = "Edge Type")+

        ggtitle(title) +
        
        theme_graph()

        if(node_color %notin% colnames(metadata)[2:ncol(metadata)]){
            message("\nset node_color parameter from the metadata\n")
            stop("the provided colum name does not exist")
        } else {
            p <- p + geom_node_point(aes(color = .data[[node_color]]), size = node_size)
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

    return(p)
}
