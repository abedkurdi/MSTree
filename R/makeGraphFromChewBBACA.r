##' @title MSTree

##' @export makeGraphFromChewBBACA

makeGraphFromChewBBACA <- function(chewbbaca_ExtractCgMLST_output, 
            max_allelic_difference = -1){
    
    # check argument
    if (missing(chewbbaca_ExtractCgMLST_output)) {
        stop("A path or dataframe of chewBBACA ExtractCgMLST is needed!")
    }

    # check input formats
    ## path or data frame
    if(class(chewbbaca_ExtractCgMLST_output) == "character"){
        alleles <- read.csv(chewbbaca_ExtractCgMLST_output, row.names = 1, 
                    check.names = FALSE, sep='\t')
    } else if(class(chewbbaca_ExtractCgMLST_output) == "data.frame") {
        alleles <- chewbbaca_ExtractCgMLST_output
    } else {
        stop("A path or dataframe of chewBBACA ExtractCgMLST is required!")
    }

    ## max_allelic_difference type
    if(!(is.numeric(max_allelic_difference) && 
            length(max_allelic_difference) == 1 && 
            max_allelic_difference == round(max_allelic_difference))){
        stop("max_allelic_difference should be an integer!")
    } else if(max_allelic_difference < -1){
        max_allelic_difference <- -1
        message("max_allelic_difference is < -1 - it is set to default (-1)")
    }

    # Replace any specific missing value codes (like -1) with NA
    alleles[alleles == -1] <- NA

    # Initialize empty matrix for pairwise allele differences
    n <- nrow(alleles)
    allele_diff_matrix <- matrix(0, n, n)
    rownames(allele_diff_matrix) <- rownames(alleles)
    colnames(allele_diff_matrix) <- rownames(alleles)

    # Compute pairwise allele differences manually
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            a <- as.numeric(alleles[i, ]) # extract the alleles for sample i
            b <- as.numeric(alleles[j, ]) # extract the alleles for sample j
            valid <- !is.na(a) & !is.na(b) # check if we have NAs
            diff_count <- sum(a[valid] != b[valid]) # calc difference
            allele_diff_matrix[i, j] <- diff_count # save the number
            allele_diff_matrix[j, i] <- diff_count # save the number
        }
    }

    # Create full graph and compute MST
    g_full <- graph_from_adjacency_matrix(allele_diff_matrix, 
                mode = "undirected", weighted = TRUE)
    mst_g <- mst(g_full)

    # Get MST edges
    mst_edges <- as_edgelist(mst_g)
    mst_edge_set <- apply(mst_edges, 1, 
                    function(x) paste(sort(x), collapse = "-"))

    # Now add edges from full graph where difference <= threshold
    edge_list <- data.frame(
        from = character(),
        to = character(),
        weight = numeric(),
        is_mst = logical(),
        stringsAsFactors = FALSE
    )

    # Add all MST edges first
    for (i in 1:nrow(mst_edges)) {
        v1 <- mst_edges[i, 1]
        v2 <- mst_edges[i, 2]
        v1_idx <- which(rownames(alleles) == v1)
        v2_idx <- which(rownames(alleles) == v2)
    
        edge_list <- rbind(edge_list, data.frame(
            from = v1,
            to = v2,
            weight = allele_diff_matrix[v1_idx, v2_idx],
            is_mst = TRUE,
            stringsAsFactors = FALSE
        ))
    }

    # SET YOUR THRESHOLD - # acineto: 9, pseudo: 12, KLB: 15, ECOL: 10
    max_difference <- max_allelic_difference

    # Add additional edges where difference <= threshold 
    if(max_difference > -1){
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                diff_val <- allele_diff_matrix[i, j]
                edge_id <- paste(sort(c(rownames(alleles)[i], 
                                rownames(alleles)[j])), collapse = "-")
    
    ## Only add if within threshold AND not already in MST
    if (diff_val <= max_difference && !(edge_id %in% mst_edge_set)) {
        edge_list <- rbind(edge_list, data.frame(
            from = rownames(alleles)[i],
            to = rownames(alleles)[j],
            weight = diff_val,
            is_mst = FALSE,
            stringsAsFactors = FALSE
                ))
            }
        }
    }
    }

    cat("Total edges:", nrow(edge_list), "\n")
    cat("MST edges:", sum(edge_list$is_mst), "\n")
    cat("Additional edges (≤",max_difference,"):",sum(!edge_list$is_mst),"\n")

    # Create final graph
    g <- graph_from_data_frame(edge_list, directed = FALSE, 
            vertices = rownames(alleles))

    
    # Label edges with allele differences
    E(g)$allele_diff <- E(g)$weight
    E(g)$is_mst <- edge_list$is_mst
    V(g)$name <- rownames(alleles)

    return(g)
}
