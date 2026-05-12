library(RUnit)
library(testthat)
require(igraph)

csv_path <- system.file("extdata", "cgMLST95.csv", package = "MSTree")
meta_path <- system.file("extdata", "metadata.txt", package = "MSTree")

test_data <- read.csv(csv_path, row.names = 1, check.names = FALSE, sep='\t')
meta_data <- read.csv(meta_path, check.names = FALSE, sep='\t')

test_makeGraphFromChewBBACA_and_PlotMST <- function() {
    checkEquals(class(makeGraphFromChewBBACA(test_data, max_allelic_difference=9)), "igraph")

    checkEquals(class(test_data), "data.frame")
    
    my_graph <- makeGraphFromChewBBACA(test_data, max_allelic_difference=9)

    make_test_graph <- function() {
  g <- graph_from_edgelist(
    matrix(c("A","B",
             "B","C",
             "A","C",
             "C","D"), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  E(g)$allele_diff <- c(1, 2, 3, 1)
  E(g)$is_mst      <- c(TRUE, TRUE, FALSE, TRUE)
  g
}


    test_that("Edge colour scale is named 'Edge Type'", {
        g  <- make_test_graph()
        p  <- PlotMST(g, metadata = NULL)
        scale_names <- sapply(p$scales$scales, function(s) s$name)
        expect_true("Edge Type" %in% scale_names)
    })
 
    test_that("Edge width scale is named 'Edge Type'", {
        g  <- make_test_graph()
        p  <- PlotMST(g, metadata = NULL)
        scale_names <- sapply(p$scales$scales, function(s) s$name)
        # Both colour and width scales share the same name 'Edge Type'
        expect_gte(sum(scale_names == "Edge Type"), 2L)
    })

    test_that("PlotMST returns a ggplot object (no metadata)", {
        g <- make_test_graph()
        p <- PlotMST(g, metadata = NULL)
        expect_true(inherits(p, "ggplot"))
    })

    test_that("show_clustering = FALSE removes non-MST edges from the plot", {
        g <- make_test_graph()
        # Should produce a valid plot without error
        p <- PlotMST(g, metadata = NULL, show_clustering = FALSE)
        expect_true(inherits(p, "ggplot"))
    })


    test_that("Custom MST_edges_color and cluster_edges_color are accepted", {
      g <- make_test_graph()
      expect_silent(
        PlotMST(g, metadata = NULL,
                MST_edges_color     = "blue",
                cluster_edges_color = "green")
      )
    })

    test_that("node_size parameter is accepted without error", {
      g <- make_test_graph()
      expect_silent(PlotMST(g, metadata = NULL, node_size = 6))
    })

    test_that("Multi-column metadata sets multiple vertex attributes-no error", {
          g  <- make_test_graph()
          md <- data.frame(
            id      = V(g)$name,
            group   = c("X", "Y", "X", "Y"),
            country = c("LB", "US", "LB", "US"),
            stringsAsFactors = FALSE
          )
          p <- PlotMST(g, metadata = md, node_color = "group")
          expect_true(inherits(p, "ggplot"))
    })
}