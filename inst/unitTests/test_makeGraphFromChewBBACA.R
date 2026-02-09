library(RUnit)
require(igraph)

csv_path <- system.file("extdata", "cgMLST95.csv", package = "MSTree")
test_data <- read.csv(csv_path, row.names = 1, check.names = FALSE, sep='\t')

test_makeGraphFromChewBBACA <- function() {
    checkEquals(class(makeGraphFromChewBBACA(test_data, max_allelic_difference=9)), "igraph")

    checkEquals(class(test_data), "data.frame")
}
