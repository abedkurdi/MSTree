MSTree
======

MSTree is a Bioconductor package that generates a minimum spanning tree based on the output of ChewBBACA pipeline.

Dependencies
-----------

The below R packages are required for installation:

+ igraph
+ ggraph
+ ggplot2

Installation
------------

1. Git clone the project directory `git clone https://github.com/abedkurdi/MSTree.git`
2. From the terminal run `R CMD build MSTree`
3. Run `R CMD check MSTree_0.99.0.tar.gz`
4. Run  `R CMD BiocCheck MSTree_0.99.0.tar.gz`
5. Install package by `R CMD INSTALL MSTree_0.99.0.tar.gz`
6. Start R and load library using `library(MSTree)`
