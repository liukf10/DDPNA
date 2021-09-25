
<!-- README.md is generated from README.Rmd. Please edit that file -->
DDPNA
=====

<!-- badges: start -->
<!-- badges: end -->
'DDPNA' is a package for Disease-Drived Differential Proteins(DFP) and Proteome-Wide Co-expression Network Associated Analysis. The goal of 'DDPNA' is offered a better methods to analyze omic data. The package is designed for proteomic data, but it is also fit for expression data in RNA-seq and metabolome. It is associated DFP and co-expression network module, and constructed a Mod-DFP network to remove lower connectivity DFP. The lower connectivity DFP is hard to get the key function in PPI and is more likely a false postive protein. The Mod-DFP network can also get DFP related proteins which is more likely a false negative protein. It provides the essential statisic analysis included t.test, ANOVA analysis to extract differential proteins. The package also provide some module analysis included PCA analysis, two enrichment analysis, Planner maximally filtered graph extraction and hub analysis. The co-expression network should constructed by other package or software.('WGCNA' package or others)

Installation
------------

You can install the developed version of DDPNA from [github](https://github.com/liukf10/DDPNA)

``` r
install.packages("devtools")
devtools::install_github("liukf10/DDPNA")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

    #> 
    #> Warning in .NAnum.proteomic_data(data, miss.value = miss.value, verbose =
    #> verbose): The sample ad059 have been removed

    #> ####### PFN Calculation commences ########
    #> permutation no.:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,
    #> Warning: package 'MEGENA' was built under R version 3.5.3
    #> Loading required package: doParallel
    #> Warning: package 'doParallel' was built under R version 3.5.3
    #> Loading required package: foreach
    #> Warning: package 'foreach' was built under R version 3.5.2
    #> Loading required package: iterators
    #> Warning: package 'iterators' was built under R version 3.5.2
    #> Loading required package: parallel
    #> Loading required package: igraph
    #> Warning: package 'igraph' was built under R version 3.5.2
    #> 
    #> Attaching package: 'igraph'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     decompose, spectrum
    #> The following object is masked from 'package:base':
    #> 
    #>     union
    #>  - # of genes: 57 
    #>  - # of hubs: 2 
    #> - generating module subnetwork figure...

References
==========

1.  Peter, L.; Steve H. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics.* **2008**, 9, 559.

2.  Won-Min Song, Bin Zhang. Multiscale Embedded Gene Co-expression Network Analysis *PLOS Computational Biology.* **2015**, 11(11), e1004574

3.  Zhenfeng Wu, et al. NormExpression: an R package to normalize gene expression data using evaluated methods. *BioRxiv* **2018**, JAN.

4.  Aravind Subramanian, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles *PNAS* **2005**, 102(43), 15545-15550.

5.  Vince V. ggbiplot. **2018**; <https://github.com/vqv/ggbiplot>.

6.  Michael C Oldham, et al. Network Methods for Describing Sample Relationships in Genomic Datasets: Application to Huntington's Disease; <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/SampleNetwork/>
