## 4D_Nucleome_Analysis_Toolbox

### Introduction

The 4D Nucleome Analysis Toolbox (4D NAT) includes functions to load Hi-C matrices, Hi-C read data, and RNA-seq from text files, normalize Hi-C data, detect TADs, plotting, translocation exploration, or. Example scripts showing examples are available at: https://github.com/laseaman/4D_Nucleome_Analysis_Toolbox

### Installation

install 4DNucleomeAanalysisToolbox.mltbx by downloading and double clicking. Run examples by downloading .m files, opening and running in matlab. Suggested order: Load_Normalize.m, Tad_methods.m, TranslocationAnalysis_100kb.m, TranslocationAnalysis_read.m, PhasePlane.m

### Help / Report Bugs

Check to make sure the toolbox is installed by looking at: Home/Add-Ons/Manage Add-Ons and verifying that the tool box is listed there. If it is not, double click on the toolbox to install. For additional support, email Laura Seaman at laseaman@umich.edu

### Examples

    Load_Normalize_Example - loads 100kb chr8 Hi-C data from a text file. Then removes unmappable regions (repetative regions like centromeres), normalizes the data using three different methods and plots the raw and normalized matrices. The normalization methods are:
      1. Spectral Identification of Topological Domains by Chen et al, Bioinformatics May 2016
      1. Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017
      1. Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome Organization by Imakaev et al, Nat Methods Oct 2012

    TAD_methods - performs 3 methods of TAD calculation and plotting on a 100kb resolution chromsome 22 Hi-C matrix from fibroblasts. Normalization methods are:
      1. Spectral Identification of Topological Domains by Chen et al, 2016
      1. Multiscale Identification of Topological Domain in Chromatin by Filippova et al, 2013
      1. Topological domains in mammalian genomes identified by analysis of chromatin interactions by Dixon et al 2012
      note: warnings "Matrix is singular to working precision." during TAD_HMM is normal.

    TranslocationAnalysis_100kb - Analyzes t(6;14) in 100 kb resolution data including identifying the site of translocation, and constructing the translocated chromosome. It also includesnormalizing Chr 6, 14, and t(6;14), and calculating TADs. Demonstrates plotting abilities including plotting Hi-C, RNA-seq, Fiedler vector, and TADs in one figure.

    TranslocationAnalysis_read -  Loads read level data from output text files produced by Homer and other software. Uses read level data for HT-29 from chromosmes 6 and 14 to identify the site of translocation at high resolution.

    PhasePlane - demonstrates loading RNA-seq data and converting raw data into binned data. Plots time series data in three dimensions with the x-y dimension showing Hi-C matrices and the z-direction showing time. Also calculate necessary data and and a phase plane which shows how structure, as measured by the Fiedler number and function as measured by the square root of the average FPKM for all genes in the chromosome (or region) for chromosome 22 fibroblast time series.

## Example data
    HT29_2D_12hr_100kb_chr8.txt - text file including 100 kb resolution chromosome 8 from a 12 hour time point of 2D grown HT-29 cells, originally published in Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017 chr22_100kb_HiC_FibTS.mat - 100 kb resolution chr 22 RNA_seq and Hi-C matrix of all time points from Functional organization of the human 4D Nucleome by Chen et al 2015.
    chr6_chr14_reads_HT29.mat - chr 6, 14 and 6-14 inter-chr read locations from the 2D 12 hr sample  published in Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017.
    readData_optional/ FibTs0_chr6_interReads.txt & FibTs0_chr12_interReads.txt - text files that contain the read pairs for all inter-chromosomal reads for chromsomes 6 and 14 respectively. Published in Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017.
    chr6_chr14_100kb_RNA_HT29.mat - chr 6 and 14 RNA-seq binned into 100 kb bins to match the Hi-C data published in Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017.
    chr6_chr14_100kb_HT29.mat - chr 6, 14 and 6-14 inter-chr matrix from the 2D 12 hr sample at 100kb resolution published in Nucleome analysis of a colorectal cancer cell line reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017.
    fib_ts_RNA_rawCounts.csv - all raw count gene expression estimates for all time points in the cell cycle synchronized proliferating  fibroblasts. Originally published in Functional organization of the human 4D Nucleome by Chen et al 2015
    GenInfo_loc.mat - gene names, brief descriptions, and locations for hg19 reference genome. Also includes the length of each chromosome in 100 kb bins.

## Function descriptions
description of key functions overview, full documentation including descriptions of all inputs and outputs are included within each function. Type help fxnName in the command window for more information.

    BlockToepNorm - performs copy number based normalization
    cbar - adds small colorbars to at the top right. Sized based on number of subplots
    Draw_TADs - plots Hi-C data with TAD boundaries overlaid
    FindGenBin - looks up a gene's location.
    FindGenLocation - Calculates which bin(s) a gene is locatated in.
    HiC_load - loads Hi-C matrices from text files (in homer format)
    HiC_plot - plots Hi-C matrix in log2 scale
    HiC_remove_cent - removes centromeres from Hi-C matrices
    HicRna_plot - plots HiC, RNA-seq, TADs, and Fiedler vector in various combinations
    HicTensorFig - plots time series data in 3D. z-dimension indicates time while each x-y plane shows the Hi-C and RNA-seq for a single time point.
    ICE - performs iterative correction and eigenvector decomposition normalization
    Load_TSV - loads tab seperated value files containing read pairs of Hi-C data.
    Mat2Val - calculate Fiedler number, Fiedler vector, Von Neumann Entropy, and more from matrices depending on provided type.
    PltFV - plots Fiedler vector bar graph with positive and negative values colored seperately
    PltPhasePlane - plots a phase plane showing a projection of structure and function over time.
    RNA2Bin - converts gene resolution RNA-seq data into bin resolution.
    TAD_DP1 - performs dynamic programming / multiscale TAD analysis
    TAD_HMM - performs Hiden Markov Model based algorithm for TAD analysis.
    TAD_Laplace - performs iterative, Fiedler number based TAD analysis
    ToepNorm - performs Toeplitz based normalization
    TranslocHiC - constructs translocated chromsomes matrices from parts

