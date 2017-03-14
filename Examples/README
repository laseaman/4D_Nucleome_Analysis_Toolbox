%% Hi-C analysis READ_ME

% Package for Hi-C data including loading, plotting, translocation exploration, normalization, and TAD calculation.

% packaged by:
% Laura Seaman
% PhD Candidate in Bioinformatics, University of Michigan
% laseaman@umich.edu

% includes functions to load Hi-C data from text files, normalization, TAD
% detection, and plotting.
% scripts showing examples are included.

% install 4DNucleomeAanalysisToolbox.mltbx by downloading and double
% clicking.

% run examples by downloading .m files, opening and running in matlab.

%% what is included

% key functions (excludes additional helper funcitons)
%   arrow - ADD DETAIL
%   BlockToepNorm - performs copy number based normalization
%   cbar - adds small colorbars to at the top right. Sized based on number
%       of subplots
%   Draw_TADs - plots Hi-C data with TAD boundaries overlaid
%   HiC_load - loads Hi-C matrices from text files (in homer format)
%   HiC_plot - plots Hi-C matrix in log2 scale
%   HiC_remove_cent - removes centromeres from Hi-C matrices
%   HicRna_plot - plots HiC, RNA-seq, TADs, and Fiedler vector in various
%       combinations
%   HicTensorFig - ADD DETAILS
%   ICE - performs iterative correction and eigenvector decomposition
%       normalization
%   Mat2Val - calculate Fiedler number, Fiedler vector, Von Neumann
%       Entropy, and more from matrices depending on provided type.
%   PltFV - plots Fiedler vector bar graph with positive and negative
%       values colored seperately
%   PltPhasePlane - ADD DETAIL
%   TAD_DP1 - performs dynamic programming / multiscale TAD analysis
%   TAD_Laplace - performs iterative, Fiedler number based TAD analysis
%   ToepNorm - performs Toeplitz based normalization
%   TranslocHiC - constructs translocated chromsomes matrices from parts

% data
%   HT29_2D_12hr_100kb_chr8.txt - text file including 100 kb resolution chromosome 
%       8 from a 12 hour time point of 2D grown HT-29 cells, originally published
%       in Nucleome analysis of a colorectal cancer cell line reveals
%       structure-function relationships by Seaman et al, Mol Canc Res xx 2017
%   chr22_100kb_RnaHiC_FibTS.mat - 100 kb resolution chr 22 RNA_seq and Hi-C matrix 
%      of all time points from Functional organization of the human 4D Nucleome 
%      by Chen et al 2015 as well as gene expression (FPKM) for all genes and 
%      time points.
%   chr6_chr14_reads_HT29.mat - chr 6, 14 and 6-14 inter-chr read locations
%       from the 2D 12 hr sample  published in Nucleome analysis of a 
%       colorectal cancer cell line reveals structure-function relationships
%       by Seaman et al, Mol Canc Res xx 2017.
%   chr6_chr14_100kb_RNA_HT29.mat - chr 6 and 14 RNA-seq binned into 100 kb
%      bins to match the Hi-C data published in Nucleome analysis of a 
%      colorectal cancer cell line reveals structure-function relationships by 
%      Seaman et al, Mol Canc Res xx 2017.
%   chr6_chr14_100kb_HT29.mat - chr 6, 14 and 6-14 inter-chr matrix from
%       the 2D 12 hr sample at 100kb resolution published in Nucleome analysis 
%       of a colorectal cancer cell line reveals structure-function relationships
%       by Seaman et al, Mol Canc Res xx 2017.

% examples
%   Load_Normalize_Example - loads Hi-C data from a text file, removes unmappable 
%       regions, normalizes the data using three different methods and plots the 
%       raw and normalized matrices.  
%   PhasePlane - plots time series data in three dimensions and a phase plane 
%       (plot of structure measured by Fiedler number vs function measured by 
%       square root of the average FPKM for all genes in the region/chromosome
%       over time) for chromsome 22 for the fibroblast times series.
%   TranslocationAnalysis - uses read level data for HT-29 from chromosmes 6 
%       and 14 to identify the site of translocation at high resolution. Also
%       identifies the site at 100 kb resolution, constructs the translocated 
%       chromosome, normalizes Chr 6, 14, and t(6;14), and calculates TADs. 
%       Demonstrates plotting abilities including plotting Hi-C, RNA-seq, Fiedler
%       vector, and TADs on one plot.
%   TAD_methods - performs 3 different methods of TAD calculation with
%       plotting on chromsome 22 from fibroblasts. note: warnings "Matrix is 
%       singular to working precision." during TAD_HMM is expected.

