clear all, close all, clc

%% Load Hi-C data

HiC_wCent = HiC_load('HT29_2D_12hr_100kb_chr8.txt');

%% remove centromere / unmappable regions

HiC = HiC_remove_cent(HiC_wCent);

%% save loaded data
% for faster running, if repeating parts. Run all once, uncomment the line
% below to save the loaded data as a matlab file. Then comment out the line
% below as well as everything above except the first line, and uncomment the 
% one under that to load the saved data.

save('chr8_1Mb_HT29.mat','HiC','HiC_wCent')
% load('chr8_1Mb_HT29.mat')

%% Toeplitz normalization
% originally described in: Spectral Identification of Topological Domains
% by Chen et al, Bioinformatics May 2016

Toep_HiC = ToepNorm(HiC);

%% copy number normalization
% originally described in: Nucleome analysis of a colorectal cancer cell line 
% reveals structure-function relationships by Seaman et al, Mol Canc Res xx 2017

% change jumpThreshold based on resolution and chromsome, set it at a
% number high enough to only select big jumps between sustained regions
% rather than just the noisy variation.

jumpThreshold = 75;
plotNorm = 1;
resolution = 1e5;
[CN_HiC,brks] = BlockToepNorm(HiC,jumpThreshold,resolution,plotNorm);
'break used during copy number normalization ',brks

%% ICE normalization
% originally described in: Iterative Correction of Hi-C Data Reveals Hallmarks 
% of Chromosome Organization by Imakaev et al, Nat Methods Oct 2012

ICE_HiC = ICE(HiC);

%% plot normalization methods

figure('Position',[100,100,800,900])

% plot raw HiC
plt=1; subplot(2,2,plt)
HiC_plot(HiC,'Raw Hi-C')
cbar(2,2,plt,[-1,floor(max(max(log2(HiC-diag(diag(HiC))+.5))))])

% plot Toeplitz normalized
plt=2; subplot(2,2,plt)
HiC_plot(Toep_HiC,'Toeplitz Normalized')
cbar(2,2,plt,[-1,floor(max(max(log2(Toep_HiC-diag(diag(Toep_HiC))+.5))))])

% plot copy number normalized
plt=3; subplot(2,2,plt)
HiC_plot(CN_HiC,'Copy Number Normalized')
cbar(2,2,plt,[-1,floor(max(max(log2(CN_HiC-diag(diag(CN_HiC))+.5))))])

% plot ICE normalized
plt=4; subplot(2,2,plt)
HiC_plot(ICE_HiC,'ICE Normalized')
cbar(2,2,plt,[-1,floor(max(max(log2(ICE_HiC-diag(diag(ICE_HiC))+.5))))])

