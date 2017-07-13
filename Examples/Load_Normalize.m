clear all, close all, clc

fname_txtInput = 'HT29_2D_12hr_100kb_chr8.txt';
fname_matFile = 'chr8_1Mb_HT29.mat';
res = 1e6; % used for saving - 1e6 = 1 Mb, 1e5 = 100kb

% set to a filename to save normalized matrices as comma seperated text
%   files. Set to 0 to not save
saveToeplitz = 0;
saveCopyNum = 'CopyNumNormHT29.txt';
saveICE = 0;

%% Load Hi-C data
% several formats are supported. Select one below based on the format.

% n x n plain text with a Homer like header. Chromosome format is chrN
HiC_wCent = HiC_load_mat(fname_txtInput);

% n x n plain text with a Homer like header. Chromosome format is N
% HiC_wCent = HiC_load_mat(fname_txtInput,0,1);

% n x n plain text no header or other format. loads matrix only, not chr bounds
% HiC_wCent = HiC_load_mat(fname_txtInput,1);

% load .cool file (binary HDF5 of sparse matrix)
% [HiC_wCent,chrBnds] = HiC_load_cool(fname_cool);

%% remove centromere / unmappable regions

[HiC,~,idx_cent] = HiC_remove_cent(HiC_wCent);

%% save loaded data
% for faster running, especially for repeating following analysis. Run this 
% and the lines above it once, with the line below uncommented to save the 
% loaded data as a matlab file. Then comment out the line below as well as
% everything above except the first line, and uncomment the line under that
% to load the saved data.

save(fname_matFile,'HiC','HiC_wCent')
% load(fname_matFile)

%% Toeplitz normalization
% originally described in: Spectral Identification of Topological Domains
% by Chen et al, Bioinformatics May 2016

Toep_HiC = ToepNorm(HiC);

% save as excel file if wanted
if ischar(saveToeplitz)
    saveHiC(saveToeplitz, Toep_HiC, res, idx_cent)
end

%% copy number normalization
% originally described in: Nucleome analysis of a colorectal cancer cell line 
% reveals structure-function relationships by Seaman et al, Mol Canc Res March 2017

jumpThreshold = 70;
% change jumpThreshold based on resolution and chromsome, set it at a
% number high enough to only select big jumps between sustained regions
% rather than just the noisy variation. For this example anything about 30
% and below 78 will lead to the same choice in breakponts (and normalized
% matrix)

plotNorm = 1; 
% plotNorm when set to 1 will produce a figure to help with selection of
% jumpThresold. The left plot shows the total number of reads per bin, the
% values after smoothing, and points that were preliminaryly selected as
% large changes. These become final changes if they are the largest in the
% region (1/8th chromosome is minimum distance between two breakpoints).
% The right plot shows the change in smoothed number of counts between sequential
% bins as well as the thresold the change must be greater than to be
% considered as a change in copy number (creating a breakpoint in the
% normalization). When set to 0, this plot is not produced.

filtMax = .01;
% this is the upper limit for the bandpass filter. Use 0.01 for 100kb
% matrices and .1 for 1 Mb matrices.

% perform normalization
[CN_HiC,brks] = BlockToepNorm(HiC,jumpThreshold,filtMax,plotNorm);
'break used during copy number normalization ',brks

% save normalized matrix as a comma seperated text file
if ischar(saveCopyNum)
    saveHiC(saveCopyNum, CN_HiC, res, idx_cent)
end

%% ICE normalization
% originally described in: Iterative Correction of Hi-C Data Reveals Hallmarks 
% of Chromosome Organization by Imakaev et al, Nat Methods Oct 2012

ICE_HiC = ICE(HiC);

% save normalized matrix as a comma seperated text file
if ischar(saveICE)
    saveHiC(saveICE, CN_HiC, res, idx_cent)
end

%% plot normalization methods

% these lines can be used to zoom into a region on the Hi-C matrix. As
% written it shows the entire chromsome. In addition to manually zooming
% through the magnifing glasses on each figure, uncommenting the alternate
% definations of zoomBounds will zoom in on: a) just the HSR or b) the
% boundary between the low and mid copy number regions c) off diagonal 
% region where the two ends of the HSR interact.

zoomBounds = [1,length(HiC), 1, length(HiC)];
%zoomBounds = [1100,length(HiC), 1100,length(HiC)]; % a) HSR
%zoomBounds = [600, 1000, 400, 480]; % b) CN boundary
%zoomBounds = [1320, 1420, 1080, 1220]; % c) HSR corner

figure('Position',[100,100,800,900])

% plot raw HiC
plt=1; subplot(2,2,plt)
HiC_plot(HiC,'Raw Hi-C');
maxVal = floor(max(max(log2(HiC-diag(diag(HiC))+.5))));
cbar(2,2,plt,[-1,maxVal])
axis(zoomBounds)

% plot Toeplitz normalized
plt=2; subplot(2,2,plt)
HiC_plot(Toep_HiC,'Toeplitz Normalized');
maxVal = floor(max(max(log2(Toep_HiC-diag(diag(Toep_HiC))+.5))));
cbar(2,2,plt,[-1,maxVal])
axis(zoomBounds)

% plot copy number normalized
plt=3; subplot(2,2,plt)
HiC_plot(CN_HiC,'Copy Number Normalized');
maxVal = floor(max(max(log2(CN_HiC-diag(diag(CN_HiC))+.5))));
cbar(2,2,plt,[-1,maxVal])
axis(zoomBounds)

% plot ICE normalized
plt=4; subplot(2,2,plt)
HiC_plot(ICE_HiC,'ICE Normalized');
maxVal = floor(max(max(log2(ICE_HiC-diag(diag(ICE_HiC))+.5))));
cbar(2,2,plt,[-1,maxVal])
axis(zoomBounds)

% VNE and FV for 1 TAD in TS dataset