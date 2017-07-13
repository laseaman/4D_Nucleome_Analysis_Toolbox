clear, close all, clc

% Load the example data
load('chr22_100kb_RnaHiC_FibTS')
% select a single time point
tp = 2; % time point to use
HiC_wCent = double(C22(:,:,tp));
[HiC,~,idx_cent] = HiC_remove_cent(HiC_wCent,0);
chrN = 22;

% to save TADs in text files, make the variables file names. Set the 
%   variables to 0 to not save TAD results.
save_iterTADs = 'fb_c22_iter_tads.csv';
save_multiTADs = 'fb_c22_multiScale_tads.csv';
save_hmmTADs = 'fb_c22_hmm_tads.csv';
res = 1e5; % resolution of Hi-C mat, used for saving, 1e5 = 100kb, 1e6 = 1Mb

%% iterative TADs

% Script for defining topological domains using the iterative method
% Based on: Spectral Identification of Topological Domains by Chen et al, 2016

TADs_iter = TAD_Laplace(HiC,.4);

% Display
Draw_TADs(HiC, TADs_iter,[0,6]);
title(['Chr ',num2str(chrN),' TADs, iterative method'])

% save to text file
if ischar(save_iterTADs)
    saveTADs(save_iterTADs,TADs_iter,res,idx_cent)
end

%% multiscale TADs

% Script for defining topological domains using the multiscale method
% Based on: Multiscale Identification of Topological Domain in Chromatin
%     by Filippova et al, 2013

% resolution value
gamma = 0.20;

% Call Algorithm
TAD_bound = TAD_DP1(HiC,gamma);

% Display
Draw_TADs(HiC, TAD_bound,[0,6]);
title(['Chr ',num2str(chrN),' TADs, multiscale method'])

% save to text file
if ischar(save_multiTADs)
    saveTADs(save_multiTADs,TAD_bound',res,idx_cent)
end

%% HMM method

% Script for defining topological domains using the HMM method.
% Based on: Topological domains in mammalian genomes identified by analysis
%     of chromatin interactions by Dixon et al 2012

% note: several warnings about "matrix is singular to working precision" are common

% Searching length
L  = 10;

% Remove unmappable region
[HiC,~,idx_cent] = HiC_remove_cent(HiC);

% rescale / normalize
% Apply a transformation
HiC = ceil(HiC);
% log transformation, and saturated by 6;
HiC = min(log(HiC),6);
% Process -inf, because log(0) = -inf
HiC(HiC == -inf) = -1;
% Shift to be positive
HiC = HiC + 1.001;

% Call Algorithm
TAD_mark=TAD_HMM(HiC,L);
 % for verbose output set, change 'verbose',0,  to 'verbose',1, in mhmm_em, line 41
TAD_boundaries = TADMark2Pos_HMM(TAD_mark);

% Display
Draw_TADs(HiC, TAD_boundaries,[0,6]);
title(['Chr ',num2str(chrN),' TADs, HMM method'])

% save to text file
if ischar(save_hmmTADs)
    saveTADs(save_hmmTADs,TAD_boundaries',res,idx_cent)
end

[length(TADs_iter), length(TAD_bound), length(TAD_boundaries)]
