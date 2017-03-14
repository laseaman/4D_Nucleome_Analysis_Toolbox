clear, close all, clc

% Load the example data
load('chr22_100kb_RnaHiC_FibTS')
% select a single time point
C22 = double(C22(:,:,2));

%% HMM method

% Script for defining topological domains using the HMM method.
% Based on: Topological domains in mammalian genomes identified by analysis
%     of chromatin interactions by Dixon et al 2012

% note: several warnings about "matrix is singular to working precision" are common

% Searching length
L  = 10;

% Remove unmappable region
[C22,~,idx_cent] = HiC_remove_cent(C22);

% rescale / normalize
% Apply a transformation
Ht = ceil(C22);
% log transformation, and saturated by 6;
Ht = min(log(Ht),6);
% Process -inf, because log(0) = -inf
Ht(Ht == -inf) = -1;
% Shift to be positive
Ht = Ht + 1.001;

% Call Algorithm
TAD_mark=TAD_HMM(Ht,L);
 % for verbose output set, change 'verbose',0,  to 'verbose',1, in mhmm_em, line 41
TAD_boundaries = TADMark2Pos_HMM(TAD_mark);

% Display
Draw_TADs(Ht, TAD_boundaries,[0,6]);
title('Chr 22 TADs, HMM method')

%% multiscale TADs

% Script for defining topological domains using the multiscale method
% Based on: Multiscale Identification of Topological Domain in Chromatin
%     by Filippova et al, 2013

% resolution value
gamma = 0.20;

% Call Algorithm
TAD_bound = TAD_DP1(Ht,gamma);

% Display
Draw_TADs(Ht, TAD_bound,[0,6]);
title('Chr 22 TADs, multiscale method')

%% iterative TADs

% Script for defining topological domains using the iterative method
% Based on: Spectral Identification of Topological Domains by Chen et al, 2016

TADs_iter = TAD_Laplace(C22,.4);

% Display
Draw_TADs(Ht, TADs_iter,[0,6]);
title('Chr 22 TADs, iterative method')

[length(TADs_iter), length(TAD_bound), length(TAD_boundaries)]
