clear all, close all, clc

f_readData = 'chr6_chr14_read_HiC_FibTs.mat';
chr1 = 6; chr2 = 14;

%% load read data from text files
% example of loading read level Hi-C data into matlab. Takes < 4 minutes to
% run on 2.4 GHz Intel Core i5, 8GB laptop. Preloaded data is included for
% convience and is loaded below this section. Only inter-chr reads are
% included although filtering steps in the next section are sufficient to
% avoid intra-chr reads if custom data is used.
% 
% When loading a whole genome's worth of data, storing each chromosome
% seperately is recommended due to time and size requirements. Homer by
% default produces a seperate file for each chromosome so loading htem each
% seperately works well.

% tic
% flocs = '../Example_readFile/';
% fnames = {'FibTs0_chr14_interReads.txt','FibTs0_chr6_interReads.txt'};
% 
% % load first file - all reads connect to chr 6
% [locA1, chrB1, locB1 ] = Load_TSV( [flocs,fnames{1}] );
% chrA1 = 6*ones(size(chrB1));
% 
% toc
% % load second file - all reads connect to chr 14
% [locA2, chrB2, locB2 ] = Load_TSV( [flocs,fnames{2}] );
% chrA2 = 6*ones(size(chrB2));
% toc
% 
% chrA = [chrA1;chrA2];
% locA = [locA1;locA2];
% chrB = [chrB1;chrB2];
% locB = [locB1;locB2];
% toc
% save(f_readData,'chrA','chrB','locA','locB')

load(f_readData)

%% find translocation in read data

% select reads that connect the two chromosomes of interest
idx1 = chrA == chr1 & chrB == chr2;
idx2 = chrA == chr2 & chrB == chr1;
locChr1 = [locA(idx1);locB(idx2)];
locChr2= [locB(idx1);locA(idx2)];

% define the region within which you want to plot all of the reads
reg1 = [132500000,133200000]; % region on chr 6 in this example
reg2 = [36000000,37000000]; % region on chr 14 in this example
% define (visually) where the translocation occurs
t1 = 132825000;
t1b = 132890000; % if there is only one break point comment this line
t2 = 36508800;

% plot the reads
figure
scatter(locChr1,locChr2,'.')
set(gca,'Ydir','reverse')
axis([reg1,reg2])
% add lines showing where the breakpoint is marked
line([t1,t1],reg2,'Color','r','LineWidth',2)
line([t1b,t1b],reg2,'Color','r','LineWidth',2) % if there is only one break point comment this line
line(reg1,[t2,t2],'Color','r','Color','r','LineWidth',2)
xlabel(['Chr ',num2str(chr1)]), ylabel(['Chr ',num2str(chr2)])
title(['Chr ',num2str(chr1),'- Chr ',num2str(chr2)])

