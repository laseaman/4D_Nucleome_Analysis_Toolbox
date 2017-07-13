close all, clc, clear all
%% load RNA from text

% fname is a comma seperated value file in which the first column is gene
%    name, all following columns are raw counts for the associated gene for 
%    different time points. Change rawSt to the first row that contains a
%    gene name based on the number of header lines included in the file.
fname = 'fib_ts_RNA_rawCounts.csv';
rowSt = 2; %row in the csv file in which the first gene name is stored.
T = importdata(fname);

RNA_gene = T.data;
gene_names = T.textdata(rowSt:end,1);
clear T

%% convert RNA to 100 kb resolution bins

% load gene annotations and chromosome lengths - hg19 reference genome. 
%   Can be used for any hg19 aligned samples.
load('GenInfo_loc.mat')

% bin at 100 kb resolution
collectGenes = 1; % set to 0 to not report back a cell in which each entry  
%    is a cell of all of the genes in that bin.
res = 1e5; % resolution, 1e5=100kb resolution matrix
[Rna_100kb,GenesInBins] = RNA2Bin(RNA_gene, gene_names, Gen_info, chrBnds, res, collectGenes);
size(Rna_100kb)

save('fib_ts_RNA_rawAndBinned','chrBnds','RNA_gene','gene_names',...
    'Rna_100kb','GenesInBins')
% run the all lines above this once to load and save the data. Then comment
% the line above out, and uncomment the one below this for faster running.

load('fib_ts_RNA_rawAndBinned')

% view genes in a bin
if collectGenes
    bn = 20;
    GenesInBins{20}
end

%% select desired chromosome

chrN = 22; % what chromosome should be selected 
RNA = Rna_100kb(chrBnds(chrN):chrBnds(chrN+1)-1,:); % contains expression 
% for the selected chromsome for all samples.

%% load HiC data

fname = 'chr22_100kb_RnaHiC_FibTS.mat';
load(fname)
HiC = C22; 

%% remove centromere

[HiC,~,idx_cent] = HiC_remove_cent(HiC);
RNA = RNA(~idx_cent,:);

%% plot all data

figure('Position',[100,20,1000,600]);
HicTensorFig( double(HiC), RNA,1);

%% customized labeling
% positioning probably needs to change for other datasets

yshift = 30; 
zshift = 15;
zshift2 = 5;

text(1.2,length(HiC)+yshift,zshift2,'RNA-seq')
text(.85,length(HiC)-yshift,-length(HiC)-zshift,'Hi-C')
text(4.4,length(HiC),zshift2,'Time','Rotation',28);

text(.9,-10,-length(HiC)-zshift,'t = 1')
text(1.9,-10,-length(HiC)-zshift,'t = 2')
text(7.9,-10,-length(HiC)-zshift,'t = 8')

%% caculate function

idxGenes = GeneLocs(:,1)==chrN;
RPKM_genes = GeneRpkm(idxGenes,:);
Function = sqrt(mean(RPKM_genes,1));

%% claculate Fiedler number from Hi-C

% remove centromeres
[HiC_noCent,~,idx_cent] = HiC_remove_cent(double(HiC),1);

Norm22 = NaN(size(HiC_noCent));
FN22 = NaN(size(HiC_noCent,3),1);
VNE22 = NaN(size(HiC_noCent,3),1);
figure('Position',[20,100,900,900])
for s = 1:size(HiC_noCent,3)
    % normalize
    Norm22(:,:,s) = ToepNorm(HiC_noCent(:,:,s));
    % Fiedler number
    FN22(s) = Mat2Val(Norm22(:,:,s),3);
    % Von Nuemann Entropy
    VNE22(s) = Mat2Val(Norm22(:,:,s),9);
    
    subplot(3,3,s)
    HiC_plot(Norm22(:,:,s),['time point ',num2str(s)]);
end

%% phase plane
% structure vs function
% Fiedler number vs mean sqrt RPKM

figure('Position',[200,400,500,500])
nSpline = 10; % calculates this number of intermediates between each time point
%      leading to a smooth curve instead of straight lines between points.
PltPhasePlane( FN22,Function,'','none',nSpline)
% for title, in previous line: replace '' with '\fontsize{16}Chr 22 Phase Plane'
set(gca, 'box', 'off')
text(FN22(1),Function(1)+.02,'t = 1','FontSize',14)
text(FN22(2)-.007,Function(2)-.006,'t = 2','FontSize',14)
text(FN22(end)+.001,Function(end)-.006,'t = 8','FontSize',14)

%% phase plane - using entropy
% structure vs function
% Fiedler number vs mean sqrt RPKM
% plots with straight arrow lines instead of spline which was shown above.

figure('Position',[200,400,500,500])
PltPhasePlane( VNE22,Function)
% for title, add to end of previous line: ,'\fontsize{16}Chr 22 Phase Plane'
xlabel({'\fontsize{15}Structure','\fontsize{12}Von Neumann Entropy'})
set(gca, 'box', 'off')
text(VNE22(1),Function(1)+.01,'t = 1','FontSize',14)
text(VNE22(2)+.001,Function(2)-.006,'t = 2','FontSize',14)

%% define TADs
% using iterative/Laplacian method on the sum of the matrix for all of the
% time points to provide a single defination for TADs for all time points.
% This allows comparing the properties of a TAD over time, since the
% boundaries do not change.

TADs_iter = TAD_Laplace(sum(HiC_noCent,3),.4);
nTAD = 28; 
res = 1e5; % resolution, 1e5=100kb, used for TAD RNA-seq

% select the TAD from the chromosome's Hi-C matrix
r = TADs_iter(nTAD):TADs_iter(nTAD+1)-1;
HiCTad = HiC_noCent(r,r,:);

%% caculate function of TAD

% find TAD bounds in original chr coordinates
TADs_wCent = AddCent(TADs_iter,idx_cent);
idxTAD = GeneLocs(:,2) < TADs_wCent (nTAD+1)*res & GeneLocs(:,3) > TADs_wCent (nTAD)*res;
% start of gene is before end of TAD & end of gene is after start of TAD
RPKM_genes = GeneRpkm(idxGenes & idxTAD,:);
% select genes on this chromosome and in this TAD

FxnTad = sqrt(mean(RPKM_genes,1)); % function of the TAD

%% claculate Fiedler number from Hi-C

FNTad = NaN(size(HiCTad,3),1);
VNETad = NaN(size(HiCTad,3),1);
for s = 1:size(HiCTad,3)
    % normalize
    NormTad = ToepNorm(HiCTad(:,:,s));
    % Fiedler number
    FNTad(s) = Mat2Val(NormTad,3);
    % Von Nuemann Entropy
    VNETad(s) = Mat2Val(NormTad,9);
end

%% TAD phase plane
% structure vs function
% Fiedler number vs mean sqrt RPKM

figure('Position',[200,400,500,500])
PltPhasePlane( FNTad,FxnTad)
% for title, add to end of previous line: ,'\fontsize{16}TAD Phase Plane'
set(gca, 'box', 'off')
text(FNTad(1),FxnTad(1)+.02,'t = 1','FontSize',14)
text(FNTad(2)-.007,FxnTad(2)-.006,'t = 2','FontSize',14)
text(FNTad(end)+.001,FxnTad(end)-.006,'t = 8','FontSize',14)

%% TAD phase plane - using entropy
% structure vs function
% Fiedler number vs mean sqrt RPKM

figure('Position',[200,400,500,500])
PltPhasePlane( VNETad,FxnTad )
% for title, add to end of previous line: ,'\fontsize{16}TAD Phase Plane'
xlabel({'\fontsize{15}Structure','\fontsize{12}Von Neumann Entropy'})
set(gca, 'box', 'off')
text(VNETad(1),FxnTad(1)+.01,'t = 1','FontSize',14)
text(VNETad(2)+.001,FxnTad(2)-.006,'t = 2','FontSize',14)

