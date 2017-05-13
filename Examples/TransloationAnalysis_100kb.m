clear all, close all, clc

%% load data and define what chrs were looking at
% define chrs being analyzed
chr1 = 6;
chr2 = 14;

% define files
f_100kb = 'Chr6_chr14_100kb_HT29.mat';
f_RNA = 'Chr6_chr14_100kb_RNA_HT29';

% load and re-name Hi-C data
load(f_100kb)
Mat_inter = C6_C14;
Mat_c1 = C6;
Mat_c2 = C14;

% laod and re-name RNA-seq data
load(f_RNA)
RNA1 = RNA6;
RNA2 = RNA14;

%% find translocation in 100 kb

% define where the translocation occurs
t1 = 1329;
t2 = 366;

% define axis - uncomment out one to select it
%ax = [0,size(Mat_inter,2),0,size(Mat_inter,1)]; % whole inter-chr matrix
ax = [t2-5,t2+5,t1-5,t1+5]; %zoom on translocation

% figure title
tit = ['Chr ',num2str(chr1),' Chr ',num2str(chr2),' interactions'];

% plot inter-chr matarix
figure
HiC_plot(Mat_inter,tit,2,0,0,0) 
% add lines marking translocation
line([t2,t2],[1,size(Mat_inter,1)],'Color','g','LineWidth',2)
line([1,size(Mat_inter,2)],[t1,t1],'Color','g','LineWidth',2)
xlabel(['Chr ',num2str(chr2)]), ylabel(['Chr ',num2str(chr1)])
cbar(1,1,1,[0,6])
axis(ax)

% print out a small region around the trans
Mat_inter(t1-2:t1+2,t2-2:t2+2)

%% reconstruct transcloated matrix

% changes based on qhich quadrant the translocation is in. 
dir = [-1,1]; 

[Mat_trans,RNA_trans] = TranslocHiC(Mat_c1,Mat_c2,Mat_inter,[t1,t2],dir,RNA1,RNA2);

% plot reconstruction
%   line locations change in dir changes
tit = ['t(',num2str(chr1),';',num2str(chr2),')']; % figure title
figure, HiC_plot(Mat_trans,tit);
line([0,length(Mat_trans)],[t1,t1],'Color','g')
line([t1,t1],[0,length(Mat_trans)],'Color','g')
cbar(1,1,1,[0,8])

%% remove centromeres

[NoCent1,~,idxCent1] = HiC_remove_cent(Mat_c1);
[NoCent2,~,idxCent2] = HiC_remove_cent(Mat_c2);
[NoCent_trans,~,idxCent_trans] = HiC_remove_cent(Mat_trans);

% remove centromeres from binned RNA-seq
Rna1_noCent   = RNA1(~idxCent1);
Rna2_noCent  = RNA2(~idxCent2);
RnaTrans_noCent = RNA_trans(~idxCent_trans);

%% normalize

plt = 0; % create plots showing selection of breakpoint - unnecessary because forcing the translocation breakpoint
filtMax = 0.01; % maximum filter frequency, standard for 100 kb matrices.
[norm1,brks1]   = BlockToepNorm(NoCent1,100,filtMax,plt);
brks1 % prints out breakpoints used during normalization
[norm2,brks2]  = BlockToepNorm(NoCent2,100,filtMax,plt);
brks2
[normTrans,brksTrans] = BlockToepNorm(NoCent_trans,0,1e5,plt,[1,t1,length(NoCent_trans)+1]);
brksTrans

%% calculate TADs

[TADs1,~,FV1] = TAD_Laplace(norm1,.6,3,brks1,1,0);
[TADs2,~,FV2] = TAD_Laplace(norm2,.6,3,brks2,1,0);
[TADTrans,~,FVTrans] = TAD_Laplace(normTrans,.6,3,brksTrans,1,0);

% the sign of the Fielder vector is arbitrary, so for plotting purposes, we
% always make the correlation with the RNA-seq positive.
if corr(FV1,Rna1_noCent) < 0, FV1 = - FV1; end
if corr(FV2,Rna2_noCent) < 0, FV2 = - FV2; end
if corr(FVTrans,RnaTrans_noCent) < 0, FVTrans = - FVTrans; end

%% analyze structure and function near translocations

% calcualte Von Neumann Entropy of Hi-C matrix and correlation between Hi-C
% Fiedler vector and RNA-seq for a small matrix centered on the site of
% translocation.

% select regions
siz = 3;
r1 = t1-siz:t1+siz;
r2 = t2-siz:t2+siz;

% calculate Von Neumann Entropy
VNE1 = Mat2Val(Mat_c1(r1,r1),9); % type 9 in von neumann entropy
VNE2 = Mat2Val(Mat_c2(r2,r2),9);
VNETrans = Mat2Val(Mat_trans(r1,r1),9);
% print out Von Neumann Entropy
['Von Neumann Entropy of chr ',num2str(chr1),', chr ',num2str(chr2),', and t(',...
    num2str(chr1),';',num2str(chr2),'): ',num2str(VNE1), ...
    ', ',num2str(VNE2),  ', ',num2str(VNETrans)]

% calculate Structure-Function Correlation
cor1 = corr(RNA1(r1),FV1(r1));
cor2 = corr(RNA2(r2),FV2(r2));
corTrans = corr(RNA_trans(r1),FVTrans(r1));
% print Structure-Function Correlation
['Structure-Function correlation of chr ',num2str(chr1),', chr ',num2str(chr2), ...
    ', and t(',num2str(chr2),';',num2str(chr2),'): ',num2str(cor1),', ', ...
    num2str(cor2),  ', ',num2str(corTrans)]

%% plot Hi-C and TADs

% define titles
tit1 = ['Chr ',num2str(chr1)]; % title: Chr 6
tit2 = ['Chr ',num2str(chr2)]; % title: Chr 14
titTrans= ['Chr t(',num2str(chr1),';',num2str(chr2),')'];

figure % Hi-C and TADs for Chr 6
HicRna_plot(1,1,tit1,norm1,TADs1,[],[])

% HiC and TADs for Chr6, 14, and t(6;14)
%    the first number is which horizontal position to put it in, the second in
%    how many will there be.
figure('Position',[5,150,1250,340])
HicRna_plot(1,3,tit1,norm1,TADs1,[],[])
HicRna_plot(2,3,tit2,norm2,TADs2,[],[])
HicRna_plot(3,3,titTrans,normTrans,TADTrans,[],[])

%% plot HiC and RNA-seq

figure %Chr 6 only
HicRna_plot(1,1,tit1,norm1,TADs1,Rna1_noCent,[])

%% plot HiC and FV

figure('Position',[5,150,1250,420])
HicRna_plot(1,3,tit1,norm1,TADs1,[],FV1)
HicRna_plot(2,3,tit2,norm2,TADs2,[],FV2)
HicRna_plot(3,3,titTrans,normTrans,TADTrans,[],FVTrans)

%% plot HiC, RNA-seq, and FV

figure
HicRna_plot(1,1,tit1,norm1,TADs1,Rna1_noCent,FV1)

figure('Position',[100,100,1100,600])
HicRna_plot(1,2,tit1,norm1,TADs1,Rna1_noCent,FV1)
HicRna_plot(2,2,tit2,norm2,TADs2,Rna2_noCent,FV2)

figure('Position',[5,100,1250,420])
HicRna_plot(1,3,tit1,norm1,TADs1,Rna1_noCent,FV1)
HicRna_plot(2,3,tit2,norm2,TADs2,Rna2_noCent,FV2)
HicRna_plot(3,3,titTrans,normTrans,TADTrans,RnaTrans_noCent,FVTrans)
