clear all, close all, clc

%% find translocation in 100 kb

load('Chr6_chr14_100kb_HT29.mat')

figure
HiC_plot(C6_C14,'Chr 6 Chr 14 interactions',2,0,0,0) % last 0
line([366,366],[1,size(C6_C14,1)],'Color','g','LineWidth',4)
line([1,size(C6_C14,2)],[1329,1329],'Color','g','LineWidth',4)
xlabel('Chr14'), ylabel('Chr6')
cbar(1,1,1,[0,6])
axis([361,371,1323,1334]) % comment to view entire inter-chr region

C6_C14(1326:1333,363:369)

%% find translocation in read data
load('chr6_chr14_reads_HT29.mat')

idx1 = chrA == 6 & chrB == 14;
idx2 = chrA == 14 & chrB == 6;
locChr6 = [locA(idx1);locB(idx2)];
locChr14= [locB(idx1);locA(idx2)];

reg6 = [132500000,133200000];
reg14 = [36000000,37000000];

figure
scatter(locChr6,locChr14,'.')
set(gca,'Ydir','reverse')
axis([reg6,reg14])
line([132825000,132825000],reg14,'Color','r','LineWidth',2)
line([132890000,132890000],reg14,'Color','r','LineWidth',2)
line(reg6,[36508800,36508800],'Color','r','Color','r','LineWidth',2)
xlabel('Chr 6'), ylabel('Chr 14'), title('Chr6-Chr14')

clear idx1 idx2 reg6 reg14 locChr6 locChr14 chrA chrB locA locB

%% reconstruct transcloated matrix

load('chr6_chr14_100kb_RNA_HT29')
brk6 = 1329;
brk14 = 366;
dir = [-1,1];

[T614,RNA614] = TranslocHiC(C6,C14,C6_C14,[brk6,brk14],dir,RNA6,RNA14);

% plot reconstruction
%   line locations change in dir changes
figure, HiC_plot(T614,'t(6;14)');
line([0,length(T614)],[brk6,brk6],'Color','g')
line([brk6,brk6],[0,length(T614)],'Color','g')
cbar(1,1,1,[0,8])

%% remove centromeres

[NoCent6,~,idxCent6]     = HiC_remove_cent(C6);
[NoCent14,~,idxCent14]   = HiC_remove_cent(C14);
[NoCent614,~,idxCent614] = HiC_remove_cent(T614);

Rna6_noCent   = RNA6(~idxCent6);
Rna14_noCent  = RNA14(~idxCent14);
Rna614_noCent = RNA614(~idxCent614);

%% normalize

plt = 0;
[norm6,brks6]   = BlockToepNorm(NoCent6,100,1e5,plt);
brks6
[norm14,brks14]  = BlockToepNorm(NoCent14,100,1e5,plt);
brks14
[norm614,brks614] = BlockToepNorm(NoCent614,0,1e5,plt,[1,brk6,length(NoCent614)+1]);
brks614

%% calculate TADs

[TADs6,~,FV6]   = TAD_Laplace(norm6,.6,3,brks6,1,0);
[TADs14,~,FV14]  = TAD_Laplace(norm14,.6,3,brks14,1,0);
[TADs614,~,FV614] = TAD_Laplace(norm614,.6,3,brks614,1,0);

if corr(FV6,Rna6_noCent) < 0, FV6 = - FV6; end
if corr(FV14,Rna14_noCent) < 0, FV14 = - FV14; end
if corr(FV614,Rna614_noCent) < 0, FV614 = - FV614; end

%% analyze structure and function near translocations

% calcualte Von Neumann Entropy of Hi-C matrix and correlation between Hi-C
% Fiedler vector and RNA-seq for a small matrix centered on the site of
% translocation.

% select regions
siz = 3;
r6 = brk6-siz:brk6+siz;
r14 = brk14-siz:brk14+siz;

% calculate Von Neumann Entropy
VNE6 = Mat2Val(C6(r6,r6),9); % type 9 in von neumann entropy
VNE14 = Mat2Val(C14(r14,r14),9);
VNE614 = Mat2Val(T614(r6,r6),9);
['Von Neumann Entropy of chr 6, chr 14, and t(6;14): ',num2str(VNE6), ...
    ', ',num2str(VNE14),  ', ',num2str(VNE614)]

% calculate Structure-Function Correlation
cor6 = corr(RNA6(r6),FV6(r6));
cor14 = corr(RNA14(r14),FV14(r14));
cor614 = corr(RNA614(r6),FV614(r6));
['Structure-Function correlation of chr 6, chr 14, and t(6;14): ',num2str(cor6), ...
    ', ',num2str(cor14),  ', ',num2str(cor614)]

%% plot Hi-C and TADs

figure
HicRna_plot(1,1,'Chr 6',norm6,TADs6,[],[])

figure('Position',[5,150,1250,340])
HicRna_plot(1,3,'Chr 6',norm6,TADs6,[],[])
HicRna_plot(2,3,'Chr 14',norm14,TADs14,[],[])
HicRna_plot(3,3,'Chr t(6;14)',norm614,TADs614,[],[])

%% plot HiC and RNA-seq

figure
HicRna_plot(1,1,'Chr 6',norm6,TADs6,Rna6_noCent,[])

%% plot HiC and FV


figure('Position',[5,150,1250,420])
HicRna_plot(1,3,'Chr 6',norm6,TADs6,[],FV6)
HicRna_plot(2,3,'Chr 14',norm14,TADs14,[],FV14)
HicRna_plot(3,3,'Chr t(6;14)',norm614,TADs614,[],FV614)

%% plot HiC, RNA-seq, and FV

figure
HicRna_plot(1,1,'Chr 6',norm6,TADs6,Rna6_noCent,FV6)

figure('Position',[100,100,1100,600])
HicRna_plot(1,2,'Chr 6',norm6,TADs6,Rna6_noCent,FV6)
HicRna_plot(2,2,'Chr 14',norm14,TADs14,Rna14_noCent,FV14)

figure('Position',[5,100,1250,420])
HicRna_plot(1,3,'Chr 6',norm6,TADs6,Rna6_noCent,FV6)
HicRna_plot(2,3,'Chr 14',norm14,TADs14,Rna14_noCent,FV14)
HicRna_plot(3,3,'Chr t(6;14)',norm614,TADs614,Rna614_noCent,FV614)
