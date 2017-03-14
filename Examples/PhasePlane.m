%% load HiC data
clear all, close all, clc
load('chr22_100kb_RnaHiC_FibTS.mat')

%% remove centromere

[C22,~,idx_cent] = HiC_remove_cent(C22);
RNA22 = RNA22(~idx_cent,:);

%% plot all data

figure('Position',[100,20,1000,600]);
HicTensorFig( double(C22), RNA22,1);
text(1.2,length(C22)+30,0,'RNA-seq')
text(.85,length(C22)-30,-length(C22)-15,'Hi-C')
text(4.2,length(C22),5,'Time','Rotation',28);

text(.9,-10,-length(C22)-15,'t = 1')
text(1.9,-10,-length(C22)-15,'t = 2')
text(7.9,-10,-length(C22)-15,'t = 8')

%% caculate function

idx22genes = GeneLocs(:,1)==22;
RPKM_22genes = GeneRpkm(idx22genes,:);
Function = sqrt(mean(RPKM_22genes,1));

%% claculate Fiedler number from Hi-C

% remove centromeres
[C22_noCent,~,idx_cent] = HiC_remove_cent(double(C22),1);

Norm22 = NaN(size(C22_noCent));
FN22 = NaN(size(C22_noCent,3),1);
VNE22 = NaN(size(C22_noCent,3),1);
figure('Position',[20,100,900,900])
for s = 1:size(C22_noCent,3)
    % normalize
    Norm22(:,:,s) = ToepNorm(C22_noCent(:,:,s));
    % Fiedler number
    FN22(s) = Mat2Val(Norm22(:,:,s),3);
    % Von Nuemann Entropy
    VNE22(s) = Mat2Val(Norm22(:,:,s),9);
    
    subplot(3,3,s)
    HiC_plot(Norm22(:,:,s));
end

%% phase plane
% structure vs function
% Fiedler number vs mean sqrt RPKM

figure('Position',[200,400,500,500])
PltPhasePlane( FN22,Function)
% for title, add to end of previous line: ,'\fontsize{16}Chr 22 Phase Plane'
set(gca, 'box', 'off')
text(FN22(1),Function(1)+.01,'t = 1','FontSize',14)
text(FN22(2)+.001,Function(2)-.006,'t = 2','FontSize',14)
text(FN22(end)+.001,Function(end)-.006,'t = 8','FontSize',14)


%% phase plane - using entropy
% structure vs function
% Fiedler number vs mean sqrt RPKM

figure('Position',[200,400,500,500])
PltPhasePlane( VNE22,Function)
% for title, add to end of previous line: ,'\fontsize{16}Chr 22 Phase Plane'
xlabel({'\fontsize{15}Structure','\fontsize{12}Von Neumann Entropy'})
set(gca, 'box', 'off')
text(VNE22(1),Function(1)+.01,'t = 1','FontSize',14)
text(VNE22(2)+.001,Function(2)-.006,'t = 2','FontSize',14)


