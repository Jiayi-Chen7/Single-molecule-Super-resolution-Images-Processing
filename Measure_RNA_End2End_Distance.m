clear all
close all
clc

%%
%initial parameters
%folder
% in_folder = 'E:/CJY/20220314/EM1/';
% mid_sur_in_folder = 'E:/CJY/20220314/EM1/MID SUR/';
% out_folder =  'E:/CJY/20220314/EM1/output/';
%the image names of zstack
AF594_channel = 'nc14 Em1 S1.lif_L5P*_Crop001_ch02.tif';
ATTO647N_channel = 'nc14 Em1 S1.lif_L5P*_Crop001_ch01.tif';
% image list
imlist_AF594 = dir(AF594_channel);
imlist_ATTO647N = dir(ATTO647N_channel);
size_in_merge_image = [1024,1024];


%% read in raw data

AF594_raw_3D = zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_AF594));
ATTO647N_raw_3D =zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_ATTO647N));

for image_I = 1:length(imlist_AF594)
    
    AF594_raw_3D(:,:,image_I) = imread([imlist_AF594(image_I).name]);
    ATTO647N_raw_3D(:,:,image_I) = imread([imlist_ATTO647N(image_I).name]);
end


%%image mask of nuclei in high magnitude image
RefineZ = (numel(imlist_AF594));



ExRNACore = 1.2%1.2;
IxRNACore = 2.2%2.2;
%%
%L5
tic  
% figure();
L5Distance_200 = [];
L5Pon_100_200 = [];
L5Pon_70_200 = [];
L5Pon_66_200 = [];
L5Pon_28_200 = [];
L5IntactRatio3End_200 = [];
L5IntactRatio5End_200 = [];
% TruncatedtL5Distance_200 = [];

L5Distance_300 = [];
L5Pon_100_300 = [];
L5Pon_70_300 = [];
L5Pon_66_300 = [];
L5Pon_28_300 = [];
L5IntactRatio3End_300 = [];
L5IntactRatio5End_300 = [];
% TruncatedtL5Distance_300 = [];
% EnhanceProcessing = 1;
BypassGauFilter = 0;
RemoveNoise = 1;
PixelSize = 17;
ManualThre = 1;
AdjDynamicRange = 1;
for RNAStackIndex = 1:RefineZ
RNAStackIndex 
SingleAF594 = AF594_raw_3D(:,:,RNAStackIndex);
SingleATTO647N = ATTO647N_raw_3D(:,:,RNAStackIndex);

VolThre = 30;
MaxIntThre = 10000;
LowMinIntThre = 380;
ThreAdj = 0;
IniGuess = [1500,300,50];
DoGSize = 15;
GaussainManualScale = 3;
 [CytoRNAMask_trueAF594,GauFilAF594,CytoRNACenTrueAF594,AdjRNARaw594] = SinglePlaneRNAMaskMeasureDistance4(SingleAF594,DoGSize,ExRNACore,IxRNACore,RNAStackIndex,AdjDynamicRange,[0 0.02],ThreAdj,[],VolThre,MaxIntThre,LowMinIntThre,GaussainManualScale);

 
VolThre = 20;
MaxIntThre = 10000;
LowMinIntThre = 350;
ThreAdj = 0;
IniGuess = [1500,250,50];
DoGSize = 12;
GaussainManualScale = 3;
[CytoRNAMask_trueATTO647N,GauFilATTO647N,CytoRNACenTrueATTO647N,AdjRNARaw647] = SinglePlaneRNAMaskMeasureDistance4(SingleATTO647N,DoGSize,ExRNACore,IxRNACore,RNAStackIndex,AdjDynamicRange,[0 0.02],ThreAdj,[],VolThre,MaxIntThre,LowMinIntThre,GaussainManualScale);
NANInd = isnan(CytoRNACenTrueAF594);
SumNANInd = sum(NANInd,2);
PassInd = find(SumNANInd<1);
CytoRNACenTrueAF594 = CytoRNACenTrueAF594(PassInd,1:2);

NANInd2 = isnan(CytoRNACenTrueATTO647N);
SumNANInd2 = sum(NANInd2,2);
PassInd2 = find(SumNANInd2<1);
CytoRNACenTrueATTO647N = CytoRNACenTrueATTO647N(PassInd2,1:2);

P = CytoRNACenTrueATTO647N ;
% T = delaunayn(P);
PQ = CytoRNACenTrueAF594;
%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance
% [kRNA2GFPGFP distRNA2GFP] = dsearchn(P,T,PQ);
[kRNA2GFPGFP distRNA2GFP] = dsearchn(P,PQ);

% [kRNA2GFPGFPUnique,ia,ic] = unique(kRNA2GFPGFP,'rows');
% distRNA2GFP2 = distRNA2GFP(ia);
%determine colocalization
%300 nm as the co-localization distance threshold
distRNA2GFPScale = distRNA2GFP.*PixelSize;

h = 2.*iqr(distRNA2GFPScale)./(length(distRNA2GFPScale).^(1/3));
n = (max(distRNA2GFPScale)-min(distRNA2GFPScale))./h;
% histogram(distRNA2GFPScale,round(n));hold on
% ColThre = ceil(200./69);
% Temp = find(distRNA2GFP < ColThre);
% TempIso = find(distRNA2GFP >= ColThre);

Pon_100_200 = length(distRNA2GFPScale(distRNA2GFPScale<=200 & distRNA2GFPScale>100))./length(distRNA2GFPScale(distRNA2GFPScale<=200));
Pon_70_200 = length(distRNA2GFPScale(distRNA2GFPScale<=200 & distRNA2GFPScale>70))./length(distRNA2GFPScale(distRNA2GFPScale<=200));
Pon_66_200 = length(distRNA2GFPScale(distRNA2GFPScale<=200 & distRNA2GFPScale>66))./length(distRNA2GFPScale(distRNA2GFPScale<=200));
Pon_28_200 = length(distRNA2GFPScale(distRNA2GFPScale<=200 & distRNA2GFPScale>28))./length(distRNA2GFPScale(distRNA2GFPScale<=200));
L5Pon_100_200 = [L5Pon_100_200;Pon_100_200];
L5Pon_70_200 = [L5Pon_70_200;Pon_70_200];
L5Pon_66_200 = [L5Pon_66_200;Pon_66_200];
L5Pon_28_200 = [L5Pon_28_200;Pon_28_200];
L5Distance_200 = [L5Distance_200;distRNA2GFPScale(distRNA2GFPScale<=200)];

%%
%300 nm threshold
Pon_100_300 = length(distRNA2GFPScale(distRNA2GFPScale<=300 & distRNA2GFPScale>100))./length(distRNA2GFPScale(distRNA2GFPScale<=300));
Pon_70_300 = length(distRNA2GFPScale(distRNA2GFPScale<=300 & distRNA2GFPScale>70))./length(distRNA2GFPScale(distRNA2GFPScale<=300));
Pon_66_300 = length(distRNA2GFPScale(distRNA2GFPScale<=300 & distRNA2GFPScale>66))./length(distRNA2GFPScale(distRNA2GFPScale<=300));
Pon_28_300 = length(distRNA2GFPScale(distRNA2GFPScale<=300 & distRNA2GFPScale>28))./length(distRNA2GFPScale(distRNA2GFPScale<=300));
L5Pon_100_300 = [L5Pon_100_300;Pon_100_300];
L5Pon_70_300 = [L5Pon_70_300;Pon_70_300];
L5Pon_66_300 = [L5Pon_66_300;Pon_66_300];
L5Pon_28_300 = [L5Pon_28_300;Pon_28_300];
L5Distance_300 = [L5Distance_300;distRNA2GFPScale(distRNA2GFPScale<=300)];

end

save L5Distance_nc14_Em1.mat L5Distance_200 L5Pon_28_200 L5Pon_66_200 L5Pon_70_200 L5Pon_100_200 L5Distance_300 L5Pon_28_300 L5Pon_66_300 L5Pon_70_300 L5Pon_100_300 
toc

function [CytoRNAMask_true,AdjGauOutims,CytoRNACenTrue,AdjRNARaw1] = SinglePlaneRNAMaskMeasureDistance4(SinglemRNA,DoGSize,ExRNACore,IxRNACore,IndexOfZstack,AdjDynamicRange,InputDynamicRange,ThreAdj,AdjVal,VolThre,MaxIntThre,LowMinIntThre,GaussainManualScale)
%%


Ex = fspecial('gaussian',10,ExRNACore);%10
Ix = fspecial('gaussian',20,IxRNACore);%20
outE = imfilter(single(SinglemRNA),Ex,'replicate'); 
outI = imfilter(single(SinglemRNA),Ix,'replicate'); 
outims = outE - outI;  

LowHighRNARaw1 = stretchlim(uint16(SinglemRNA))
if AdjDynamicRange == 0
AdjRNARaw1 = imadjust(uint16(SinglemRNA),LowHighRNARaw1,[]);
else
   AdjRNARaw1 = imadjust(uint16(SinglemRNA),InputDynamicRange,[]);
end 


% figure();
% imshow(AdjRNARaw1)
% figure();
% imshow(outims,[]);

Gauoutims = imgaussfilt(single(outims),1);
% figure();
% imshow(Gauoutims);
LowHighRNA = stretchlim(uint16(Gauoutims));
AdjGauOutims = imadjust(uint16(Gauoutims),LowHighRNA,[]);
% figure();
% imshow(AdjGauOutims)



[T EM] = graythresh(AdjGauOutims);
if ThreAdj == 1
   T = T+AdjVal;
   RawRNAMask = imbinarize(AdjGauOutims,T);
else
    RawRNAMask = imbinarize(AdjGauOutims,T);
end
% figure();
% imshow(RawRNAMask)

% bw2 = ~bwareaopen(~RawRNAMask,3);
% D = bwdist(~bw2);
% D = -D;
% %local minimum mark
% LocalMinMask = imextendedmin(D,2);
% %varified local minima location 
% % imshowpair(D,LocalMinMask,'blend');
% % Modify the distance transform so it only has minima at the desired locations,
% D2 = imimposemin(D,LocalMinMask);
% L = watershed(D2);
% %the logical ture of the watershed segmentation located in the region of
% %the raw RNA mask
% ZerosCross = RawRNAMask & (L==0);
% %assign the position find above to 0 in the mask image 
% RawRNAMask(ZerosCross == 1) = 0;


% [B,L] = bwboundaries(RawRNAMask);
% Bdr_fig=zeros(2048);
% for k = 1:length(B)
%    boundary = B{k};
%    for j=1:length(boundary)
%     Bdr_fig(boundary(j,1),boundary(j,2)) =1;
%    end
% end
%  
% figure();
% Com = imoverlay(uint16(SinglemRNA),Bdr_fig,'red');
% imshow(Com);



%remove the single mRNA in the nuclei region and transcription foci

% RawRNAMaskDilate = imdilate(RawRNAMask,strel("disk",1));
RawRNAMaskOpen = imopen(RawRNAMask,strel("disk",1));
CytoRNAMask = RawRNAMaskOpen;
% figure();
% imshow(CytoRNAMask)
% class(CytoRNAMask)
% class(AdjOutims)
%%
%plot
% % RawRNAAfterMask = double(SinglemRNA).*CytoRNAMask;
% % LowHighRNARaw = stretchlim(uint16(RawRNAAfterMask));
% % AdjRNARaw = imadjust(uint16(RawRNAAfterMask),LowHighRNARaw,[]);
% % figure();
% % imshow(AdjRNARaw(1:300,1:300));

%%
CytoRNALL = bwlabel(CytoRNAMask);
CytoRNALLsts = regionprops3(CytoRNALL,SinglemRNA,'WeightedCentroid','Volume','MeanIntensity','MaxIntensity');
CytoRNACen = [CytoRNALLsts.WeightedCentroid];
CytoRNAVol = [CytoRNALLsts.Volume];
CytoRNAMeanInt = [CytoRNALLsts.MeanIntensity];
CytoRNAMaxInt = [CytoRNALLsts.MaxIntensity];
% figure();
% histogram(CytoRNAVol);
figure();
histogram(CytoRNAMeanInt);

%Determination of the initial guess for automatic fitting 
[Val2, Edges2] = histcounts(CytoRNAMeanInt);
[M I] = max(Val2);

para02 = [max(Val2),max(Edges2(I)),std(Val2)];
ThreeTermGau = @(x,xdata) x(1)*exp(-((xdata - x(2)).^2./(2.*x(3).^2)));
lb2 = [0,0,0];
ub2 = [inf,inf,inf];
edges2 = Edges2(2:end) - (Edges2(2)-Edges2(1))/2;
options = optimset('Display','off');
[paraoutputThreeTerm,resnormoutputThreeTerm,~,exitflagoutputThreeTerm,~] = lsqcurvefit(ThreeTermGau,para02,edges2,Val2,lb2,ub2,options);

%plot
[x, y] = prepareCurveData(edges2,Val2);

GauTerm1 =  paraoutputThreeTerm(1)*exp(-((x - paraoutputThreeTerm(2)).^2./(2.*paraoutputThreeTerm(3).^2)));
figure();
plot(x,GauTerm1,'LineWidth',3);
hold on
scatter(x,y,30,"black","filled");
hold off
xlabel("Mean Intensity (a.u.)");
set (gca,'linewidth',2,'fontsize',12,'FontName','Arial');
ylabel("The number of spots");
set (gca,'linewidth',2,'fontsize',12,'FontName','Arial');
box on;
axis square
set(gcf,'position',[100 100 500 500]);
lgd = legend('Fitted','Raw data points');
lgd.FontSize = 12;
lgd.FontName = 'Arial';

MinIntThrePre = paraoutputThreeTerm(2)+GaussainManualScale.*paraoutputThreeTerm(3)
if (MinIntThrePre<LowMinIntThre)
    MinIntThre = LowMinIntThre;
else
    MinIntThre = MinIntThrePre;
end

% Q = quantile(Val2,0.95)
% MinIntThre = Q;
IndEnhance = find(CytoRNAMeanInt >= MinIntThre & CytoRNAVol>VolThre & CytoRNAMeanInt < MaxIntThre);

    CytoRNAMask_true = ismember(CytoRNALL,IndEnhance);
   CytoRNATrueLL = bwlabel(CytoRNAMask_true);
CytoRNATrueLLsts = regionprops3(CytoRNATrueLL,SinglemRNA,'WeightedCentroid','Volume','MeanIntensity','MaxIntensity');
CytoRNACenTrue = [CytoRNATrueLLsts.WeightedCentroid];


figure();
imshow(AdjRNARaw1);

figure();
imshow(AdjRNARaw1);
hold on
scatter(CytoRNACenTrue(:,1),CytoRNACenTrue(:,2),15,'m','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);


end