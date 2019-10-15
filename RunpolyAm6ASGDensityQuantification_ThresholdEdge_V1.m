%%This code is used to quantify SG/PB and the distribution of polyA and m6A signal in SG/PB with features of single cell cytoplasm selection
% Ye Fu
% yefu01@fas.harvard.edu
% Harvard University, 
% 2019
% for four channels

%% Define path to analysis codes
matlabStormPath = 'N:\MATLAB_Analysis\';
%% Define paths and paramters for images to be analyzed
mainDirPath = 'Q:\STORM2_Data\20190829_U2OS+NaAsO2+BlueLight_siRNA_pYF127_128_BFP-405_DCP1A-488_G3BP1-561_Y3-680_Y1_BFP-750\A4_siC_pYF128_750-BFP\';

% Image parameters

analysisPath = mainDirPath;

%% Add analysis path
disp('Adding analysis code path...');
addpath(matlabStormPath, '-begin');
disp(['Adding matlab-storm: ' matlabStormPath])
disp('------------------------------------------------------------------');
subfolderPaths = genpath(matlabStormPath); 
disp('Adding subfolder Paths');
addpath(subfolderPaths);
disp('------------------------------------------------------------------');



% %% Step 1: Threshold determination %Comment this section for SG identification
% analysisSavePath = [mainDirPath, 'Overlapanalysis_ThresholdEdge\'];
% disp(['Set analysis folder as: ' analysisSavePath]);
% if ~exist(analysisSavePath)
%     mkdir(analysisSavePath)
% end
% 
% 
% 
% %%Run threshold
% thresholdAll = [];
% 
% for i = 0 : 9
% try threshold = FuncCalThreshold_Edge ('561Conv_', ['000', num2str(i)], mainDirPath);
% 
%     thresholdAll = [thresholdAll; threshold];
% end
% end
%  
% for i = 10 : 83
% try threshold = FuncCalThreshold_Edge ('561Conv_00', [num2str(i)], mainDirPath);
%     thresholdAll = [thresholdAll; threshold];
% end
% end
% 
% threshold_median = median(thresholdAll);
% xlswrite([analysisSavePath 'threshold.csv'],[threshold_median; thresholdAll]);
% 
%  %find threshold_median
% msgbox('Threshold calculation done');
% pause;
% % 
% 
% threshold_use = 100; %Median of all samples
% 
% %% To evaluate the edge threshold before Cell Selection 
% 
% for i = 0 : 9
% try FuncThreshold_Edge ('561Conv_', ['000', num2str(i)], threshold_use, mainDirPath);
% %4 color
% end
% end
% 
% 
% for i = 10 : 50
% try FuncThreshold_Edge ('561Conv_00', [num2str(i)], threshold_use, mainDirPath);
% end
% end
% 
% 
% for i = 50 : 90
% try FuncThreshold_Edge ('561Conv_00', [num2str(i)], threshold_use, mainDirPath);
% end
% end
% msgbox('Threshold done');
% 
% %% Step 2: Segment cell and calculate SG
% threshold_use = 200;
% FuncpolyAm6ASGDensityQuantification_ThresholdEdge ('00', threshold_use, mainDirPath); % Run one first to check if there is any error
% for i = 1 : 9
% try FuncpolyAm6ASGDensityQuantification_ThresholdEdge (['0', num2str(i)], threshold_use, mainDirPath);
% end
% end
% 
% 
% for i = 10 : 50
% try FuncpolyAm6ASGDensityQuantification_ThresholdEdge ( ['',num2str(i)], threshold_use, mainDirPath);
% end
% end
% 
% for i = 50 : 90
% try FuncpolyAm6ASGDensityQuantification_ThresholdEdge ( ['',num2str(i)], threshold_use, mainDirPath);
% end
% end

% %% Step 2: Segment cell and calculate SG
% threshold_use = 200;
%  FuncpolyAm6ASGDensityQuantification_ThresholdEdge_RGB ('00', threshold_use, mainDirPath); % Run one first to check if there is any error
% for i = 1 : 9
% try FuncpolyAm6ASGDensityQuantification_ThresholdEdge_RGB (['0', num2str(i)], threshold_use, mainDirPath); %No 750
% end
% end
% 
% 
% for i = 10 : 50
% try FuncpolyAm6ASGDensityQuantification_ThresholdEdge_RGB ( ['',num2str(i)], threshold_use, mainDirPath);
% end
% end

%% Step 2: Segment cell and calculate SG
threshold_use = 100;
%FuncSGStatSNAP_ThresholdEdge_NoPB ('0001', threshold_use, mainDirPath); % Run one first to check if there is any error
for i = 0 : 9
try FuncSGStatSNAP_ThresholdEdge_NoPB (['000', num2str(i)], threshold_use, mainDirPath); %for SNAP
end
end

for i = 10 : 50
try FuncSGStatSNAP_ThresholdEdge_NoPB ( ['00',num2str(i)], threshold_use, mainDirPath);
end
end

for i = 50 : 90
try FuncSGStatSNAP_ThresholdEdge_NoPB ( ['00',num2str(i)], threshold_use, mainDirPath);
end
end

% 
% %% Step 2: Segment cell and calculate SG
% threshold_use = 100;
% FuncSGStatCry2O_ThresholdEdge ('01', threshold_use, mainDirPath); % Run one first to check if there is any error
% for i = 0 : 9
% try FuncSGStatCry2O_ThresholdEdge(['0', num2str(i)], threshold_use, mainDirPath); %for SNAP
% end
% end
% 
% for i = 10 : 50
% try FuncSGStatCry2O_ThresholdEdge( num2str(i), threshold_use, mainDirPath);
% end
% end
% msgbox('Analysis done');


function threshold = FuncCalThreshold_Edge (SGname, SerNo, mainDirPath)

SGImageName = [SGname, SerNo];

%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);

%% Process data
SGimage(88,405) = SGimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background
[~,threshold] = edge(SGimage1,'sobel');
end



function FuncThreshold_Edge (SGname, SerNo, threshold, mainDirPath)
%% Define path to analysis codes


fudgeFactor = 2.5; %Edge threshold multiplier
size_threshold = 10; %bwopen size threshold
SGImageName = [SGname SerNo];

% Image parameters
analysisSavePath = [mainDirPath, ['ThresholdEdge_Evaluation_' num2str(threshold) '_' num2str(fudgeFactor) '_' num2str(size_threshold) '\']];

%% Set up analysis folder

if ~exist(analysisSavePath)
    mkdir(analysisSavePath)
end




%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);


%% Process data
SGimage(88,405) = SGimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background

BWs = edge(SGimage1,'sobel',threshold * fudgeFactor);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', 2);
BWer1 = imerode(BWnobord,seD);
xbw = bwareaopen(BWer1, size_threshold);

g1 = imshow(SGimage, [0, max(max(SGimage))]);
saveas (g1, [analysisSavePath SGImageName '_SGImage.png']);
close(gcf);

g2 = imshow(xbw);
saveas (g2, [analysisSavePath SGImageName '_SGMask.png']);
close(gcf);
end

function FuncpolyAm6ASGDensityQuantification_ThresholdEdge (SerNo, threshold, mainDirPath)
%Use threshold, threshold = 0, and threshold = 5 here. 

m6AImageName = ['2_647Conv_' SerNo];
PAImageName = ['1_750Conv_' SerNo]; %PolyA image name 
SGImageName = ['3_561Conv_' SerNo];
PBImageName = ['4_488Conv_' SerNo];

% m6AImageName = ['647Conv_' SerNo];
% PAImageName = ['750Conv_' SerNo]; %PolyA image name 
% SGImageName = ['561Conv_' SerNo];
% PBImageName = ['488Conv_' SerNo];

% Image parameters
pixelSize = 153; % nm
background750 = 150; %sCMOS background
background647 = 150;
fudgeFactor = 2.5;
size_threshold = 10;
analysisSavePath = [mainDirPath, ['ThresholdEdge_Evaluation_' num2str(threshold) '_' num2str(fudgeFactor) '_' num2str(size_threshold) '\']];
cd(analysisSavePath);

%% Further define image parameters
% parameters for visualizing conventional images
pixelNum = 600;


%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[PBimage, ~] = ReadDaxL([mainDirPath PBImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[PAimage, ~] = ReadDaxL([mainDirPath PAImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[m6Aimage, ~] = ReadDaxL([mainDirPath m6AImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);



%% Process data
SGimage(88,405) = SGimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background
BWs = edge(SGimage1,'sobel',threshold * fudgeFactor);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', 2);
BWer1 = imerode(BWnobord,seD);
xbw = bwareaopen(BWer1, size_threshold); %Only select SG larger than 10 pixel after dilation
xbw0 = bwareaopen(BWer1, 2);
xbw1 = bwareaopen(BWer1, 5);

PBbg = imopen(PBimage,strel('disk',10));
PBimage1 = imsubtract(PBimage, PBbg);
PBimage2 = double(PBimage1)/max(max(double(PBimage1))); %Convert to double and 0 to 1 range.
ybw = imbinarize(PBimage2,graythresh(PBimage2)); %Pbody mask

PAbg = uint16(background750*ones(size(PAimage)));
PAimage = imsubtract(PAimage,PAbg);
m6Abg = uint16(background647*ones(size(m6Aimage)));
m6Aimage = imsubtract(m6Aimage,m6Abg);

%% Map
SGPortion = double(PAimage).*xbw;
NonSGPortion = double(PAimage).*(~xbw);

SGPortion0 = double(PAimage).*xbw0;
NonSGPortion0 = double(PAimage).*(~xbw0);

SGPortion1 = double(PAimage).*xbw1;
NonSGPortion1 = double(PAimage).*(~xbw1);

PBPortion = double(PAimage).*ybw;
NonPBPortion = double(PAimage).*(~ybw);

m6ASGPortion = double(m6Aimage).*xbw;
m6ANonSGPortion = double(m6Aimage).*(~xbw);

m6ASGPortion0 = double(m6Aimage).*xbw0;
m6ANonSGPortion0 = double(m6Aimage).*(~xbw0);

m6ASGPortion1 = double(m6Aimage).*xbw1;
m6ANonSGPortion1 = double(m6Aimage).*(~xbw1);

m6APBPortion = double(m6Aimage).*ybw;
m6ANonPBPortion = double(m6Aimage).*(~ybw);

%% XY projection image for ROI determination
% Select ROI
g = figure; imshow(PAimage, [0, max(max(PAimage))]); colormap hot;
WhetherROI = questdlg('Do you want to select single cell?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
wait(h);
i = 1;
CellMask{i} = createMask(h);
WhetherNucl = questdlg('Do you want to select nucleus?'); %ask question
    if(strcmp(WhetherNucl, 'Yes'))
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    i = 1;
    NuclMask{i} = createMask(h);
    end
moreROI = questdlg('Do you want to select more cells?'); %ask question
%keep looping till the user select 'No'. 
while(strcmp(moreROI, 'Yes'))
    i = i + 1;
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    CellMask{i} = createMask(h);
    WhetherNucl = questdlg('Do you want to select nucleus?'); %ask question
    if(strcmp(WhetherNucl, 'Yes'))
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    NuclMask{i} = createMask(h);
    end
    moreROI = questdlg('Do you want to select more cells?');
end
% else
%     CellMask{1} = ones(size(xbw)); %If commented, no output. Otherwise
%     use full image to analyze.
%     NuclMask{1} = zeros(size(xbw));
end

saveas (g, [analysisSavePath PAImageName '_polyA.png']);
close(g);
g1 = imshow(SGimage, [0, max(max(SGimage))]);
saveas (g1, [analysisSavePath PAImageName '_SGImage.png']);
close(gcf);
g2 = imshow(xbw);
saveas (g2, [analysisSavePath PAImageName '_SGMask.png']);
close(gcf);

%% Quantify PA Intensity Ratio in SG and PB 
for i = 1:length(CellMask)
xbwROI{i} = xbw.*CellMask{i};
xbw0ROI{i} = xbw0.*CellMask{i};
xbw1ROI{i} = xbw1.*CellMask{i};
ybwROI{i} = ybw.*CellMask{i};
%% feature extraction - size distribution (area, pixels)
SG{i} = bwconncomp(xbwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats{i} = regionprops(SG{i});
As{i} = [SGstats{i}.Area];


%statistical measurements
if isempty(As{i})
    A(i).SGmean = 0;
    A(i).SGstd = 0;
    A(i).SGmedian = 0;
else
    A(i).SGmean = mean(As{i});
    A(i).SGstd = std(As{i});
    A(i).SGmedian = median(As{i});
end

A(i).SGtotalarea = sum(As{i});
A(i).SGNo = length(As{i});
A(i).SGtotalIntensity = sum(sum(SGimage1.*xbwROI{i}));
A(i).G3BP1totalIntensity = sum(sum(SGimage1.*(CellMask{i})));
A(i).G3BP1SGRatio = A(i).SGtotalIntensity/A(i).G3BP1totalIntensity;    
A(i).G3BP1mean = A(i).G3BP1totalIntensity/sum(sum(SGimage1.*CellMask{i}>0));

SG0{i} = bwconncomp(xbw0ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats0{i} = regionprops(SG0{i});
As0{i} = [SGstats0{i}.Area];

%statistical measurements
if isempty(As0{i})
    A(i).SG0mean = 0;
    A(i).SG0std = 0;
    A(i).SG0median = 0;
else
    A(i).SG0mean = mean(As0{i});
    A(i).SG0std = std(As0{i});
    A(i).SG0median = median(As0{i});
end
A(i).SG0totalarea = sum(As0{i});
A(i).SG0No = length(As0{i});
A(i).SG0totalIntensity = sum(sum(SGimage1.*xbw0ROI{i}));
A(i).G3BP10SGRatio = A(i).SG0totalIntensity/A(i).G3BP1totalIntensity;   

SG1{i} = bwconncomp(xbw1ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats1{i} = regionprops(SG1{i});
As1{i} = [SGstats1{i}.Area];
%statistical measurements
if isempty(As1{i})
    A(i).SG1mean = 0;
    A(i).SG1std = 0;
    A(i).SG1median = 0;
else
    A(i).SG1mean = mean(As1{i});
    A(i).SG1std = std(As1{i});
    A(i).SG1median = median(As1{i});
end
A(i).SG1totalarea = sum(As1{i});
A(i).SG1No = length(As1{i});
A(i).SG1totalIntensity = sum(sum(SGimage1.*xbw1ROI{i}));
A(i).G3BP11SGRatio = A(i).SG1totalIntensity/A(i).G3BP1totalIntensity;   

PASG{i} = SGPortion.*CellMask{i}.*(~NuclMask{i});
A(i).PASGtotal = sum(sum(PASG{i}));
A(i).SGsize = sum(sum(PASG{i}~=0));
A(i).PASGmean = A(i).PASGtotal/A(i).SGsize;
PANonSG{i} = NonSGPortion.*CellMask{i}.*(~NuclMask{i});
A(i).PANonSGtotal = sum(sum(PANonSG{i}));
A(i).NonSGsize = sum(sum(PANonSG{i}~=0));
A(i).PANonSGmean = A(i).PANonSGtotal/A(i).NonSGsize;
A(i).PASGRatio = A(i).PASGtotal/(A(i).PASGtotal+A(i).PANonSGtotal);

PASG0{i} = SGPortion0.*CellMask{i}.*(~NuclMask{i});
A(i).PASG0total = sum(sum(PASG0{i}));
A(i).SG0size = sum(sum(PASG0{i}~=0));
A(i).PASG0mean = A(i).PASG0total/A(i).SG0size;
PANonSG0{i} = NonSGPortion0.*CellMask{i}.*(~NuclMask{i});
A(i).PANonSG0total = sum(sum(PANonSG0{i}));
A(i).NonSG0size = sum(sum(PANonSG0{i}~=0));
A(i).PANonSG0mean = A(i).PANonSG0total/A(i).NonSG0size;
A(i).PASG0Ratio = A(i).PASG0total/(A(i).PASG0total+A(i).PANonSG0total);

PASG1{i} = SGPortion1.*CellMask{i}.*(~NuclMask{i});
A(i).PASG1total = sum(sum(PASG1{i}));
A(i).SG1size = sum(sum(PASG1{i}~=0));
A(i).PASG1mean = A(i).PASG1total/A(i).SG1size;
PANonSG1{i} = NonSGPortion1.*CellMask{i}.*(~NuclMask{i});
A(i).PANonSG1total = sum(sum(PANonSG1{i}));
A(i).NonSG1size = sum(sum(PANonSG1{i}~=0));
A(i).PANonSG1mean = A(i).PANonSG1total/A(i).NonSG1size;
A(i).PASG1Ratio = A(i).PASG1total/(A(i).PASG1total+A(i).PANonSG1total);


PB{i} = PBPortion.*CellMask{i}.*(~NuclMask{i});
A(i).PBtotal = sum(sum(PB{i}));
A(i).PBsize = sum(sum(PB{i}~=0));
A(i).PBmean = A(i).PBtotal/A(i).PBsize;
NonPB{i} = NonPBPortion.*CellMask{i}.*(~NuclMask{i});
A(i).NonPBtotal = sum(sum(NonPB{i}));
A(i).NonPBsize = sum(sum(NonPB{i}~=0));
A(i).NonPBmean = A(i).NonPBtotal/A(i).NonPBsize;
A(i).PBRatio = A(i).PBtotal/(A(i).PBtotal+A(i).NonPBtotal);

m6ASG{i} = m6ASGPortion.*CellMask{i}.*(~NuclMask{i});
A(i).m6ASGtotal = sum(sum(m6ASG{i}));
A(i).m6ASGsize = sum(sum(m6ASG{i}~=0));
A(i).m6ASGmean = A(i).m6ASGtotal/A(i).m6ASGsize;
m6ANonSG{i} = m6ANonSGPortion.*CellMask{i}.*(~NuclMask{i});
A(i).m6ANonSGtotal = sum(sum(m6ANonSG{i}));
A(i).m6ANonSGsize = sum(sum(m6ANonSG{i}~=0));
A(i).m6ANonSGmean = A(i).m6ANonSGtotal/A(i).m6ANonSGsize;
A(i).m6ASGRatio = A(i).m6ASGtotal/(A(i).m6ASGtotal+A(i).m6ANonSGtotal);

m6ASG0{i} = m6ASGPortion0.*CellMask{i}.*(~NuclMask{i});
A(i).m6ASG0total = sum(sum(m6ASG0{i}));
A(i).m6ASG0size = sum(sum(m6ASG0{i}~=0));
A(i).m6ASG0mean = A(i).m6ASG0total/A(i).m6ASG0size;
m6ANonSG0{i} = m6ANonSGPortion0.*CellMask{i}.*(~NuclMask{i});
A(i).m6ANonSG0total = sum(sum(m6ANonSG0{i}));
A(i).m6ANonSG0size = sum(sum(m6ANonSG0{i}~=0));
A(i).m6ANonSG0mean = A(i).m6ANonSG0total/A(i).m6ANonSG0size;
A(i).m6ASG0Ratio = A(i).m6ASG0total/(A(i).m6ASG0total+A(i).m6ANonSG0total);

m6ASG1{i} = m6ASGPortion1.*CellMask{i}.*(~NuclMask{i});
A(i).m6ASG1total = sum(sum(m6ASG1{i}));
A(i).m6ASG1size = sum(sum(m6ASG1{i}~=0));
A(i).m6ASG1mean = A(i).m6ASG1total/A(i).m6ASG1size;
m6ANonSG1{i} = m6ANonSGPortion1.*CellMask{i}.*(~NuclMask{i});
A(i).m6ANonSG1total = sum(sum(m6ANonSG1{i}));
A(i).m6ANonSG1size = sum(sum(m6ANonSG1{i}~=0));
A(i).m6ANonSG1mean = A(i).m6ANonSG1total/A(i).m6ANonSG1size;
A(i).m6ASG1Ratio = A(i).m6ASG1total/(A(i).m6ASG1total+A(i).m6ANonSG1total);

m6APB{i} = m6APBPortion.*CellMask{i}.*(~NuclMask{i});
A(i).m6APBtotal = sum(sum(m6APB{i}));
A(i).m6APBsize = sum(sum(m6APB{i}~=0));
A(i).m6APBmean = A(i).m6APBtotal/A(i).m6APBsize;
m6ANonPB{i} = m6ANonPBPortion.*CellMask{i}.*(~NuclMask{i});
A(i).m6ANonPBtotal = sum(sum(m6ANonPB{i}));
A(i).m6ANonPBsize = sum(sum(m6ANonPB{i}~=0));
A(i).m6ANonPBmean = A(i).m6ANonPBtotal/A(i).m6ANonPBsize;
A(i).m6APBRatio = A(i).m6APBtotal/(A(i).m6APBtotal+A(i).m6ANonPBtotal);

A(i).PASGEnrich = A(i).PASGmean/A(i).PANonSGmean; %PolyA
A(i).m6ASGEnrich = A(i).m6ASGmean/A(i).m6ANonSGmean; %m6A
A(i).m6AvsASG = A(i).m6ASGEnrich/A(i).PASGEnrich;

A(i).PASG0Enrich = A(i).PASG0mean/A(i).PANonSG0mean; %PolyA
A(i).m6ASG0Enrich = A(i).m6ASG0mean/A(i).m6ANonSG0mean; %m6A
A(i).m6AvsASG0 = A(i).m6ASG0Enrich/A(i).PASG0Enrich;

A(i).PASG1Enrich = A(i).PASG1mean/A(i).PANonSG1mean; %PolyA
A(i).m6ASG1Enrich = A(i).m6ASG1mean/A(i).m6ANonSG1mean; %m6A
A(i).m6AvsASG1 = A(i).m6ASG1Enrich/A(i).PASG1Enrich;

A(i).PBEnrich = A(i).PBmean/A(i).NonPBmean;
A(i).m6APBEnrich = A(i).m6APBmean/A(i).m6ANonPBmean;
A(i).m6AvsAPB = A(i).m6APBEnrich/A(i).PBEnrich;
end
%%Save results

filename = [analysisSavePath PAImageName '_polyA_m6A_SG_PB_Ratio.xlsx'];
InfoTable = struct2table(A);
writetable(InfoTable, filename);
end


   
     

function FuncpolyAm6ASGDensityQuantification_ThresholdEdge_RGB (SerNo, threshold, mainDirPath)
%Use threshold, threshold = 0, and threshold = 5 here. No Nuclear selection

m6AImageName = ['2_647Conv_' SerNo];
%PAImageName = ['1_750Conv_' SerNo]; %PolyA image name 
SGImageName = ['3_561Conv_' SerNo];
PBImageName = ['4_488Conv_' SerNo];

% m6AImageName = ['647Conv_' SerNo];
% PAImageName = ['750Conv_' SerNo]; %PolyA image name 
% SGImageName = ['561Conv_' SerNo];
% PBImageName = ['488Conv_' SerNo];

% Image parameters
pixelSize = 153; % nm
%background750 = 150; %sCMOS background
background647 = 150;
fudgeFactor = 2.5;
size_threshold = 10;
analysisSavePath = [mainDirPath, ['ThresholdEdge_Evaluation_' num2str(threshold) '_' num2str(fudgeFactor) '_' num2str(size_threshold) '\']];
cd(analysisSavePath);

%% Further define image parameters
% parameters for visualizing conventional images
pixelNum = 600;


%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[PBimage, ~] = ReadDaxL([mainDirPath PBImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
%[PAimage, ~] = ReadDaxL([mainDirPath PAImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[m6Aimage, ~] = ReadDaxL([mainDirPath m6AImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);



%% Process data
SGimage(88,405) = SGimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background
BWs = edge(SGimage1,'sobel',threshold * fudgeFactor);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', 2);
BWer1 = imerode(BWnobord,seD);
xbw = bwareaopen(BWer1, size_threshold); %Only select SG larger than 10 pixel after dilation
xbw0 = bwareaopen(BWer1, 2);
xbw1 = bwareaopen(BWer1, 5);

PBbg = imopen(PBimage,strel('disk',10));
PBimage1 = imsubtract(PBimage, PBbg);
PBimage2 = double(PBimage1)/max(max(double(PBimage1))); %Convert to double and 0 to 1 range.
ybw = imbinarize(PBimage2,graythresh(PBimage2)); %Pbody mask

%PAbg = uint16(background750*ones(size(PAimage)));
%PAimage = imsubtract(PAimage,PAbg);
m6Abg = uint16(background647*ones(size(m6Aimage)));
m6Aimage = imsubtract(m6Aimage,m6Abg);

% %% Map
% SGPortion = double(PAimage).*xbw;
% NonSGPortion = double(PAimage).*(~xbw);
% 
% SGPortion0 = double(PAimage).*xbw0;
% NonSGPortion0 = double(PAimage).*(~xbw0);
% 
% SGPortion1 = double(PAimage).*xbw1;
% NonSGPortion1 = double(PAimage).*(~xbw1);
% 
% PBPortion = double(PAimage).*ybw;
% NonPBPortion = double(PAimage).*(~ybw);

m6ASGPortion = double(m6Aimage).*xbw;
m6ANonSGPortion = double(m6Aimage).*(~xbw);

m6ASGPortion0 = double(m6Aimage).*xbw0;
m6ANonSGPortion0 = double(m6Aimage).*(~xbw0);

m6ASGPortion1 = double(m6Aimage).*xbw1;
m6ANonSGPortion1 = double(m6Aimage).*(~xbw1);

m6APBPortion = double(m6Aimage).*ybw;
m6ANonPBPortion = double(m6Aimage).*(~ybw);

%% XY projection image for ROI determination
% Select ROI
g = figure; imshow(SGimage, [0, max(max(SGimage))]); colormap hot;
WhetherROI = questdlg('Do you want to select single cell?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
wait(h);
i = 1;
CellMask{i} = createMask(h);
% WhetherNucl = questdlg('Do you want to select nucleus?'); %ask question
%     if(strcmp(WhetherNucl, 'Yes'))
%     h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
%     wait(h);
%     i = 1;
%     NuclMask{i} = createMask(h);
%     end
moreROI = questdlg('Do you want to select more cells?'); %ask question
%keep looping till the user select 'No'. 
while(strcmp(moreROI, 'Yes'))
    i = i + 1;
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    CellMask{i} = createMask(h);
%     WhetherNucl = questdlg('Do you want to select nucleus?'); %ask question
%     if(strcmp(WhetherNucl, 'Yes'))
%     h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
%     wait(h);
%     NuclMask{i} = createMask(h);
%     end
    moreROI = questdlg('Do you want to select more cells?');
end
% else
%     CellMask{1} = ones(size(xbw)); %If commented, no output. Otherwise
%     use full image to analyze.
%     NuclMask{1} = zeros(size(xbw));
end

% saveas (g, [analysisSavePath PAImageName '_polyA.png']);
% close(g);
% g1 = imshow(SGimage, [0, max(max(SGimage))]);
saveas (g, [analysisSavePath SGImageName '_SGImage.png']);
close(gcf);
g2 = imshow(xbw);
saveas (g2, [analysisSavePath SGImageName '_SGMask.png']);
close(gcf);

%% Quantify PA Intensity Ratio in SG and PB 
for i = 1:length(CellMask)
xbwROI{i} = xbw.*CellMask{i};
xbw0ROI{i} = xbw0.*CellMask{i};
xbw1ROI{i} = xbw1.*CellMask{i};
ybwROI{i} = ybw.*CellMask{i};
%% feature extraction - size distribution (area, pixels)
SG{i} = bwconncomp(xbwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats{i} = regionprops(SG{i});
As{i} = [SGstats{i}.Area];


%statistical measurements
if isempty(As{i})
    A(i).SGmean = 0;
    A(i).SGstd = 0;
    A(i).SGmedian = 0;
else
    A(i).SGmean = mean(As{i});
    A(i).SGstd = std(As{i});
    A(i).SGmedian = median(As{i});
end

A(i).SGtotalarea = sum(As{i});
A(i).SGNo = length(As{i});
A(i).SGtotalIntensity = sum(sum(SGimage1.*xbwROI{i}));
A(i).G3BP1totalIntensity = sum(sum(SGimage1.*(CellMask{i})));
A(i).G3BP1SGRatio = A(i).SGtotalIntensity/A(i).G3BP1totalIntensity;    
A(i).G3BP1mean = A(i).G3BP1totalIntensity/sum(sum(SGimage1.*CellMask{i}>0));

SG0{i} = bwconncomp(xbw0ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats0{i} = regionprops(SG0{i});
As0{i} = [SGstats0{i}.Area];

%statistical measurements
if isempty(As0{i})
    A(i).SG0mean = 0;
    A(i).SG0std = 0;
    A(i).SG0median = 0;
else
    A(i).SG0mean = mean(As0{i});
    A(i).SG0std = std(As0{i});
    A(i).SG0median = median(As0{i});
end
A(i).SG0totalarea = sum(As0{i});
A(i).SG0No = length(As0{i});
A(i).SG0totalIntensity = sum(sum(SGimage1.*xbw0ROI{i}));
A(i).G3BP10SGRatio = A(i).SG0totalIntensity/A(i).G3BP1totalIntensity;   

SG1{i} = bwconncomp(xbw1ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats1{i} = regionprops(SG1{i});
As1{i} = [SGstats1{i}.Area];
%statistical measurements
if isempty(As1{i})
    A(i).SG1mean = 0;
    A(i).SG1std = 0;
    A(i).SG1median = 0;
else
    A(i).SG1mean = mean(As1{i});
    A(i).SG1std = std(As1{i});
    A(i).SG1median = median(As1{i});
end
A(i).SG1totalarea = sum(As1{i});
A(i).SG1No = length(As1{i});
A(i).SG1totalIntensity = sum(sum(SGimage1.*xbw1ROI{i}));
A(i).G3BP11SGRatio = A(i).SG1totalIntensity/A(i).G3BP1totalIntensity;   
% 
% PASG{i} = SGPortion.*CellMask{i}.*(~NuclMask{i});
% A(i).PASGtotal = sum(sum(PASG{i}));
% A(i).SGsize = sum(sum(PASG{i}~=0));
% A(i).PASGmean = A(i).PASGtotal/A(i).SGsize;
% PANonSG{i} = NonSGPortion.*CellMask{i}.*(~NuclMask{i});
% A(i).PANonSGtotal = sum(sum(PANonSG{i}));
% A(i).NonSGsize = sum(sum(PANonSG{i}~=0));
% A(i).PANonSGmean = A(i).PANonSGtotal/A(i).NonSGsize;
% A(i).PASGRatio = A(i).PASGtotal/(A(i).PASGtotal+A(i).PANonSGtotal);
% 
% PASG0{i} = SGPortion0.*CellMask{i}.*(~NuclMask{i});
% A(i).PASG0total = sum(sum(PASG0{i}));
% A(i).SG0size = sum(sum(PASG0{i}~=0));
% A(i).PASG0mean = A(i).PASG0total/A(i).SG0size;
% PANonSG0{i} = NonSGPortion0.*CellMask{i}.*(~NuclMask{i});
% A(i).PANonSG0total = sum(sum(PANonSG0{i}));
% A(i).NonSG0size = sum(sum(PANonSG0{i}~=0));
% A(i).PANonSG0mean = A(i).PANonSG0total/A(i).NonSG0size;
% A(i).PASG0Ratio = A(i).PASG0total/(A(i).PASG0total+A(i).PANonSG0total);
% 
% PASG1{i} = SGPortion1.*CellMask{i}.*(~NuclMask{i});
% A(i).PASG1total = sum(sum(PASG1{i}));
% A(i).SG1size = sum(sum(PASG1{i}~=0));
% A(i).PASG1mean = A(i).PASG1total/A(i).SG1size;
% PANonSG1{i} = NonSGPortion1.*CellMask{i}.*(~NuclMask{i});
% A(i).PANonSG1total = sum(sum(PANonSG1{i}));
% A(i).NonSG1size = sum(sum(PANonSG1{i}~=0));
% A(i).PANonSG1mean = A(i).PANonSG1total/A(i).NonSG1size;
% A(i).PASG1Ratio = A(i).PASG1total/(A(i).PASG1total+A(i).PANonSG1total);

% 
% PB{i} = PBPortion.*CellMask{i}.*(~NuclMask{i});
% A(i).PBtotal = sum(sum(PB{i}));
% A(i).PBsize = sum(sum(PB{i}~=0));
% A(i).PBmean = A(i).PBtotal/A(i).PBsize;
% NonPB{i} = NonPBPortion.*CellMask{i}.*(~NuclMask{i});
% A(i).NonPBtotal = sum(sum(NonPB{i}));
% A(i).NonPBsize = sum(sum(NonPB{i}~=0));
% A(i).NonPBmean = A(i).NonPBtotal/A(i).NonPBsize;
% A(i).PBRatio = A(i).PBtotal/(A(i).PBtotal+A(i).NonPBtotal);

m6ASG{i} = m6ASGPortion.*CellMask{i};
A(i).m6ASGtotal = sum(sum(m6ASG{i}));
A(i).m6ASGsize = sum(sum(m6ASG{i}~=0));
A(i).m6ASGmean = A(i).m6ASGtotal/A(i).m6ASGsize;
m6ANonSG{i} = m6ANonSGPortion.*CellMask{i};
A(i).m6ANonSGtotal = sum(sum(m6ANonSG{i}));
A(i).m6ANonSGsize = sum(sum(m6ANonSG{i}~=0));
A(i).m6ANonSGmean = A(i).m6ANonSGtotal/A(i).m6ANonSGsize;
A(i).m6ASGRatio = A(i).m6ASGtotal/(A(i).m6ASGtotal+A(i).m6ANonSGtotal);

m6ASG0{i} = m6ASGPortion0.*CellMask{i};
A(i).m6ASG0total = sum(sum(m6ASG0{i}));
A(i).m6ASG0size = sum(sum(m6ASG0{i}~=0));
A(i).m6ASG0mean = A(i).m6ASG0total/A(i).m6ASG0size;
m6ANonSG0{i} = m6ANonSGPortion0.*CellMask{i};
A(i).m6ANonSG0total = sum(sum(m6ANonSG0{i}));
A(i).m6ANonSG0size = sum(sum(m6ANonSG0{i}~=0));
A(i).m6ANonSG0mean = A(i).m6ANonSG0total/A(i).m6ANonSG0size;
A(i).m6ASG0Ratio = A(i).m6ASG0total/(A(i).m6ASG0total+A(i).m6ANonSG0total);

m6ASG1{i} = m6ASGPortion1.*CellMask{i};
A(i).m6ASG1total = sum(sum(m6ASG1{i}));
A(i).m6ASG1size = sum(sum(m6ASG1{i}~=0));
A(i).m6ASG1mean = A(i).m6ASG1total/A(i).m6ASG1size;
m6ANonSG1{i} = m6ANonSGPortion1.*CellMask{i};
A(i).m6ANonSG1total = sum(sum(m6ANonSG1{i}));
A(i).m6ANonSG1size = sum(sum(m6ANonSG1{i}~=0));
A(i).m6ANonSG1mean = A(i).m6ANonSG1total/A(i).m6ANonSG1size;
A(i).m6ASG1Ratio = A(i).m6ASG1total/(A(i).m6ASG1total+A(i).m6ANonSG1total);

m6APB{i} = m6APBPortion.*CellMask{i};
A(i).m6APBtotal = sum(sum(m6APB{i}));
A(i).m6APBsize = sum(sum(m6APB{i}~=0));
A(i).m6APBmean = A(i).m6APBtotal/A(i).m6APBsize;
m6ANonPB{i} = m6ANonPBPortion.*CellMask{i};
A(i).m6ANonPBtotal = sum(sum(m6ANonPB{i}));
A(i).m6ANonPBsize = sum(sum(m6ANonPB{i}~=0));
A(i).m6ANonPBmean = A(i).m6ANonPBtotal/A(i).m6ANonPBsize;
A(i).m6APBRatio = A(i).m6APBtotal/(A(i).m6APBtotal+A(i).m6ANonPBtotal);

%A(i).PASGEnrich = A(i).PASGmean/A(i).PANonSGmean; %PolyA
A(i).m6ASGEnrich = A(i).m6ASGmean/A(i).m6ANonSGmean; %m6A
%A(i).m6AvsASG = A(i).m6ASGEnrich/A(i).PASGEnrich;

%A(i).PASG0Enrich = A(i).PASG0mean/A(i).PANonSG0mean; %PolyA
A(i).m6ASG0Enrich = A(i).m6ASG0mean/A(i).m6ANonSG0mean; %m6A
%A(i).m6AvsASG0 = A(i).m6ASG0Enrich/A(i).PASG0Enrich;

%A(i).PASG1Enrich = A(i).PASG1mean/A(i).PANonSG1mean; %PolyA
A(i).m6ASG1Enrich = A(i).m6ASG1mean/A(i).m6ANonSG1mean; %m6A
%A(i).m6AvsASG1 = A(i).m6ASG1Enrich/A(i).PASG1Enrich;

%A(i).PBEnrich = A(i).PBmean/A(i).NonPBmean;
A(i).m6APBEnrich = A(i).m6APBmean/A(i).m6ANonPBmean;
%A(i).m6AvsAPB = A(i).m6APBEnrich/A(i).PBEnrich;
end
%%Save results

filename = [analysisSavePath SGImageName '_SG_PB_Ratio.xlsx'];
InfoTable = struct2table(A);
writetable(InfoTable, filename);
end     
     
     
function FuncSGStatSNAP_ThresholdEdge_NoPB (SerNo, threshold, mainDirPath)
%%This is a function used to quantify SG/PB Number and Size distribution in single cells.

%% Define paths and paramters for images to be analyzed
SGImageName = ['561Conv_' SerNo];
BFPImageName = ['488Conv_' SerNo];
SNAPImageName = ['750Conv_' SerNo];
Y3ImageName = ['647Conv_' SerNo]; %Can be any signals, Y1 staining for siCtr or knockdown.
% SGImageName = ['488Conv_' SerNo];
% Y3ImageName = ['647Conv_' SerNo]; %Y3
% SNAPImageName = ['561Conv_' SerNo];
% BFPImageName = ['750Conv_' SerNo]; %Can be any signals, Y1 staining for siCtr or knockdown.
% SGImageName = ['647Conv_' SerNo]; %Mannual file name 
% BFPImageName = ['488Conv_' SerNo];
% SNAPImageName = ['561Conv_' SerNo];
% Y3ImageName = ['750Conv_' SerNo];

% SGImageName = ['2_647Conv_' SerNo]; %Mannual file name 
% BFPImageName = ['4_488Conv_' SerNo];
% SNAPImageName = ['3_561Conv_' SerNo];
% Y3ImageName = ['1_750Conv_' SerNo];
% Image parameters
BFPbg = 120;
SNAPbg = 120;
%PBthreshold = 1000;
SNAPThreshold = 500; %For select SNAP positive cells
pixelSize = 153; % nm
fudgeFactor = 2.5;
size_threshold = 10;
analysisSavePath = [mainDirPath, ['ThresholdEdge_Evaluation_' num2str(threshold) '_' num2str(fudgeFactor) '_' num2str(size_threshold) '\']];
%% Set up analysis folder

if ~exist(analysisSavePath)
    mkdir(analysisSavePath)
end
cd(analysisSavePath);


%% Further define image parameters
% parameters for visualizing conventional images
pixelNum = 600;

%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[BFPimage, ~] = ReadDaxL([mainDirPath BFPImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[Y3image, ~] = ReadDaxL([mainDirPath Y3ImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[SNAPimage, ~] = ReadDaxL([mainDirPath SNAPImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);


%% Process data
SGimage(88,405) = SGimage(88,404);
BFPimage(88,405) = BFPimage(88,404);
SNAPimage(88,405) = SNAPimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background
BWs = edge(SGimage1,'sobel',threshold * fudgeFactor);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', 2);
BWer1 = imerode(BWnobord,seD);
xbw = bwareaopen(BWer1, size_threshold); %Only select SG larger than 10 pixel after dilation
xbw0 = bwareaopen(BWer1, 2);
xbw1 = bwareaopen(BWer1, 5);

BFPbgm = uint16(BFPbg*ones(size(BFPimage)));
BFPimage1 = imsubtract(BFPimage, BFPbgm);
BFPimage2 = double(BFPimage1);
BFPimage3 = double(BFPimage1)/max(max(double(BFPimage1)));

SNAPbgm = uint16(SNAPbg*ones(size(SNAPimage)));
SNAPimage1 = imsubtract(SNAPimage, SNAPbgm);
SNAPimage2 = double(SNAPimage1);
SNAPimage3 = double(SNAPimage1)/max(max(double(SNAPimage1)));
% 
% PBbg = imopen(PBimage,strel('disk',10));
% PBimage1 = double(imsubtract(PBimage, PBbg));
% %PBimage2 = double(PBimage1)/max(max(double(PBimage1))); %Convert to double and 0 to 1 range.
% ybw = imbinarize(PBimage1,PBthreshold); %Pbody mask

if max(max(SNAPimage)) > SNAPThreshold

%% XY projection image for ROI determination
% Select ROI
g = figure; imshow(SNAPimage, [0, max(max(SNAPimage))]); colormap hot;
WhetherROI = questdlg('Do you want to select single cell?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
wait(h);
i = 1;
CellMask{i} = createMask(h);
moreROI = questdlg('Do you want to select more cells?'); %ask question
%keep looping till the user select 'No'. 
while(strcmp(moreROI, 'Yes'))
    i = i + 1;
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    CellMask{i} = createMask(h);
    moreROI = questdlg('Do you want to select more cells?');
end


saveas (g, [analysisSavePath SGImageName '_Seg.png']);
close(g);
g1 = imshow(SGimage, [0, max(max(SGimage))]);
saveas (g1, [analysisSavePath SGImageName '_SGImage.png']);
g2 = imshow(xbw);
saveas (g2, [analysisSavePath SGImageName '_SGMask.png']);
close(gcf);
g3 = imshow(BFPimage3);
saveas (g3, [analysisSavePath SGImageName '_BFP.png']);
close(gcf);
g5 = imshow(SNAPimage3);
saveas (g5, [analysisSavePath SGImageName '_SNAP.png']);
close(gcf);

else
   close(g); 
end
for i = 1:length(CellMask)
xbwROI{i} = xbw.*CellMask{i};
xbw0ROI{i} = xbw0.*CellMask{i};
xbw1ROI{i} = xbw1.*CellMask{i};
% ybwROI0{i} = logical(ybw.*CellMask{i});
%% feature extraction - size distribution (area, pixels)
SG{i} = bwconncomp(xbwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats{i} = regionprops(SG{i});
A{i} = [SGstats{i}.Area];
%statistical measurements
if isempty(A{i})
    a(i).SGmean = 0;
    a(i).SGstd = 0;
    a(i).SGmedian = 0;
else
    a(i).SGmean = mean(A{i});
    a(i).SGstd = std(A{i});
    a(i).SGmedian = median(A{i});
end
a(i).SGtotalarea = sum(A{i});
%a(i).SGNo = length(A{i});
a(i).SGNo = length(A{i});
a(i).SGtotalIntensity = sum(sum(SGimage1.*xbwROI{i}));
a(i).G3BP1totalIntensity = sum(sum(SGimage1.*CellMask{i}));
a(i).G3BP1SGRatio = a(i).SGtotalIntensity/a(i).G3BP1totalIntensity;


SG0{i} = bwconncomp(xbw0ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats0{i} = regionprops(SG0{i});
As0{i} = [SGstats0{i}.Area];

%statistical measurements
if isempty(As0{i})
    a(i).SG0mean = 0;
    a(i).SG0std = 0;
    a(i).SG0median = 0;
else
    a(i).SG0mean = mean(As0{i});
    a(i).SG0std = std(As0{i});
    a(i).SG0median = median(As0{i});
end
a(i).SG0totalarea = sum(As0{i});
a(i).SG0No = length(As0{i});
a(i).SG0totalIntensity = sum(sum(SGimage1.*xbw0ROI{i}));
a(i).G3BP10SGRatio = a(i).SG0totalIntensity/a(i).G3BP1totalIntensity;   

SG1{i} = bwconncomp(xbw1ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats1{i} = regionprops(SG1{i});
As1{i} = [SGstats1{i}.Area];
%statistical measurements
if isempty(As1{i})
    a(i).SG1mean = 0;
    a(i).SG1std = 0;
    a(i).SG1median = 0;
else
    a(i).SG1mean = mean(As1{i});
    a(i).SG1std = std(As1{i});
    a(i).SG1median = median(As1{i});
end
a(i).SG1totalarea = sum(As1{i});
a(i).SG1No = length(As1{i});
a(i).SG1totalIntensity = sum(sum(SGimage1.*xbw1ROI{i}));
a(i).G3BP11SGRatio = a(i).SG1totalIntensity/a(i).G3BP1totalIntensity;  


% 
% ybwROI{i} = bwareafilt(ybwROI0{i},[4 300]); %Filter SG pixel number
% PB{i} = bwconncomp(ybwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
% PBstats{i} = regionprops(PB{i});
% B{i} = [PBstats{i}.Area];
% %statistical measurements
% a(i).PBmean = mean(B{i});
% a(i).PBstd = std(B{i});
% a(i).PBmedian = median(B{i});
% a(i).PBtotalarea = sum(B{i});
% a(i).PBNo = length(B{i});
% a(i).PBtotalIntensity = sum(sum(PBimage1.*ybwROI{i}));
% a(i).DCP1totalIntensity = sum(sum(PBimage1.*CellMask{i}));
% a(i).DCP1PBRatio = a(i).PBtotalIntensity/a(i).DCP1totalIntensity;
Y3{i} = double(Y3image).*CellMask{i};
a(i).Y3total = sum(sum(Y3{i}));
a(i).Y3size = sum(sum(Y3{i}~=0));
a(i).Y3mean = a(i).Y3total/a(i).Y3size;

BFP{i} = BFPimage2.*CellMask{i};
a(i).Y1total = sum(sum(BFP{i}));
a(i).Y1size = sum(sum(BFP{i}~=0));
a(i).Y1mean = a(i).Y1total/a(i).Y1size;

SNAP{i} = SNAPimage2.*CellMask{i};
a(i).SNAPtotal = sum(sum(SNAP{i}));
a(i).SNAPsize = sum(sum(SNAP{i}~=0));
a(i).SNAPmean = a(i).SNAPtotal/a(i).SNAPsize;
end
%%Save results

filename = [analysisSavePath SGImageName '_SG_PB_BFP_SNAP_Info.xlsx'];
InfoTable = struct2table(a);
writetable(InfoTable, filename);

else
end
end  


     
function FuncSGStatCry2O_ThresholdEdge (SerNo, threshold, mainDirPath)
%%This is a function used to quantify SG/PB Number and Size distribution in single cells.

%% Define paths and paramters for images to be analyzed
% SGImageName = ['561Conv_' SerNo];
% BFPImageName = ['488Conv_' SerNo];
% SNAPImageName = ['750Conv_' SerNo];
% Y3ImageName = ['647Conv_' SerNo]; %Can be any signals, Y1 staining for siCtr or knockdown.
% SGImageName = ['488Conv_' SerNo];
% Y3ImageName = ['647Conv_' SerNo]; %Y3
% SNAPImageName = ['561Conv_' SerNo];
% BFPImageName = ['750Conv_' SerNo]; %Can be any signals, Y1 staining for siCtr or knockdown.
% SGImageName = ['647Conv_' SerNo]; %Mannual file name 
% BFPImageName = ['488Conv_' SerNo];
% SNAPImageName = ['561Conv_' SerNo];
% Y3ImageName = ['750Conv_' SerNo];

SGImageName = ['2_647Conv_' SerNo]; %Mannual file name 
%BFPImageName = ['4_488Conv_' SerNo]; %Null, space filler
SNAPImageName = ['3_561Conv_' SerNo]; %Actually BFP image
%Y3ImageName = ['1_750Conv_' SerNo]; %Null, space filler
% Image parameters
BFPbg = 120;
SNAPbg = 120;
%PBthreshold = 1000;
SNAPThreshold = 15000; %For select SNAP positive cells
pixelSize = 153; % nm
fudgeFactor = 2.5;
size_threshold = 10;
analysisSavePath = [mainDirPath, ['ThresholdEdge_Evaluation_' num2str(threshold) '_' num2str(fudgeFactor) '_' num2str(size_threshold) '\']];
%% Set up analysis folder

if ~exist(analysisSavePath)
    mkdir(analysisSavePath)
end
cd(analysisSavePath);


%% Further define image parameters
% parameters for visualizing conventional images
pixelNum = 600;

%% Load data
[SGimage, ~] = ReadDaxL([mainDirPath SGImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
%[BFPimage, ~] = ReadDaxL([mainDirPath BFPImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
%[Y3image, ~] = ReadDaxL([mainDirPath Y3ImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);
[SNAPimage, ~] = ReadDaxL([mainDirPath SNAPImageName '.dax'], 'startFrame', 2, 'endFrame', 2, 'verbose', false);


%% Process data
SGimage(88,405) = SGimage(88,404);
% BFPimage(88,405) = BFPimage(88,404);
SNAPimage(88,405) = SNAPimage(88,404);
%background estimation (non uniform illumination)
SGbg = imopen(SGimage,strel('disk',10));
%background removal (flatten background level)
SGimage1 = double(imsubtract(SGimage,SGbg));
%segment grains from background
BWs = edge(SGimage1,'sobel',threshold * fudgeFactor);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
BWsdil = imdilate(BWs,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', 2);
BWer1 = imerode(BWnobord,seD);
xbw = bwareaopen(BWer1, size_threshold); %Only select SG larger than 10 pixel after dilation
xbw0 = bwareaopen(BWer1, 2);
xbw1 = bwareaopen(BWer1, 5);

% BFPbgm = uint16(BFPbg*ones(size(BFPimage)));
% BFPimage1 = imsubtract(BFPimage, BFPbgm);
% BFPimage2 = double(BFPimage1);
% BFPimage3 = double(BFPimage1)/max(max(double(BFPimage1)));

SNAPbgm = uint16(SNAPbg*ones(size(SNAPimage)));
SNAPimage1 = imsubtract(SNAPimage, SNAPbgm);
SNAPimage2 = double(SNAPimage1);
SNAPimage3 = double(SNAPimage1)/max(max(double(SNAPimage1)));
% 
% PBbg = imopen(PBimage,strel('disk',10));
% PBimage1 = double(imsubtract(PBimage, PBbg));
% %PBimage2 = double(PBimage1)/max(max(double(PBimage1))); %Convert to double and 0 to 1 range.
% ybw = imbinarize(PBimage1,PBthreshold); %Pbody mask

if max(max(SNAPimage)) > SNAPThreshold

%% XY projection image for ROI determination
% Select ROI
g = figure; imshow(SNAPimage, [0, max(max(SNAPimage))]); colormap hot;
WhetherROI = questdlg('Do you want to select single cell?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
wait(h);
i = 1;
CellMask{i} = createMask(h);
moreROI = questdlg('Do you want to select more cells?'); %ask question
%keep looping till the user select 'No'. 
while(strcmp(moreROI, 'Yes'))
    i = i + 1;
    h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
    wait(h);
    CellMask{i} = createMask(h);
    moreROI = questdlg('Do you want to select more cells?');
end


saveas (g, [analysisSavePath SGImageName '_Seg.png']);
close(g);
g1 = imshow(SGimage, [0, max(max(SGimage))]);
saveas (g1, [analysisSavePath SGImageName '_SGImage.png']);
g2 = imshow(xbw);
saveas (g2, [analysisSavePath SGImageName '_SGMask.png']);
close(gcf);
% g3 = imshow(BFPimage3);
% saveas (g3, [analysisSavePath SGImageName '_BFP.png']);
% close(gcf);
g5 = imshow(SNAPimage3);
saveas (g5, [analysisSavePath SGImageName '_BFP.png']);
close(gcf);

else
   close(g); 
end
for i = 1:length(CellMask)
xbwROI{i} = xbw.*CellMask{i};
xbw0ROI{i} = xbw0.*CellMask{i};
xbw1ROI{i} = xbw1.*CellMask{i};
% ybwROI0{i} = logical(ybw.*CellMask{i});
%% feature extraction - size distribution (area, pixels)
SG{i} = bwconncomp(xbwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats{i} = regionprops(SG{i});
A{i} = [SGstats{i}.Area];
%statistical measurements
if isempty(A{i})
    a(i).SGmean = 0;
    a(i).SGstd = 0;
    a(i).SGmedian = 0;
else
    a(i).SGmean = mean(A{i});
    a(i).SGstd = std(A{i});
    a(i).SGmedian = median(A{i});
end
a(i).SGtotalarea = sum(A{i});
%a(i).SGNo = length(A{i});
a(i).SGNo = length(A{i});
a(i).SGtotalIntensity = sum(sum(SGimage1.*xbwROI{i}));
a(i).G3BP1totalIntensity = sum(sum(SGimage1.*CellMask{i}));
a(i).G3BP1SGRatio = a(i).SGtotalIntensity/a(i).G3BP1totalIntensity;


SG0{i} = bwconncomp(xbw0ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats0{i} = regionprops(SG0{i});
As0{i} = [SGstats0{i}.Area];

%statistical measurements
if isempty(As0{i})
    a(i).SG0mean = 0;
    a(i).SG0std = 0;
    a(i).SG0median = 0;
else
    a(i).SG0mean = mean(As0{i});
    a(i).SG0std = std(As0{i});
    a(i).SG0median = median(As0{i});
end
a(i).SG0totalarea = sum(As0{i});
a(i).SG0No = length(As0{i});
a(i).SG0totalIntensity = sum(sum(SGimage1.*xbw0ROI{i}));
a(i).G3BP10SGRatio = a(i).SG0totalIntensity/a(i).G3BP1totalIntensity;   

SG1{i} = bwconncomp(xbw1ROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
SGstats1{i} = regionprops(SG1{i});
As1{i} = [SGstats1{i}.Area];
%statistical measurements
if isempty(As1{i})
    a(i).SG1mean = 0;
    a(i).SG1std = 0;
    a(i).SG1median = 0;
else
    a(i).SG1mean = mean(As1{i});
    a(i).SG1std = std(As1{i});
    a(i).SG1median = median(As1{i});
end
a(i).SG1totalarea = sum(As1{i});
a(i).SG1No = length(As1{i});
a(i).SG1totalIntensity = sum(sum(SGimage1.*xbw1ROI{i}));
a(i).G3BP11SGRatio = a(i).SG1totalIntensity/a(i).G3BP1totalIntensity;  


% 
% ybwROI{i} = bwareafilt(ybwROI0{i},[4 300]); %Filter SG pixel number
% PB{i} = bwconncomp(ybwROI{i}); %SG{i} = bwconncomp(xbwROI{i}) is more memory efficient
% PBstats{i} = regionprops(PB{i});
% B{i} = [PBstats{i}.Area];
% %statistical measurements
% a(i).PBmean = mean(B{i});
% a(i).PBstd = std(B{i});
% a(i).PBmedian = median(B{i});
% a(i).PBtotalarea = sum(B{i});
% a(i).PBNo = length(B{i});
% a(i).PBtotalIntensity = sum(sum(PBimage1.*ybwROI{i}));
% a(i).DCP1totalIntensity = sum(sum(PBimage1.*CellMask{i}));
% a(i).DCP1PBRatio = a(i).PBtotalIntensity/a(i).DCP1totalIntensity;

% Y3{i} = double(Y3image).*CellMask{i};
% a(i).Y3total = sum(sum(Y3{i}));
% a(i).Y3size = sum(sum(Y3{i}~=0));
% a(i).Y3mean = a(i).Y3total/a(i).Y3size;
% 
% BFP{i} = BFPimage2.*CellMask{i};
% a(i).BFPtotal = sum(sum(BFP{i}));
% a(i).BFPsize = sum(sum(BFP{i}~=0));
% a(i).BFPmean = a(i).Y1total/a(i).BFPsize;

SNAP{i} = SNAPimage2.*CellMask{i};
a(i).BFPtotal = sum(sum(SNAP{i}));
a(i).BFPsize = sum(sum(SNAP{i}~=0));
a(i).BFPmean = a(i).BFPtotal/a(i).BFPsize;
end
%%Save results

filename = [analysisSavePath SGImageName '_SG_BFP_Info.xlsx'];
InfoTable = struct2table(a);
writetable(InfoTable, filename);

else
end
end 