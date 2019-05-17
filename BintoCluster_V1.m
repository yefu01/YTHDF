clear all;
close all;
a1=20;
b=2;
mainDirPath = 'P:\STORM1_Data\'; %Datapath for .bin file
addpath('N:\MATLAB_Analysis\matlab-storm\Functions\Plotting and Display');
addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');

filename = dir(fullfile(strcat([mainDirPath '\'],'*.bin')));
for ii=1:length(filename)
    
    mList = ReadMasterMoleculeList([mainDirPath filename(ii).name],'ZScale',167);
    data = importdata([mainDirPath filename(ii).name(1:end-8) 'drift.txt']);
    x_drift = data(:,2); %pixel
    y_drift = data(:,3);
    z_drift = data(:,4);
    for i=1:length(mList.frame)
        mList.xc(i) = mList.x(i) - x_drift(mList.frame(i));
        mList.yc(i) = mList.y(i) - y_drift(mList.frame(i));
        mList.zc(i) = mList.z(i) - z_drift(mList.frame(i));
    end
    ROI=[ min(mList.yc (mList.c==1)) max(mList.yc(mList.c==1));min(mList.xc(mList.c==1)) max(mList.xc(mList.c==1))];

    renderedStack_1 = RenderMList([mList.xc(mList.c==1) mList.yc(mList.c==1)],  ...
    'gaussianWidth', 0.1, ...
    'ROI', ROI, ...
    'imageScale', a1);
    min_length=size(renderedStack_1,1);
    min_width=size(renderedStack_1,2);
    
    %Background subtraction (optional) 20190121 Does not change much
    %though.
    backgroundIm_1 = imopen(renderedStack_1, strel('disk', 10));
    renderedStack_1 = renderedStack_1 - backgroundIm_1;
    
    fighandle1 = imshow(renderedStack_1);
    saveas(fighandle1, [mainDirPath filename(ii).name(1:end-8) 'RenderImage.png']);
%     renderedStack_1Binary = double(renderedStack_1)/max(max(double(renderedStack_1)));
%     xbw = imbinarize(renderedStack_1Binary,graythresh(renderedStack_1Binary));
    %xbw = imbinarize(renderedStack_1,2*graythresh(renderedStack_1));
    xbw = imbinarize(renderedStack_1,graythresh(renderedStack_1));
    figure();
    fighandle2 =  imshow(xbw);
    saveas(fighandle2, [mainDirPath filename(ii).name(1:end-8) 'BWImage.png']);
    RenderLabeled_1 = bwlabel(xbw);
    Renderstats_1 = regionprops(RenderLabeled_1);
    for i = 1:length(Renderstats_1)
    Renderstats_1(i).Area_in_nm = Renderstats_1(i).Area*(167/a1)^2; %Area in a unit of (pixel/a1)^2
    end
    Cluster_Info = struct2table(Renderstats_1);
    xlsfilename = [mainDirPath filename(ii).name(1:end-8) 'RenderCluster.xls'];
    writetable(Cluster_Info, xlsfilename);
    close all;
end