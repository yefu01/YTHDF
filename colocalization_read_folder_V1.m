clear all;
close all;
a1=20;
a2=20;
b=2;
ipath = 'D:\STORM1_Data\20181003_U2OS_G3BP1-647-680\A3_NaAsO2_G3BP1-RbmAb-647_MsmAb-680\split\SG_regions\Medium\';
addpath('N:\MATLAB_Analysis\matlab-storm\Functions\Plotting and Display');
addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');

filename = dir(fullfile(strcat([ipath '\'],'*.bin')));
for ii=1:length(filename)
    [mList, memoryMap] = ReadMasterMoleculeList([ipath '\' filename(ii).name],'fieldsToLoad',{'xc','yc','z', 'c'},'ZScale',167);
    if sum(mList.c==2) && sum(mList.c==3) ~= 0 
    ROI=[ min(mList.yc (mList.c==2)) max(mList.yc(mList.c==2));min(mList.xc(mList.c==2)) max(mList.xc(mList.c==2))];
    renderedStack_1 = RenderMList([mList.xc(mList.c==2) mList.yc(mList.c==2)],  ...
    'gaussianWidth', 0.1, ...
    'ROI', ROI, ...
    'imageScale', a1);
%     a2=a1*length(mList.xc(mList.c==3))/length(mList.xc(mList.c==2));
    %renderedStack_2 = RenderMList([mList.xc(mList.c==1) mList.yc(mList.c==1)],  ...
    %renderedStack_2 = RenderMList([mList.xc(mList.c==3) mList.yc(mList.c==3)],  ...
    renderedStack_2 = RenderMList([mList.xc(mList.c==3|mList.c==1) mList.yc(mList.c==3|mList.c==1)],  ...
    'gaussianWidth', 0.1, ...
    'ROI', ROI, ...
     'imageScale', a2);
%     renderedStack_2=renderedStack_2/max(renderedStack_2(:))*max(renderedStack_1(:));
%     renderedStack_1=renderedStack_1/max(renderedStack_1(:));
    renderedStack_2=imresize(renderedStack_2,a1/a2);
    min_length=min(size(renderedStack_1,1),size(renderedStack_2,1));
    min_width=min(size(renderedStack_1,2),size(renderedStack_2,2));
    
    backgroundIm_1 = imopen(renderedStack_1, strel('disk', 10));
    renderedStack_1 = renderedStack_1 - backgroundIm_1;
    backgroundIm_2 = imopen(renderedStack_2, strel('disk', 10));
    renderedStack_2 = renderedStack_2 - backgroundIm_2;
    XC_2D(ii) = corr2 (renderedStack_1,renderedStack_2);
    imbw_1=imbinarize(renderedStack_1(1:min_length,1:min_width),graythresh(renderedStack_1)*b);%,'adaptive','Sensitivity',0.4);
    imbw_2=imbinarize(renderedStack_2(1:min_length,1:min_width),graythresh(renderedStack_2)*b);%,'adaptive','Sensitivity',0.4);
    index=find(imbw_1==1|imbw_2==1);
    XC(ii) = corr2 (renderedStack_1(index),renderedStack_2(index));
    XCBW(ii) = corr2 (imbw_1(index),imbw_2(index));
%         XC(ii) = corr2 (imbw_1,imbw_2);
%     figure(1);imshow(imbw_1);
%     figure(2);imshow(imbw_2);
    ratio_1(ii)=length(find(imbw_1==1&imbw_2==1))/length(find(imbw_1==1));
    ratio_2(ii)=length(find(imbw_1==1&imbw_2==1))/length(find(imbw_2==1));

    renderedStack=cat(3,renderedStack_1(1:min_length,1:min_width),renderedStack_2(1:min_length,1:min_width));
    renderedStack=cat(3,renderedStack,zeros(size(renderedStack,1),size(renderedStack,2)));
%     figure(3);imshow(renderedStack, []);
%     imwrite(coloredSTORM,[ipath '\' filename(ii).name(1:end-4) '.png'])
    else
        XC_2D(ii) = NaN;
        XC(ii) = NaN;
        XCBW(ii) = NaN;
        ratio_1(ii) = NaN;
        ratio_2(ii) = NaN;
    end
end




%   'type1 XC'
%   aamplitude_random=XC;%(randperm(length(aamplitude)));
%   replicate=[nanmean(aamplitude_random(1:round(length(aamplitude_random)/3))),nanmean(aamplitude_random(round(length(aamplitude_random)/3)+1:round(2*length(aamplitude_random)/3))),nanmean(aamplitude_random(round(2*length(aamplitude_random)/3)+1:end))];
%   ste=std(replicate)/sqrt(3)
%   average=nanmean(replicate)
%   replicate
  samplesize = sum(~isnan(XC));
  XClength = length(XC);
  XC_2D(XClength+1) = nanmean (XC_2D(1:XClength))
  XC_2D(XClength+2) = nanstd(XC_2D(1:XClength))./sqrt(samplesize)  
  XC(XClength+1) = nanmean (XC(1:XClength))
  XC(XClength+2) = nanstd(XC(1:XClength))./sqrt(samplesize)
  XCBW(XClength+1) = nanmean (XCBW(1:XClength))
  XCBW(XClength+2) = nanstd(XCBW(1:XClength))./sqrt(samplesize)
  ratio_1(XClength+1) = nanmean (ratio_1(1:XClength))
  ratio_1(XClength+2) = nanstd(ratio_1(1:XClength))./sqrt(samplesize)
  ratio_2(XClength+1) = nanmean (ratio_2(1:XClength))
  ratio_2(XClength+2) = nanstd(ratio_2(1:XClength))./sqrt(samplesize)
  Table = [XC_2D',XC',XCBW', ratio_1',ratio_2'];
  xlswrite([ipath 'XC-c13_a' num2str(a1) '_b' num2str(b) '.xls'], Table)