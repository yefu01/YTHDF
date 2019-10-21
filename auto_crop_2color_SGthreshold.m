% how to use it:
% 1 Open STORM-L.dax file, and load two-color STORM-L.txt molecule list.
% 2 Append STORM-R.txt molecule list, save the combined molecule list as .bin file in the same directory as .dax file.
% 3 Load 2d-680_647_View_Crop_G3BP680.ini configuration file to Insight3.
% 4 Save images as _SG.png file from Insight3 in green color in the same directory as the .dax file
% 5 Run this program
% 
% spectial note: one .bin file MUST correspond to one *_SG.png. the
% program are ONLY designed when png and bin files are matched in the same
% order

% Ye Fu, modified from Boran Han's code


clear all;
close all;
addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');

auto_crop_2color_SG('D:\STORM1_Data\20181003_U2OS_G3BP1-647-680\A3_NaAsO2_G3BP1-RbmAb-647_MsmAb-680\split\');
function auto_crop_2color_SG(binpath)
mkdir([binpath 'SG_regions']);
ipath=[binpath 'SG_regions\'];
mkdir([ipath 'Large']);
mkdir([ipath 'Medium']);
filename_bin = dir(fullfile(strcat([binpath '\'],'*647*_Split-L_H1.bin')));
for ii=1:length(filename_bin)
    
B = regexp(filename_bin(ii).name,'\d*','Match');
for i= 1:length(B)
  if ~isempty(B{i})
      Num(i,1)=str2double(B{i});
  else
      Num(i,1)=NaN;
  end
end

image_file=['STORM_647_',B{2}, '_SG.png'];
binname = filename_bin(ii).name;
if ~exist([binpath image_file], 'file')
   continue;
end

% data = importdata([binpath binname]);
[mList, memoryMap] = ReadMasterMoleculeList([binpath '\' binname]);
Mlist=[];
for i=1:length(mList.x)
    %,'fieldsToLoad',{'xc','yc','z'},'ZScale',167);
    Mlist(i).x=mList.x(i);
    Mlist(i).y=mList.y(i);
    Mlist(i).xc=mList.xc(i);
    Mlist(i).yc=mList.yc(i);
    Mlist(i).h=mList.h(i);
    Mlist(i).a=mList.a(i);
    Mlist(i).w=mList.w(i);
    Mlist(i).phi=mList.phi(i);
    Mlist(i).ax=mList.ax(i);
    Mlist(i).bg=mList.bg(i);
    Mlist(i).i=mList.i(i);
    Mlist(i).c=mList.c(i);
    Mlist(i).density=mList.density(i);
    Mlist(i).frame=mList.frame(i);
    Mlist(i).length=mList.length(i);
    Mlist(i).link=mList.link(i);
    Mlist(i).z=mList.z(i);
    Mlist(i).zc=mList.zc(i);
end
% data = importdata([binpath binname(1:end-8) 'drift.txt']);
%   x_drift = data(:,2); %pixel -> nm
%   y_drift = data(:,3);
%   z_drift = data(:,4);
% % [x_drift,y_drift] = XcorrDriftCorrect(mList);
% for i=1:length(Mlist)
% Mlist(i).xc = Mlist(i).x - x_drift(Mlist(i).frame);
% Mlist(i).yc = Mlist(i).y - y_drift(Mlist(i).frame);
% Mlist(i).zc = Mlist(i).z - z_drift(Mlist(i).frame);
% end

im=double(imread([binpath image_file]));
im_length = length(im);
im = im(:,:,2)./max(max(im(:,:,2)));
%im = im./max(max(im));%for not green 
%im = rgb2gray(im); %for not green
im = imgaussfilt(im,1.2);

%xbw = imbinarize(im,2*graythresh(im));
xbw = imbinarize(im,graythresh(im)/1.2); %2* for crop using G3BP1-680
xbw(1:floor(16*im_length/256), :) = 0;
xbw(floor(240*im_length/256):im_length, :) = 0;
xbw(:, 1:floor(16*im_length/256)) = 0;
xbw(:, floor(240*im_length/256):im_length) = 0;
% 
% [crop_index_row,crop_index_col]=find(im(:,:,1)==255&im(:,:,2)==0);
% for i=1:length(crop_index_row)
% crop_image(crop_index_row(i),crop_index_col(i))=1;
% end
% crop_image = imfill(crop_image,'holes');
crop_image=bwareaopen(xbw, floor(2^2*(im_length/256)^2)); % select cluster larger than 2 pixel*2 pixel (167 nm diameter) 4*4 for m6A
%Discard data at edge to eliminate Edge effect

g2 = imshow(crop_image);
saveas (g2, [ipath binname(1:end-4) '_SGMark.png']);
close(gcf);
% figure;
% imshow(crop_image);
% imbinarize(crop_image)



crop_image_label=bwlabel(crop_image);
% 
% imshow(crop_image_label);
data_crop=[];
%% 
for j=1:length(Mlist)
Mlist(j).boxnumber=0;
end


for i=1:max(max(crop_image_label))
    [crop_index_row,crop_index_col]=find(crop_image_label==i);
    index=find([Mlist.xc]*im_length/256>min(crop_index_col)&[Mlist.xc]*im_length/256<max(crop_index_col)&[Mlist.yc]*im_length/256<max(crop_index_row)&[Mlist.yc]*im_length/256>min(crop_index_row));
    for j=index
        if ismember([round(Mlist(j).xc*im_length/256),round(Mlist(j).yc*im_length/256)],[crop_index_col,crop_index_row],'rows')
          Mlist(j).boxnumber=i;
%         else Mlist(j).boxnumber=0;
        end
    end
end
%% 

for i=1:max(max(crop_image_label))
    data_crop=Mlist(find([Mlist.boxnumber]==i));
%     figure; scatter([data_crop.x]*im_length/256,[data_crop.y]*im_length/256,'.')
    [crop_index_row,crop_index_col]=find(crop_image_label==i);
    crop_center_x=mean(crop_index_row);
    crop_center_y=mean(crop_index_col);
%     distance=pdist2([crop_center_x,crop_center_y],[crop_index_row,crop_index_col]);
%     [sorteddis,sortingindices]=sort(distance, 'descend');
%     corner_x=crop_index_row(sortingindices(1:10));
%     corner_y=crop_index_col(sortingindices(1:10));

    corner_x=max(crop_index_row);
    corner_y=min(crop_index_col(find(crop_index_row==corner_x)));    
    theta_1=atan((corner_y-crop_center_y)./(corner_x-crop_center_x))/pi*180;
    if (corner_y-crop_center_y)<0 & (corner_x-crop_center_x)<0
        theta_1=theta_1+180;
    end
%     angle=unique(round(theta/20));
%     if length(angle)<2
%         continue
%     end
%     theta_1=mean(theta(find(round(theta/20)==angle(1))));
%     theta_2=mean(theta(find(round(theta/20)==angle(2))));
    corner_y=max(crop_index_col);
    corner_x=max(crop_index_row(find(crop_index_col==corner_y))); 
    theta_2=atan((corner_y-crop_center_y)./(corner_x-crop_center_x))/pi*180;
    
    if abs(theta_1-theta_2)>=90
        rotate_angle=-(theta_1+theta_2)/2;
    else rotate_angle=90-(theta_1+theta_2)/2;
    end
    rotate_angle=-deg2rad(rotate_angle);
    R=[cos(rotate_angle),-sin(rotate_angle);sin(rotate_angle), cos(rotate_angle)];
    coordinate_new=R*[[data_crop.xc];[data_crop.yc]];
    for jj=1:length(coordinate_new)
    data_crop(jj).xc=coordinate_new(1,jj)-min(coordinate_new(1,:))+5;
    data_crop(jj).yc=coordinate_new(2,jj)-min(coordinate_new(2,:))+5;
    end
    SGarea_ind = sum(sum(crop_image_label==i));
    %if SGarea_ind < 2^2*(im_length/256)^2
        %WriteMoleculeList(data_crop, [ipath 'Small\' binname(1:end-4) '_' num2str(i) '.bin']) %Discard less than 4 pixels
    if (SGarea_ind  >= 2^2*(im_length/256)^2) && (SGarea_ind < 4^2*(im_length/256)^2) %Medium SG  16 pixels
        WriteMoleculeList(data_crop, [ipath 'Medium\' binname(1:end-4) '_' num2str(i) '.bin'])
    else 
        WriteMoleculeList(data_crop, [ipath 'Large\' binname(1:end-4) '_' num2str(i) '.bin'])
    end
%     clear data_crop
    
end
end
end