%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
% this code can split channel, fit molecules, and then assign the molecules
% to two channels

addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');
addpath('M:\Matlab\Dual_View_Boran');
dataPath = ''; %Datapath here
analysisPath = dataPath;
pixelSize = 167; 
dataFiles = dir(fullfile(strcat([dataPath],'*storm*.dax')));
%Main folder is dataPath, with data in subfolders
if ~isempty(dataFiles)

dataFiles = dir(fullfile(strcat(dataPath,'*storm*.dax')));

for ii=1:length(dataFiles)
    try
    FuncMasterRunFile_insight_V1(dataPath, 'split\', dataFiles(ii).name(1:end-4),'split\');
    end
end


else
    subFolders = dir(fullfile(strcat([dataPath],'*_*'))); 
for i = 1:length(subFolders) 
    dataPathSub= [dataPath,subFolders(i).name,'\'];
    dataFiles = dir(fullfile(strcat(dataPathSub,'*storm*.dax'))); %Main folder is dataPath, with data in subfolders


for ii=1:length(dataFiles)
    try
    FuncMasterRunFile_insight_V1(dataPathSub, 'split\', dataFiles(ii).name(1:end-4),'split\');
    end
end
 
end
end

