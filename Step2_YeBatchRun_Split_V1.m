%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
% this code can split channel, fit molecules, and then assign the molecules
% to two channels

addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');
addpath('M:\Matlab\Dual_View_Boran');

global daoSTORMexe;
daoSTORMexe='path=C:\Python27\;C:\Users\Hazen\storm-analysis\windows_dll\; && set PYTHONPATH=%PYTHONPATH%;C:\Users\Hazen\storm-analysis\;  && python.exe C:\Users\Hazen\storm-analysis\3d_daostorm\mufit_analysis.py';
dataPath = 'P:\STORM1_Data\20181003_U2OS_G3BP1-647-680\';
analysisPath = dataPath;
pixelSize = 167; 
threshold=200;
dataFiles = dir(fullfile(strcat([dataPath],'*storm*.dax')));
%Main folder is dataPath, with data in subfolders
if ~isempty(dataFiles)

for ii=1:length(dataFiles)
    function_split_two_channel(dataPath, dataFiles(ii).name)
end

else
    subFolders = dir(fullfile(strcat([dataPath],'*_*'))); %Added by Ye 20180404, to process subfolders. Folder name contains '_'
for i = 1:length(subFolders) %modified to 2
    dataPathSub= [dataPath,subFolders(i).name,'\'];
    dataFiles = dir(fullfile(strcat([dataPathSub],'*storm*.dax'))); %Main folder is dataPath, with data in subfolders



for ii=1:length(dataFiles)
    function_split_two_channel(dataPathSub, dataFiles(ii).name)
end

end
end

