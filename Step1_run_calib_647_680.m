%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
% this code can split channel, fit molecules, and then assign the molecules
% to two channels

addpath('C:\Users\Ye\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Ye\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Ye\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Ye\matlab-storm\Functions\IO\');
addpath('C:\Users\Ye\matlab-storm\Functions\Misc\');
%c:\MinGW64\mingw\bin\;c:\MinGW64\mingw\lib\;c:\Users\Hazen\lib\;
global daoSTORMexe;
daoSTORMexe='path=C:\Python27\;C:\Users\Hazen\storm-analysis\windows_dll\; && set PYTHONPATH=%PYTHONPATH%;C:\Users\Hazen\storm-analysis\;  && python.exe C:\Users\Hazen\storm-analysis\3d_daostorm\mufit_analysis.py';
dataPath = ''; %dataPath for Beads image here
analysisPath = dataPath;
pixelSize = 167; 
threshold=50;

dataFiles = dir(fullfile(strcat([dataPath],'*Conv*.dax')));

for ii=1:length(dataFiles)
    function_split_two_channel(dataPath, dataFiles(ii).name)
end

dataFiles = dir(fullfile(strcat([dataPath  'split\'],'*Conv*.dax')));

for ii=1:length(dataFiles)
    analysisLabel = '';
    desiredFrame = 1;
    
        parametersName = ['daoSTORMParameters_' analysisLabel '.xml'];
        
        WriteDaoSTORMParameters([analysisPath parametersName], ...
        'iterations', 20, ... % Do not fit overlapping molecules
        'threshold', threshold, ...
        'sigma', 1.6, 'frame_step',500,... % The best guess for the PSF size in terms of pixels
        'pixel_size', pixelSize, 'drift_correction', true);

%% Run analysis of all fov

        daoSTORM([dataPath 'split\' dataFiles(ii).name], ... % Path to dax to analyze
        [analysisPath parametersName], ... % Path to daoSTORM configuration file
        'overwrite', true, ... % Overwrite any existing analysis
        'numParallel', 1, ... % Only process one file at a time
        'savePath', [analysisPath 'split\'], ... % Location of the mlist files
        'mListType', [analysisLabel 'mlist'], ... % A label to mark the analysis for each file
        'outputInMatlab', true); % Display the output of daoSTORM in matlab
  
%         mListFiles = dir(fullfile([dataPath dataFiles(ii).name(1:end-4) '_mlist.bin']));
       
  
end

for ii= 10%1:length(dataFiles)
    run_calib_647_680([dataPath, 'split\'], dataFiles(ii).name(1:end-6));
end

