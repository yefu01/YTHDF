% clear all
% close all

% calibrate in 3D the two channels.
% generate the warping matrix from Cy5 channel to the Cy3 channel.

function function_split_two_channel(dataPath, FileName)
output_folder='split\';
[Datafiles, InfoFile] = ReadDax([dataPath FileName],'startFrame', 1, 'endFrame', 1000);
FrameNumber = InfoFile.number_of_frames;
for j = 0: floor(FrameNumber/1000)
[Datafiles, InfoFile] = ReadDax([dataPath FileName],'startFrame', j*1000+1, 'endFrame', (j+1)*1000);
Datafiles=double(Datafiles);

for i = 1:size(Datafiles, 3)
    LeftDatafiles(:,:,i) = Datafiles(1:256,1:256,i);
    RightDatafiles(:,:,i) = Datafiles(1:256,257:512,i);
end
 %figure;imagesc(LeftDatafiles(:,:,1))
    LeftDatafiles=uint16(LeftDatafiles);
    RightDatafiles=uint16(RightDatafiles);
if j>0
    WriteDAXFiles_append(LeftDatafiles, newInfFileL);
    WriteDAXFiles_append(RightDatafiles, newInfFileR);
else
    InfoFile.vstart=1;
    InfoFile.hstart=1;
    InfoFile.frame_dimensions=[256 256];
    InfoFile.frame_size = InfoFile.frame_size/2;
    newInfFileL=InfoFile;
    newInfFileR=InfoFile;
    newInfFileL.localName=[InfoFile.localName(1:end-4) '-L.inf'];
    filteredDataPath = [dataPath output_folder];
    newInfFileL.localPath = filteredDataPath;
    WriteDAXFiles(LeftDatafiles, newInfFileL);
    newInfFileR.localName=[InfoFile.localName(1:end-4) '-R.inf'];
    newInfFileR.localPath = filteredDataPath;
    WriteDAXFiles(RightDatafiles, newInfFileR);
end
end
