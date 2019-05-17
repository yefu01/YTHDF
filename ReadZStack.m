function [ImageStack, InfoFile] = ReadZStack(FileName, NumImage, StepInterval)

[MovieFP, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', NumImage);
for i = 1:floor((NumImage)/StepInterval)-1
    ImageStack(:,:,i) = mean(MovieFP(:,:, (i-1)*StepInterval+2:i*StepInterval+1),3);
end
