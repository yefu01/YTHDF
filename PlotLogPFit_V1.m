%%Fit cluster list with cluster area to -Log(P) = aR^3+bR^2+c;

close all;
clear all;
dataPath = 'P:\';
datafiles = dir(fullfile(strcat(dataPath, '*RenderCluster*.xls')));


Area = {};
Localization = {};
for i = 1:length(datafiles)
filename = [dataPath datafiles(i).name];
Data = xlsread(filename);
Area{i} = Data(:,8);
end

Total.Area = [];
Total.Valid = [];
for i = 1:length(Area)
    Total.Area = [Total.Area; Area{i}];
end
Total.R = sqrt(Total.Area./pi);
Total.V = 4/3*pi.*(Total.R).^3; 
for i = 1:length(Total.Area)
    Total.Valid(i) = (Total.R(i) > 50);
end
Total.Valid = Total.Valid';
%Threshold R = 50 nm, A = 7854 nm^2, V = 523600
% xlswrite([dataPath 'CluserSum.xls'],[Total.Localization, Total.R, Total.Area, Total.Valid]);
A = Total.V(Total.Valid==1); %If size is too big, write valid R only
Radius = (A*3/(4*pi)).^(1/3);
mkdir([dataPath 'Rendered_Analysis\Start80R50_1E5_320\']);
dataPath = [dataPath 'Rendered_Analysis\Start80R50_1E5_320\'];
%xlswrite([dataPath 'RenderedCluserSum.xls'],A); %Can't write too large dataset
save([dataPath 'RenderedCluserSum.mat'], 'A','Radius');
%% Fit (-logP)-V
binsize = 100000;
edges = [4/3*pi*50^3:binsize:4/3*pi*350^3]; 
[N,edges] = histcounts(Total.V(Total.Valid==1),edges);
TotalNo=sum(N);
P = N./TotalNo;
NegLogP = -log(P);
A = [((edges(2:end)-binsize/2).*3/(4.*pi)).^(1/3)',NegLogP'];
ANonZero = A(isfinite(A(:,2)),:);

StartBin = 1;
EndBinStart = length(ANonZero);


for i = 1:EndBinStart
     %for polyfit end
     EndBin2 = EndBinStart - i;
    [fitresult, gof] = poly320Fit(ANonZero(StartBin:EndBin2,1), ANonZero(StartBin:EndBin2,2)); %, ANonZero(StartBin:EndBin2,4));
    Rc = -2*fitresult.p2/(3*fitresult.p1);
    if ANonZero(EndBin2,1) < (Rc * 0.8)
        break
    else
        close gcf;
    end
    end
ANonZero(:,3) = ANonZero(:,2) - fitresult.p4;
x = [10:2:Rc];
y = fitresult.p1.*x.^3+fitresult.p2.*x.^2+fitresult.p4;
fitcurve = [x',y'];

hold on;
gof.Rc = Rc;
p = [fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4];
ci = confint(fitresult);
SD = ci-p; %Standard Deviation
SD = SD(2,:)./1.96;
SDRc = (fitresult.p2*SD(1)-fitresult.p1*SD(2))/fitresult.p1^2;
gof.SDRc = SDRc;
X = A(:,1);
Y = A(:,2);
h = scatter(X, Y);
xlabel('Radius/nm');
ylabel('-log(P)');
txt = {['a=',num2str(fitresult.p2), '+-',num2str(SD(2))],['b=',num2str(fitresult.p1), '+-',num2str(SD(1))],['c=',num2str(fitresult.p4), '+-',num2str(SD(4))],['Rc=',num2str(gof.Rc),'+-',num2str(gof.SDRc)],['r^2=',num2str(gof.rsquare)]};
text(90,4,txt);

mkdir(dataPath);
savefig([dataPath 'CluserSumFit.fig']);
saveas(gcf, [dataPath 'CluserSumFit.png']);

save([dataPath 'fitresults.mat'],'fitresult','gof','X','Y', 'ANonZero','fitcurve');