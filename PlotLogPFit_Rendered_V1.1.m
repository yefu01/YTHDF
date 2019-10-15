%%For Classical Nucleation Model data fitting
%Ye Fu
%Harvard University
%yefu01@fas.harvard.edu

close all;
clear all;
dataPath = 'D:\STORM1_Data\20190108_U2OS_siY1Y3_G3BP1-647\B1_siCtr_None\'; %Folder name here
datafiles = dir(fullfile(strcat(dataPath, '*RenderCluster*.xls'))); %Cluster file name

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
    Total.Valid(i) = (Total.R(i) > 50); %Lower Threshold R = 50 nm, A = 7854 nm^2, V = 523600
end
Total.Valid = Total.Valid';
% xlswrite([dataPath 'CluserSum.xls'],[Total.Localization, Total.R,
% Total.Area, Total.Valid]); %optional
A = Total.V(Total.Valid==1); %If size is too big, write valid R only
Radius = (A*3/(4*pi)).^(1/3);
mkdir([dataPath 'Rendered_Analysis\Start80R50_1E5\']);
dataPath = [dataPath 'Rendered_Analysis\Start80R50_1E5\'];
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
SeqInf = find(isinf(A(:,2)),10,'first');
EndBin = SeqInf(10);
% SeqInf = find(isinf(A(:,2)),5,'first'); %To discard after 5th empty bin
% EndBin = SeqInf(5);
for i = 2:(length(SeqInf)-1)
    if (SeqInf(i)-SeqInf(i-1))==1 && (SeqInf(i+1)-SeqInf(i)==1)
        EndBin = SeqInf(i-1);
    end
end
%EndBin = find(isinf(A(:,2)),1,'first'); %This line only for finding 1st inf and discard
AShort = A(1:EndBin,:);
ANonZero = AShort(isfinite(AShort(:,2)),:);
EndBinStart = length(ANonZero);
StartBin = 1;

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
%plot(fitresult, ANonZero(StartBin:EndBin2,1), ANonZero(StartBin:EndBin2,2), 'residuals')
fitx = ANonZero (StartBin:EndBin2,1);
fity = fitresult.p1.*fitx.^3+fitresult.p2.*fitx.^2+fitresult.p4;
chi2 = chi_squared(ANonZero(StartBin:EndBin2,2),fity,3); %Pearson chi-squared goodness of fit test
p_value_chi2 = chi2cdf(chi2,gof.dfe);
gof.p = p_value_chi2;
fitresidual = [fitx, (ANonZero(StartBin:EndBin2,2) - fity)];
ANonZero(:,3) = ANonZero(:,2) - fitresult.p4;
x = [10:2:Rc];
y = fitresult.p1.*x.^3+fitresult.p2.*x.^2+fitresult.p4;
fitcurve = [x',y'];
hold on;
p = [fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4];
ci = confint(fitresult);
CI = ci-p; %95 confidence interval
SD = CI(2,:)./1.96; %Standard error. %norminv((1-0.95)/2) = 1.96 Normal distribution
tStat_at_zero = abs(p./SD);
pvalue = 2*tcdf(tStat_at_zero, gof.dfe, 'upper'); %P value from zero
X = A(:,1);
Y = A(:,2);
h = scatter(X, Y, 'g.');
h = scatter(ANonZero(:,1), ANonZero(:,2), 'r.');
h = scatter(ANonZero(StartBin:EndBin2,1), ANonZero(StartBin:EndBin2,2), 'k.');
set(gca,'xlim',[50, 320]);
plot(fitresult, 'predobs');
grid off;
legend(gca,'off');
xlabel('Radius/nm');
ylabel('-log(P)');
txt = {['a = ',num2str(fitresult.p2,'%.2s'), ', CI = ',num2str(CI(2),'%.1s'), ', p = ', num2str(pvalue(2),'%.2s')],['b = ',num2str(fitresult.p1,'%.2s'), ', CI = ',num2str(CI(1),'%.1s'), ', p = ', num2str(pvalue(1),'%.2s')],['c = ',num2str(fitresult.p4,'%.2s'), ', CI = ',num2str(CI(4),'%.1s'), ', p = ', num2str(pvalue(4),'%.2s')],['Rc = ',num2str(Rc, '%.1f')],['Adjusted R^2 = ',num2str(gof.adjrsquare,'%.4f')], ['SSE = ',num2str(gof.sse,'%.1f')], ['pchi2 = ', num2str(gof.p, '%.2s')]};
text(110,4,txt,'FontSize', 12);
set(gca, 'FontSize', 12);
set(gca,'xlim',[50, 375], 'ylim', [0, 15]); %13 for G3BP1
mkdir(dataPath);
savefig([dataPath 'CluserSumFit.fig']);
saveas(gcf, [dataPath 'CluserSumFit.png']);

figure()
scatter(fitresidual(:,1),fitresidual(:,2), 'k.');
hline = refline([0 0]);
hline.Color = 'r';
set(gca,'ylim',[-2.5, 2.5]);
saveas(gcf, [dataPath 'Fitresiduals.png']);
saveas(gcf, [dataPath 'Fitresiduals.fig']);
save([dataPath 'fitresults.mat'],'fitresult','gof','X','Y', 'ANonZero','fitcurve');


function chi2 = chi_squared(y,fit,P,eb)
% returns *reduced* chi^2 value for use in data modelling
% "y" is a vector of data, "fit" is a vector of model values (size(fit)=size(y)), P is the number of
% parameters fit in the model, and eb is a vector of error bars (1-to-1 correspondnce with y)
% Ref: John R. Taylor, "An Introduction to Error Analysis", (2nd ed., 1997)
% 11/11/01 Mike Scarpulla.  Please direct questions or comments to scarps@uclink.berkeley.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    error('Wrong number of arguments passed to "chi_squared"')
end
% if error bars are not availible, evaluate chi^2 by normalizing deviation^2 by magnitude of data.
% This assumes that the STDEV of a value scales as SQRT(value).  USE WITH THIS CAVEAT IN MIND
if nargin==3
    N = max(size(y));
    terms = ((y-fit).^2)./abs(y);
    chi2 = 1/(N-P)*sum(terms);
end
%if error bars are availible, normalize the deviation to the expectred error
if nargin==4
    N = max(size(y));
    terms = ((y-fit)./eb).^2;
    chi2 = 1/(N-P)*sum(terms);
end
end