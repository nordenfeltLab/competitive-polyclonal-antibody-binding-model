%% Variables
clear all
n = 50; %number of conc points
%conc = 1000;
conc = logspace(1,8,n) ; % concentration antibodies
%conc = [ 0 0.0001 0.3 1 3 10 30 100 300 1000 3000 10000];
%conc = [0 0.00671 0.02013 0.06710 0.2013 0.67100 2.013 6.71000 6.71000 13.42000 26.84000 53.68000 67.10000];
%n = length(conc);
noAntibodies = 30;
AntibodyLength = 7;
BacProteinLength = 33; 
Range = 11.154; %Fab affinity distribution MAGNITUDE log10
MeanValue = 2.5; %Fab affinity MeanValue MAGNITUDE log10

%Fc properties for  bacterial protein
noFcSites= 1;
%KdFc = [0.8068]; %in micrograms/ml
%KdFc = [0.51]; %in micrograms/ml MPROTEIN
KdFc = [0.08];
FcBindingRegion = [ 25 ]; %nonlinear protein, regions set arbitrarily, test differences?
k = 0;
[TotalBinding, FabBinding, FcBinding] = setCalcConc(n, conc, noAntibodies, AntibodyLength, BacProteinLength,...
    MeanValue, Range, noFcSites, KdFc, FcBindingRegion, k);

