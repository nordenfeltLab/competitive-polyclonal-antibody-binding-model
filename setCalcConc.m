%Calculate and plot binding curves: Total, Fab and Fc.
function [TotalBinding FabBinding FcBinding ] = setCalcConc(n, conc, noAntibodies, AntibodyLength, BacProteinLength,...
    MeanValue, Range, noFcSites, KdFc, FcBindingRegion, k)

    abDist = preparePolyclonal(noAntibodies,AntibodyLength, BacProteinLength,FcBindingRegion,MeanValue,Range,KdFc, k);

%% Calculating Binding probability for the different concentrations
A = zeros(length(conc).*2,BacProteinLength);
for c =1:n
    [A([c.*2-1 c.*2],:)] = equilibriumSearch(n, conc(c), noAntibodies, AntibodyLength, BacProteinLength, abDist, noFcSites, KdFc, FcBindingRegion, @calcMultipleIgG, k);
    %fprintf('count = %d\n',count)
end

FabBinding = sum(A(1:2:end,:),2)./BacProteinLength;
FcBinding = sum(A(2:2:end,:),2)./BacProteinLength;
TotalBinding = FcBinding+FabBinding;
%[TotalBinding FabBinding FcBinding ]
% % Plotting Binding Probability Curve
% 
% hold on
% x = 10.^x;
%  
% %plot(x, AbConc)
% set(gca,'xscale','log')
% stem(x,iBindingProbability(1,:)./(BacProteinLength));
% set(gca,'xscale','log')
% stem(x,iBindingProbability(2,:)./(BacProteinLength));
% 
% %plot(x, AbConc,x,iBindingProbability(2,:)./BacProteinLength.*AbConc, x,iBindingProbability(1,:)./BacProteinLength.*AbConc)
% set(gca,'xscale','log')
% legend('Fc','Fab')
% %legend('Distribution','Fc', 'Fab')
% xlabel('KdFab')
% ylabel('Probability of binding')
plotResults(TotalBinding, FabBinding, FcBinding, conc)

end
