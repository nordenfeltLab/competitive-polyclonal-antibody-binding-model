function [BindingProbability,expectationValue,iBindingProbability, x,y, AbConc] = calcMultipleIgG(conc, noAntibodies, AntibodyLength, BacProteinLength,FcBindingRegion,abDist,k)
AbConc = abDist.y .* conc;

y = abDist.y;
x = log10(1./(abDist.KdFab.^2) - 1);
[BindingProbability,expectationValue,iBindingProbability] = bindingCalc(AbConc, noAntibodies, AntibodyLength, BacProteinLength, FcBindingRegion, abDist.FabBindingRegion, abDist.KdFc, abDist.KdFab, k);
end
