%Calculate and plot binding curves: Total, Fab and Fc.
function [BindingProbability, expectationValue, iBindingProbability, x,y, AbConc] = equilibriumSearch(n, c, noAntibodies, AntibodyLength, BacProteinLength,...
    abDist, noFcSites, KdFc, FcBindingRegion, bindingCalculationFun, k)
    
    concFree(2) = c;
    count = 0;
    finished = 0;
    while ~finished
        concFree(1) = concFree(2);
        [BindingProbability, expectationValue, iBindingProbability, x,y, AbConc] = bindingCalculationFun(concFree(2),noAntibodies,AntibodyLength,BacProteinLength,FcBindingRegion, abDist, k);
        
        %concBound= (sum(A(c.*2-1)) + sum(A(c.*2)));
        %concFree(2) = conc(c) - concBound; %P(s)_i
        concFree(2) = c - expectationValue; 
        
        count = count+1;
        finished = ~(abs(concFree(2)-concFree(1))./concFree(1) > 0.0001);
    end

end
