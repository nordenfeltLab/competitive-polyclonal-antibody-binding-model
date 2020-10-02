function [abDist] = preparePolyclonal(noAntibodies,AntibodyLength, BacProteinLength,FcBindingRegion,MeanValue,Range,KdFc, k)
rng(k);
min = MeanValue - Range./2;
max = MeanValue + Range./2;
abDist = struct();

x = linspace(min,max,noAntibodies);

if noAntibodies == 1
    abDist.y = 1;
else
    abDist.y = normpdf(x,MeanValue,Range/4*0.95).*(Range./(noAntibodies-1)); %sum should be conc, hence scaling to range/noAnts
end


%Dissociation constants expressed as probability; input Kd is dissociation constant, Kd in transfer matrix is equillibrium constant.
%% Creating normal distribution of affinity magnitude
abDist.KdFab = (1./(1+10.^x)).^(1/2);
%KdFab =(1./(1+x)).^(1/2); %square root because there are two fab sites per antibody.
abDist.KdFc = 1./(1+KdFc);

abDist.FabBindingRegion = randomBindingSites(AntibodyLength, BacProteinLength, noAntibodies);
end

function bindingSites = randomBindingSites(AntibodyLength, BacProteinLength, noAntibodies)
% give random binding sites/epitopes
bindingSites = ceil(randi([ceil(AntibodyLength./2),BacProteinLength-ceil(AntibodyLength./2)], 1, noAntibodies));
end
