% %%
function [BindingProbability,expectationValue,iBindingProbability] = bindingCalc(AbConc, noAntibodies, AntibodyLength, BacProteinLength, FcBindingRegion, FabBindingRegion, KdFc, KdFab, k)

noFcSites = length(FcBindingRegion);

nopStates = 3.*noAntibodies.*AntibodyLength + 1; %Number of possible states
TransferMatrix = cell(1, BacProteinLength); % an array with all transfer matrices

%% Free Concentration 
% solution to calculating free concentration of antibodies. 
% Holds only for lambda = 1

%BacProteinConc = 1; %Constant
%K = KdFab.*prod(KdFc);
%for n = 1:length(TotalConc)
%BoundConc = (sqrt((BacProteinConc-TotalConc(n)).^2.*K.^2 + 2.*K.*(BacProteinConc+TotalConc(n)+1))+K.*(BacProteinConc+TotalConc(n))+1)./2.*K;
%end

%BoundConc = (-sqrt((BacProteinConc-TotalConc).^2.*K.^2 + 2.*K.*(BacProteinConc+TotalConc+1))+K.*(BacProteinConc+TotalConc)+1)./2.*K
%BoundConc = - (sqrt((BacProteinConc-TotalConc
%FreeConc = TotalConc-BoundConc;
%FreeConc = AbConc;
%% Concentration

%% Fab Binding 

bindingProbability = [repmat(KdFab.*AbConc,2,1);ones(1,noAntibodies).*KdFc.*(AbConc)];
%bindingProbability = [repmat(KdFab.*AbConc,2,1);zeros(1,noAntibodies)];
%bindingProbability(3,5) = KdFc.*sum(AbConc);
bindingProbability = bindingProbability(:);

bindingRegion = [repmat(FabBindingRegion,2,1);ones(1,noAntibodies).*FcBindingRegion];
bindingRegion = bindingRegion(:);

diagonalVector = ones(1,nopStates-1);
diagonalVector((floor(AntibodyLength./2))+1:AntibodyLength:nopStates-1) = 0;
prob = sparse(nopStates,nopStates);
prob([1 2:AntibodyLength:nopStates],1:AntibodyLength:nopStates) = 1;

prob(2:nopStates+1:end) = diagonalVector;

for n = 1:BacProteinLength
    %roi = computeMatrix(n,AntibodyLength,nopStates,bindingRegion(:));
    TransferMatrix{n} = prob;
    %TransferMatrix{n}(roi,:) = prob(roi,:);
end
for i = 1:length(bindingRegion)
    index = ceil(AntibodyLength./2)+(i-1).*AntibodyLength;
    %assert(2+index <= nopStates)
    TransferMatrix{bindingRegion(i)}(1+index,index) = bindingProbability(i);
end
%% Fc Binding
% FOR MULTIPLE FC-Sites
% 
% %noFcSites = 2;
% for k = 1:noFcSites
%     m= FcBindingRegion(k)-floor(AntibodyLength./2):FcBindingRegion(k)+floor(AntibodyLength./2);
%     %t = m(m>1 & m<BacProteinLength);
%     for x = 1:AntibodyLength;
%         if m(x) > 0 && m(x) < BacProteinLength+1
%             for i=ceil(2.*AntibodyLength):(3.*AntibodyLength):nopStates-(x+1)
%              TransferMatrix{m(x)}(i+x+1,i+x) = 1;
%              TransferMatrix{m(1)}(i+2,1:AntibodyLength:nopStates) = 1;
%             end
%         end
%     end
%     
%     %TransferMatrix{m(1)}([1 2:AntibodyLength:nopStates],1:AntibodyLength:nopStates) = 1
%     locations = floor(2.5.*AntibodyLength):(3.*AntibodyLength):nopStates-2;
%     for l = locations    
%         value = KdFc(k).*sum(AbConc);
%         TransferMatrix{FcBindingRegion(k)}(l+3,l+2) = value;
%         %m= FcBindingRegion(k)-floor(AntibodyLength./2):FcBindingRegion(k)+floor(AntibodyLength./2);
%     end
% end 

%% Calculating P(i)_s for fab and fc and entire partition function
%Partition function Z(i)

FirstColumnVector =  zeros(nopStates,1)';
FirstColumnVector([1:AntibodyLength:nopStates]) = 1;

LastColumnVector = zeros(nopStates,1);
LastColumnVector([1 2:AntibodyLength:nopStates]) = 1;


Bound{1} = FirstColumnVector;
Bound{2} = LastColumnVector;

MatrixCumProduct = cell(2,BacProteinLength);
[MatrixCumProduct, PartitionFunction] = MultiplyMatrices(FirstColumnVector,LastColumnVector,MatrixCumProduct, nopStates,BacProteinLength,TransferMatrix);

%% Projectionoperators for fab and fc boltzmann weigthed possible states
 
%ProjectionOperatorFab = diag([ 0 repmat([ones(1,2.*AntibodyLength) zeros(1,AntibodyLength)],1,noAntibodies)]);
ProjectionOperatorFab = sparse(1:nopStates,1:nopStates,[ 0 repmat([ones(1,2.*AntibodyLength) zeros(1,AntibodyLength)],1,noAntibodies)]);
%ProjectionOperatorFab = diag(ones(1,nopStates));
%ProjectionOperatorFab = sparse(ProjectionOperatorFab);
BoltzmannWeigthedZFab = ComputeBoltzmann(MatrixCumProduct,ProjectionOperatorFab,BacProteinLength,nopStates);
%BoltzmannWeigthedZFab = CompBoltz(TransferMatrix,ProjectionOperatorFab,BacProteinLength,Bound);
ProbabilityFabBinding = BoltzmannWeigthedZFab./PartitionFunction;
 
%ProjectionOperatorFc = diag([ 0 repmat([zeros(1,2.*AntibodyLength) ones(1,AntibodyLength)],1,noAntibodies)]);
ProjectionOperatorFc = sparse(1:nopStates,1:nopStates,[ 0 repmat([zeros(1,2.*AntibodyLength) ones(1,AntibodyLength)],1,noAntibodies)]);
%ProjectionOperatorFc = diag(ones(1,nopStates));
%ProjectionOperatorFc = sparse(ProjectionOperatorFc);
BoltzmannWeigthedZFc = ComputeBoltzmann(MatrixCumProduct,ProjectionOperatorFc,BacProteinLength,nopStates);
%BoltzmannWeigthedZFc = CompBoltz(TransferMatrix,ProjectionOperatorFc,BacProteinLength,Bound);
ProbabilityFcBinding = BoltzmannWeigthedZFc./PartitionFunction;

% %% Identifying individual antibodies
% iProjectionOperatorFc = cell(1,noAntibodies);
% iProjectionOperatorFab = cell(1,noAntibodies);
% iBoltzmannWeigthedZFab = zeros(noAntibodies,BacProteinLength);
% iBoltzmannWeigthedZFc = zeros(noAntibodies,BacProteinLength);
% 
% 
% for a = 1:noAntibodies
%     fab = zeros(1, nopStates);
%     fc = zeros(1, nopStates);
%     fab((2+3.*AntibodyLength.*(a-1)):(1 + 3.*AntibodyLength.*(a-1)+2.*AntibodyLength)) = 1;
%     fc((2+2.*AntibodyLength+3.*AntibodyLength.*(a-1)):(1 + 3.*AntibodyLength.*(a-1)+3.*AntibodyLength)) = 1;
%     iProjectionOperatorFab{1,a} = sparse(diag(fab));
%     iProjectionOperatorFc{1,a} = sparse(diag(fc));
%     iBoltzmannWeigthedZFab(a,:) = ComputeBoltzmann(MatrixCumProduct,iProjectionOperatorFab{1,a},BacProteinLength,nopStates);
%     iBoltzmannWeigthedZFc(a,:) = ComputeBoltzmann(MatrixCumProduct,iProjectionOperatorFc{1,a},BacProteinLength,nopStates);
% end
% iBoltzmannWeigthedZFab = sum(iBoltzmannWeigthedZFab');
% iBoltzmannWeigthedZFc = sum(iBoltzmannWeigthedZFc');
% iProbabilityFcBinding = iBoltzmannWeigthedZFc./PartitionFunction;
% iProbabilityFabBinding = iBoltzmannWeigthedZFab./PartitionFunction;
% iBindingProbability =[ iProbabilityFcBinding; iProbabilityFabBinding ];
%%
iBindingProbability = zeros(2, noAntibodies);

 %% Projectionoperator for each binding weight and expectationvalue
 sumN = 0;
 for k =1:noAntibodies.*3;
 PO = sparse(nopStates,nopStates);
 Element =ceil(AntibodyLength./2)+1+AntibodyLength.*(k-1);
 PO(ceil(AntibodyLength./2)+1+AntibodyLength.*(k-1),ceil(AntibodyLength./2)+1+AntibodyLength.*(k-1)) = 1;
 %PO = sparse(PO);
 sumN = sumN + sum(ComputeBoltzmann(MatrixCumProduct,PO,BacProteinLength,nopStates));
 end
 expectationValue = (sumN./PartitionFunction)/BacProteinLength;
 BindingProbability = [ProbabilityFabBinding;ProbabilityFcBinding];
end
