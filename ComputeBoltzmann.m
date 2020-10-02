function BoltzmannWeightedZ = ComputeBoltzmann( MatrixCumProduct, ProjectionOperator, BacProteinLength, nopStates )
%Computes allowed states using MatrixCumproduct and the projection
%operators for Fab and Fc states.
 BoltzmannWeightedZ = zeros(1,BacProteinLength);
 for n =1:BacProteinLength
    BoltzmannWeightedZ(n) = MatrixCumProduct{1,n}*ProjectionOperator*MatrixCumProduct{2,n};
 end
end
