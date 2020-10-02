function [ MatrixCumProduct, PartitionFunction ] = MultiplyMatrices(FirstColumnVector,LastColumnVector, MatrixCumProduct, nopStates,BacProteinLength,TransferMatrix)
%Performs matrix multiplication "forward" and "backwards", saves results in
%MatrixCumProduct, calculates partition function.

 for n = 1:BacProteinLength
     MatrixCumProduct{1,n}= FirstColumnVector;
     FirstColumnVector = FirstColumnVector*TransferMatrix{n};
 end

 PartitionFunction = FirstColumnVector*LastColumnVector;
 
 for n = BacProteinLength:-1:1
     LastColumnVector = TransferMatrix{n}*LastColumnVector;
     MatrixCumProduct{2,n}= LastColumnVector;
 end
end

