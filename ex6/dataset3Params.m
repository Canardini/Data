function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 0.01;
sigma = 0.01;

bestIndexi=1;
bestIndexj=1;
minerror=1e20;
for i=1:3
  Cp=C*10^(i-1);
  for j=1:3  
    sigmap=sigma*10^(j-1);  
    model= svmTrain(X, y,Cp, @(x1, x2) gaussianKernel(x1, x2, sigmap));  
    predictions=svmPredict(model,Xval);
    err=mean(double(predictions~=yval));
    if(err<minerror)
      bestIndexi=i;
      bestIndexj=j;
      minerror=err;
    end
  endfor
endfor


C=C*10^(bestIndexi-1);
sigma=sigma*10^(bestIndexj-1); 


% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%







% =========================================================================

end
