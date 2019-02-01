function [J, grad] = linearRegCostFunction(X, y, theta, lambda)
%LINEARREGCOSTFUNCTION Compute cost and gradient for regularized linear 
%regression with multiple variables
%   [J, grad] = LINEARREGCOSTFUNCTION(X, y, theta, lambda) computes the 
%   cost of using theta as the parameter for linear regression to fit the 
%   data points in X and y. Returns the cost in J and the gradient in grad

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
n=size(theta);
grad = zeros(n);

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost and gradient of regularized linear 
%               regression for a particular choice of theta.
%
%               You should set J to the cost and grad to the gradient.
%

reg=0;
for i=1:m
  p=0;
  for j=1:n
    p=p+theta(j)*X(i,j);
  endfor
  p=p-y(i,1);
  
  for j=1:n
    grad(j,1)=grad(j,1)+p* X(i,j);
  endfor
  p=p*p;
  J=J+p;
endfor

grad(1,1)=grad(1,1)/m;
 for j=2:n
   reg=reg+theta(j)*theta(j);
   grad(j,1)=grad(j,1)+lambda*theta(j);
   grad(j,1)=grad(j,1)/m;
 endfor
 reg=reg*lambda*0.5/m;
 
 J=J*0.5/m;
 
 J=J+reg;











% =========================================================================

grad = grad(:);

end
