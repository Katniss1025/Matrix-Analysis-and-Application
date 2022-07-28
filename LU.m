function [L,U] = LU(A)
% Aiming to find L and U such that A = LU,
% where L is lower triangular matrix
%       U is upper tirangular matrix
% created by Xiaodan Yin, 11/18/2021
L = eye(size(A));
U = eye(size(A));
[rows,~] = size(A);
 for i = 1: rows-1
     L(i+1:rows,i) = A(i+1:rows,i) ./ A(i,i);
     A(i+1:rows,:) = A(i+1:rows,:) - A(i,:) .* A(i+1:rows,i) ./ A(i,i);
 end
 U = A;
end