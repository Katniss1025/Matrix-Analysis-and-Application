function [Q,R] = QR(A)
% Using Gram-Schmidt method to find Q,R such that A = QR
% where R is upper triangular matrix
% A:m*n, Q:m*n,R:n*n
% created by Xiaodan Yin,11/18/2021
    [m,n] = size(A);
    Q = zeros(size(A));
    R = zeros(n,n);
    R(1,1) = sqrt(sum(A(:,1).^2)); %r11 = ||a1||
    Q(:,1) = A(:,1)/R(1,1); % q1 = a1/r11;
    for i = 2:n
        for j = 1:i-1
            R(j,i) = Q(:,j)' * A(:,i);
            Q(:,i) = Q(:,i) + R(j,i)*Q(:,j);
        end
        Q(:,i) = A(:,i) - Q(:,i);
        R(i,i) = sqrt(sum(Q(:,i).^2));
        Q(:,i) = Q(:,i)/R(i,i);
    end
end