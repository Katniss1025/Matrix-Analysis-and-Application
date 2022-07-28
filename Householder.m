function [Q,R] = Householder(A)
    [m,n] = size(A);
    k = min(m,n);
    R = zeros(m,m,k);
    for i = 1:k
        R(:,:,i) = eye(m,m);
    end
    A_1 = A;
    for i = 1:m-1
        x = A_1(i:end,i);
        e1 = zeros(size(x));
        e1(1) = 1;
        if isreal(x(1))
            mu = 1;
        else
            mu = x(1)/norm(x(1));
        end
        u = x-mu*norm(x)*e1;
        [d,~] = size(u);
        R(i:end,i:end,i) = eye(d) - 2*u*u'/(u'*u);
        A_1 = R(:,:,i)*A_1;
    end
    Q = 1;
    for i = k:-1:1
        Q = Q * R(:,:,i);
    end
    Q = Q';
    R = A_1;
end