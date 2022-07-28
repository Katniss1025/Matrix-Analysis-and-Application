function [Q,R]=Givens(A)
    [m,n] = size(A);
    R = A;
    Q = eye(m);
    for i = 1:n-1 
        for j = i+1:m 
            x = R(:,i);
            rt = GivensTrans(x,i,j);%Jæÿ’Û
            %r = blkdiag(eye(i-1),rt)
            Q = Q*rt';
            R = rt*R;
        end
    end
end

function [R,y] = GivensTrans(x,i,j)
    xi = x(i);          
    xj = x(j);
    r = sqrt(xi^2+xj^2);
    c = xi/r;
    s = xj/r;
    R = eye(length(x));
    R(i,i) = c;
    R(i,j) = s;
    R(j,i) = -s;
    R(j,j) = c;
    y = x(:);
    y([i,j]) = [r,0];
end

