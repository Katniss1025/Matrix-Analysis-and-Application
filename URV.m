function [U,R,V] = URV(A)
    U1 = orth(A);
    U2 = null(A');
    V1 = orth(A');
    V2 = null(A);
    U = [U1 U2];
    V = [V1 V2];
    R = U'*A*V;
end
