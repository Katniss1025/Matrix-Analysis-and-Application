clear all;
%% main
Main();
function Main
    A = shuruA();
    calculate(A);
end

%% function
function calculate(A)
    para = shuru_para();
    switch para
        case 'LU'
            if ~isSquare(A)
                disp('A is not a square matrix, please try again!');
                return;
            end
            if sequential(A)==0
                disp('A can not be decomposed, please try again！');
                return;
            end
            [L,U] = LU(A);%求A的LU分解
            disp('matrix L=');
            disp(L);
            disp('matrix U=');
            disp(U);
            disp('det(A)=');
            disp(LU_det(L,U));%计算A的行列式
            b = calculateX();
            if b ~= 1e10 %yes
                x = Ax_b_LU(L,U,b);
                disp('Ax=b solution x=');
                disp(x);
            end
            return;
        case 'QR'
            flag = isLinearIndependent(A);
            if flag == 0
                disp('A with liner dependent columns can not be decomposed as A = QR, please try again!');
                return;
            else
                [Q,R] = QR(A);
                disp('matrix Q=');
                disp(Q);
                disp('matrix R=');
                disp(R);
                if isSquare(A)
                    disp('det(A)=');
                    detA = QR_det(Q,R);
                    disp(detA);
                else
                    disp('A is not square matrix,we can not calculate det(A)')
                end
                b = calculateX();
                if b ~= 1e10 %yes                   
                    x = Ax_b_QR(Q,R,b);
                    disp('Ax=b solution x=');
                    disp(x);
                end
                return;
            end
        case 'Householder'
            [Q,R] = Householder(A);
            disp('matrix Q=');
            disp(Q);
            disp('matrix R=');
            disp(R);
            if isSquare(A)
                disp('det(A)=');
                detA = QR_det(Q,R);
                disp(detA);
            else
                disp('A is not square matrix,we can not calculate det(A)')
            end
            b = calculateX();
            if rank(A)==size(A,2)
                x = Ax_b_QR(Q,R,b);
                disp('Ax=b solution x=');
                disp(x);
            elseif rank(A)~=rank([A,b])
                disp('Ax=b has no solution.');
            else
                disp('特解为:');
                disp(A\b);
                disp('通解为：');
                disp(null(A,'r'));
            end
            return;
        case 'Givens'
            [Q,R] = Givens(A);
            disp('matrix Q=');
            disp(Q);
            disp('matrix R=');
            disp(R);
            if isSquare(A)
                disp('det(A)=');
                detA = QR_det(Q,R);
                disp(detA);
            else
                disp('A is not square matrix,we can not calculate det(A)')
            end
            b = calculateX();
            if rank(A)==size(A,2)
                x = Ax_b_QR(Q,R,b);
                disp('Ax=b solution x=');
                disp(x);
            elseif rank(A)~=rank([A,b])
                disp('Ax=b has no solution.');
            else
                disp('特解为:');
                disp(A\b);
                disp('通解为：');
                disp(null(A,'r'));
            end
            return;
        case 'URV'
            [U,R,V] = URV(A);
            disp('matrix U=');
            disp(U);
            disp('matrix R=');
            disp(R);
            disp('matrix V=');
            disp(V);
            if isSquare(A)
                disp('det(A)=');
                detA = URV_det(U,R,V);
                disp(detA);
            else
                disp('A is not square matrix,we can not calculate det(A)')
            end
            b = calculateX();
            if b ~= 1e10 %yes
                [pinvA,x] = URV_b_LU(U,R,V,b);
                disp('Moore-Penrose Pseudoinverse of A=');
                disp(pinvA);
                disp('Ax=b solution x=');
                disp(x);
            end
            
            return;
       
        otherwise
            disp('The parameter is not valid, please try again！');
            calculate(A);
    end
end

function A = shuruA
    A = '1.请输入即将分解的矩阵A：';
    try
        A = input(A);
    catch
        warning('A is not valid, please check the input!');
    end
end
function para = shuru_para
    para = '2.请输入矩阵分解参数（如LU）:';
    para = input(para,'s');
end

function flag = sequential(A) %square A 是否存在det=0的顺序主子式
    flag = 1;
    n = size(A);
    for i = 1:n
        mat = A(1:i,1:i);
        if det(mat)==0
            flag = 0;
            break;
        end
    end
end

function flag = isSquare(A)
    [m,n] = size(A);
    flag = 0;
    if m==n
        flag = 1;
    end
end

function detA = LU_det(L,U)
    [n1,~] = size(L);
    [n2,~] = size(U);
    detL = 1;
    detU = 1;
    for i = 1:1:n1
        detL = detL * L(i,i);
    end
    for i = 1:1:n2
        detU = detU * U(i,i);
    end
    detA = detL * detU;
end

function b = calculateX
    para_calculateX = 'Calculate the solution of Ax=b?(yes/no)';
    para_calculateX = input(para_calculateX,'s');
    if isequal(para_calculateX,'yes')
        b = 'please enter the vetor b:';
        try
            b = input(b);
        catch
            warning('b is not valid');
        end
    elseif isequal(para_calculateX ,'no')
        b = 1e10;
    else 
        disp('Not Valid!');
        return;
    end
end

function x = Ax_b_LU(L,U,b)
    % L:lower triangular
    % U:Upper triangular
    % LUx = b,Ly = b
    x = zeros(length(b),1);
    y = zeros(length(b),1);
    for i = 1:length(b) %回代法求y
        for j = 1:1:i-1
            y(i) = y(i) + L(i,j)*y(j);
        end
        y(i) = (b(i)-y(i))/L(i,i);
    end
    % Ux = y,回代法求x
    for i = length(b):-1:1
        for j = i+1:length(b)
            x(i) = x(i) + U(i,j)*x(j);
        end
        x(i) = (y(i)-x(i))/U(i,i);
    end
end

function flag = isLinearIndependent(A)
    r = rank(A);
    [~,n] = size(A);
    if r == n
        flag = 1;
    else 
        flag = 0;
    end
end

function x = Ax_b_QR(Q,R,b)
    % Ax = b, A = QR
    % QRx = b
    % Rx = Q'b = y
    % R: upper triangular
    x = zeros(size(Q*R,2),1);
    try
        y = Q'* b;
    catch
        warning("Error,b's dimension is not correct!");
    end
    % Rx = y,回代法求x
    for i = length(x):-1:1
        for j = i+1:length(x)
            x(i) = x(i) + R(i,j)*x(j);
        end
        x(i) = (y(i)-x(i))/R(i,i);
    end
end
function detA = QR_det(Q,R)
    % Q is n*n orthonormal matrix
    % |Q| = (-1)^((n-trace(A))/2)
    detQ = (-1)^((length(Q)-trace(Q))/2);
    % R is upper triangular matrix
    detR = 1;
    for i = 1:length(R)
        detR = detR * R(i,i);
    end
    detA = det(Q) * detR;
end
function detA = URV_det(U,R,V)
    detA = det(U)*det(R)*det(V');
end
function [pseu_inv,x] = URV_b_LU(U,R,V,b)
    m = size(V,2);
    n = size(U',1);
    R1 = zeros(m,n);
    R1(1:rank(R),1:rank(R)) = inv(R(1:rank(R),1:rank(R)));
    pseu_inv = V*R1*U';
    x = pseu_inv*b;
end