function [x,rho,res,iter] = jacobi(A,b,x0,nmax,prec)
    %A = input('Enter the A matrix in brackets this way [x; x; x; x;]    ')
    %b = input('Enter the b vector in brackets this way [x; x; x; x;]    ')
    load('matrix200.mat')
    x0 = input('Enter the first initial solution this way [x; x; x; x;]     ')
    nmax = input('Enter the number of iterations you want me to do      ')
    prec = input('Enter the tolerance of this system         ')
    % ALTERNATIVELY:
    %%% load('matriu100.mat')
    %%% load('matriu150.mat')
    %%% load('matriu200.mat')
    rho = max(abs(eig(A)))
    if(rho >= 1)
        error('Error: the method is not convergent for this system')
    end
    n = length(A);
    [m,c] = size(A);
    if m ~= c 
        error('Error: the system introduced is not nxn')
    end
    for i=1:m
       if(A(i,i) < prec)
           error('Error: the system introduced is singular')
       end
    end
    D = diag(diag(A));
    B = eye(n) - (inv(D)*A);
    c = inv(D)*b;
    suma = 0;
    x = x0;
    for k = 2: nmax
        %%%Convergence criteria
        if (norm(b-A*x,2)/norm(b-A*x0,2)<prec)
            break
        end
        iter = k
        x = B*x+c
    end
    res = x;
    return
end
    function [x,rho,res,iter] = gaussS(A,b,x0,nmax,prec)
A = input('Enter the A matrix in brackets this way [x; x; x; x;]    ')
b = input('Enter the b vector in brackets this way [x; x; x; x;]    ')
x0 = input('Enter the first initial solution this way [x; x; x; x;]     ')
nmax = input('Enter the number of iterations you want me to do      ')
prec = input('Enter the tolerance of this system         ')
%%% ALTERNATIVELY: xlsread(for .xls), load(for ascii or .DAT), ...%%%
%%% ALTERNATIVELY:


%%% load('matriu100.mat')
%%% load('matriu150.mat')
%%% load('matriu200.mat')

rho = max(abs(eig(A)))

if(rho >= 1)
    error('Error: the method is not convergent for this system')
end
n = length(A);

[m,c] = size(A);
if m ~= c 
    error('Error: the system introduced is not nxn')
end
for i=1:m
   if(A(i,i) < prec)
       error('Error: the system introduced is singular')
   end
end
D = diag(diag(A));
L = tril(A,-1);
B =(D+w.*L)*((1-w).*D-w.*U);
c = w.*inv(L+D)*b;
suma = 0;
x = x0;

for k = 2: nmax
    if (norm(b-A*x,2)/norm(b-A*x0,2)<prec)
            break
        end
    iter = k
    x = B*x+c
end
res = x;
return
end
function [x,rho,res,iter] = overRelaxation(A,b,x0,w,nmax,prec)
A = input('Enter the A matrix in brackets this way [x; x; x; x;]    ')
b = input('Enter the b vector in brackets this way [x; x; x; x;]    ')
x0 = input('Enter the first initial solution this way [x; x; x; x;]     ')
nmax = input('Enter the number of iterations you want me to do      ')
prec = input('Enter the tolerance of this system         ')
w = input('Enter the omega of this system         ')
%%% ALTERNATIVELY: xlsread(for .xls), load(for ascii or .DAT), ...%%%
%%% ALTERNATIVELY:


%%% load('matriu100.mat')
%%% load('matriu150.mat')
%%% load('matriu200.mat')

rho = max(abs(eig(A)))

if(rho >= 1)
    error('Error: the method is not convergent for this system')
end
n = length(A);

[m,c] = size(A);
if m ~= c 
    error('Error: the system introduced is not nxn')
end
for i=1:m
   if(A(i,i) < prec)
       error('Error: the system introduced is singular')
   end
end
D = diag(diag(A));
L = tril(A,-1);
B = eye(n) - (inv(D)*A);
c = inv(D)*b;
suma = 0;
x = x0;

for k = 2: nmax
    if (norm(b-A*x,2)/norm(b-A*x0,2)<prec)
            break
        end
    iter = k
    x = B*x+c
end
res = x;
return
end


%for i = 1:nmax
    %for j = 1:nmax
     %   if(i != j) 
        %    suma = suma + A(i,j)*x(j);
    %    end
  %  end
  %  x(i) = suma/A(i,i);
%end

%for i = 1:nmax
    %for j = 1:nmax
     %   if(i != j) 
        %    suma = suma + A(i,j)*x(j);
    %    end
  %  end
  %  x(i) = suma/A(i,i);
%end


    %for i = 1:nmax
        %for j = 1:nmax
         %   if(i != j) 
            %    suma = suma + A(i,j)*x(j);
        %    end
      %  end
      %  x(i) = suma/A(i,i);
    %end