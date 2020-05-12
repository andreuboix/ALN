
function [coefs, norm2Res] = polminquad(x,y,grau,plt)
    format long;
    %%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%
    grau = input("Enter the degree of the quadratic model   ")
    x = input("Enter the x vector in brackets this way [x; x; x; x;]    ")
    y = input("Enter the y=f(x) vector in brackets this way [x; x; x; x;]    ")
    plt = input("Draw or not? 1 or 0   ")
    
    %%%%%%%%%%%%%%%% ALTERNATIVELY, x,y = xlsread('XXXX.xls') %%%%%%%%%%%%%%%%%%%%%
    m = length(x);
    q = grau;
    
    %%%%%%%% ERRORS %%%%%%%%%%%%%%%%%%%%%%%%
    if(m ~= length(y))
        error('Error: different lengths of x and y')
    end
    %%%%%%%%%%%%%%%%%  PREPARATION OF THE OVERDETERMINED SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:q
       A(:,i) = x.^(q+1-i);
    end
    A(:,q+1) = ones(1,m);
   
    [m,n] = size(A);
   
    %%%%%%%%      QR DECOMPOSITION       %%%%%%%%%%%%%%%
    Q = A;
    for k = 1:n
        R(k,k) = norm(Q(:,k),2);
        Q(:,k) = Q(:,k)/R(k,k);
        R(k, k+1:n) = Q(:,k)'*Q(:,k+1:n);
        Q(:,k+1:n) = Q(:,k+1:n)-Q(:,k)*R(k,k+1:n);
    end
    
    
    %%%%%%%% RESOLUTION OF THE OVERDETERMINED SYSTEM %%%%%%%%%%%%%
    coefs = R\(Q'*y);
    norm2Res = norm(A*coefs - y,2);
    identity = eye(q+1);
    falseidentity = (transpose(Q))*Q;
    norm_inf_qq = norm(falseidentity - identity,inf);
    
    disp('Job done! The coefficients of the a vector are:');
    disp(coefs);
    disp('The norm of the residual vector is: ');
    disp(norm2Res);
    disp('The norm of Q^TQ - Id is:');
    disp(norm_inf_qq);
    %%%%%%%%%%%%%%%%%%%%%%%%%%% WE DRAW IFF PLT == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plt == 1;
    figure()
    scatter(x,y,'k','filled');
    xx = x(1):0.01: x(end);
    yy = polyval(coefs,xx);
    hold on
    plot(xx,yy,'-g');
    xlabel('x');
    ylabel('y');
    title(['Quadratic model of the polynomial aproximation, degree ',int2str(grau)]);
    hold off
    end
end






