

%%% Mètode de la potència Estàndard. Suposem que té un VAP dominant.
%%% Computem el valor absolut d'aquest VAP dominant!
function [lambda] = powermethod(A,x0,nmax)
    format long
    A = input('Enter the A matrix this way [...; ...; ...;]     ')
    x0 = input('Enter an approximation of the dominant eigenvector this way [.;.;.;]     ')
    nmax = input('Enter the maximum number of iterations you would want         ')
    x = x0;
    %%% RAYLEIGH QUOFICIENT %%%
    for i=1:nmax
        x = A*x;
        x = x/max(x);
    end
    %%% DISPLAYING RESULTS %%%
    lambda1 = A*x/x;
    disp('So, the approximation for the dominant eigenvalue is:')
    [a,ii] = max(abs(lambda1));
    lambda = max(max(abs(lambda1(ii,:))));
    disp(lambda)
    disp('And the approximation for the dominant eigenvector is:')
    vep = abs(x)/abs((x0(1))^nmax);
    disp(vep)
    disp('The error commited is:')
    disp(norm(max(abs(eig(A))) - lambda,2))
end



%%% Mètode de la potència inversa estàndard. Suposem que té un VAP
%%% dominant.
function [lambda] = invpowermethod(A,x0,nmax)
    format long
    A = input('Enter the A matrix this way [...; ...; ...;]     ')
    x0 = input('Enter an approximation of the dominant eigenvector this way [.;.;.;]     ')
    nmax = input('Enter the maximum number of iterations you would want         ')
    x = x0;
    A1 = inv(A);
    %%% RAYLEIGH QUOFICIENT %%%
    for i=1:nmax
        x = A1*x;
        x = x/max(x);
    end
    %%% DISPLAYING RESULTS %%%
    lambda1 = A1*x/x;
    disp('So, the approximation for the dominant eigenvalue is:')
    [a,ii] = max(abs(lambda1));
    lambda = max(max(abs(lambda1(ii,:))));
    disp(lambda)
    disp('And the approximation for the dominant eigenvector is:')
    vep = abs(x)/abs((x0(1))^nmax);
    disp(vep)
    disp('The error commited is:')
    disp(norm(max(abs(eig(A1))) - lambda,2))
end




%%% Mètode de la potència desplaçada. Suposem que té un VAP dominant.
%%% Si lambda VAP de A de VEP v, $\implies$ lambda - q VAP de A-qId de VEP
%%% v, q $\in \mathbb{R}$
function [lambda] = desppowermethod(A,x0,nmax)
    format long
    A = input('Enter the A matrix this way [...; ...; ...;]     ')
    x0 = input('Enter an approximation of the dominant eigenvector this way [.;.;.;]     ')
    nmax = input('Enter the maximum number of iterations you would want         ')
    q = input('Enter the translation you would want         ')
    x = x0;
    A2 = A - q*(eye(length(A)));
    %%% RAYLEIGH QUOFICIENT %%%
    for i=1:nmax
        x = A2*x;
        x = x/max(x);
    end
    %%% DISPLAYING RESULTS %%%
    lambda1 = A2*x/x;
    disp('So, the approximation for the dominant eigenvalue is:')
    [a,ii] = max(abs(lambda1));
    lambda = max(max(abs(lambda1(ii,:))));
    disp(lambda)
    disp('And the approximation for the dominant eigenvector is:')
    vep = abs(x)/abs((x0(1))^nmax);
    disp(vep)
    disp('The error commited is:')
    disp(norm(max(abs(eig(A2))) - lambda,2))
end

function [lambda] = despinvpowermethod(A,x0,nmax)
    format long
    A = input('Enter the A matrix this way [...; ...; ...;]     ')
    x0 = input('Enter an approximation of the dominant eigenvector this way [.;.;.;]     ')
    nmax = input('Enter the maximum number of iterations you would want         ')
    q = input('Enter the translation you would want         ')
    x = x0;
    A3 = inv(A - q*(eye(length(A))));
    %%% RAYLEIGH QUOFICIENT %%%
    for i=1:nmax
        x = A3*x;
        x = x/max(x);
    end
    %%% DISPLAYING RESULTS %%%
    lambda1 = A3*x/x;
    disp('So, the approximation for the dominant eigenvalue is:')
    [a,ii] = max(abs(lambda1));
    lambda = max(max(abs(lambda1(ii,:))));
    disp(lambda)
    disp('And the approximation for the dominant eigenvector is:')
    vep = abs(x)/abs((x0(1))^nmax);
    disp(vep)
    disp('The error commited is:')
    disp(norm(max(abs(eig(A3))) - lambda,2))
end



