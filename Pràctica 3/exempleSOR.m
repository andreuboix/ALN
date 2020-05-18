A = readmatrix('matriu-A.dat.txt');
n = sqrt(length(A));
A= reshape(A,n,n);


D = diag(diag(A));
L = tril(A,-1);
U = triu(A,+1);

%find() per trobar max i min.
for i=1:2000
    omega(i)= 0.001*i;
    Bomega= inv(D+omega(i)*L)*((1-omega(i))*D-omega(i)*U);
    rho = max(abs(eig(Bomega)));
    scatter(omega(i),rho,'b','filled');
    hold on;
end
hold off;