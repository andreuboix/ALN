A = input("Enter the A Matrix in brackets this way [xxx; xxx; xxx; xxx;]")
function [Q,R] = QRDECOMPOSITION(A)
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
copiaA = A;
for i = 1:n
    R(i,i) = norm(Q(:,i),2);
    Q(:,i) = Q(:,i)\R(i,i);
    for j = i:n
        R(i,j) = dot(A(:,i),A(:,j));
        Q(:,j) = Q(:,j) - R(i,j) * Q(:,i);
    end
end
return Q,R;
end

% Now solving Rx = Q^(t)y;
function [x] = solveutri(Q,R,y)
x = R\(transpose(Q)*y);
return x;
end
