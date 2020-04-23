A = input("Enter the A Matrix in brackets this way [xxx; xxx; xxx; xxx;]")
%A = readfile('')
format long
function [Q,R] = qr(A)
[m,n] = size(A);
Q = A;
copiaA = A;
for i = 1:n
    R(i,i) = norm(Q(:,i),2);
    Q(:,i) = Q(:,i)/R(i,i);
    for j = i:n
        R(i,j) = dot(A(:,i),A(:,j));
        Q(:,j) = Q(:,j) - R(i,j) * Q(:,i);
    end
end
return
end