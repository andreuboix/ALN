function [maxs] = twonorm(A)
char_eqn = charpoly(A);
vaps = roots(char_eqn);
maxs = vaps(1);
n = length(A);
for i=1:n
    if(maxs < vaps(i)) maxs = vaps(i);
    end
end
end