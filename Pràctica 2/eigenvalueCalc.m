function [Avaps] = vap(A)
char_eqn = charpoly(A);
vaps = roots(char_eqn);
Avaps = diag(vaps);
return
end
