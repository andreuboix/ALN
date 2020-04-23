% Now solving Rx = Q^(t)y;
format long
function [x] = solveutri(Q,R,y)
y = z;
x = R\(transpose(Q)*y);
return
end
