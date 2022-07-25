function dxdt = sis(t,x)

global D A

dxdt = (-D + (eye(length(x)) - diag(x))*A)*x;

end

