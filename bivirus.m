function dxdt = bivirus(t,x)
%BIVIRUS Differential equation of bivirus dynamics
%   The global variables are the parameters of the bivirus model. Also, dfx
%   gives the Jacobian of the bivirus system. The vector x = [x^1', x^2']'
%   is a column vector where x^i is the n-vector of infection probability
%   for all nodes, for virus i


global D1 D2 B1 B2 n dfx

x1 = x(1:n); x2 = x(n+1:end);   %The two vectors of virus infection prob
alpha1 = 10;
alpha2 = 10;
dx1dt = alpha1.*(-D1 + (eye(n) - diag(x1) - diag(x2))*B1)*x1;

dx2dt = alpha2.*(-D2 + (eye(n) - diag(x1) - diag(x2))*B2)*x2;

dxdt = [dx1dt', dx2dt']';
end

