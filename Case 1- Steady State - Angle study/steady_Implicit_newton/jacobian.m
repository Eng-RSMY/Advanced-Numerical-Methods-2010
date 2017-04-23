function g = jacobian(f,x,h,varargin) 
%% this functin is for calculating the Jacobian matrix.
if nargin < 3, h = 1e-4; end
h2 = 2*h; N = length(x); x = x(:).'; I = eye(N);
for n = 1:N
g(:,n) = (feval(f,x + I(n,:)*h,varargin{:}) ...
-feval(f,x - I(n,:)*h,varargin{:}))'/h2;
end