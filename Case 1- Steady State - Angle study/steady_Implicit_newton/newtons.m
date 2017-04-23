function [x,fx,xx] = newtons(f,x0,TolX,MaxIter,varargin)
%% Defining the initial condition 
h = 1e-2; TolFun = 10^-4; EPS = 1e-4;
fx = feval(f,x0,varargin{:});
Nf = length(fx); Nx = length(x0);
if Nf ~= Nx, error('Incompatible dimensions of f and x0!'); end
if nargin < 4, MaxIter = 100; end
if nargin < 3, TolX = EPS; end
xx(1,:) = x0(:).'; 
%% this part is for handling Jocobina matrix
for k = 1: MaxIter
    
dx = -jacobian(f,xx(k,:),h,varargin{:})\fx(:);
xx(k + 1,:) = xx(k,:) + dx.';
fx = feval(f,xx(k + 1,:),varargin{:}); fxn = norm(fx);

if fxn < TolFun | norm(dx) < TolX, break; end

end
x = xx(k + 1,:);
%if k == MaxIter, fprintf('The best in %d iterations\n',MaxIter), end