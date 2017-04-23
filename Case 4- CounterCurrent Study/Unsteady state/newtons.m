function [x,fx,xx] = newtons(f,x0,TolX,MaxIter,varargin)
%% Defining the initial condition 
h = 1e-2; TolFun = 1e-6; EPS = 1e-6;
fx = feval(f,x0,varargin{:});
Nf = length(fx); Nx = length(x0);
if Nf ~= Nx, error('Incompatible dimensions of f and x0!'); end
if nargin < 4, MaxIter = 100; end
if nargin < 3, TolX = EPS; end
xx(1,:) = x0(:).'; 
House_Holder_iter = 10;
IterMod = 1;

%% this part is for handling Jocobina matrix
for k = 1: MaxIter
    %% these lines are for Broyden method
    if k == 1, S = 1; end

    H = -jacobian(f,xx(k+IterMod-1,:),h,varargin{:})^-1;
    dx = H * fx(:);
    
    %% these lines are for House holder method
    for itr=1:House_Holder_iter 
        
        if k > 1 && S > 0.3

            eta = norm(feval(f,xx(k+IterMod-1,:),varargin{:}))/norm(feval(f,xx(k+IterMod-2,:),varargin{:}));
            S = ((1+6*eta)^0.5-1)/(3*eta);

        end
        
        xx(k + IterMod,:) = xx(k + IterMod-1,:) + S * dx.';
        fx_k = fx;
        fx = feval(f,xx(k + IterMod,:),varargin{:});
        Y_K = fx - fx_k;
        H = H - ((H*Y_K' + S*dx)*dx'*H)/(dx'*H*Y_K');   
        dx = H * fx(:);
        IterMod = IterMod + 1; 
    end
    IterMod = IterMod -1;
    
    if norm(fx) < TolFun || norm(dx) < TolX , break; end
end
x = xx(k + IterMod,:);


function g = jacobian(f,x,h,varargin) 
%% this functin is for calculating the Jacobian matrix.
g=zeros(size(x,2));
if nargin < 3, h = 1e-4; end
h2 = 2*h; N = length(x); x = x(:).'; I = eye(N);
for n = 1:N
    g(:,n) = (feval(f,x + I(n,:)*h,varargin{:}) ...
    -feval(f,x - I(n,:)*h,varargin{:}))'/h2;
end
