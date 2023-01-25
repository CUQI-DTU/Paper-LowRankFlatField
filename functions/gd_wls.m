function [u] = gd_wls(A,F,Y,options)
% GD_WLS   Gradient descent method for X-ray tomographic
% reconstructions based on the model
%
%   minimize    J(u) + lambda*hTV(u,delta)
%   subject to  u >= 0 
%
% where u is the attenuation image, lambda is a nonnegative regularization
% parameter.   
%
% The function J(u)  is the weighted least squares approximation.
%
%       J(u) = 0.5*(A*x-b)*inv(Sigma)*(A*x-b)
%
%       with b = -log(Y(:)./repmat(mean(F,2),p))
%       Sigma  = diag(1./y)
%
% The function hTV denotes discrete total variation with smoothing,
% i.e.,   
%
%   hTV(u;delta) = sum_i Huber(||Di*u||_2; delta)
% 
% where Huber(t,delta) denotes the Huber norm with parameter delta.
%   
% The inputs A, Y, and F are required: A (m-by-n) is the system matrix
% (or a Spot operator), Y (r-by-p) is the matrix of measurements
% where each column corresponds to a projection, and F (r-by-s) is
% a matrix with s flat-field samples. The system matrix should map
% the vector of attenuation coefficients u into the vectorized
% sinogram, ie., b = A*x if B = reshape(b,r,p) is the sinogram.
%
% The optional input 'options' is a struct with one or more fields:
%
%   'maxiters'   Maximum number of iterations  (default: 200)
%   'tolf'       Tolerance, rel. obj. value    (default:1e-8)
%   'rho'        Gradient step multiplier      (default: 1.5)
%   'lambda'     TV-regularization parameter   (default: 0.0)
%   'tau'        TV-smoothing parameter        (default: 1e-2)
%   'u0'         Starting point                (default: zeros(n,1))
%   'uhold'      List of indices of intermediate iterates to return
%
% The return value u is the vectorized reconstruction of the
% attenuation coefficients, and iterinfo is a struct with
% information pertaining to the iterates and the reconstruction:
%
%   'ngrad'      Normalized gradient image (vectorized)
%   'optu'       Attenuation image optimality cond. (vectorized)
%   'uhold'      Intermediate iterates (if options.uhold is specified)
%
% Reference:
%   Hari Om Aggrawal, Martin S. Andersen, Sean Rose, and Emil Sidky,
%   "A Convex Reconstruction Model for X-ray tomographic Imaging with
%   Uncertain Flat-fields", submitted to IEEE Transactions on
%   Computational Imaging, 2017. 
%
% License: 
%   GPL-3
%
% Authors: 
%   Hari Om Aggrawal (hariom85@gmail.com) 
%   Martin S. Andersen (mskan@dtu.dk)
% 
% Date: 
%   April 19, 2017
% 
% Modified by Katrine O. Bangsgaard, October 2021.

    
% Check/extract problem dimensions
[r,p] = size(Y);
n = size(A,2);
N = sqrt(n);
s = size(F,2);
assert(size(A,1) == r*p)
assert(size(F,1) == r)

% Attenuation hyperparameter (TV-regularization)
if isfield(options,'lambda')
    lambda = options.lambda;
    tv = 0.0;
else
    lambda = 0.0;
    tv = 0.0;
end

% Maximum numer of iterations
if isfield(options,'maxiters')
    maxiters = options.maxiters;
else
    maxiters = 200;
end

% Stopping tolerance for relative change in objective function
if isfield(options,'tolf')
    tolf = options.tolf;
else
    tolf = 1e-8;
end

% Save intermediate iterations
if isfield(options,'uhold')
    uhold_idx = options.uhold;
    uhold = zeros(n,length(uhold_idx));
end

% Step-size multiplier
if isfield(options,'rho')
    rho = options.rho;
    assert(rho > 0)
    if rho >= 2.0
        warning('rho >= 2.0: convergence not guaranteed')
    end
else
    rho = 1.5;
end

% Total variation smoothing
if isfield(options,'tau')
    tau = options.tau;
    assert(tau > 0)
else
    tau = 1e-2;
end

% Verbose
if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 0;
end

% Flat-field ML estimate
vh = mean(F,2);

% Check model and set up objective function
y = Y(:);
if all(y>0)

% Preprocessing of sinogram
if isfield(options,'preproc')
    yy = reshape(y./repmat(vh,p,1),r,p)';
    preproc = options.preproc;
    if preproc == 1
        b = (xRemoveStripesVertical(-log(yy),3,'db5',0.9))';
        if mod(N,2) ~= 0 
            b = b(1:N,:);
        end
    else
        b = tsino_nonlocal_medfilt(-log(yy'), 31);
    end

    b = b(:);
else
    b = -log(y./repmat(vh,p,1));
end


J = @(u) Jwls(u,A,b,y);


% Starting point
if isfield(options,'u0')
    u = options.u0;
else
    u = zeros(n,1);
end

% Circular mask
if isfield(options,'mask')
    mask = options.mask;
else
    mask = false;
end
if mask    
    [XX,YY] = meshgrid(linspace(-1,1,N),linspace(-1,1,N));
    M = (XX(:).^2 + YY(:).^2 > 1);
    clear('XX','YY');
end    
    
% Finite-difference matrix with Neumann boundary conditions
Dn = spdiags([ones(N-1,1),-ones(N-1,1);0,-1],[0,1],N,N);
Dy = kron(speye(N),Dn);
Dx = kron(Dn,speye(N));

% Compute Lipschitz constant (upper bound)
if ismatrix(A) && ~isa(A,'opFoG')
    if verbose
        fprintf(1,'Computing Lipschitz constant\n');
    end
    L = normest(spdiags(sqrt(y),0,r*p,r*p)*A)^2;
elseif isa(A,'opFoG')
    % Power iteration
    if verbose
        fprintf(1,'Estimating Lipschitz constant via power iteration\n');
    end
    tmp = 1/sqrt(n)*ones(n,1);
    for it = 1:20
        tmp = A'*(y.*(A*tmp/norm(tmp)));
    end
    L = norm(tmp);
    clear('tmp');
end
L = L + lambda*8/tau;
t = rho/L;    

% Gradient descent method
if verbose
    fprintf(1,'%4s %12s %8s %8s\n','it.','objval','opt','step'); 
end
fhist = zeros(maxiters+1,1);
fold = 0;
for it = 1:maxiters+1
    if exist('uhold_idx','var')==1 && any(uhold_idx==it)
        uhold(:,uhold_idx==it) = u;
    end
    % Compute objective value and gradient
    [val,grad] = J(u);
    fhist(it) = val;    
    if lambda > 0
        [tv,tvgrad] = huber_tv(Dx,Dy,u,tau);
        grad = grad + lambda*tvgrad;
        fhist(it) = fhist(it) + lambda*tv;
    end
    if verbose && mod(it-1,50) == 0
        nopt = t*abs((u>0).*grad + (u<=0).*min(grad,0));
        if mask
            nopt(M) = 0;
        end
        fprintf(1,'%4i %12.6e %8.2e %8.2e\n',it-1,fhist(it),norm(nopt),t*norm(grad));
    end
    
    % Stopping criteria
    STOP(1) = abs(fhist(it)-fold)  <= tolf*(1+abs(fold));
    STOP(2) = (it >= maxiters + 1);

    if STOP(1) || STOP(2)
        break
    end
    fold = fhist(it); % update function value
    
    % Update u
    u = u-t*grad;
    %u =max(0,u-t*grad);
    if mask
        u(M) = 0;
    end
    
end
else
    u = zeros(n,1);
end

% Return info struct
% iterinfo = struct('fhist',fhist);
% 
% if exist('uhold','var') == 1
%     iterinfo.uhold = uhold;
% end
% 
% iterinfo.objval = [val, tv];
% iterinfo.ngrad = grad*t;
% if mask
%     iterinfo.ngrad(M) = 0;
% end
% iterinfo.Au = A*u;
% iterinfo.optu = abs((u>0).*iterinfo.ngrad + (u<=0).*min(iterinfo.ngrad,0));
% 
% iterinfo.options = options;
% iterinfo.maxiter = it-1;
% iterinfo.relobjval = abs(fhist(it)-fold)/(1+abs(fold));

end

function [tv, grad] = huber_tv(Dx,Dy,u,tau)
% Huber total variation
g = [Dx*u Dy*u];
ng = sqrt(sum(g.^2,2));
tv = sum((ng <= tau).*(ng.^2/(2*tau)) + (ng > tau).*(ng - tau/2));
grad = Dx'*(g(:,1)./max(tau,ng)) + Dy'*(g(:,2)./max(tau,ng));
end

function [val, grad] = Jwls(u,A,b,y)
% Weighted least-squares approximation
r = A*u-b;
grad = A'*(y.*r);
val = 0.5*r'*(y.*r);
end

