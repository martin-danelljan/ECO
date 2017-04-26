function [x,flag,relres,iter,resvec,p,rho,r_old] = pcg_ccot(A,b,opts,M1,M2,x0,p,rho,r_old,varargin)

% This is a modified version of Matlab's pcg function, that performs 
% preconditioned conjugate gradient.

if (nargin < 2)
    error(message('MATLAB:pcg:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);

m = sum(cellfun(@numel, b));
n = m;

tol = opts.tol;
maxit = opts.maxit;

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol <= eps
    warning(message('MATLAB:pcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:pcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 3) || isempty(maxit)
    maxit = min(n,20);
end

if opts.debug
    n2b = norm_cdcf(b); % Norm of rhs vector, b
end

if ((nargin >= 4) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:pcg:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 5) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:pcg:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

x = x0;

% Set up for the method
flag = 1;

r = cellfun(@minus, b, iterapp('mtimes',afun,atype,afcnstr,x,varargin{:}), 'uniformoutput', false);

if opts.debug
    normr = norm_cdcf(r);                   % Norm of residual
    normr_act = normr;
end

% Set old r
if nargin < 9 || opts.CG_use_FR
    r_old = [];
end

if opts.debug
    resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
    resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
else
    resvec = [];
    relres = [];
end

if nargin < 8 || isempty(rho)
    rho = 1;
end

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    if existM1
        y = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
    else % no preconditioner
        y = r;
    end
    
    if existM2
        z = iterapp('mldivide',m2fun,m2type,m2fcnstr,y,varargin{:});
    else % no preconditioner
        z = y;
    end
    
    rho1 = rho;
    rho = inner_product_cdcf(r, z);
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1 && (nargin < 7 || isempty(p)))
        p = z;
    else
        if opts.CG_use_FR
            % Use Fletcher-Reeves
            beta = rho / rho1;
        else
            % Use Polak-Ribiere
            rho2 = inner_product_cdcf(r_old, z);
            beta = (rho - rho2) / rho1;
        end
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        beta = max(0, beta);
        p = cellfun(@(z,p) z + beta * p, z, p, 'uniformoutput', false);
    end
    q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:});
    pq = inner_product_cdcf(p, q);
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        if opts.CG_standard_alpha
            alpha = rho / pq;
        else
            alpha = inner_product_cdcf(p, r) / pq;
        end
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Save old r if not using FR formula for beta
    if ~opts.CG_use_FR
        r_old = r;
    end
    
    x = cellfun(@(x,p) x + alpha * p, x, p, 'uniformoutput', false);             % form new iterate
    
    if ii < maxit || opts.debug
        r = cellfun(@(r,q) r - alpha * q, r, q, 'uniformoutput', false);
    end
    
    if opts.debug
        normr = norm_cdcf(r);
        normr_act = normr;
        resvec(ii+1,1) = normr;
    end
end                                % for ii = 1 : maxit

iter = ii;
if opts.debug
    relres = normr_act / n2b;
end

% truncate the zeros from resvec
if opts.debug
    if ((flag <= 1) || (flag == 3))
        resvec = resvec(1:ii+1,:);
    else
        resvec = resvec(1:ii,:);
    end
end
