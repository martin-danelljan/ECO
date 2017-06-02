function [x,resvec,state] = pcg_ccot(A,b,opts,M1,M2,ip,x0,state)

% This is a modified version of Matlab's pcg function, that performs 
% preconditioned conjugate gradient.


% tol = opts.tol;
maxit = opts.maxit;

if ~isfield(opts, 'init_forget_factor')
    opts.init_forget_factor = 1;
end

if opts.debug
    n2b = sqrt(ip(b,b)); % Norm of rhs vector, b
end

existM1 = ((nargin >= 4) && ~isempty(M1));
existM2 = ((nargin >= 5) && ~isempty(M2));

x = x0;

% Load the CG state
p = [];
rho = 1;
r_prev = [];
load_state = nargin > 7 && ~isempty(state) && opts.init_forget_factor > 0;
if load_state
    if isfield(state, 'p')
        p = state.p;
    end
    if isfield(state, 'rho') && ~isempty(state.rho)
        rho = state.rho / opts.init_forget_factor;
    end
    if isfield(state, 'r_prev') && ~opts.CG_use_FR
        r_prev = state.r_prev;
    end
end

% Set up for the method
state.flag = 1;

% r = cellfun(@minus, b, iterapp('mtimes',afun,atype,afcnstr,x,varargin{:}), 'uniformoutput', false);
r = cellfun(@minus, b, A(x), 'uniformoutput', false);

if opts.debug
    normr = sqrt(ip(r,r));                   % Norm of residual
    normr_act = normr;
end


if opts.debug
    resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
    resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
else
    resvec = [];
    relres = [];
end

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    if existM1
        y = M1(r);
    else % no preconditioner
        y = r;
    end
    
    if existM2
        z = M2(y);
    else % no preconditioner
        z = y;
    end
    
    rho1 = rho;
    rho = ip(r, z);
    if ((rho == 0) || isinf(rho))
        state.flag = 4;
        break
    end
    
    if (ii == 1 && isempty(p))
        p = z;
    else
        if opts.CG_use_FR
            % Use Fletcher-Reeves
            beta = rho / rho1;
        else
            % Use Polak-Ribiere
            rho2 = ip(r_prev, z);
            beta = (rho - rho2) / rho1;
        end
        if ((beta == 0) || isinf(beta))
            state.flag = 4;
            break
        end
        beta = max(0, beta);
        p = cellfun(@(z,p) z + beta * p, z, p, 'uniformoutput', false);
    end
    
    q = A(p);
    pq = ip(p, q);
    if ((pq <= 0) || isinf(pq))
        state.flag = 4;
        break
    else
        if opts.CG_standard_alpha
            alpha = rho / pq;
        else
            alpha = ip(p, r) / pq;
        end
    end
    if isinf(alpha)
        state.flag = 4;
        break
    end
    
    % Save old r if not using FR formula for beta
    if ~opts.CG_use_FR
        r_prev = r;
    end
    
    % form new iterate
    x = cellfun(@(x,p) x + alpha * p, x, p, 'uniformoutput', false);
    
    if ii < maxit || opts.debug
        r = cellfun(@(r,q) r - alpha * q, r, q, 'uniformoutput', false);
    end
    
    if opts.debug
        normr = sqrt(ip(r,r));
        normr_act = normr;
        resvec(ii+1,1) = normr;
    end
end

iter = ii;
if opts.debug
    relres = normr_act / n2b;
end

% truncate the zeros from resvec
if opts.debug
    if ((state.flag <= 1) || (state.flag == 3))
        resvec = resvec(1:ii+1,:);
    else
        resvec = resvec(1:ii,:);
    end
end

% Save the state
if nargout > 2
    state.p = p;
    state.rho = rho;
    if ~opts.CG_use_FR
        state.r_prev = r_prev;
    end
end