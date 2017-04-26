function y = iterapp(op,afun,atype,afcnstr,x,varargin)
%ITERAPP   Apply matrix operator to vector and error gracefully.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X) applies matrix operator AFUN to vector
%   X. If ATYPE is 'matrix, then AFUN is a matrix and the OP is applied
%   directly. OP is either 'mtimes' or 'mldivide'.
%   ATYPE and AFCNSTR are used in case of error.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X,P1,P2,...) allows extra arguments to
%   AFUN(X,P1,P2,...) although this usage is now discouraged in favor of
%   using anonymous functions.
%   AFUN(X,P1,P2,...,PN,TFLAG) should accept a TFLAG as its final input if
%   the calling function is BICG, LSQR or QMR. TFLAG is either 'transp' or
%   'notransp' depending on whether A' OP X or A OP X is required.
%   ITERAPP is designed for use by iterative methods like PCG which
%   require matrix operators AFUN representing matrices A to operate on
%   vectors X and return A*X and may also have operators MFUN representing
%   preconditioning matrices M operate on vectors X and return M\X.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2013 The MathWorks, Inc.

if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' * x;
            else
                y = afun * x;
            end
        case 'mldivide'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' \ x;
            else
                y = afun \ x;
            end
        otherwise
            error(message('MATLAB:iterapp:InvalidOp'))
    end
else
%     try
        if (nargin >= 6) && isequal(varargin{end},'notransp')
            % A request for A*x coming from BICG, LSQR and QMR
            try
                % New syntax: we now request afun(x,P1,P2,...,PN,'notransp')
                y = afun(x,varargin{:});
            catch
                % Old syntax: we used to accept afun(x,P1,P2,...,PN)
                y = afun(x,varargin{1:end-1});
            end
        else
            % A request for A*x
            % coming from BICGSTAB, CGS, GMRES, MINRES, PCG or SYMMLQ
            % with the call always afun(P1,P2,...,PN)
            % or a request for A'*x coming from
            % BICG, LSQR and QMR in the afun(x,P1,P2,...,PN,'transp') case
            y = afun(x,varargin{:});
        end
%     catch ME
%         error(message('MATLAB:iterapp:InvalidInput', atype,afcnstr, ME.message));
%     end

%     if ~iscolumn(y)
%         error(message('MATLAB:iterapp:MustReturnColumn', atype, afcnstr));
%     end
end
