function c = cumprod(a,varargin)

% CUMPROD   overloaded method for LOGICAL
%
% CUMPROD(X,DIM)
%

objwarn(mfilename)

c = cumprod( double(a) , varargin{:} );
