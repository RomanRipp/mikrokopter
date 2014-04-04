function c = cumprod(a,varargin)

% CUMPROD   overloaded method for UINT8
%
% CUMPROD(X,DIM)
%

c = uint8( cumprod( double(a) , varargin{:} ));
