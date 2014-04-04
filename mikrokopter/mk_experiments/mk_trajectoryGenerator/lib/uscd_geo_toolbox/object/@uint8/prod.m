function c = prod(a,varargin)

% PROD   overloaded method for UINT8
%
% PROD(X,DIM)
%

c = uint8( prod( double(a) , varargin{:} ) );
