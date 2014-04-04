function c = diff(a,varargin)

% DIFF   overloaded method for UINT8
%
% DIFF(X,N,DIM)
%

c = uint8( diff( double(a) , varargin{:} ) );
