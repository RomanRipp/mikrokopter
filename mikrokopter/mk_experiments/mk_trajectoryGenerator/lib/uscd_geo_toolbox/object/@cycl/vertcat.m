function c = vertcat(varargin)

% VERTCAT  Vertical concatenation of CYCL-arrays
%
% B = VERTCAT( A1 , A2 , A3 , ... );
%

c = cat(1,varargin{:});
