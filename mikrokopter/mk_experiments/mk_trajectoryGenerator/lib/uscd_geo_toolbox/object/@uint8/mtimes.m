function c = mtimes(a,b)

% MTIMES   Matrix multiply, overloaded method for UINT8
%


c = uint8( mtimes( double(a) , double(b) ) );
