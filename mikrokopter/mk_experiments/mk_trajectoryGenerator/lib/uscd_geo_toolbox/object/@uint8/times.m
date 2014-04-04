function c = times(a,b)

% TIMES   Array multiply, overloaded method for UINT8
%


c = uint8( times( double(a) , double(b) ) );
