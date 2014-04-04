function  [ h , form ] = epsstr( h );

% EPSSTR  Transform Number exact into String,
%          
% using Matlab's floating point relative accuracy.
%
%  Form = EPSSTR;   
%    returns Format for using with SPRINTF
%
%  [ String , Form ] = EPSSTR( Number ); 
%    returns exact String for Number
%
%  Form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 )
%


form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

if nargin < 1

  h = form;

  return

end


if ~isnumeric(h)  |  ( prod(size(h)) > 1 )
 error('Handle must be a single Numeric.');
end


  h = sprintf(form,h);
