function txt = displ(v)

% DISPL  returns String to Display of Variable
%
% String = DISPL( V )
%

txt = '';

if nargin == 0
   return
end

% Format for single double
form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

% Maximum Number of Characters of a String to display
nstr_max = 24;


si = size(v);
nv = si;

if isnumeric(v)
   cv = 'numeric';
else
   cv = class(v);
end

switch cv

  %------------------------------------------------
  case 'numeric'

     if nv == 0

        txt = '[]';

     elseif nv == 1

        txt = sprintf(form,double(v));

     end

  %------------------------------------------------
  case 'char'

     if nv == 0

        txt = '''''';

     elseif chkstr(v,1)  &  ( nv <= nstr_max )

        txt = sprintf('''%s''',v);

     end

  %------------------------------------------------
  case 'cell' 

     if nv == 0

        txt = '{}';

     elseif nv == 1

        txt = sprintf('{%s}',displ(v{1}));

     end

end
% switch

if ~isempty(txt)
   return
end


% Build String from Size and class

ff = [ '[]' ; '{}' ];

ff = ff(1+strcmp(cv,'cell'),:);

txt = sprintf('%.0fx',si);

txt = cat(2,ff(1),txt(1:end-1),' ',cv,ff(2));


%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)

% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

