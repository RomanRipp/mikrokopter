function [str,val,uni] = splithead(str,sep,ref);

% SPLITHEAD  Splits a HeaderString in KeyWord / Value
% 
% [STR,VAL,UNI] = SPLITHEAD( STRING , SEP , REF )
%
% example: STRING =  'STR (UNI{1}) = VAL [UNI{2}]'
%
% default: SEP = '='
%          REF = ['|/\([{'    % Start for PAR or UNI
%                 '|/\)]}']   % End   for PAR or UNI
%
% Check for Parenthesis-Stements in Values like for Units
%
% see also: SPLITUNIT
%
 
Nout = nargout;

val = '';
uni = {'' ''};

if isempty(str)
   return
end

if nargin < 2
   sep = '=';
end

if nargin < 3
   ref = cat( 1 , '|/\([{' , ...
                  '|/\)]}'       );
end

ok = ~isempty(sep);

if ok
   if prod(size(sep)) == 1
      ii = find( str == sep );
   elseif size(str,2) > size(sep,2)
      ii = findstr( str , sep );
   else
      ii = [];
   end
   ok = ~isempty(ii);
   if ok
      ii = ii(1);
      val = str( ii+size(sep,2) : end );
      str = str( 1 : ii-1 );
   end
end

str = rmblank(str,2);

if ~isempty(str) & ~isempty(ref)
    [str,uni{1}] = splitunit(str,ref);
end

if ~isempty(val) & ( Nout >= 2 )
    val = rmblank(val,2);
    if ~isempty(val) & ~isempty(ref)
        [val,uni{2}] = splitunit(val,ref);
    end
end

