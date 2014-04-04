function [Msg,c,t] = recpar(h,typ);

% RECPAR  returns ParentHistory of Handle
%
% [Msg,HandleHist,TypeHist] = RECPAR( Handle , StopType )
%
%  recurse trough the Parents unto ( ParentType == StopType )
%
%  StopType starts with UpperCase, returns History excluding
%     Handle with ( ParentType == StopType )
%
%  default: StopType = 'Root'  (recurse unto 'figure')
%
%  HandleHist(end) == Handle
%    TypeHist(end) == HandleType
%

Msg = '';
 c  = zeros(0,1);
 t  =  cell(0,1);

if nargin < 1
   Msg = 'Input Handle is missing.';
   return
end

if nargin < 2
   typ = 'Root';
end

%-----------------------------------------------

if isempty(h)
   return
end

ok = ( isnumeric(h) &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   Msg = 'First Input must be a Single Handle.';
   return
end

if ~( ischar(typ) & ~isempty(typ) & ...
      ( prod(size(typ)) == size(typ,2) ) )
   Msg = 'Type must be a String';
   return
end

%-----------------------------------------------

c = h;
t = { get(h,'type') };

z = 1;

t0 = lower(typ);

while ~( c(1) == 0 )  &  ( ~strcmp(t{1},t0) | ( z == 1 ) )

   z = z + 1;

   c = cat( 1 ,         get(c(1),'parent')  , c );

   t = cat( 1 , { lower(get(c(1),'type')) } , t );

end


if strcmp( t{1} , t0 )

  n = 1 + strcmp( typ(1) , upper(typ(1)) );

  c = c(n:z);
  t = t(n:z);

else

   Msg = [ 'Handle has no Parents with Type '''  typ '''.' ]; 

end

%-----------------------------------------------
   
