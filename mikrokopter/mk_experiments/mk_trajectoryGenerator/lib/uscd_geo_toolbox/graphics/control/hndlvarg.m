function [h,v,nr] = hndlvarg(v,n,typ,tag);

% HNDLVARG  Checks Input-varargin for Handle
%
%  [H,V,InputNr] = HANDLVARG( V , n , Type , Tag )
%
%  Search in V{n(1)} .. V{n(2)} for Handle of Type and Tag
%
%  If no Handle found and Type =  'figure' | 'axes', 
%    GCF | GCA is returned 
%
%  see also: CHKHNDL
%

Nin = nargin;

h  = [];
nr = 0;

if Nin == 0
   v = cell(1,0);
   return
end


if ~iscell(v)
   error('Input must be a CellArray.');
end

nv = prod(size(v));

if Nin < 2
   n = [];
end

if Nin < 3
   typ = 'figure';
end

chkin = { typ };

if Nin == 4
   chkin = cat( 2 , chkin , {tag} );
end

%------------------------------------------------------------

if isempty(n)
   n = [ 1  nv ];
else
   n = n(:)';
   n = n( 1 : max(1,size(n,2)) );
   n = min( n , nv );
   if size(n,2) == 1
      n = [ 1  n ];
   end
end


%-------------------------------------------------------------

ok = 0;

for nr = n(1) : n(2)

   ok = chkhndl(v{nr},chkin{:});

   if ok
      break
   end

end

%-------------------------------------------------------------

if ok

  h = v{nr};

  v(nr) = [];   % Reduce v

  return

end

nr = 0;

if ~any( strcmp( typ , { 'figure' 'axes' } ) )
   return
end

h = get( 0 , 'currentfigure' );
if isempty(h)
   return
end

if strcmp( typ , 'axes' )
   h = get( h , 'currentaxes' );
end

ok = chkhndl(h,chkin{:});

if ~ok
   h = [];
end



%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,typ,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%
%  Tag   CharArray or CellStringArray which Tags, the Handle has to be.
%          The Wildcard '*' can be used in the Strings for Tag
%

Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

is_typ = chkstr(typ,0);

if is_typ

  ok = isempty(typ);
  if ~ok
      ok = strcmp( get(h,'type') , typ  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = typ;
  
end


%-------------------------------------------------------------------------
% Check Tag

[ok,tag] = chkcstr(tag,0);

if ~ok
    msg = 'Input Tag must be a CharArray or CellStringArray.';
    return
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = strwcmp(t,tag);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,cmp,wc,nn,cc] = strwcmp(str,cmp,wc,nn,cc)

% STRWCMP   Compare Strings, including WildCards
%
% OK = STRWCMP( STR , CMP , WildCard )
%
% STR: CharacterArray or CellStringArray to compare with Comp
%
% CMP: CharacterArray or CellStringArray with strings to compare with STR
%         strings can contains WildCards
%
% WildCard: specify WildCard to use, default: '*'
%
% OK : logical Array with same size of STR, contains 1 if 
%        any strings of CMP match string of STR
%
% Special Output:  [ok,cmp,wc,nn,cc] = STRWCMP( STR , CMP , [WildCard] )
%
%  to use in follwing statements: ok = STRWCMP( STR , cmp , wc , nn , cc );
%
%  which makes it a bit faster.
%
% see also: STRCMP, FINDSTR
%

Nin  = nargin;
Nout = nargout;

Msg = '';
 nl = char(10);

%***************************************************************
% Check Inputs

%---------------------------------------------------------------
% String

if Nin < 1
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
           'First Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% CompareString

if Nin < 2
   cmp = cell(0,0);
end

if ~( iscell(cmp) & isempty(cmp) )
   [ok,cmp] = chkcstr(cmp,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Second Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% WildCard

if Nin < 3
   wc = '*';
elseif ~( ischar(wc) & ( prod(size(wc)) == 1 ) )
   Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
         'WildCard must be a Single Character.' ];
end
  
%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%***************************************************************

si = size(str);

if ( isempty(str) | isempty(cmp) ) & ( Nout <= 1 ) 
   ok = zeros(si);
   return
end

cmp = cmp(:);

if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
   ok = ones(si);
   return
end

%***************************************************************
% Analyze CompareStrings

nc = size(cmp,1);

ok = ( Nin == 5 );

if ok
   ok = ( isequal(size(cc),[nc 1]) & iscell(cc) & ...
          isequal(size(nn),[nc 2]) & isnumeric(nn)    );
   if ok
      try
         ok = cat( 1 , cc{:} ); 
         ok = ( size(ok,2) == 3 );
      catch
         ok = 0;
      end
   end
end


%--------------------------------------------------
if ~ok
%--------------------------------------------------

  cc    = cell(nc,1);
  cc(:) = { zeros(0,3) };  % { [ Start End N  ] }

  nn = zeros(nc,2);        %   [ Ncmp  sum(N) ] 

  for ii = 1 : nc

    if ~isempty(cmp{ii})

       iwc = ( double(cmp{ii}) == double(wc) );

       if ~any( iwc )

          nn(ii,:) = size(cmp{ii},2);
          cc{ii}   = [ 1  nn(ii,:) ];

       else

         %--------------------------------------------
         % Remove Duplicate WildCards

         iwc = find( iwc );
         if ~( prod(size(iwc)) == 1 )
            jj = find( diff(iwc) == 1 );
            if ~isempty(jj)
               cmp{ii}(iwc(jj+1)) = [];
            end
         end

         %--------------------------------------------
         % Get Start End
     
         iwc = ( double(cmp{ii}) == double(wc) );
  
          n  = size(iwc,2);

          if ( n == 1 ) & ( iwc == 1 ) & ( Nout <= 1 ) % Wildcard only
             ok = ones(si);
             return
          end

          i0 = ~iwc(1);
          i1 = ~iwc(n);

         iwc = cat( 2 , ones(1,i0) , iwc , ones(1,i1) );

         iwc = find( iwc );

         iwc = iwc(:);

           n = size(iwc,1) - 1;

         cc{ii} = zeros(n,3);

         if n > 0      

            cc{ii}(:,[1 2]) = cat( 2 , iwc(1:n)+1 , iwc((1:n)+1)-1 ) - i0;

            cc{ii}(:,3) = cc{ii}(:,2) - cc{ii}(:,1) + 1;
 
         end

         nn(ii,:) = cat( 2 , size(cmp{ii},2) , sum(cc{ii}(:,3),1) );

       end

    end

  end

%--------------------------------------------------
end
%--------------------------------------------------

if ( Nout > 1 )

  if ( isempty(str) | isempty(cmp) )
     ok = zeros(si);
     return
  end

  if any(strcmp(cmp,wc))  % Wildcard only
     ok = ones(si);
     return
  end

end

%***************************************************************
% Compare

ok = zeros(si);

for ii = 1 : prod(si)
 
    s2 = size(str{ii},2);

    for jj = 1 : nc

        ok(ii) = ( ( s2 == 0 ) & ( nn(jj,1) == 0 ) );
       
        if ok(ii)
           break
        end
       
        if ( s2 >= nn(jj,2) ) & ~( nn(jj,1) == 0 )

           ok(ii) = compare(str{ii},cmp{jj},cc{jj});

           if ok(ii)
              break
           end
            
        end

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = compare( str , cmp , cc )

sc = size(cmp,2);

ok = 1;

for ii = 1 : size(cc,1) 

    s2 = size(str,2);

    ok = ( ok &  ( s2 >= cc(ii,3) ) );

    if ~ok
       break
    end

    ic  = ( cc(ii,1) : cc(ii,2) );

    i01 = ( cc(ii,[1 2]) == [ 1  sc ] );

    i0  = ( 1 : cc(ii,3) );
    i1  = s2 - cc(ii,3) + i0;

    i2  = cat( 2 , findstr( str , cmp(ic) ) , 0 );
    i2 = i2(1);
    i2 = ( i2 + cc(ii,3) - 1 ) * ~( i2 == 0 );

    i2  = s2       * (  i01(1) &  i01(2) & strcmp(str    ,cmp    ) ) + ...
          cc(ii,3) * (  i01(1) & ~i01(2) & strcmp(str(i0),cmp(ic)) ) + ...
          s2       * ( ~i01(1) &  i01(2) & strcmp(str(i1),cmp(ic)) ) + ...
          i2       * ( ~i01(1) & ~i01(2) );

    ok = ( ok & ~( i2 == 0 ) );

    if ~ok
       break
    end

    str = str( i2+1 : s2 );

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   str = cellstr(str);
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

