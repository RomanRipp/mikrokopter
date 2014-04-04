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

if ~isempty(cmp)
    if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
       ok = ones(si);
       return
    end
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


