function [msg,id] = chkident(id,dlm,fid)

% CHKIDENT  Checks for valid Identifer
%
% [Msg,IDENT] = CHKIDENT( IDENT , Delimiter , FileIdent )
%
% IDENT is the identifer for DataFiles, a String or CellArray of Strings.
% Each String can be a individual FileName of a DataFile or 
%  a coded Identifer with the Separator by DBGET/DBSET,
%  default is the Colon-Character ":".
%
%   Coded Ident: [Base]:[Name]:[Type]:[Numbers]
%                [Base]:[Name]:[Type]:[Pattern]
% 
% Base, Name, Type and Numbers or Pattern can be Elements of a CellArray,
%  which can have up to 4 Columns. Use this form for numeric Values of Numbers.
%
%  Example:  { '[Base]:[Name]:[Type]'  [Numbers] }       2 Columns
%
%            { '[Base]' '[Name]' '[Type]'  [Numbers] }   4 Columns
%

msg = '';


Nin = nargin;

if Nin < 2
   dlm = dbget('Separator');
end

if Nin < 3
   fid = {filesep '.'} % Multiple allowed  !!!
end

%---------------------------------------------------------

s = size(id); 
p = prod(s); 
v = ( max(s) == p );  % True for Vector

%---------------------------------------------------------

[ok,str] = chkcstr(id);

if ~( ( ok | iscell(id) ) &  ( s(1) * s(2) == p ) )
    msg = 'Identifer must be a 2D-CellArray.';
    return
end

%---------------------------------------------------------
% Check for individual IdentCellArray

chk = ( iscell(id) & any( s(2) == [ 2  3  4 ] ) );

if ~( ok | chk )
    msg = 'Identifer must be a CharacterArray or CellArray of Strings.';
    return
end

%---------------------------------------------------------

if ok
   id = str;
   % Vector and evt. separated ==> Check individual Strings
   if  chk & v
       ok = any( cat(2,id{:}) == dlm );
       if ~ok
           % Check for Files only
           ok = 1;
           for f = id(:)'
               k = 0;
               for c = fid
                   k = ( k | any( f{1} == c{1} ) );
               end
               ok = ( ok & k );
           end
       end
   elseif ~chk
       id = str;
   end    
end

if ok
   % Append a Column of empty Elements 
   id = id(:);
   id = id(:,[1 1]);
   id(:,2) = {[]};    % !!! Empty numeric, NOT String {''} !!!
   return
end

%---------------------------------------------------------
% Check Elements for Empty, Strings, Numeric or Invalid (NaN)

ini = zeros(s);

for ii = 1 : s(1)
    for jj = 1 : s(2)
        v = id{ii,jj};
        ok = chkstr(v,1);                     % NonEmpty String
        if ~ok & isnumeric(v) & ~isempty(v)   % NonEmpty Numeric
            v = v(:);
            ok = 2 * all( isfinite(v) & ( v >= 0 ) & ( mod(v,1) == 0 ) );
        end
        ok = double(ok);
        if ( ok == 0 )
           if ~isempty(v)
               ok = NaN;
           else
               id{ii,jj} = '';  % Empty String !!!
           end
        end
        ini(ii,jj) = ok;
    end
end

%---------------------------------------------------------
% Find bad Rows

i0 = ( ini == 0 );  % Empty Elements
i1 = ( ini == 1 );  % True for String 
i2 = ( ini == 2 );  % True for Numeric

j2 = double(~i2);               % True for NonNumeric
j2 = ~( cumprod(j2,2) | i2 );   % True for Elements following Numerics in Row
j2 = any( j2 & i1 , 2 );        % True for Strings  following Numerics in Row

s2 = sum( i2 , 2 );  % Number of Numeric Elements in Row

i1 = any( i1 , 2 );  % True for String  in Row
i2 = any( i2 , 2 );  % True for Numeric in Row

i0 = all( i0 , 2 );  % True for Empty Row


bd = cat( 2 , any(isnan(ini),2) , ...        % 1: Invalid Elements
              ( i2 & ~i1 )      , ...        % 2: Numeric only
              ( s2 >   1  )     , ...        % 3: Multiple Numerics
                j2              , ...        % 5: Strings follows Numerics
                i0                    );     % 6: Only Empty Elements

msg = cell(0,1);

if any( bd(:,1) )
   msg = cat(1,msg,{'CellElements of Identifer must be a String or positive Integers.'});
end

if any( bd(:,2) )
   msg = cat(1,msg,{'A single IdentiferRow must have a StringElement.'});
end

if any( bd(:,3) )
   msg = cat(1,msg,{'A single IdentiferRow can have only one numeric Element.'});
end

if any( bd(:,4) )
   msg = cat(1,msg,{'No Strings can follow a numeric Element.'});
end


if ~isempty(msg) 
    msg = sprintf('%s\n',msg{:});
    msg = msg(1:(end-1));
    return
end

msg = '';

%---------------------------------------------------------
% Remove all bad Rows

ok = ~any(bd,2);

if ~all(ok)
    if ~any(ok)
        id = cell(0,2);
        return
    end
    ok  = find(ok);
    id  =   id(ok,:);
    ini =  ini(ok,:);
end

%---------------------------------------------------------
% Concatinate Strings with Separator (Delimiter)
%  Let Numbers unmodified


if 0
   % Numerics ---> Strings
   i2 = ( ini == 2 );  % True for Numerics
   if any(i2(:))
      i2 = find(i2);
      for ii = i2(:)'
           id{ii} = cat(2,'[',sprintf(' %.0f ',id{ii}),']');
          ini(ii) = 1;
      end
   end
end

i1 = ( ( ini == 1 ) | ( ini == 0 ) );  % True for String 

j1 = double(~i1);               % True for NonString
j1 = ( cumprod(j1,2) | i1 );    % True for Elements before and incl. Strings

n1 = sum( j1 , 2 );             % Number of Elements

ok = ( n1 >= 1 );

if ~any(ok)
    id = cell(0,2);
    return
end

if ~all(ok)
    ok = find(ok);
    id = id(ok,:);
    n1 = n1(ok);
end

if ~all( n1 == 1 )
    i1 = find( n1 > 1 );
    for ii = i1(:)'
        v = id([ii ii],1:n1(ii));
        v(2,:) = {dlm};
        v(2,n1(ii)) = {''};
       id{ii,1} = sprintf('%s',v{:});
    end
end

% Numerics ---> 2. Column

i2 = double( ~( ini == 2 ) );   % Not True  for Numerics
i2 = sum(cumprod(i2,2),2) + 1;  % Column of Numeric
i2 = i2 .* ( i2 <= s(2) );

nr    = cell(size(id,1),1);
nr(:) = {[]};

if any(i2)

   i1 = find(i2);
   i2 = i2(i1);

   ii = i1 + (i2-1)*s(1);

   nr(i1) = id(ii);

end

id = cat( 2 , id(:,1) , nr );
