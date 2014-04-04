function  [n,c] = colspec(v,x)

% COLSPEC   Returns Matlab's or XRGB color specifier 
%
% Get List of Colors:
%
% [ Name , RGB ] = COLSPEC      List of ColorSpecifier
% 
% [ Name , RGB ] = COLSPEC(1)   List of ColorSpecifier and XRGB-Colors
% 
% Find a Color:
%
%   RGB  = COLSPEC( Name )       Color corresponds with Name, Matlab's only
%   RGB  = COLSPEC( Name , 1 )   Color corresponds with Name, incl. XRGB-ColorNames
%
%   Name = COLSPEC( RGB )        Find Name corresponds with Color, Matlab's only
%
%   Name = COLSPEC( RGB , 1+R*i) Find Name corresponds with Color, incl. XRGB
%
%   R is the ToleranceRadius to search for a XRGB-Color, default: 0.03
%
% see also: XRGB
%

Nin = nargin;

nv = 0;

if Nin == 0
   nv = 1;        % Not V
   v  = [];
   x  = 0;
elseif Nin < 2
   nv = ( isnumeric(v) & ( prod(size(v)) == 1 ) );
   if nv
      x = v;
   else
      x = 0;
   end
end

rad = imag(x);       % ToleranceRadius

x   = isequal(real(x),1);  % Use XRGB

%------------------------------------------------------

n = {   'r'        [ 1  0  0 ]
        'g'        [ 0  1  0 ]
        'b'        [ 0  0  1 ]
        'c'        [ 0  1  1 ]
        'm'        [ 1  0  1 ]
        'y'        [ 1  1  0 ]
        'k'        [ 0  0  0 ]
        'w'        [ 1  1  1 ]
        'red'      [ 1  0  0 ]
        'green'    [ 0  1  0 ]
        'blue'     [ 0  0  1 ]
        'cyan'     [ 0  1  1 ]
        'magenta'  [ 1  0  1 ]
        'yellow'   [ 1  1  0 ]
        'black'    [ 0  0  0 ]
        'white'    [ 1  1  1 ]   };

c = cat(1,n{:,2});
n = n(:,1);

%------------------------------------------------------

if nv
   m = size(n,1) / 2;         % Unique Colors
   c = c(1:m,:);
   n = n(1:m);
   if x
      try
         [nx,cx] = xrgb(1);      % Unique Colors
         n  = cat( 1 , n , nx );
         c  = cat( 1 , c , cx );
      catch
         warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
      end
   end
   return
end

%------------------------------------------------------

s = size(v);
p = prod(s);

%----------------------------------------------
if isnumeric(v) & isequal(s,[1 3])
%----------------------------------------------

   d = c - v(ones(1,size(c,1)),:);
   d = all(d==0,2);

   if  any(d)

       d = find(d);
       d = d(1);
       n = n{d};
       c = c(d,:);

   elseif x

       try

          rad = rad + 0.03 * ( rad == 0 );

          [n,c] = xrgb(1);  % Unique Colors

           m = ones(size(c,1),1);

           d = ( c - v(m,:) );

           d = sum( d.^2 , 2 );

           ok = ( d <= rad^2 );

           if any(ok)
              nk =  sum(ok);
              ok = find(ok);
              if nk > 1
                 [d,ii] = min(d(ok));
                  ok = ok(ii);
              end
              n = n{ok};
              c = c(ok,:);
           else
              n = '';
              c = [];
           end

        catch

           warn(sprintf('Can''t get Color by XRGB.\n%s',lasterr))
           
           n = '';
           c = [];

        end

   else
       n = '';
       c = [];
   end

%----------------------------------------------
elseif ischar(v) & ~isempty(v) & ( p == s(2) )
%----------------------------------------------

   v = lower(v);

   for ii = 1 : size(n,1)

       nn = n{ii}( 1 : min(p,size(n{ii},2)) );

       ok = strcmp( nn , v );      

       if ok
          break
       end

   end


   if  ok
       n = n{ii};
       c = c(ii,:);
   elseif x
       try
          [c,h,n] = xrgb(v);
       catch
          warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
          n = '';
          c = [];
       end
   else
       n = '';
       c = [];
   end

   v = c;
   c = n;
   n = v;

%----------------------------------------------
else
%----------------------------------------------

   error('Input must be an RGB-Tripel or String.');

%----------------------------------------------
end
%----------------------------------------------

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function warn(wrn)

% Display Warning-Message with Beep

ww = warnstat;

if strcmp(ww,'off') | isempty(wrn)
   return
end

warning('on');
  
fprintf(1,'\n%s',char(7));

warning(wrn);

fprintf(1,'\n%s',char(7));

warning(ww);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

