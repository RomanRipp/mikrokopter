function [ok,varargout] = mkgrid(varargin);

% MKGRID  Generation of common Grid from Inputs
%
% [ Ok, Y1 , Y2 , Y3 , ... ] = MKGRID( X1 , X2 , X3 , ... )
%
%   Ok =  0   if MatrixDimensions are not agree
%             [ Y1 , Y2 , ... ] == [ X1 , X2 , ... ]
%
%         1   if MatrixDimensions are agree
%             [ Y1 , Y2 , ... ] are Matrices with equal Size
%
%        NaN  if missing Input
%
%
%   MatriceInputs must have an equal Size.
% 
%   VectorInputs are checked successive if they match the Length
%    of the previous Vectors and all MatriceInputs, 
%     or the Length in the DimensionOrder of InputOrder. 
%
%
% MKGRID( '#<MemoryLimit>' , ... )  sets a Limit 
%
%   of the Number of OutputElements to: MemoryLimit * 2^20 / 8
%
%   The Unit of MemoryLimit is [Mbytes] in case of class "double".
%
%   In case of exeeding MemoryLimit, VectorInputs returns reshaped
%    in correct DimensionOrder but not gridded !
%
%   default: MemoryLimit = 256 [Mbytes] ==  2 x [ 4096 x 4096 ] double
%                                       ==  8 x [ 2048 x 2048 ] double
%                                       == 32 x [ 1024 x 1024 ] double
%
%   example: MKGRID( '#512' , ... )   % sets the MemoryLimit to 512 Mbyte
%
%            The Maximum Number of Elements in Output is:  512 * 2^20 / 8
%
% 
% See also: NDGRID, MESHGRID
%

%-------------------------------------------------

mb = 256;     % default MemoryLimit [Mbytes] !!!!!!

bp = char(7); % WarningBeep

%-------------------------------------------------

 in =  nargin;
out = nargout - 1;

out = max(0,out);

varargout = cell(1,out);

%-------------------------------------------------
% Check for 1. MemoryInput

if in > 0

   v = varargin{1};
   s = size(v);
   p = prod(s);

   ok = ( ischar(v) & ( s(2) == p ) & ( s(2) > 1 ) );
   if ok
      ok = ( v(1) == '#' );
   end

   if ok

      try 
         m = eval(v(2:s(2)));
      catch
        ok = 0;
         m = 0;
      end

      wrn = '';

      if ~ok
          wrn = sprintf('Cann''t evaluate Expression for MemoryLimit "%s".',v);
      else
         ok = ( isnumeric(m) & ( prod(size(m)) == 1 ) );
         if ok
            ok = ( m > 0 );
         end
         if ~ok
             wrn = sprintf('MemoryLimit in "%s" must be a single numeric larger ZERO.',v);
         end
      end

      if ~ok
          ww = warnstat;
          if ~strcmp(ww,'off')
              warning('on')
              warning([ bp  wrn ]);
              warning(ww);
          end
      else
          varargin = varargin(2:in);
          in = in - 1;
          mb = m;
      end

   end

end

%-------------------------------------------------
         
if in < 1
   ok = NaN;
   return
end

out = min(in,out);

%-------------------------------------------------
% Check NDim

d  = zeros(in,1);

for ii = 1 : in
    d(ii) = ndims(varargin{ii});
end

nd = max(d);  % NDIM

%-------------------------------------------------
% Get Sizes

sz = ones(in,nd);

for ii = 1 : in
    sz(ii,1:d(ii)) = size(varargin{ii});
end

%-------------------------------------------------
% Check all Equal 

ok = sz - sz(ones(1,in),:);
ok = ( ok == 0 );
ok = all(ok(:));

if ok
   for ii = 1 : out
       varargout{ii} = varargin{ii};
   end
   return
end

%-------------------------------------------------
% Check for Vectors

mv = max(sz,[],2);         % Maximum per Matrice

v = ( mv == prod(sz,2) );  % True for Vector

isv = any(v);

if isv

   vd = double( mv(:,ones(1,nd)) == sz );

   vd = sum(cumprod(1-vd,2),2) + 1;

   vd = vd .* v;            % Orig VectorDimension

   wd = vd;                 % New  VectorDimension

   %-----------------------------------
   % Sort VectorDimension by InputOrder

   sv         = ones(in,max(in,nd));
   sv(:,1:nd) = sz;

      iv    = find(v);
   sv(iv,:) = 1;

   sv(iv+(iv-1)*in) = mv(iv);

else

   iv = [];

   vd = zeros(in,1);
   wd = zeros(in,1);
   sv = sz;

end

%-------------------------------------------------
% Check Length of Dimensions

ok = chkdimsize(sz);

%-------------------------------------------------
% Check Vectors successive

if ~ok & isv

    %-------------------------------------------------
    % Loop over Vectors in InputOrder
    % Check with previous Vectors and all Matrices

    for ii = iv(:)'

        w       = zeros(in,1);
        w(1:ii) = v(1:ii);      % Vector until II
        
        iw = ( ~v | w );        % Vector until II or Matrice
        nw =  sum(iw);
        iw = find(iw);

        %------------------------------------------------------
        % Check Size Including Vector

        ok = chkdimsize(sz(iw,:));

        if ~ok & ( nw > 1 )

            %------------------------------------------------------
            % Check Length of Vector with Size exluding Vector

            iw = iw( 1 : (nw-1) );

            d  = max(sz(iw,:),[],1);  % Maximum per Dimension

            jj = ( d == mv(ii) );     % Length of Vector match Maximum

            ok = any(jj);

            if ok

                     wd(ii)  = sum(cumprod(~jj,2),2) + 1;
               sz(ii,  :   ) = ones(1,nd);
               sz(ii,wd(ii)) = mv(ii);

            elseif nd < ii                     % Expand NDIM

               sz = cat(2,sz,ones(in,ii-nd));
               nd = ii;

                     wd(ii)  = ii;
               sz(ii,  :   ) = ones(1,nd);
               sz(ii,wd(ii)) = mv(ii);

            end

        end

    end

    %-------------------------------------------------
    % Check Length of Dimensions again

    ok = chkdimsize(sz);

    %-------------------------------------------------
    % Check Length of Dimensions, 
    %  using Input-Order sorted VectorDimensions

    if ~ok

        ok = chkdimsize(sv);

        if ok
           wd(iv) = iv;
               sz = sv;
               nd = size(sv,2);
        end

    end
    %-------------------------------------------------
 
end

%-------------------------------------------------
% Check for Changed VectorDimensions

if ok & ( out < in ) & isv

   jj = ( (out+1) : in );

   if ~all( vd(jj) == wd(jj) ) 
       ww = warnstat;
       if ~strcmp(ww,'off')
           warning('on')
           warning([ bp  'VectorDimension not agree.' ]);
           warning(ww);
       end
   end

end

%-------------------------------------------------
% Check for Return

if ~ok | ( out == 0 )
    for ii = 1 : out
        varargout{ii} = varargin{ii};
    end
    return
end

%-------------------------------------------------
% Check for VectorPermutation

jj = ( 1 : out );

jj = ~( vd(jj) == wd(jj) );

if any(jj)

   jj = find(jj);

   for ii = jj(:)'

       varargin{ii} = varargin{ii}(:);

       prm = ( 1 : nd ) - wd(ii) + 1;  % 1 at wd

       prm = prm - nd * ( ceil(prm/nd) - 1 );
 
       varargin{ii} = permute(varargin{ii},prm);

   end

end

%-------------------------------------------------
% Check Memory

d = max(sz,[],1);

m = out * prod(d,2) * 8 / 2^20;

if m > mb

   ww = warnstat;
   if ~strcmp(ww,'off')
       warning('on')
       warning(sprintf('%sMemoryLimit of %.3g Mbytes exeeded.',bp,mb));
       warning(ww);
   end

    for ii = 1 : out
        varargout{ii} = varargin{ii};
    end

    return

end

%-------------------------------------------------
% Build Matrices

o = cell(1,nd);
k = cell(1,nd);

for ii = 1 : nd
    o{ii} = ones(1,d(ii));
    k{ii} = ( 1 : d(ii) );
end

for ii = 1 : out

     q = k;

    jj = ( sz(ii,:) == 1 );

    if any(jj)
           jj  = find(jj);
         q(jj) = o(jj);
    end

    varargout{ii} = varargin{ii}(q{:});

end


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkdimsize(s);

% Check if all DimLength's are to equal Maximum or ONE
%

d = max(s,[],1);                   % Maximum per Dimension

ok = s - d(ones(1,size(s,1)),:);

ok = ( ( ok == 0 ) | ( s == 1 ) );

ok = all(ok(:));


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

