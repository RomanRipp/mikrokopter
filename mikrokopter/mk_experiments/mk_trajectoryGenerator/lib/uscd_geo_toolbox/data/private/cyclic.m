function x = cyclic(x,p,s,d);

% CYCLIC  Creates a monotonic Vector of cyclic Values
%
%   X = CYCLIC( X , Period , Level , DIM )
%
%  Level  Level to check a PeriodOverrunning
%          normalized to Period, optional
%          default: 0.5
%
%  DIM    Dimension to work along, optional
%          default: first nonsingleton Dimension
%
%  Take care, that X contains NO NaN's!
% 
%  Example: use to montonize CompassData
%

Nin = nargin;

%***********************************************
% Check Inputs

msg = cell(0,1);

%-----------------------------------------------
% Basic

if Nin < 2
   error('Inputs X and Period are missing.')
end

if ~isnumeric(x)
    msg = cat(1,msg,{'X must be numeric'});
end

sx = size(x);

%-----------------------------------------------
% Period

if Nin < 1
   msg = cat(1,msg,{'Period is missing.'});
else
   ok = ( ( prod(size(p)) == 1 ) & isnumeric(p) );
   if ok
      ok = ( p > 0 );
   end
   if ~ok
       msg = cat(1,msg,{'Period must be a single positive Numeric.'});
   end
end

%-----------------------------------------------
% Schwelle

if Nin < 3
   s = [];
end

if isempty(s)
   s = 0.5;
else
   ok = ( ( prod(size(s)) == 1 ) & isnumeric(s) );
   if ok
      ok = ( s > 0 );
   end
   if ~ok
       msg = cat(1,msg,{'Level must be a single positive Numeric.'});
   end
end

%-----------------------------------------------
% DIM

if Nin < 4
   d = [];
end

if isempty(d)
   d = sum( cumprod( double( sx <= 1 ) ) );
   d = d + ( d < size(sx,2) );
else
   ok = ( ( prod(size(d)) == 1 ) & isnumeric(d) );
   if ok
      ok = ( ( d > 0 ) & ( mod(d,1) == 0 ) );
   end
   if ~ok
       msg = cat(1,msg,{'DIM must be a single positive Integer.'});
   end
end

%-----------------------------------------------

if ~isempty(msg)
    msg = sprintf('   %s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg);
end

%***********************************************

if isempty(x) | ( sx(d) == 1 )
   return
end

x = x - p * floor(x/p);  % [ 0 .. p )

z    = sx;
z(d) = 1;

if any(isnan(x))

   ii = find(~isnan(x));

   if isempty(ii) 
      return
   end

   dx = cat( d , zeros(z) , diff(x(ii),1,d) );

   dx = ( abs(dx) > s*p ) .* (-sign(dx));

   x(ii) = x(ii) + cumsum( p*dx , d );

else

   dx = cat( d , zeros(z) , diff(x,1,d) );

   dx = ( abs(dx) > s*p ) .* (-sign(dx));

   x = x + cumsum( p*dx , d );

end
