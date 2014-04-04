function c = bright(c,f,z)

% BRIGHT   Brighten or Darken of ColorValues
%
% C = BRIGHT( C , F , Z )
%
%  C  ColorMatrice with Values between 0 and 1
%     CellArray with String as first Element
%       for using FEVAL(C{:}), returns ColorMatrice
%
%  F  Factor to bright/dark ColorValues linear
%       0 ..  1  brighten
%       0 .. -1  darken
%       0        no change
%
%  Z = CENTER + EXP * i  Defines ExponentialFunction for F 
%
%      The Factor will multiplied by an ExponentialFunction,
%       along the first Dimension of C.
%
%      The Function is ONE at Center and decrease to ZERO with EXP.
% 
%      CENTER and EXP corresponding to first Dimension of C
%       Absolut Values <= 1 will used normalized,
%                       > 1 will used absolut
%        to Length of C in first Dimension  
%


Nin = nargin;

if Nin < 2
   error('Not enough InputArguments.');
end

if ~( isnumeric(f) & ( prod(size(f)) == 1 ) );
    error('Factor must be a single numeric.');
end

if ~( abs(f) <= 1 )
    error('Factor must be between -1 and 1.');
end

ok = isnumeric(c);

if ok
   if isempty(c)
      return
   end
else
   ok = iscell(c);
   if ok
      if isempty(c)
         c = [];
         return
      end
      ok = ( strcmp( class(c{1}) , 'char' )  & ~isempty(c{1})  & ...
            ( prod(size(c{1})) == size(c{1},2) ) );
      if ok
         try
            c = feval(c{:});
         catch
            ok = 0;
         end
      end
   end
end
 
if ~ok
    error('C must be a numeric or CellArray for using FEVAL(C{:}).');
end

if f == 0
   return
end

if any( abs(c(:)-0.5) > 0.5 )
   warning('C should have Values between 0 and 1.');
end


s = sign(f);

if Nin > 2

   if ~( isnumeric(z) & ( prod(size(z)) == 1 ) );
       error('Z must be a single numeric.');
   end

   n = size(c,1);
   m = n - 1;

   v = ( 0 : m )' / m;

   p = abs(imag(z));
   z = real(z);

   if p == 0
      p = 1 / ( pi * exp(1) );
   end

   if abs(z) > 1
      z = ( z - 1 ) / m;
   end
   if ( p > 1 )
      p = ( p - 1 ) / m;
   end

   v = abs( v - z );

   v = exp(-v/p);

   f = f * v;

   nd = ndims(c);
   ind = cell(1,nd);
   ind{1} = ( 1 : n );
   for ii = 2 : nd
       ind{ii} = ones(1,size(c,ii));
   end
  
   f = f(ind{:});

end

if s > 0
   c = ( ( 1 - f ) .* c + f );
else
   c = ( f + 1 ) .* c;
end
