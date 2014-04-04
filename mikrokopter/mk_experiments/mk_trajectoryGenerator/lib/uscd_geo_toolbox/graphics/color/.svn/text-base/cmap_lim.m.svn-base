function c = cmap_lim(fcn,nc,lim)

% CMAP_LIM  returns bright-limitated Colormap
%
% CMAP_LIM( ColorMap    , ColorNumber , Limit )
% CMAP_LIM( ColorMapFcn , ColorNumber , Limit )
%
% Options:
%
%  imag(ColorNumber) ~= 0  ==>  Invert ColorMap
%     
%  real(ColorNumber) >  0  ==>   Bright --> Dark  
%                                  Blue --> Red  
%
%  real(ColorNumber) <  0  ==>     Dark --> Bright  |
%                                   Red --> Blue
%                

if nargin < 1
   fcn = 'hot';
end

if nargin < 2
   nc = 64;
end

if nargin < 3
   lim = [ 0.5  2.5 ];
end

%------------------------------------------------

if real(nc) == 0
  c = zeros(0,3);
  return
end

invert = ~( imag(nc) == 0 );
flip   =  ( real(nc) <  0 );

nc = abs(real(nc));

%------------------------------------------------

ok = ( ischar(fcn) & ~isempty(fcn) & ...
       ( prod(size(fcn)) == size(fcn,2)   ) ); 

if ok

  try  
    c = feval(fcn,2*max(nc,256));
  catch
    ok = 0;
  end

else

  ok = ( isnumeric(fcn) & ~isempty(fcn) & ...
         ( prod(size(fcn)) == size(fcn,1)*size(fcn,2) ) & ...
         ( size(fcn,2) == 3 ) );

  if ok
     ok = ( all(isfinite(fcn)) & all( fcn(:) >= 0 ) & all( fcn(:) <= 1 ) );
  end

  if ok
     c = fcn;
  end

end

if ~ok
   c = hot(2*max(nc,256));
end

%-----------------------------------------------------

c = 1*invert + (1-2*invert) * c ;

%-----------------------------------------------------
   
cs = sum( c , 2 );
ok = find( ( lim(1) <= cs )  &  ( cs <= lim(2) ) );


if prod(size(ok)) > 1
  c = c(ok,:);
end

ncm = size(c,1);

if ncm == 1

  c = c( ones(nc,1) , : );

else

  c = interp1( (1:ncm)' , c , linspace(1,ncm,nc)' );

end

f = 1 - 2*flip;

if f*sum(c(1,:),2) < f*sum(c(nc,:),2)
% Last == Darkest

   c = c(nc:-1:1,:);

else
% Last == Nearest to Red

  h = rgb2hsv(c);
  h = abs(h(:,1)-0.5);

  if f*h(1,1) > f*h(nc,1)

   c = c(nc:-1:1,:);

  end

end


