function c = wavemap(clm,scl,n);

% WAVEMAP  ColorMap with Cosine-modified HSV-Developing
%
% ColorMap = WAVEMAP( ColorLimit , ColorScale , ColorNumber )
%
% ColorLimit = [ RGB1 ; RGB2 ] | [ Hue1  Hue2 ]
%             
% ColorScale = [ HueScale SatScale ValScale ]
%
%              Scale HSV linear developing between Colors by Cosine
%
%     Scale = [ +Amplitude + Deviation*i ]  ¯\/¯
%             [ -Amplitude + Deviation*i ]  _/\_
%
%                Amplitude = [   0  .. 1 ]
%                Deviation = [ -Inf .. 1 )
%
% HSV_Scale =  0 | NaN  for no Scale
%

Nin = nargin;

if Nin < 1
   clm = [];
end

if Nin < 2
   scl = [];
end

if Nin < 3
   n = [];
end

%**************************************************************
% Defaults

if isempty(clm)
   clm = [ 3/6  5/6 ];  % Cyan --> Magenta
end

def = [ 1  1  NaN ];

if isempty(scl)
   scl = def;
end

if isempty(n)
   fig = get( 0 , 'currentfigure' );
   if isempty(fig)
      c = get( 0 , 'defaultfigurecolormap' );
   else
      c = get( fig , 'colormap' );
   end
   n = size(c,1);
end

%**************************************************************
% Check Inputs

msg = cell(0,1);

%-------------------------------------------------------------
% ColorNumber

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n > 1 ) & ( mod(n,1) == 0 ) );
end

if ~ok
    msg = cat(1,msg,{'ColorNumber must be a positive Integer larger 1.'});
end

%-------------------------------------------------------------
% ColorLimit

ok = isnumeric(clm);
if ok
   ok = all(isfinite(clm(:)));
end

if ~ok
    msg = cat(1,msg,{'ColorLimit must have finite Numerics.'});
else
    sz = size(clm);
    pz = prod(sz);
    if pz == 2
       clm = cat( 2 , clm(:) , ones(2,2) );
    elseif isequal(sz,[2 3])
       if ~all( abs(clm-0.5) <= 0.5 )
           msg = cat(1,msg,{'RGB-ColorLimits must be in Range [ 0 .. 1 ].'});
       else
           clm = rgb2hsv(clm);
       end
    else
       msg = cat(1,msg,{'Invalid ColorLimit.'});
    end
end

%-------------------------------------------------------------
% Factors

pz = prod(size(scl));

ok = ( isnumeric(scl) & ( pz <= 3 ) );
if ok
   ok = all( isfinite(scl(:)) | isnan(scl(:)) );
end

if ~ok
    msg = cat(1,msg,{'ColorScale must be max. 3-element finite Numerics or NaN.'});
elseif pz < 3
    scl = cat( 2 , permute(scl(:),[2 1]) , def(pz+1:end) );
else
    scl = permute(scl(:),[2 1]);
end

%-------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%**************************************************************
% Build ColorMap

dev = imag(scl);
amp = real(scl);

sgn = sign(amp);
amp = abs(amp);

isn = isnan(amp);
if any(isn)
   isn = find(isn);
   dev(isn) = NaN;
end

amp = min( amp , 1 );


%-------------------------------------------------------------
% BaseMap

%%% clm(2,1) = clm(2,1) + ( clm(2,1) <= clm(1,1) );

x = linspace(0,1,n)';

c = interp1([0 1]',clm,x);

ok = ~( isnan(dev) | ( amp == 0 ) );

if ~any(ok)
   c = hsv2rgb(c);
   return
end

%-------------------------------------------------------------
% Scale Hue / Sat / Val

for ii = find(ok)

    p = 1 / ( 1 - dev(ii) + ( dev(ii) == 1 ) );

    y = x + dev(ii);

    y = ( 1 + sgn(ii) * cos( 2*pi * y/p ) ) / 2;

    y = y / ( y(1) + ( y(1) == 0 ) );

    y = amp(ii) * y  + ( 1 - amp(ii) );

    if ii == 1  % Hue

       y = cumsum(y) - y(1);

       y = y / y(n);

       c(:,ii) = clm(1,1) + y * ( clm(2,1) - clm(1,1) );

    else

      c(:,ii) = y .* c(:,ii);

    end

end

c(:,1) = c(:,1) - floor(c(:,1));

c = hsv2rgb(c);





