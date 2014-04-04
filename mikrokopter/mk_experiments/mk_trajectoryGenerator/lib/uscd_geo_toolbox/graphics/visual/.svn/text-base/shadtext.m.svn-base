function h = shadow(varargin)

% SHADTEXT  Draws Text with Shadow
%
% SHADTEXT( ... , Property , Value , ... )
%
% SHADTEXT returns the created TextHandles
%
% For Syntax and Property-Value-Pairs see TEXT.
% 
% Special Properties for SHADTEXT, beginning with "@":
%
% '@Width'  , WDT           Width of Shadow in [points]
%                            use a NonZero imaginary Part
%                            for FontSize-Normalized Units
%
% '@Number' , NR            Number of surrounding ShadowText's
%                            NR+1 TextHandles will created
%
% '@Shift'  , [ RAD  DIR ]  Shift Shadow by Radius [points] into
%                            Direction [deg], use a NonZero imaginary
%                            Part in RAD for Width-Normalized Units
%
% '@Color'  , [ R  G  B ]   ShadowColor
%              ColorSpec
%                 
%
% defaults:  WDT   = 0.1+i     (1/10 of FontSize)
%            NR    = 120
%            Shift = [ 0  0 ]
%            Color = inverted TextColor
%

wdt = 0.1+i;
nr  = 120;
col = [];
dev = [0 0];

Nin = nargin;
vin = varargin;

Nout = nargout;

prp = { '@width' '@number' '@shift' '@color' };

%*****************************************************************
% Check Inputs for Parameter

%--------------------------------------------------
if Nin > 1
%--------------------------------------------------

msg = cell(0,1);

ok = zeros(1,Nin);

for ii = 1 : Nin-1

    p = vin{ii};
    s = size(p);
 
  %----------------------------------------------------------
  if ~ok(ii)
  %----------------------------------------------------------

    ok(ii) = ( ischar(p) & ( prod(s) == s(2) ) & ( s(2) >= 2 ) );
    if ok(ii)
        p = lower(p);
       ok(ii) = any(strcmp(prp,p));
    end

    if ok(ii)
      
       v  = vin{ii+1};
       pv = prod(size(v));

       switch p

         %------------------------------------------------------
         case '@width'

           ok(ii) = ( isnumeric(v) & ( pv == 1 ) );

           if ok(ii)
              ok(ii) = isfinite(v); 
           end

           if ok(ii)
              wdt = v;
           else
              msg = cat(1,msg,{'Width must be numeric.'});
           end

         %------------------------------------------------------
         case '@number'

           ok(ii) = ( isnumeric(v) & ( pv == 1 ) );

           if ok(ii)
              ok(ii) = ( isfinite(v) & ( mod(v,1) == 0 ) ); 
           end

           if ok(ii)
              nr = max(v,0);
           else
              msg = cat(1,msg,{'Number must be a Integer.'});
           end
 
         %------------------------------------------------------
         case '@shift'

           ok(ii) = ( isnumeric(v) & ( pv == 2 ) );

           if ok(ii)
              ok(ii) = all(isfinite(v));
           end

           if ok(ii)
              dev = v;
           else
               msg = cat(1,msg,{'Deviation must have 2 finite Elements: [ Rad Dir ].'});
           end
   
         %------------------------------------------------------
         case '@color'

           ok(ii) = ( ischar(v) & ~isempty(v) );

           if ok(ii)
              v = lower(v(1));
              ok(ii) = any( v == 'rygcbmkw' );
           else
              ok(ii) = ( isnumeric(v) & ( pv == 3 ) );
              if ok(ii)
                 v = v(:)';
                 ok(ii) = all( ( 0 <= v ) & ( v <= 1 ) );
              end
           end

           if ok(ii)
              col = v;
           else
              msg = cat(1,msg,{'Color must be a ColorString or RGB-Tripel.'});
           end

       end

       ok(ii+1) = ok(ii);

    end

  %----------------------------------------------------------
  end  % ~ok(ii)
  %----------------------------------------------------------

end

if ~isempty(msg)
    error(sprintf('%s\n',msg));
end

vin(find(ok)) = [];

%--------------------------------------------------
end
%--------------------------------------------------

%*****************************************************************
% Create First Text

h = NaN * zeros(nr+1,1);

try
   h(1) = text(vin{:});
catch
   error(sprintf('Error call TEXT.\n%s',lasterr));
end

if nr == 0
   return
end

%*****************************************************************
% Create Surrounding Text's

axe = get(h(1),'parent');

uni = get(h(1),'units');
      set(h(1),'units','points');

for ii = 1 : nr
    h(ii+1) = copyobj(h(1),axe);
end

%*****************************************************************
% Parameter

if isempty(col)
   col = 1 - get(h(1),'color');
end

if ~( imag(wdt) == 0 )
   wdt = real(wdt) * get(h(1),'fontsize');
end

rad = dev(1);
dir = dev(2) * pi/180;

if ~( imag(rad) == 0 )
    rad = real(rad) * wdt;
end

phi = 2*pi * ( 0 : nr )' / (nr+1);

%*****************************************************************

pos = get(h(1),'position');

pos = pos(ones(1,nr+1),:);

pos(:,1) = pos(:,1) + wdt * cos(phi) + rad * cos(dir);
pos(:,2) = pos(:,2) + wdt * sin(phi) + rad * sin(dir);

for ii = 1 : nr
    set(h(ii),'position',pos(ii,:),'color',col);
end

set(h,'units',uni);

h = h( nr+1 : -1 : 1 );

if Nout == 0
   clear h
end
