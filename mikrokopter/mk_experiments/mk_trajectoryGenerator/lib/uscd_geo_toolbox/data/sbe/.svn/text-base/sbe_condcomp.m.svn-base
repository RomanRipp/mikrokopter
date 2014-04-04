function c = sbe_condcomp(c,pref,p,t);

% SBE_CONDCOMP Compressibility Compensation of Sea-Bird Conductivity Sensors
%
% Correct conductivity of Sea-Bird-Instruments (MicroCat/SeaCat) 
%  without build-in pressure sensor.
%
% C = SBE_CONDCOMP( Cond , Pref , Pres )
%
% C = Cond * ( 1 + CPCOR  * PREF ) / 
%            ( 1 + CPCOR  * Pres )
%
% C = SBE_CONDCOMP( Cond , PREF , Pres , Temp )
%
% C = Cond * ( 1 + CPCOR  * PREF + CTCOR * Temp ) / 
%            ( 1 + CPCOR  * Pres + CTCOR * Temp )
%
%
% defaults:  CPCOR = -9.57e-08
%            CTCOR =  3.25e-06
%
% use a 3-Element Vector as second Input to set CPCOR and CTCOR:
%
% C = SBE_CONDCOMP( C , [ PREF CPCOR CTCOR ]  , ... )
%

Nin = nargin;

if Nin < 2
   error('Not enough Input Arguments.')
end

cp = -9.57e-08;
ct =  3.25e-06; 

%*******************************************************************
% Check Inputs

msg = {};

%-------------------------------------------------------------------
% Pref

sz = size(pref); pz = prod(sz);

ok = ( isnumeric(pref) & ( 1 <= pz ) & ( pz <= 3 ) );
if ok
   ok = all(isfinite(pref));
end

if ~ok
    msg = cat(1,msg,{'Second Input must be [ Pref CPCOR CTCOR ].'}); 
elseif pz > 1
    cp = pref(2);
    if pz > 2
       ct = pref(3);
    end
    pref = pref(1);
end

%-------------------------------------------------------------------
% Pres and Temp

sc = size(c); pc = prod(sc);

if Nin > 2

   sz = size(p); pz = prod(sz);
   if ~( isequal(sc,sz) | ( pc == 1 ) | ( pz == 1 ) )
         msg = cat(1,msg,{'Size of Cond and Pres must be agree.'}); 
   end

   if ~isnumeric(p)
         msg = cat(1,msg,{'Pres must be numeric.'}); 
   end

end

if Nin > 3

   sz = size(t); pz = prod(sz);
   if ~( isequal(sc,sz) | ( pc == 1 ) | ( pz == 1 ) )
         msg = cat(1,msg,{'Size of Cond and Temp must be agree.'}); 
   end

   if ~isnumeric(t)
         msg = cat(1,msg,{'Temp must be numeric.'}); 
   end

   si = size(p); pi = prod(si);
   if ~( isequal(si,sz) | ( pi == 1 ) | ( pz == 1 ) )
         msg = cat(1,msg,{'Size of Pres and Temp must be agree.'}); 
   end

end

%-------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg);
end

%*******************************************************************

switch Nin

  case 3   % SBE_CONDCOMP( Cond , Pref , Pres )

     c = c .* ( 1 + cp*pref ) ./ ( 1 + cp*p );

  case 4   % SBE_CONDCOMP( Cond , PREF , Pres , Temp )

     c = c .* ( 1 + cp*pref + ct*t ) ./ ( 1 + cp*p + ct*t );

end
