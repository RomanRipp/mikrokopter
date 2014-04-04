function [dist,ang1,ang2]=sodano(lon1,lat1,lon2,lat2)

% SODANO Distances on an ellipsoidal body
%
% [Dist,Ang1,Ang2] = SODANO(lon1,lat1,lon2,lat2)
%
% calculates the distance between two positions in m
% for an ellipsoidal body
%
% second position may be array of positions
%

% translated from FORTRAN, IFREMER-Brest / EPSHOM-Brest
% Gerd Krahmann, IfM Kiel, Jun 1993
% last change 14. Jun 1994 G. Krahmann


%
%  values of half large and short axis of the ellipsoid WGS72: dgax72, dpax72
%  new values WGS84: ea84, eb84
%

pi19=pi;3.141592653589793238;
      %  3.141592653589793116

rterre = 6371229.0;
dgax72 = 6378135.00;
dpax72 = 6356750.52;
ea84 = 6378137.000;
eb84 = 6356752.314;
    
degrad = pi19/180.0;
a = 6370.0e03;
b = a;

rlon1 = lon1*degrad;
rlat1 = lat1*degrad;

rlat2 = lat2*degrad;
rlon2 = lon2*degrad;
a = ea84;
b = eb84;

if length(rlon2)>1
  bad=find( (rlon2==rlon1) & (rlat2==rlat1) );
  if ~isempty(bad)
    rlon2(bad)=9000*ones(1,length(bad));
  end
end

%=======================================================================
%
% sodan8 (rlat1,rlon1,rlat2,rlon2,dist,a,b,azab,azba)
%
%=======================================================================
%
%  Version CDC.1.01            28 mars 1989     source EPSHOM
%  -------
%
%  Modifications
%  -------------
%  Correction calcul des azimuths                            T.Terre
%  Correction latitudes identiques              ../01/88     F.Jaulgey
%
%  Objet :  Calcul de la distance geodesique (formule de sodano) et des
%  -----    azimuths.
%
%  Entree:    rlat1, rlon1 : latitude et longitude du point 1 en radians
%  ------     rlat2, rlon2 : idem pour le point 2
%             a, b         : demi-grand axe et demi petit axe de l'ellip-
%                            soide de reference
%=======================================================================
%
flat=1.-b/a;
flat2=flat.*flat;
f1=flat2*1.25;
f2=flat2*0.50;
f3=flat2*0.25;
f4=flat2*0.125;
f5=flat2*0.0625;
f6=flat2+flat;
f7=f6+1.0;
f8=f6*0.5;
beta1  = atan((1.-flat)*sin(rlat1)/cos(rlat1));
sbeta1 = sin(beta1);
cbeta1 = cos(beta1);
beta2  = atan((1.-flat).*sin(rlat2)./cos(rlat2));
sbeta2 = sin(beta2);
cbeta2 = cos(beta2);
dlat = rlat1 - rlat2;
dlon=rlon1-rlon2;
adell=abs(dlon);
if(adell>=pi19)
  adell=2*pi19-adell;
end
sidel=sin(adell);
codel=cos(adell);
a1=sbeta1*sbeta2;
b1=cbeta1*cbeta2;
cophi=a1+b1.*codel;
tmp = (sbeta2*cbeta1-sbeta1*cbeta2.*codel);
tmp0 = (sidel.*cbeta2);
siphi=sqrt(tmp0.*tmp0+tmp.*tmp);
if (siphi==0) 
  dist =0;
else
  c=b1.*sidel./siphi;
  em=1.-c.*c;
  phi=asin(siphi);
  if(cophi<0.)
    phi=pi-phi;
  end
  phisq=phi.*phi;
  csphi=1../siphi;
  ctphi=cophi./siphi;
  psyco=siphi.*cophi;
  term1=f7*phi;
  term2=a1.*(f6*siphi-f2.*phisq.*csphi);
  term3=em.*(f2.*phisq.*ctphi-f8.*(phi+psyco));
  term4=a1.*a1.*f2.*psyco;
  term5=em.*em.*(f5.*(phi+psyco)-f2.*phisq.*ctphi-...
      f4.*psyco.*cophi.*cophi);
  term6=a1.*em.*f2.*(phisq.*csphi+psyco.*cophi);
  dist=b.*(term1+term2+term3-term4+term5+term6);
end
ctaz1=((sbeta2./cbeta2).*cbeta1-sbeta1.*codel);
ctaz2=(sbeta2.*codel-(sbeta1./cbeta1).*cbeta2);
azab=atan2(sidel,ctaz1);
azba=atan2(sidel,ctaz2);
ang1=azab;
ang2=azba;

if length(rlon2)>1
  if ~isempty(bad)
    dist(bad)=zeros(1,length(bad));
    ang1(bad)=zeros(1,length(bad));
    ang2(bad)=zeros(1,length(bad));
  end
end
