function ost=oxsat(s,t)
% OXSAT Oxygen saturation.
%
%  OST = OXSAT(S,T) return the oxygen saturation in [ml/l] 
%  as a function of salinity S (PSS-78) and in situ temperature T (øC).
%
%  Conversion between mg and ml:  OST[ml/l] = OST[mg/l] * 0.7 [ml/mg]
%
%  References:
%  Weiss, R. F., The solubility of nitrogen, oxygen and argon in water
%    and seawater. Deep-Sea Res., 17, 721-735, 1970.

%  Christian Mertens, IfM Kiel
%  $Revision: 1.0 $ $Date: 1996/01/05 13:20:22 $


a = [-173.4292 249.6339 143.3483 -21.8492];
b = [-0.033096 0.014259 -0.0017]; 

t = 0.01*(t + 273.15);
ost = a(1) + a(2)./t + a(3)*log(t) + a(4)*t + s.*(b(1) + (b(2) + b(3)*t).*t);
ost = exp(ost);

