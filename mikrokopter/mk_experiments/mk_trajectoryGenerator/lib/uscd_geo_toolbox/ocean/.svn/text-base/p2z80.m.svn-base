function z=p2z80(p,l);

% P2Z80 Computes depth given the pressure at some latitude
%          Z=Z2P80(P,LAT) gives the depth D (m) at a pressure P (dbars)
%          at some latitude LAT (degrees).
%
%          This probably works best in mid-latitude oceans, if anywhere!
%
%          Ref: Saunders, "Practical Conversion of Pressure to Depth",
%              J. Phys. Oceanog., April 1981.
% 

%         I copied this directly from the UNESCO algorithms.


% CHECK VALUE: P80=7500.004 DBARS;FOR LAT=30 DEG., DEPTH=7321.45 METERS

if nargin < 2
   l = [];
end

if isempty(p)
   z = [];
   return
elseif isempty(l)
   l = 45;
   warning('Use Latitude 45.')
end

      l = 5.92e-3 + 5.25e-3 * sin( abs(l*pi/180) ).^2;
      p = ( 1 - l  - p*4.42e-6 ) .^ 2;
      z = ( (1-l).^2 - p  ) / 8.84e-6;
