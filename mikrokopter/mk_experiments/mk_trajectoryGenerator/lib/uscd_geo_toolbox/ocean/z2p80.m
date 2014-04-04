function p = z2p80(z,l);

% Z2P80 Computes pressure given the depth at some latitude
%
%          P = Z2P80(D,LAT) gives the pressure P (dbars) at a depth D (m)
%
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

if isempty(z)
   p = [];
   return
elseif isempty(l)
   l = 45;
   warning('Use Latitude 45.')
end

      l = 5.92e-3 + 5.25e-3 * sin( abs(l*pi/180) ).^2;
      z = sqrt( (1-l).^2 - 8.84e-6 * z  );
      p = (  1 - l - z ) / 4.42e-6;
