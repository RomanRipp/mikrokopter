function [c,t] = bbody(n,scl)

% BBODY BlackBody Color Spectra
%
% [ RGB , T ]  = BBODY( N )
%
% Returns the RGB-Values and Temperature for 
%  an N-point BB-ColorSpectra between 500K <= T <= 12500K.
%
% [ RGB , T ] = BBODY( [ T0 T1 ] ); with  N = (T1-T0)/100
% [ RGB , T ] = BBODY( [ T0 T1 N ] )
%
% Returns the RGB-Values and Temperature for 
%  an N-point BB-ColorSpectra between T0 <= T <= T1.
%
% [ RGB , T ] = BBODY( T )
%
%  Returns the RGB-Values for the Temperatures defined by T,
%    use BBODY(-W) for a single TemperatureValue, 
%    use a Column-Vector for two or three elements in T.
%
% A second Input defines the Gamma-Adjustment:
%
%   SPECTRA( ... , Gamma ), default: Gamma = 1.0
%
%--------------------------------------------------------
%
% This approximation is based on data provided 
%   by Mitchell Charity
%
%  http://www.vendian.org/mncharity/dir3/blackbody
%
%--------------------------------------------------------
%
% see also: SPECTRA
%

if nargin < 1
   n = [];
end

lim = [ 500  12500 ];

pot = 1.0;  % GAMMA ADJUST 

t = [];

p = prod(size(n));

if     p == 2
   lim = n;
     n = [];
elseif p == 3
   lim = n([1 2]);
     n = n(3);
elseif p > 3
     t = n(:);
     n = [];
end

if ~isempty(n)
    ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
    if ok
       if n < 0
          t = -n;
          n = 1;
       else
          ok = ( ( n > 0 ) & ( mod(n,1) == 0 ) );
       end
    end
    if ~ok
        error('N must be an Integer larger Zero.');
    end
end


if nargin == 2
   if ~( isnumeric(scl) & ( prod(size(scl)) == 1 ) )
       error('Second Input must be a 1-element Numeric.');
   end
   pot = scl(1);
end

if isempty(t)
   if isempty(n)
     n = ceil(abs(diff(lim)/100));
     n = max(1,n);
   end
   t = linspace(lim(1),lim(2),n)';
end


n = size(t,1);

c = zeros(n,3);

 
c(:,1) =   56100000 * t .^ (-3 / 2) + 148;

c(:,2) = ( 35200000 * t .^ (-3 / 2) + 184 ) .* ( t > 6500 ) + ...
         ( 100.04 * log(t) - 623.6 ) .* ( t <= 6500 );
         
c(:,3) =   194.18 * log(t) - 1448.6;

c = min(max(c/255,0),1) .^ pot;

% c = ( c .* fak(:,[1 1 1]) ) .^ pot;
 
