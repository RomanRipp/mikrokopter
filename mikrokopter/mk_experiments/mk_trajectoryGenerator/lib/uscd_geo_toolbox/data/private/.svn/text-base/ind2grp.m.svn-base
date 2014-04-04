function [i0,l] = ind2grp(ii);

% IND2GRP  Built StartIndex and Length from IndexVector
%
% [ StartIndex , GroupLength ] = IND2GRP( Index )
%
% See also: GRP2IND, GRPMEAN, CUMSUM, CUMPROD
%

Nout = nargout;

i0 = zeros(0,1);
l  = zeros(0,1);

if isempty(ii);
   return
end

ii = ii(:);
n  = size(ii,1);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , n+1 );

l  = diff(i0,1,1);

i0 = ii(i0(1:end-1));
