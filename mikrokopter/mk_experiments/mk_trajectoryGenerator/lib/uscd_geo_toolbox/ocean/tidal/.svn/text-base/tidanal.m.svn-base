function [off,am,ph,fit,A,co,mest]=tidanal(time,px,sg,lim)
 
% TIDANAL  Tidal Analysis using SingularValueDecomposition
%
% [off,am,ph,Y,A,C,mest] = TIDANAL(T,X,Sigma,Limit)
%
% Tidal analysis of time-series X with corresponding time-vector T,
%  using angular frequencies sigma ( 360 / Period[hour] )
%    
%  T      time-vector [days]  1 Column   [L by 1 ]
%  X    N time-series         N Columns  [L by N ]
%  
%  Sigma  M circular frequencies, 360 * freq, 360 / P, [1/hour]
%
%  Limit  relative cutoff for Eigenvalues of SVD-Analysis, default: 1e-3
%
%  Limit = NaN gets a solution solution in the least squares sense to the
%           under- or overdetermined system of equations A*C = X.
%           where A = [ 1 cos(sigma1) sin(sigma1) cos(sigma2) sin(sigma2) ... ] 
% Outputs:
%
%  off     Offset for each time-series, [ 1 by N ]
%  am      Amplitude, [ M by N ]
%  ph      Phases,    [ M by N ]
%
%  Y       Tidal fit for each time-series [ L by N ]
%
%  A       EigenVectors, [ L by 1+2*M ]
%           [ Offset cos(sigma1) sin(sigma1) cos(sigma2) sin(sigma2) ... ]
%  C       Coefficients, [ M by N ], Y = A * C
%
%  mest    SVD Matrice, mest = V * L * U'
%
%           [U,L,V]=svd(A,0);  
%           exclude EigenValues where diag(L)/max(diag(L)) < Limit
%           mest = V * L * U'
%           C = mest * X
%
% Calculation of the tidal fit:
%
%      1.    Y = A * C
%
%  or  2.    Y = off + am * cos( 2*pi * (sigma*time*24-ph) / 360 )
%
%
% see also: TIDFIT, SVD, MLDIVIDE
%

Nin = nargin;
if Nin < 4
   lim = 1e-3;
end

is_svd = isfinite(lim);

time = time(:) * 24;   % [hour]


% T  = 360 ./ sigma
  
freq = sg(:)' / 360;   % [1/hour]

n = size(time,1);
m = size(freq,2);


m = 2 * m + 1;

ics = ( 2 : 2 : (m-1) );  % Index for Cosine
isn = ics + 1;            % Index for Sine

A = ones(n,m);

A(:,ics) = cos(2*pi*time*freq);
A(:,isn) = sin(2*pi*time*freq);


if is_svd

   % SVD    Singular value decomposition.
   % [U,S,V] = SVD(X) produces a diagonal matrix S, of the same 
   % dimension as X and with nonnegative diagonal elements in
   % decreasing order, and unitary matrices U and V so that
   % X = U*S*V'.

   fprintf(1,'Run SVD Analysis ...')

   [U,L,V]=svd(A,0);

   fprintf(1,'done');

   si    = min(size(L));

   ld = diag(L(1:si,1:si));

   jj = ( ld./max(ld) >= lim );

   if any(~jj)
      fprintf(', exclude %.0f of %.0f EigenValues (Treshold %.1e)\n', ...
              sum(~jj),si,lim)
   end

   jj = find(jj);

   mest = V(:,jj) * L(jj,jj)^(-1) * U(:,jj)';

   co = mest*px;
 
else

   mest = [];

   fprintf(1,'Divide Matrices ...')

   co = A \ px;

   fprintf(1,'done');

end

fprintf(1,'\n');

am = sqrt( co(isn,:).^2 + co(ics,:).^2 );

ph = atan2( co(isn,:) , co(ics,:) ) * 180/pi;

off = co(1,:);

if nargout >= 4
   fit = A * co;                            % fitted curve
else
   fit = [];
end
