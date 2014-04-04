function sg = tidfit(time,off,am,ph,sg)

% TIDFIT  Tidal fit using results from TIDANAL
%
%  Y = TIDFIT(time,off,am,ph,Sigma)
%
%  time    time-vector [days]  1 Column   [L by 1 ]
%  off     Offset for each time-series, [ 1 by N ]
%  am      Amplitude, [ M by N ]
%  ph      Phases,    [ M by N ]
%
%  Sigma   M circular frequencies, 360 * freq, 360 / P, [1/hour]
%
%  Y       Tidal fit for N time-series, N Columns  [L by N ]
%  
% Calculation of the tidal fit:
%
%  Y = off + am * cos( 2*pi * (sigma*time*24-ph) / 360 )
%
%
% see also: TIDANAL
%



 time = time(:);   % [ L by 1 ]
 off  =  off(:)';  % [ 1 by N ]
 sg   =   sg(:)';  % [ 1 by M ]

 n = size(off,2);
 m = size( sg,2);

 if ~( isequal(size(am),[m n]) & isequal(size(ph),[m n]) )
     error('Size of Amplitude and Phase must match Size of Offset and Sigma.')
 end

 am = permute(am,[3 2 1]);  % [ 1 by N by M ];
 ph = permute(ph,[3 2 1]);  % [ 1 by N by M ];
 
 n1 = ones(1,n);
 t1 = ones(1,size(time,1));

 sg = 2 * pi / 360 * sg * 24;
 ph = 2 * pi / 360 * ph;

 sg = time * sg;                  % [ L by M ]
 sg = permute(sg(:,:,n1),[1 3 2]); % [ L by N by M ]

 sg = am(t1,:,:) .* cos( sg - ph(t1,:,:) );

 sg = sum(sg,3);                  % [ L by N ]

 sg = sg + off(t1,:);             % Add Offset

