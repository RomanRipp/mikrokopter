% TIDTEST  test and check for TIDANAL and TIDFIT
%
% Displayed values should be close to zero.
%
 
load sigma99

ind = ( 11 : 11 : 66 );  % Select Frequencies to test

sg = sigma(ind);
sg = sg(:);

sz = [ size(sg,1)  3 ];


off = round( 100 * randn([1 sz(2)]) );  % Test Offsets

am = round( 100 * rand(sz) );  % Test Amplitudes
ph = round( 360 * rand(sz) );  % Test Phases


t = ( 0 : 1/24/4 : 360 )';        % Time-Vector for Test, [days]

o1 = ones(size(t));
o2 = ones(size(off));

%--------------------------------------------------------
% Built TimeSeries

x = o1 * off;

for ii = 1 : sz(1)
 
   x = x +  am(ii*o1,:) .* cos( 2*pi * (24*(sg(ii)*t)*o2-ph(ii*o1,:)) / 360 );

end

%--------------------------------------------------------
% Tidal Analysis

[offs,am1,ph1,y,A,c,mest] = tidanal(t,x,sigma);

%--------------------------------------------------------
% Tidal Fit

z = tidfit(t,offs,am1,ph1,sigma);

%--------------------------------------------------------
% Check Differences

Offset_Diff = [ off - offs ]

Ampltitude_Diff = [ am-am1(ind,:)  ]

Phase_Diff = [ ph-(ph1(ind,:)+360*(ph1(ind,:)<0)) ]

Fit_Deviation = [ max(abs(z-x),[],1) ; max(abs(z-y),[],1) ]

am1(ind,:) = 0;

Res_Amplitude_Max = max(am1,[],1)
