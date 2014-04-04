% Another Test for TIDANAL

T = [  0.5 1 14 28 180 ]; % Period [day]
A = [  10  5  1  2 0.5 ]; % Amplitude  

hw = [ 2000 01 01 12 00 00 
       2000 01 01 01 00 00
       2000 01 12 03 00 00
       2000 01 28 04 00 00
       2000 03 08 23 00 00 ];  % Date of HighWater

hw = datenum( hw(:,1) , hw(:,2) , hw(:,3) , ...
              hw(:,4) , hw(:,5) , hw(:,6)       )';  % Day of HighWater


day = ( datenum(1999,01,01) : 1/24 : datenum(2001,01,01) )';

ph = 360 * hw./T;              % Phase of HighWater
ph = ph - 360*floor(ph/360);


f = 1./T;  % Frequency [1/day]

n = size(day,1);

u = A(ones(1,n),:) .* cos( 2*pi*day*f - ph(ones(1,n),:)*pi/180 );

figure, hold on

plot(day,u);

plot(day,sum(u,2),'k-','linewidth',3);

timeaxis(gca,'zoom');


sg = 360 ./ (24*T);  % Sigma

[off,am,ph,fit] = tidanal(day,sum(u,2),sg');


% T = 360/sg / 24;  % day
%
% cos( 2*pi*hw*f - ph*pi/180 ) == 1
%      2*pi*hw*f - ph*pi/180   == n*2*pi
%
% hw*f - ph/360 == n
%
% hw/T - ph/360 == n; mod(n,1) == 0 
%
% hw = T * ( n + ph/360 ) 
%
