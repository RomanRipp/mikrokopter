function [gmt,clk,dt] = utc;

% UTC returns actual UTC-Time on a UNIX-System
%
% UTC_Time = UTC
%
% UTC_Time = [ YYYY MM DD hh mm ss ] 
%
% use UNIX: date --utc +"%Y %m %d %H %M %S"
%
% A second Output returns the Vector of Matlab's CLOCK
%
% A third Output returns the Deviation: Local - UTC
%
% see also: CLOCK
%

w = [];

if isunix
   [s,w] = unix('date --utc +"%Y %m %d %H %M %S"');
     clk = clock;
else
   warning('UNIX required.')
   clk = clock;
end


ok = ~isempty(w);
if ok
   w  = eval(['[' w ']'],'NaN');
   ok = ( isnumeric(w) & isequal(size(w),[1 6]) );
end

if ok
   gmt = w;  
else
   gmt = NaN * ones(1,6);
end


clk = floor(clk);


if ( nargout == 0 ) & ~isnan(gmt(1))

   fprintf(1,'%s\n',datestr(gmt,0))

   clear gmt

elseif nargout == 3

   dt = cat( 1 , clk([3 4 5 6]) , gmt([3 4 5 6]) );

   dt = dt(1,:) - dt(2,:);  % CLK - GMT

   if abs(dt(1)) > 1
          dt(1)  = -sign(dt(1));
   end

   dt(4) = dt(4) .* ( abs(dt(4)) > 1 );  % Seconds
   
   dt = sum( dt .* [ 1  1/24  1/24/60  1/24/3600 ] ); 

end

  
