function time=str2day(ss,ini)

% STR2DAY  converts string with any TimeNotation into decimal days
%
% dd = str2day( TimeString , INI ) 
%
%   takes the First Number in the TimeString as
%
%    Day  when INI == 1,
%    Hour      INI == 2,
%    Min       INI == 3, 
%    Sec       INI == 4.
%
% If there is any Error during evaluating the string, 
%  or the string contains an empty value, STR2DAY returns NaN.
%


if nargin < 2
 ini = 1;
end

if isempty(ss) | ~ischar(ss)
 ss = ' ';
end


ss = [ ' '  ss(1,:)  ' ' ];


%------------------------------------------
% Accept only [ 43 45 46 48 .. 57 ]
% "+-.0123456789"
 
ii = find( ss < 43 | ss > 57 | ss == 47 | ss == 44);

ss(ii) = 32+0*ii;


%------------------------------------------
% Accept "." only if Number before or after


ind = [ 2 : length(ss)-1 ];

ii = find( [ss(ind  ) == 46 ]  & ... 
          ~[ss(ind-1) >= 48 ]  & ...
          ~[ss(ind+1) >= 48 ]        );

ss(ii+1) = 32+0*ii;


%------------------------------------------
% Accept "+-" only if any Number follows

jj = find( [ ss == 43 ]  |  [ ss == 45 ] );

ii = find( jj > max(find(ss>=48)) );

ss(jj(ii)) = 32+0*jj(ii);




val = eval( [ '[' char(ss) ']' ] , '[]' );


if isempty(val)

 time = nan;

else

  time = [0 0 0 0];

  ende = ini-1+length(val);
  ende = ende+(4-ende)*[ende>4];

  time( ini : ende ) = val((ini:ende)-ini+1);

if 0
  fak = -1+2*[ time(ini) >=0 ]; 
  if ini < ende
   time(ini+1:ende) = time(ini+1:ende)*fak;
  end
end

  time = sum( time .* [ 1  1/24  1/24/60  1/24/3600 ] );

end
