function par_edt(par,group,action,clb,varargin)

% CLB_EDT CallBack for EDIT-Buttons


ud = get(par,'userdata');

ch = getfield( ud.Children , group );

he = getfield( ch , action                 );  % Edit
hs = getfield( ch , ['Slider' action(end) ]);  % Slider


str  = get( he , 'string'    );
form = get( hs , 'userdata' );

[val,str,msg1,msg2] = val2str(str,form);

if isempty(val) | ~isempty(msg1) | ~isempty(msg2)
   val = NaN;
end
 
str0 = get( he , 'userdata' );

if isequal( val , get(hs,'value') )
   set( he , 'string' , str0 );
   return
end


if ( get(hs,'min') <= val )  &  ( val  <= get(hs,'max') );


   set( hs   , 'value'    , val );       
   set( he   , 'string'   , str , ...
               'userdata' , str      );

  if clb
   
    % Call other Functions

  end 

else

   set( he , 'string' , str0 );


    [val1,str1] = val2str(get(hs,'min'),form);
    [val2,str2] = val2str(get(hs,'max'),form);

    warndlg([ 'Value must between ' ...
               str1 ' and ' str2  '.'] ,'Invalid Input','warn');

end                                                              




%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [val,str,msg1,msg2]=val2str(val,form)

% VAL2STR  Converts between Valus and Formated Strings
%
% [ Value , String , Message1 , Message2 ] = VAL2STR( Argument , Format )
%
% Argument could be a Value or a String, Format gives the Format for the
%  String.
%
%  Valid Formats are:
%                                                   '#'   String
%  '%##.##'   for using by SPRINTF
%
%  'geo#'     for Geographic Coordinates  Longitude: 1     ##°##.##' <E/W>
%                                         Latitude : 2     ##°##.##' <N/S>
%
%  'time#'    for TimeConversation, first Value Day: 1    DD HH:MM
%                                              Hour: 2       HH:MM
%                                              Hour: 3       HH:MM:SS
%                                            Minute: 4          MM:SS
%
% the Outputs are the String in the given Format 
%  and the correspondending Value 
%    (this Value may be different (rounded)  from the InputValue,
%      depending on the Format )
%
% require functions STR2DAY, DAY2DHMS
%


str  = '';
msg1 = '';
msg2 = '';


if isempty(val) 
 val = [];
 str = '';
 return
end

if isempty(form)
  if ischar(val)
   str = val;
   val = [];
  end
  return
end


if ischar(val)
 if all(val==' ')
  str = val;
  val = [];
  return
 end
else
 if ~isnumeric(val)
  msg1 = ' VAL2STR: Inputs must be Numeric or String (Char).';
  val = [];
  str = '';
  return
 end
end




form = lower(form(:)');

ini = sum( cumsum( strcmp( form(end) , {'1' '2' '3' '4'} ) ) == 0 ) + 1;

ini = ini - 4 * ( ini == 5 );

is_geo  = ( ~isempty( strmatch('geo' ,form) )           );
is_time = ( ~isempty( strmatch('time',form) ) & ~is_geo );



ww = warnstat;

warning('off');


if ischar(val)

 %*******************************************************
 % STRING --> VALUE
 %-------------------------------------------------------

 

 if     is_geo

  ws =  ( ( any(lower(val)=='w')  &  ( ini == 1 ) ) | ...
          ( any(lower(val)=='s')  &  ( ini == 2 ) )       );

  val=str2day(val,2)*24*(1-2*ws);

 elseif is_time

  fak = [ 1  24  24  24*60 ];
  Ini = [ 1   2   2   3    ];
  val = str2day(val,Ini(ini)) * fak(ini);

 else
 
  ok = 1;
  eval( [ 'val = ['  val  '];' ] , 'ok=0;' )

  if ~ok
   msg1 = [ ' VAL2STR : ' lasterr ];
   val = [];
  end

 end

end
% isstr(val) STRING --> VALUE



if isempty(msg1)

%*************************************************
% VALUE --> STRING --> VALUE  using given FORMAT
%-------------------------------------------------

 if ~isempty( findstr(form,'%') )

  str = sprintf(form,val);

  ok = 1;

  eval( [ 'val = ['  str  '];' ] , 'ok=0;' )

  if ~ok
   msg2 = [' VAL2STR: ' lasterr ];
  end

 
 elseif is_geo  |  is_time

    fak = [ 1+23*is_geo   24   24  24*60 ];
   
    [DD,hh,mm,ss] = day2dhms(val/fak(ini) );

   MM = round(mm+ss/60);
   HH = hh + MM/60*( abs(MM) == 60 );
   MM = MM*(~( abs(MM) == 60 ));

   if is_geo

        ws = [ val < 0 ];
        WS = char( (69+18*ws)*(ini==1) + (78+5*ws)*(ini==2) );

        hh = hh+24*DD;

        gf = sprintf('%.0f',4);

       str = sprintf(['%'  gf  '.0f' char(176) ' %5.2f'' ' WS ],abs([hh mm+ss/60]));  

       val = hh + (mm+ss/60)/60;

   elseif is_time

      switch ini
       case 1
         str = sprintf('%3.0f  %2.2d:%2.2d',[DD HH MM]);  
         val = DD + HH/24 + MM/24/60 ;
       case 2  
         str = sprintf('%2.2d:%2.2d',[HH+DD*24 MM]);
         val = HH+DD*24 + MM/60 ; 
       case 3
          ss = round(ss);
          mm = mm + ss/60*[ abs(ss) == 60 ];
          ss = ss*(~( abs(ss) == 60 ));

         str = sprintf('%2.2d:%2.2d:%2.2d',[hh+DD*24  mm ss]);
         val = hh*60+DD*24*60 + mm + ss/60 ;
       case 4
          ss = round(ss);
          mm = mm + ss/60*[ abs(ss) == 60 ];
          ss = ss*(~( abs(ss) == 60 ));

         str = sprintf('%2.2d:%2.2d',[hh*60+DD*24*60+mm ss]);
         val = hh*60+DD*24*60 + mm + ss/60 ;
      end
      % ini
   end
   % form
  end
  % '%'

end


warning(ww);


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [dd,hh,mm,ss] = day2dhms(day)

% DAY2DHMS   converts  dezimal day  into  Day Hour Mine Sec 
%
%  [Day,Hour,Min,Sec] = day2dhms(day)
%

dd = fix( day );
hh = fix( 24 * (day-fix(day)) );
mm = fix( 60 * ( 24 * (day-dd) - hh ));
ss  = round( 24*3600*((day-dd)-hh/24-mm/(24*60)) );

mm_ss = fix(ss/60);
   ss = ss - mm_ss * 60;
   mm = mm + mm_ss;

hh_mm = fix(mm/60);
   mm = mm - hh_mm * 60;
   hh = hh + hh_mm;

dd_hh = fix(hh/24);
   hh = hh - dd_hh * 24;
   dd = dd + dd_hh;




