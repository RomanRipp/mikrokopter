function  Msg = lim_cdf(infile,outfile,buffer)

% LIM_CDF  Remove UNLIMITED Dimension from an NetCDF-File
%
%  Msg = LIM_CDF( InFile , OutFile , Buffer )
%
%-------------------------------------------------------------------
% Inputs:
%
%  InFile   FileName of NetCDF-File, contains UNLIMITED Dimensions
%
% OutFile   FileName of new NetCDF-File, contains same Data as InFile,
%            but no UNLIMITED Dimensions
%
%  Buffer   Intervall along the UNLIMITED Dimension to Read and Write 
%             the UNLIMITED Variables, optional
%
%-------------------------------------------------------------------
% Output:
%
%  Msg      String contains detected Errors
%
%
%-------------------------------------------------------------------
% Required m-files
%
%   LOOK_CDF,   READ_CDF,   WRITE_CDF (CREATE_CDF)
%
%-------------------------------------------------------------------
%
%  see also: LOOK_CDF
%


Msg = '';
nl  = char(10);


Nin = nargin;


%*******************************************************
% Check Inputs

%-------------------------------------------------------
if Nin < 2
 Msg = 'Minimum InFile and OutFile requested.';
 return
end


%-------------------------------------------------------

ok = ( ischar(infile)  &  ~isempty(infile) &  ...
       ( prod(size(infile)) == size(infile,2) )   );
if ~ok
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'InFile must be a nonempty String.' ];
end


%-------------------------------------------------------

ok = ( ischar(outfile)  &  ~isempty(outfile) &  ...
       ( prod(size(outfile)) == size(outfile,2) )   );
if ~ok
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'OutFile must be a nonempty String.' ];
end


%-------------------------------------------------------

if isequal(infile,outfile)
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'InFile and OutFile must be different.' ];
end


%-------------------------------------------------------

if Nin < 3
  buffer = [];
end

ok = ( isnumeric(buffer)  &  ( prod(size(buffer)) <= 1 ) );
if ok
  if ~isempty(buffer)
    ok = (  ( buffer > 0 )  &  ( mod(buffer,1) == 0 )  ); 
  end
end

if ~ok
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'Buffer must be a single Integer larger ZERO or EMPTY.' ];
end


%-------------------------------------------------------

if ~isempty(Msg)
  return
end

 
%******************************************************************
% Get Informations about NetCDF-File

fprintf(nl)
fprintf([ '  LOOK_CDF( '   infile   ' )'   nl nl ]);
 

[Msg,dim,var,att] = look_cdf(infile);

% dim = { DimName DimLength  ''          }
%                            'UNLIMITED'
%
% var = { VarName VarType Ndim [dim] Nattr }
%
% att = Attributes
%

if ~isempty(Msg)
  return
end


%-----------------------------------------------------------------
% Get UNLIMITED Dimension

% Index of UNLIMITED Dimension

ud = find( strcmp( dim(:,3) , 'UNLIMITED' ) );


if isempty(ud)
  Msg = [ 'No UNLIMITED Dimension in NetCDF-File: ' infile ];
  return
end


fprintf([ ' UNLIMITED Dimension: '   dim{ud,1}   ...
            sprintf(' = %.0f',dim{ud,2})  nl nl ]);


%-----------------------------------------------------------------
% Get UNLIMITED Variables

nv = size(var,1);

uv = zeros(nv,1);

for ii = 1 : nv

  if var{ii,3}
  % NDim > 0
    % Check 1. Dimension of Variable with  "ud" 
    uv(ii) = ( var{ii,4}(1) == ud-1 );
  end

end


if ~any(uv)
  Msg = [ 'No Variable with UNLIMITED Dimension in NetCDF-File: ' infile ];
  return
end


% Index of Variables without UNLIMITED Dimensions

cv = find( ~uv );


% Index UNLIMITED Variables

uv = find( uv );


fprintf([ ' UNLIMITED Variables: '   ...
           strhcat(var(uv,1),'  ')   nl nl ]);


%******************************************************************
% Prepare new NetCDF-File

%-----------------------------------------------------------------
% Expand "var"  with Value and Status
%
% var = { VarName VarType Ndim [dim] Nattr Value Status }
%
 
  var = cat(2,var,cell(nv,2));

  var(:,7) = { [] };  % Value Empty
  var(:,8) = {  0 };


%-----------------------------------------------------------------
% Read  Variables without UNLIMITED Dimensions

fprintf([ '  READ_CDF( '   infile   ' ), '  ...
          'Variables without UNLIMITED Dimensions'   nl nl ]);


% Fill Variables without UNLIMITED Dimensions into var1 
%  incl.  their Value

[Msg,dim1,var(cv,:)] = read_cdf(infile,'var',var(cv,1));

% var = { VarName VarType Ndim [dim] Nattr Value Status }
%


if ~isempty(Msg)
  return
end

% Check Status 
read_err = find( cat(1,var{:,8}) == -1 );
if ~isempty(read_err)
  Msg = [ 'Error read Variables, NCMEX  Status -1 : '  ...
           nl   strhcat(var(read_err,1),'  ')  ];
  return
end 


%******************************************************************
% Write new NetCDF-File

%  UNLIMITED Variables still empty
%

fprintf([ ' WRITE_CDF( '   outfile   ' ), '  ...
          'fill Data of Variables without UNLIMITED Dimensions'   nl nl ]);

 [CreateMsg,WriteMsg] = write_cdf(outfile,dim(:,[1 2]),var(:,[1 2 3 4 7]),att);

 Msg = [ CreateMsg  nl(1:(end*(~isempty(CreateMsg))))  WriteMsg ];

 if ~isempty(Msg)
   return
 end



%******************************************************************
% Fill Data of UNLIMITED Variables

%-----------------------------------------------------------------
% Open NetCDF-Files

fid1 = ncmex( 'open' , infile  , 'nowrite' );
fid2 = ncmex( 'open' , outfile , 'write'   );

if fid1 == -1
  Msg = [ 'Error using NCMEX( ''OPEN'' , '''  infile ''' , ''NoWrite'' )' ];
end

if fid2 == -1
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'Error using NCMEX( ''OPEN'' , '''  outfile ''' , ''Write'' )' ];
end

if ~isempty(Msg)
   return
end


%-----------------------------------------------------------------
% Loop

nn = 0;  % Start for UNLIMITED

nu = dim{ ud , 2 };  % Actual Value of UNLIMITED Dimension

if isempty(buffer)
  buffer = nu;
end

N  = ceil( nu / buffer );  % Max. Number of Loop's

   % Initialize StartTime
   t0 = looptime;      
   t1 = [];
   zz = 0;       % Counter

fprintf([ ' NCMEX, read and write Data of UNLIMITED Variables'  nl ...
          '        Buffer: '  sprintf('%.0f',buffer)  nl nl ]);


while nn <  nu

  zz = zz + 1;

  % Display Counter
 
  fprintf( nl );
  fprintf( '   %3.0f of %3.0f,   %6.0f',[ zz  N  nn ] );  
  fprintf( nl );


  n = buffer + ( ( nu - nn ) - buffer ) * ...
               ( ( nu - nn ) < buffer );


  % Read and Write UNLIMITED Variables
                
  for vv = uv(:)'

    count    = cat( 2 , dim{ var{vv,4} , 2 } );
    count(1) = n;

    start    = 0 * var{vv,4};
    start(1) = nn;

    [ val , stat ]  = ncmex( 'vargetg' , fid1 , vv-1 , start , count , [ 1 1] );

     if stat == -1
        Msg = [ 'Error using NCMEX( ''VARGET'' )' ];
        break
     end

      stat          = ncmex( 'varput' , fid2 , vv-1 , start , count , val );

     if stat == -1
        Msg = [ 'Error using NCMEX( ''VARPUT'' )' ];
        break
     end
   
  end
  % vv

  if ~isempty(Msg)
    break
  end


  nn = nn + buffer;

  if N > 1
  % Get LoopTime and Display with 12 Blanks
 
    t1 = looptime(ii,N,t0,t1,12);
  
  end

end
% while


%-----------------------------------------------------------------
% Close NetCDF-Files

stat1 = ncmex('close',fid1);
stat2 = ncmex('close',fid2);


if stat1 == -1
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'Error using NCMEX( ''CLOSE'' , '''  infile ''' )' ];
end

if stat2 == -1
  Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
          'Error using NCMEX( ''CLOSE'' , '''  outfile ''' )' ];
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ~( ischar(str)  |  iscellstr(str) )
   error('StringArray must be a CharArray or CellStringArray.');
end

if iscellstr(str)
  str = char(str);
end

str = double(str);

    jj    = find( sum( ( str == 32 ) , 2 ) == size(str,2) );
str(jj,:) = [];
 
str = cellstr(char(str));


str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function [d1,dt,mdt,et,de,txt] = looptime(ii,N,d00,d0,bl)

% LOOPTIME   calculates elapsed and estmated time of a Loop 
%
% [ AT, LT, MLT, ET, ED, TXT ] = LOOPTIME( LN, LL, ST, LLT, BLN );
%
%----------------------------------------------------------------
%  Inputs:
%
%    LN     LoopNumber
%    LL     LoopLength
%    ST        StartTime   [day]
%   LLT     LastLoopTime   [day]
%   BLN      BlankNumber  
%
%  If  BlankNumber is given ( 0 | 1 | .. ) , a TimeInformation will
%   displayed in Matlab's CommandWindow, with the Blanks at
%   the beginning of each Line.
%  If BlankNumber is empty or not given, no Information will displayed.
%
%
%----------------------------------------------------------------
%  Outputs:
%
%    AT       ActualTime   [day]  ( LastLoopTime for next Loop )
%    LT         LoopTime   [sec]   Duration of last loop
%   MLT     MeanLoopTime   [sec]   Mean Duration of all Loops
%    ET    EstimatedTime   [sec]
%    ED          EndTime   [day]  
%   TXT    TimeInfomations  (CellStringArray)
%  
% 
%----------------------------------------------------------------
%  Initialize the Loop before with LoopNumber == 0 :
%
%       StartTime = looptime;
%    LastLoopTime = [];                    
%
%
%----------------------------------------------------------------
% Note:  LOOPTIME, running on a PC, Pentium II, 300MHz, Linux, 
%          Matlab 5.3 (R11), free Workspace and RAM needs about:  
%           0.007 sec without DisplayOutput ( 0.7 sec per 100 loops), 
%           0.01  sec incl.   DisplayOutput (  1  sec per 100 loops),
%            incl. DisplayOutput the Time depends on the SytemGrafics.
%
%
%----------------------------------------------------------------
% Example for a loop using INTERP1:
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@qqqq
%
%  ll = 1e5;                  % Length of Vectors
%  XI = linspace(0,2*pi,ll);
%  
%   N = 10;                   % LoopLength
%
%  % Initialize StartTime
%  t0 = looptime;      
%  t1 = [];
%
%  % Start Loop
%  for ii = 1 : N
%
%     % Display Counter
%
%       fprintf(char(10));
%       fprintf( '   %3.0f of %3.0f ',[ ii  N ] );  
%       fprintf(char(10));
%
%     %-------------------------------
%     % Begin LoopBlock
%
%       X = 2*pi*rand(ll,1);
%       Y = interp1(XI,sin(XI),X);
%
%     % End LoopBlock
%     %-------------------------------
% 
%     % Get LoopTime and Display with 12 Blanks
%
%       t1 = looptime(ii,N,t0,t1,12);
%
%  end
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@qqqq
%
% see also: DATENUM, DATEVEC
%
%


Nin = nargin;

 % Get Actual Time

 t1  = clock;  % [ YY MM DD hh mm ss ]
 d1  = datenum(t1(1),t1(2),t1(3),t1(4),t1(5),t1(6));


if Nin < 1
 ii = 0;
end
if Nin < 2
 N = 1;
end
if Nin < 3
 d00 = [];
end
if Nin < 4
 d0 = [];
end
if Nin < 5
 bl = [];
end

if isempty(d00)
  d00 = d1;
end

if isempty(d0)
  d0 = d00;
end

 
 % All Units in  [ day ]

  dt = ( d1 - d0  );                   % Last Loop
 mdt = ( d1 - d00 ) / (ii+(ii==0));    % Mean of Loop's

 et = mdt*(N-ii);            % Estimated Time
 de = d1 + et;               % End Time
 dl = d1 - d00;              % Elapsed Time


 [dd,hh,mm,ss]=day2dhms(dl);
 del = [ hh+24*dd mm ss ];   % Elapsed TimeVector

 [dd,hh,mm,ss]=day2dhms(et);
 det = [ hh+24*dd mm ss ];   % Estimated TimeVector


 % Day ==> Seconds
 
  dt =  dt * 24 * 3600;
 mdt = mdt * 24 * 3600;
  et =  et * 24 * 3600;

 % TimeInformation

 dform = '%3.0f:%2.2d:%2.2d';

 txt = {   '     Loop Time: '   sprintf('%6.2f sec', dt) 
           'Mean Loop Time: '   sprintf('%6.2f sec',mdt)
           '  Elapsed Time: '   sprintf(dform,del)
           'Estimated Time: '   sprintf(dform,det)  
           '       EndDate: '   datestr(de,0)       };

 % Display Times

 if ~isempty(bl)
 
   nl  = char(10);  % NewLine
 
   str      = cell(size(txt,1),2);
   str(:,1) = { char(32*ones(1,bl)) };   % Blanks  at Begin of Line
   str(:,2) = { nl };                    % Newline at End   of Line

   str      = [ str(:,1)  txt  str(:,2) ]';

   fprintf([ nl cat(2,str{:}) nl ])

 end



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
