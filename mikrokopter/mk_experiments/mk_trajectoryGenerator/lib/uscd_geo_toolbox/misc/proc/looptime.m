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
%   BLN      BlankNumber   before each Line of TimeInformation
%
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
%   TXT  TimeInformations  (Char)
%  
% If the last Output TXT is not requested,
%   and BLN is not EMPTY, the TimeInformations will
%   displayed in Matlab's CommandWindow.
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


Nin  = nargin;
Nout = nargout;

nl = char(10);

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


% Display in CommandWindow

txt_out = ( ( Nout < 6 )  & ~isempty(bl) );

 
if isempty(d00)
  d00 = d1;
end

if isempty(d0)
  d0 = d00;
end

if isempty(bl);
   bl = 0;
end

bl = char(32*ones(1,bl));

 
 % All Units in  [ day ]

  dt = ( d1 - d0  );                   % Last Loop
 mdt = ( d1 - d00 ) / (ii+(ii==0));    % Mean of Loop's


 et = mdt - dt;
 et = dt + et * ( abs(et) < 0.1*mdt );

 et = et * ( N - ii );       % Estimated Time
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

 txt = cell(5,3);
 txt(:,1) = { bl };
 txt(:,3) = { nl };

 txt(:,2) = {  [ '     Loop Time: '   sprintf('%6.2f sec', dt) ]
               [ 'Mean Loop Time: '   sprintf('%6.2f sec',mdt) ]
               [ '  Elapsed Time: '   sprintf(dform,del) ]
               [ 'Estimated Time: '   sprintf(dform,det) ] 
               [ '       EndDate: '   datestr(de,0)      ] };

 txt = permute( txt , [ 2  1 ] );

 txt = cat( 2 , nl , txt{:} , nl );

 % Display Times

 if txt_out
 % No TextOutput requested

   fprintf( txt );

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
