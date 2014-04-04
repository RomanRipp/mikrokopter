function dat = gebco_wvs(form,area);

% GEBCO_WVS  Read's GEBCO-WorldVectorShoreLine
%
%  XY = GEBCO_WVS( GebcoFileFormat , Area );
%
%   GebcoFileFormat:  Format, using SPRINTF, to build FileNames
%
%   Area:  Area, for wich Data are to Extract:
%          [ LonMin LonMax  LatMin  LatMax ]
%
%   XY:   ShoreLine: [ Lon  Lat ] 
%         (2 Columns, Lines separated by NaN)
%
% Possible GebcoFileFormats:
%
%  ShortFormat        FileFormat                 Scale         Resolution at Eq.
%
%  'k100' |          'k100_%3.3d_%3.3d.bin' ==  1 :    100,000     100m
%  'k250' |          'k250_%3.3d_%3.3d.bin' ==  1 :    250,000     250m
%  'm1'   |            'm1_%3.3d_%3.3d.bin' ==  1 :  1,000,000    1000m
%  'm3'   |            'm3_%3.3d_%3.3d.bin' ==  1 :  3,000,000    3000m
%  'm10'  | 'm12'  |  'm12_%3.3d_%3.3d.bin' ==  1 : 10,000,000   10000m
%  'm30'  | 'm43'  |  'm43_%3.3d_%3.3d.bin' ==  1 : 30,000,000   30000m
%
% Other Formats may be used, if they exist in a 10-degree Resolution.
%
% The FileNotation is: sprintf( GebcoFileFormat , [ Lat Lon ] )
% where Lat and Lon defines the lower left Corner of an 10-degree Rectangular.
%
%  Lat = [ 0 : 10 : 350 ] == [ -180 .. 170 ]
%  Lon = [ 0 : 10 : 170 ] == [  -90 ..  80 ]
%
% Other Formats may be used, if they exist in the same 10-degree Notation.
%
% Note: DataFiles must exist on Disk and Directories should be
%         included in Matlab's SearchPath
%       
% see also: SHORE_GMT, READ_GMT
%
% Example:  >> addpath D:\data\gebco\wvs\1to250k
%
%           >> xy = gebco_wvs('k250',[ 8 32  53 67 ]);
%           >> % xy = gebco_wvs('k250_%3.3d_%3.3d.bin',[ 8 32  53 67 ]);
%
%           >> figure, plot( xy(:,1) , xy(:,2) );
%
%

if nargin < 1
   error('Not enough InputArguments.')
end

if ~( ischar(form)  &  ~isempty(form)  &  ( prod(size(form)) == size(form,2) ) )
   error('GebcoFileFormat must be a String.');
end

form = lower(form);

switch form
       case 'm10'
            form = 'm12';
       case 'm30'
            form = 'm43';
end

if any( strcmp( form , { 'k100' 'k250' 'm1' 'm3' 'm12' 'm43' } ) )
   form = cat( 2 , form , '_%3.3d_%3.3d.bin' );
end

if nargin < 2
   area = [ -180  180  -90  90 ];
else
   area  = area(:)';
   ok = ( isnumeric(area) & ( size(area,2) >= 4 ) );
   if ok
      area = area(1:4);
      ok = all(isfinite(area));
   end
   if ~ok
      error('Area must be a 4-Element Vector.')
   end
end

%----------------------------------------------------------

machineformat = 'ieee-le';  % MachineFormat

prec = 'long';              % PRECISION for FREAD
  rl = 2;                   % RecordLength

bsi =  10;  % BinSize

%----------------------------------------------------------

[bx,by,xl,yl,x0] = area2bin(area,bsi);

nx = size(bx,2);
ny = size(by,2);

%----------------------------------------------------------

dat = cell(nx*ny,1);

dat(:) = { ones(0,2) };


for ix = 1 : nx
  for iy = 1 : ny

    name = sprintf(form,by(iy),bx(ix));

    fid = fopen(name,'r',machineformat);

    if fid ~= -1

      % [ NumberOfCodes  0            ]
      % [ Code           RecordOffset ]
      %  Data

      fread(fid,2*rl,prec);  

      cc = ny * (ix-1) + iy;

      dat{cc} = fread(fid,prec);

      fclose(fid);

      dat{cc} = reshape(dat{cc},rl,size(dat{cc},1)/rl)';

      % Flip Columns  [ Lat  Lon ] --> [ Lon  Lat ] 
      dat{cc} = dat{cc}(:,[2 1]);
 
      % Extract Segments, separate by NaN's
      % No Data defined by Coords == -1

      % Start of Segment at zz with  [ Coords Altitude]

      zz = 1;
       c = 1;
      while ( zz < size(dat{cc},1) )  &  ( c > 0 )
        c = dat{cc}(zz,1);   % Coords
        dat{cc}(zz,:) = NaN;
        zz = zz+c+1;
      end

      % Data till last Segment, Scale with 1e-6
      
      dat{cc}= dat{cc}(1:min(size(dat{cc},1),zz-1),:) / 1e6;

      dat{cc}(:,1) = dat{cc}(:,1) + x0(ix);

      % Find Data inside Area

        ok = ( ( dat{cc}(:,1) >= xl(1) ) & ...
               ( dat{cc}(:,1) <= xl(2) ) & ...
               ( dat{cc}(:,2) >= yl(1) ) & ...
               ( dat{cc}(:,2) <= yl(2) )       );

      % Take 1 DataSet before and after too
 
        ok( find( diff(ok) ==  1 ) + 0 ) = 1;
        ok( find( diff(ok) == -1 ) + 1 ) = 1;


      % Segment Separator (NaN) is also ok

        ok = ( ok | isnan(dat{cc}(:,1)) );


      % Set all bad Data to NaN

        dat{cc}(find(~ok),:) = NaN;


      % Remove duplicate NaN-Row's

         jj = find(isnan(dat{cc}(:,1)));
         jj = jj( find( diff(jj) == 1 ) + 1 );

         dat{cc}(jj,:) = [];

         if isnan(dat{cc}(1,1))
            dat{cc} = dat{cc}(2:end,:);
         end

         if ~isempty(dat{cc})
             if ~isnan(dat{cc}(end,1))
                 dat{cc} = cat(1,dat{cc},[NaN NaN]);
             end
         end

    end
    % fid ~= -1
  end
  % i2
end    
% i1  

dat = cat(1,dat{:});

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bx,by,xl,yl,x0] = area2bin(area,bsi)

% GMT-Area: ( 0 : bsi : 360 ), ( 90 : -bsi : -90 )

area = area(:)';

xl = area([1 2]);
yl = area([3 4]);

px = 360;
nx = px / bsi;

py = 180;
ny = py / bsi;

%*********************************************
% Transform Longitude

xb = xl(1);   % Base

xl = xl - px * floor((xl+px/2)/px);  % [ -px/2 .. px/2 )

xl(2) = xl(2) + px * ( ( xl(2) == -px/2 ) | ( xl(2) <= xl(1) ) );

bx = xl + px/2;    % [ -p/2 .. p/2 ) -- GEBCO_WVS -->  [ 0 .. p )

bx = [ floor(bx(1)/bsi)  ceil(bx(2)/bsi) ];

bx(2) = bx(2) + nx * ( bx(2) <= bx(1) );

bx = ( bx(1) : bx(2)-1 );

bx = bx + nx * ( ( bx < 0 ) - ( bx >= nx ) );

bx = bsi * bx;

%*********************************************
% Transform Latitude

yl = yl + py/2;                            % [ -py/2 .. py/2 ] --> [ 0 .. py ]

yl = yl - 2*py * floor( yl / (2*py) );     %  [ 0 .. 2*py ]

yl = yl + ( 2*py - 2*yl ) .* ( yl > py );  %  /\: [ 0 .. py ]

yl = yl - py/2;                            % [ -py/2  py/2 ]

yl = yl( [ 1  2 ] + [ 1 -1 ] * ( yl(2) < yl(1) ) );

by = yl + py/2;  % [ -p/2 .. p/2 ) -- GEBCO_WVS -->  [ 0 .. p )

by = bsi * ( floor(by(1)/bsi) :  ceil(by(2)/bsi)-1 );

%*********************************************
% Get BinOffset: SW-Corner of Bin

x0 = bx - px/2;

dx = diff(x0,1,2);
lg = ( x0(1) > xl(1) );

if any( dx < 0 )

   i0 = find( dx < 0 );
   i0 = i0(1);          % Last Value before ZERO

   ok       = 0*x0;
   ok(1:i0) = 1;

   x0 = x0 - px * ( ok * lg - ~ok * ~lg );

elseif lg

   x0 = x0 - px;

end

% y0 = by - py/2;

%*********************************************
% Set to XBase

x0 = x0 + ( xb - xl(1) );
xl = xl + ( xb - xl(1) );

x0 = x0 - ( bx - px/2 );  % XOffset
