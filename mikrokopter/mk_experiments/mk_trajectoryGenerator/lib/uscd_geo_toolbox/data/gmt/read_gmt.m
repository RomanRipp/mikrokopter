function [xy,lv] = read_gmt(file,varargin)

% READ_GMT   Read GMT-Data from binned NetCDF-Files
%
% XY = READ_GMT( File , Area )
%
%   Area:  Area, for wich Data are to Extract:
%          [ LonMin LonMax  LatMin  LatMax ]
%
%   XY:   Lines: [ Lon  Lat ] 
%         (2 Columns, Lines separated by NaN)
%
% [ XY , Level ] = READ_GMT( File , Area , Levels )
%
%  Returns only Lines with specified Levels, see below
%
% To use a ShortName for the File, use up to 3 Letters, specifies one
%  of the DataTypes, the Resolution and the Level, descripted below.
%
% Example:  File = 'sh'  == Shore, High Resolution == binned_shore_h.cdf  
%           File = 'sh2' == Shore, High Resolution, Level 2
%
% Special Feature of SquareKilometers of Feature for ShoreLines:
%
% [ XY , Level+i*SquareKilometer ] = READ_GMT( File , ... , +MinSkm*i )
% [ XY , Level+i*SquareKilometer ] = READ_GMT( File , ... , -MaxSkm*i )
%
% [ XY , Level+i*SquareKilometer ] = READ_GMT( File , ... , [Min Max]Skm*i )
%
%  Returns only Lines with a Minimum/Maximum Value of 0.1*SquareKilometers,
%    only valid for  Shore-(Coast)-LineData
%  Note: A single ZERO-Number is also accepted, 
%        the Area of the parent Segment is returned.
%
%
% see also:  LOOK_CDF (required), SHORE_GMT, WRT_GMT, GEBCO_WVS
%
%**********************************************************************
% Description of GMT-DataBase
%
%----------------------------------------------------------------------
% DataTypes:
%                                  File
%  s = Shore                    == binned_shore_#.cdf
%  c = Coast (equal to shore)
%  r = River                    == binned_river_#.cdf
%  b = Boundaries               == binned_border_#.cdf
%  t = Tectonics                == binned_tectonic_#.cdf
%
%  
%----------------------------------------------------------------------
% Resolutions "#" for DataTypes:
%                                 Scale         Resolution at Eq.
%  f = Full                  ==  1 :    100,000     100m
%  h = High                  ==  1 :  1,000,000    1000m
%  i = Intermediate          ==  1 :  3,000,000    3000m
%  l = Low                   ==  1 : 10,000,000   10000m
%  c = Crude                 ==  1 : 30,000,000   30000m
%
%
%----------------------------------------------------------------------
% River Levels:
%
%  1 = Permanent major rivers
%  2 = Additional major rivers
%  3 = Additional rivers
%  4 = Minor rivers
%  5 = Intermittent rivers - major
%  6 = Intermittent rivers - additional
%  7 = Intermittent rivers - minor
%  8 = Major canals
%  9 = Minor canals
% 10 = Irrigation canals
%
%  You can also use the following Letters in the ShortName for DataType:
%
%  r = All permanent rivers (1-4)
%  i = All intermittent rivers (5-7)
%  c = All canals (8-10)
%
%  Example:  File = 'rhc' == River, High Resolution, All canals
%            File = 'rh1' == River, High Resolution, Permanent major rivers
%            File = 'rh0' == River, High Resolution, Irrigation canals
%
%
%----------------------------------------------------------------------
% CoastLine Levels:
%
%   1 = WVS   land/ocean
%   2 = WDBII land/lake
%   3 = WDBII lake/island-in-lake
%   4 = WDBII island-in-lake/lake-in-island-in-lake
%
%
%----------------------------------------------------------------------
% Boundary Levels:
%
%   1 = National boundaries
%   2 = State boundaries within the Americas
%   3 = Marine boundaries
%
%
%----------------------------------------------------------------------
% Tectonic Levels:
%
%    1 = Ridges
%    2 = Trenches
%    3 = Transform zones
%    4 = Fracture zones
%    5 = Magnetic lineations
%    6 = Hotspots locations
%
%


xy = zeros(0,2);
lv = zeros(0,1);

if nargin == 0
   return
end

sc = 65535;

[msg,file,area,level,mskm] = checkin(file,varargin);

if ~isempty(msg)
    error(msg)
end

is_level = ( nargout > 1 );

%***************************************************************
% Check for VariableNames, get ID's

%  { Required  ShortName  {NetCDF-VarName(s)} }

ini = { 1 'bsi'  { 'Bin size in minutes'             }
        1 'nbf'  { 'N bins in file'                  }
        1 'nsf'  { 'N segments in file'              }
        1 'npf'  { 'N points in file'                }
        1 'snr'  { 'N segments in a bin'             }
        1 'sid'  { 'Id of first segment in a bin'    }
        1 'pid'  { 'Id of first point in a segment'  }
        1 'lon'  { 'Relative longitude from SE corner of bin' 'Relative longitude from SW corner of bin' }
        1 'lat'  { 'Relative latitude from SE corner of bin'  'Relative latitude from SW corner of bin'  }
        0 'lev'  { 'Hierarchial level of a segment'  }
        0 'emb'  { 'Embedded <npts, levels, exit, entry> for a segment' 'Embedded npts levels exit entry for a segment' } 
        0 'skm'  { 'Area in 0.1km^2 of the parent polygon of a segment' 'Ten_times_the_km_squared_area_of_the_parent_polygon_of_a_segment'}  };

%       'snr'  'N segments in a bin'
%       'pnr'  'N points for a segment'

[msg,vnr] = chk_gmt(file,ini);

if ~isempty(msg)
    error(msg)
end

embed_ok = ( ~isempty(vnr.emb) & isempty(vnr.lev) );
if embed_ok
   vnr.lev = vnr.emb;
end

level_ok =   ~isempty(vnr.lev);
  skm_ok = ( ~isempty(vnr.skm) & ~all(isnan(mskm)) );

%***************************************************************
% Open NetCDF-File

[fid,stat] = ncmex('open',file,'nowrite');

if ~( ( fid > 0 ) & ( stat == 0 ) )
  file1 = which(file);
  [fid,stat] = ncmex('open',file1,'nowrite');
  if ~( ( fid > 0 ) & ( stat == 0 ) )
     error(['Can''t open NetCDF-File: ' file ]);
     return
  end
end

%***************************************************************
% Get BinSize and BinNumbers, Offset

[bsi,status] = ncmex('vargetg',fid,vnr.bsi,0,1,1);

bsi = bsi / 60;    % [min] --> [deg]

[bnr,x0,y0,xl,yl] = area2bin(area,bsi);

%***************************************************************
% Get Segments

n = prod(size(bnr));

snr = zeros(n,1);   % N segments in a bin
sid = zeros(n,1);   % Id of first segment in a bin

pnr = cell(n,1);    % N points for a segment
pid = cell(n,1);    % Id of first point in a segment

lev = cell(n,1);    % Hierarchial level of a segment
lok = cell(n,1);    % LevelOk of a segment

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if 0 % embed_ok
   i01 = cell(n,1); % [ Entry ; Exit ] of a segment
                    % 0 From South Side of bin
                    % 1 From East  Side of bin
                    % 2 From North Side of bin
                    % 3 From West  Side of bin
                    % 4 Closed Polygon
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

np  = zeros(n,1);   % Total Number of points for a segment
ns  = zeros(n,1);   % Total Number of valid segments

[nbf,status] = ncmex('vargetg',fid,vnr.nbf,0,1,1);
[nsf,status] = ncmex('vargetg',fid,vnr.nsf,0,1,1);
[npf,status] = ncmex('vargetg',fid,vnr.npf,0,1,1);

for ii = 1 : n

  %  [sid(ii),status] = ncmex('vargetg',fid,vnr.sid,bnr(ii)-1,1,1);
  %  if bnr(ii) < nbf
  %     [next_id,status] = ncmex('vargetg',fid,vnr.sid,bnr(ii)-0,1,1);
  %  else
  %      next_id = nsf;
  %  end
  %  snr(ii) = next_id - sid(ii);

    [snr(ii),status] = ncmex('vargetg',fid,vnr.snr,bnr(ii)-1,1,1);
    [sid(ii),status] = ncmex('vargetg',fid,vnr.sid,bnr(ii)-1,1,1);

    if ~( snr(ii) == 0 )

        [pid{ii},status] = ncmex('vargetg',fid,vnr.pid,sid(ii),snr(ii),1);

        if sid(ii)+snr(ii) < nsf
          [next_id,status] = ncmex('vargetg',fid,vnr.pid,sid(ii)+snr(ii),1,1);
        else
           next_id = npf;
        end

        pnr{ii} = diff(cat(2,pid{ii},next_id),1,2);

%       [pnr{ii},status] = ncmex('vargetg',fid,vnr.pnr,sid(ii),snr(ii),1);
        
        %********************************************************************
        if ~level_ok

            lev{ii} = NaN * pid{ii};
            lok{ii} = ones(size(pid{ii}));

        %********************************************************************
        else

           [lev{ii},status] = ncmex('vargetg',fid,vnr.lev,sid(ii),snr(ii),1);

           %------------------------------------------------------------------
           if embed_ok  % 'Embedded <npts, levels, exit, entry> for a segment'

              % GMT_SHORE.C
              %
              %   seg[s].level = (seg_info[s] >> 6) & 7;
              %   seg[s].n     = (seg_info[s] >> 9);
              %   seg[s].entry = (seg_info[s] >> 3) & 7;
              %   seg[s].exit  =  seg_info[s] & 7;
              %
              %  ">>" == bitshift right (divide)
              %  "&"  == bitand 
       
             nr = bitshift( lev{ii} , -9 );
             i0 = bitand( bitshift( lev{ii} , -3 ) , 7 );
             i1 = bitand( lev{ii} , 7 );

              lev{ii} = bitand( bitshift( lev{ii} , -6 ) , 7 );

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%              i01{ii} = cat( 1 , bitand( bitshift( lev{ii} , -3 ) , 7 ) , ...
%                                 bitand( lev{ii} , 7 )
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           end
           %------------------------------------------------------------------

           if isempty(level)
              lok{ii} =  ones(1,snr(ii));
           else
              lok{ii} = zeros(1,snr(ii));
                 for ll = level
                     lok{ii} = ( lok{ii} | ( lev{ii} == ll ) );
              end
           end

        end
        % level_ok

        %********************************************************************
        if skm_ok

           [skm,status] = ncmex('vargetg',fid,vnr.skm,sid(ii),snr(ii),1);

            lev{ii} = lev{ii} + i*skm;

            if ~isnan(mskm(1))
                lok{ii} = ( lok{ii} & ( skm >=  mskm(1) ) );
            end
            if ~isnan(mskm(2))
                lok{ii} = ( lok{ii} & ( skm <=  mskm(2) ) );
            end
 
        end

        %********************************************************************

         np(ii) = sum(pnr{ii}.*lok{ii});
         ns(ii) = sum(lok{ii});

    end

end

if sum(np) == 0 
   return
end

%***************************************************************
% Get Points of Bins

xy = NaN * zeros( 2+is_level , sum(np)+sum(ns) );

i0 = cat(1,0,cumsum(np+ns,1));  % StartIndex, ZERO-Based

for ii = 1 : n

    if ~( ns(ii) == 0 )

        %---------------------------------------------------
        if all(lok{ii})

           % StartIndex in xy

           i1 = cumsum( cat( 2 , 1 , pnr{ii}(1:(snr(ii)-1))+1 ) , 2 );

           ind = grp2ind(i1,pnr{ii}) + i0(ii);

           if snr(ii) == 1
              ind = ind(:)';
           end

           [xy(1,ind),status] =  ncmex('vargetg',fid,vnr.lon,pid{ii}(1),np(ii),1);
           [xy(2,ind),status] =  ncmex('vargetg',fid,vnr.lat,pid{ii}(1),np(ii),1);

           if is_level

              jj = ones(size(ind));
              di = diff(ind);
              kk = find( di == 1 );
              jj(kk+1) = 0;
              jj = cumsum(jj,2);

              xy(3,ind) = lev{ii}(jj);

              xy(3,ind(end)+1) = lev{ii}(jj(end));

              kk = find( di > 1 );
              jj = jj(kk);
              kk = ind(kk) + 1;
              xy(3,kk) = lev{ii}(jj);

           end

        %---------------------------------------------------
        else

           ok = find(lok{ii});

           pn = pnr{ii}(ok);
           pd = pid{ii}(ok);

           % StartIndex in xy

           i1 = cumsum( cat( 2 , 0 , pn(1:(ns(ii)-1))+1 ) );

           for jj = 1 : ns(ii)

               ind = ( 1 : pn(jj) ) + i1(jj) + i0(ii);

               [xy(1,ind),status] =  ncmex('vargetg',fid,vnr.lon,pd(jj),pn(jj),1);
               [xy(2,ind),status] =  ncmex('vargetg',fid,vnr.lat,pd(jj),pn(jj),1);

               if is_level
                  xy(3,ind) = lev{ii}(ok(jj));
                  xy(3,ind(pn(jj))+1) = lev{ii}(ok(jj));
               end

           end

           ind = grp2ind(i1+1,pn) + i0(ii);

        end

        %---------------------------------------------------
        % Values of xy :   -32768 .. -1  |   0 .. 32767

        xy([1 2],ind) = ( xy([1 2],ind) + (sc+1) * ( xy([1 2],ind) < 0 ) ) / sc;       

        xy([1 2],ind) = xy([1 2],ind) * bsi;

        xy(1,ind) = xy(1,ind) + x0(ii);
        xy(2,ind) = xy(2,ind) + y0(ii);


    end

end

status = ncmex('close',fid);


%***************************************************************
% Extract Data inside Area

xy = permute(xy,[2 1]);

ok = ( ( xy(:,1) >= xl(1) ) & ...
       ( xy(:,1) <= xl(2) ) & ...
       ( xy(:,2) >= yl(1) ) & ...
       ( xy(:,2) <= yl(2) )       );

% Take 1 DataSet before and after too
ok( find( diff(ok) ==  1 ) + 0 ) = 1;
ok( find( diff(ok) == -1 ) + 1 ) = 1;

% Segment-Separator (NaN) is also ok

ok = ( ok | isnan(xy(:,1)) );

xy(find(~ok),[1 2]) = NaN;

%***************************************************************
% Remove following NaN-Values

jj = find(isnan(xy(:,1)));
jj = jj( find(diff(jj)==1) + 1 );

xy(jj,:) = [];

if isnan(xy(1,1))
   xy(1,:) = [];
end


if is_level
   lv = xy(:,3);
   xy = xy(:,[1 2]);
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,area,level,mskm] = checkin(file,vin);

msg = cell(0,1);

area  = [];
level = [];
mskm  = NaN * zeros(1,2);
 
%********************************************************************
% Check File

if ~chkstr(file,1)
    msg = cat( 1 , msg , { 'FileName must be a nonempty String.' } );
end

ok = ( exist(file,'file') == 2 );

nf = size(file,2);

if ~ok | ( nf <= 3 )

    form = 'binned_%s_%s.cdf';

    switch lower(file(1))
      case 't'
           typ = 'tectonic';
      case 'r'
           typ = 'river';
      case 'b'
           typ = 'border';
      case { 's' 'c' }
           typ = 'shore';
      otherwise
           typ = '';
    end

    ok = ~isempty(typ);

    if ~ok

       msg = cat( 1 , msg , { 'Invalid Identifer for Datatype.' } );

    else

       % Get Resolution
       if nf >= 2
          res = lower(file(2));
       else
          res = 'l';
       end

       % Get Levels

       if nf == 3
          lev = lower(file(3));
          if ( double('0') <= double(lev) ) & ( double(lev) <= double('9') )
             level = eval(lev);
             level = level + 10 * ( level == 0 );
          elseif strcmp(typ,'river')
             switch lower(file(3))
               case 'r'
                    level = ( 1 : 4 );
               case 'i'
                    level = ( 5 : 7 );
               case 'c'
                    level = ( 8 : 10 );
             end
          end
       end

       file1 = sprintf(form,typ,res);

       ok = ( exist(file1,'file') == 2 );
       if ok
            file = file1;
       else
           msg = cat( 1 , msg , { [ 'File not exists: ' file1 ] } );
       end

    end

elseif ~ok

   msg = cat( 1 , msg , { [ 'File not exists: ' file ] } );

end 

%********************************************************************
% Check other Inputs

n = prod(size(vin));

ok = zeros(n,1);

for ii = 1 : n
    
    v = vin{ii};

    ok(ii) = isnumeric(v);
    if ok(ii)
       ok(ii) = all(isfinite(v(:)));
       if ok(ii)
          if isequal(size(v),[1 4]) & isempty(area)
             area = v;
          elseif ( prod(size(v)) == 1 ) & all( real(v) == 0 )
             v = imag(v);
             mskm(1+(v<0)) = abs(v);
          elseif ( prod(size(v)) == 2 ) & all( real(v) == 0 )
             v = imag(v);
             mskm = v(:)';
          else
             level = cat( 2 , level , v(:)' );
          end
       end
    end

end

if isempty(area)
   area = [ -180 180 -90 90 ];
end

if ~isempty(level)
    level = sort(level);
    level( find(diff(level)==0) + 1 ) = [];
end

if ~any(isnan(mskm)) & ( mskm(2) < mskm(1) )
    msg = cat( 1 , msg , { 'Intervall for SquareKilometers must increasing.' } );
end
    
if any(~ok)
   msg = cat( 1 , msg , { 'Following Inputs must be finite Numerics.' } );
end

if isempty(msg)
   msg = '';
   return
end

sep = char([10 32 32 32]);

msg = sprintf('Invalid Inputs.%s%s',sep,strhcat(msg,sep));

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bnr,x0,y0,xl,yl] = area2bin(area,bsi)

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

dx = xl(1) - xl(2);

xl(2) = xl(2) + px * ( floor(dx/px) + 1 ) * ( dx >= 0 );

bx = [ floor(xl(1)/bsi)  ceil(xl(2)/bsi) ];

bx(2) = bx(2) + nx * ( bx(2) <= bx(1) );

bx = ( bx(1) : bx(2)-1 ) + 1;

bx = bx - nx * ( ceil( bx / nx ) - 1 );

%*********************************************
% Transform Latitude

yl = yl + py/2;                            % [ -py/2 .. py/2 ] --> [ 0 .. py ]

yl = yl - 2*py * floor( yl / (2*py) );     %  [ 0 .. 2*py ]

yl = yl + ( 2*py - 2*yl ) .* ( yl > py );  %  /\: [ 0 .. py ]

yl = yl - py/2;                            % [ -py/2  py/2 ]

by = -yl + py/2;

by = by( [ 1  2 ] + [ 1 -1 ] * ( by(2) < by(1) ) );

by = ( floor(by(1)/bsi) :  ceil(by(2)/bsi)-1 )' + 1;

yl = yl( [ 1  2 ] + [ 1 -1 ] * ( yl(2) < yl(1) ) );

%*********************************************
% Get BinNumber

sb = [ size(by,1) size(bx,2) ];

bnr = bx(ones(sb(1),1),:) + ( by(:,ones(1,sb(2))) - 1 ) * nx;

%*********************************************
% Get BinOffset: SW-Corner of Bin

dx = ( diff(bx,1,2) < 0 );

if any(dx)

   dx     = find(dx) + 1;
   ok     = 0 * bx;
   ok(dx) = 1;

   bx = bx + nx * cumsum(ok);

end

x0 = bsi * ( bx - 1 );

dx = x0(1) - xl(1);

x0 = x0 - px * ( ceil(dx/px)  );

y0 = py/2 - by * bsi;

x0 = x0(ones(sb(1),1),:);
y0 = y0(:,ones(1,sb(2)));

%*********************************************

bnr = permute(bnr,[2 1]);
x0  = permute(x0,[2 1]);
y0  = permute(y0,[2 1]);

bnr = bnr(:);

x0  = x0(:);
y0  = y0(:);

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,c] = chk_gmt(file,ini)

c      = ini(:,[2 2])';
c(2,:) = { {[]} };

c = struct(c{:});

[msg,d,v] = look_cdf(file);

if ~isempty(msg)
    return
end

n  = size(ini,1);
ok = zeros(n,1);

for ii = 1 : n

    ini{ii,3} = cat( 2 , ini{ii,3} ,  strrep(ini{ii,3},' ','_') );

    for nn = ini{ii,3}
  
        jj = strcmp( v(:,1) , nn{1} );

        ok(ii) = any(jj);

        if ok(ii)
           break
        end

    end

    if ok(ii)
 
       c = setfield( c , ini{ii,2} , find(jj)-1 );
 
    else

      ini{ii,3} = ini{ii,3}{1};

    end

end

ok = ( ok | ~cat(1,ini{:,1}) );

if all(ok)
   return
end

ok = find(~ok);

msg = 'Doesn''t found following Variables in NetCDF-File: %s%s%s';

sep = char([10 32 32 32]);

msg = sprintf(msg,file,sep,strhcat(ini(ok,3),sep));

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);


%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

