function [xyd,lev,skm,inf] = shore_gmt(file,varargin)

% SHORE_GMT_DB   Read GMT-ShoreLine-Data from binned NetCDF-Files
%
%                !!! Debugging Mode, old Version of AREA2BIN !!!
%
% XY = SHORE_GMT_DB( File , Area )
%
%   Area:  Area, for wich Data are to Extract:
%          [ LonMin LonMax  LatMin  LatMax ]
%
%   XY:   CellArray of closed Lines:  { [ Lon ;  Lat ] } 
%         (2 Rows per CellElement)
%
% [ XY , Level ] = SHORE_GMT( File , Area , Levels )
%
%  Returns only Lines with specified Levels, see below.
%
% [ XY , Level , SquareKilometer ] = SHORE_GMT( File , ... , +MinSkm*i )
% [ XY , Level , SquareKilometer ] = SHORE_GMT( File , ... , -MaxSkm*i )
%
% [ XY , Level , SquareKilometer ] = SHORE_GMT( File , ... , [Min Max]Skm*i )
%
%  Returns only Lines with a Minimum/Maximum Value of 0.1*SquareKilometers
%
% To use a ShortName for the File, use up to 2 Letters,
%  specifies the Resolution and the Level, descripted below.
%
% Example:  File = 'h'  == High Resolution == binned_shore_h.cdf  
%           File = 'h2' == High Resolution, Level 2
%
%
% Note: SHORE_GMT defines the global variables "cfg" "cnf" "bin" "inf" "hst"
%
% see also:  LOOK_CDF (required), READ_GMT, WRT_GMT, GEBCO_WVS
%
%**********************************************************************
% Description of GMT-ShoreLine-DataBase
%
%----------------------------------------------------------------------
% Resolutions for DataTypes:
%                                 Scale         Resolution at Eq.
%  f = Full                  ==  1 :    100,000     100m
%  h = High                  ==  1 :  1,000,000    1000m
%  i = Intermediate          ==  1 :  3,000,000    3000m
%  l = Low                   ==  1 : 10,000,000   10000m
%  c = Crude                 ==  1 : 30,000,000   30000m
%
%
%----------------------------------------------------------------------
% Levels:
%
%   1 = WVS   land/ocean
%   2 = WDBII land/lake
%   3 = WDBII lake/island-in-lake
%   4 = WDBII island-in-lake/lake-in-island-in-lake
%
%

Show = 0;   % Debugging: Show Development of Segments !!!

xyd = cell(0,1);
lev = zeros(0,1);
skm = zeros(0,1);

if nargin == 0
   return
end

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
        1 'emb'  { 'Embedded <npts, levels, exit, entry> for a segment' 'Embedded npts levels exit entry for a segment' }
        1 'skm'  { 'Area in 0.1km^2 of the parent polygon of a segment' 'Ten_times_the_km_squared_area_of_the_parent_polygon_of_a_segment'}  };


[msg,vnr] = chk_gmt(file,ini);

if ~isempty(msg)
    error(msg)
end

%***************************************************************
% Open NetCDF-File

fid = ncmex('open',file,'nowrite');

if fid == -1
  file1 = which(file);
  fid = ncmex('open',file1,'nowrite');
  if fid == -1
     error(['Can''t open NetCDF-File: ' file ]);
     return
  end
end

%***************************************************************
% Define global Variables

global cfg cnf bin inf hst

cfg = struct( 'fid' , { fid } , ...
              'vnr' , { vnr } , ...
              'nbf' , { [] } , ...
              'nsf' , { [] } , ...
              'npf' , { [] } , ...
              'bsi' , { [] } , ...
              'sc'  , { 65535 } , ...
              'nx'  , { [] }         );

cnf = struct( 'xl'    , { [] } , ...
              'yl'    , { [] } , ...
              'level' , { level } , ...
              'mskm'  , { mskm  } , ...
              'rcl'   , { 100 }       );
 
%***************************************************************
% Get BinSize and BinNumbers, Offset

[cfg.bsi,status] = ncmex('vargetg',fid,vnr.bsi,0,1,1);

cfg.bsi = cfg.bsi / 60;    % [min] --> [deg]

[bnr,x0,y0,cnf.xl,cnf.yl,cfg.nx,cfg.ny] = area2bin(area,cfg.bsi);

bin = cat(2,bnr,x0,y0);

%***************************************************************
% Get Segments of bins

n = size(bin,1);

np  = zeros(n,1);
bok = zeros(n,1);

[cfg.nbf,status] = ncmex('vargetg',fid,vnr.nbf,0,1,1);
[cfg.nsf,status] = ncmex('vargetg',fid,vnr.nsf,0,1,1);
[cfg.npf,status] = ncmex('vargetg',fid,vnr.npf,0,1,1);

inf = struct( 'lev' , { [] } , ...     % Hierarchial level of a segment
              'skm' , { [] } , ...     % Square kilometers of a segment
              'i01' , { [] } , ...     % [ Entry ; Exit ] ID of a segment at the border of bin
              'x01' , { [] } , ...     % [ First ; Last ] Longitude of a segment, crossing a border of bin
              'y01' , { [] } , ...     % [ First ; Last ] Latitude  of a segment, crossing a border of bin
              'cls' , { [] } , ...     % True for closed in bin
              'ins' , { [] } , ...     % True for Data inside Area
              'chk' , { [] } , ...     % CheckValue for connection
              'hst' , { [] } , ...
              'xyd' , { [] }       );  % Data for border crossing segments, inside Area

dat = struct( 'snr' , { [] } , ...     % N segments in a bin
              'pid' , { [] } , ...     % N points for a segment
              'pnr' , { [] } , ...     % Id of first point in a segment
              'xyd' , { [] }        ); % XY-Data of points

dat = dat(ones(n,1),:);
inf = inf(ones(n,1),:);

for ii = 1 : n

    %----------------------------------------------------------------
    % SegmentInfo

    [ dat(ii).snr , dat(ii).pid , dat(ii).pnr , inf(ii).lev , ...
      inf(ii).skm , inf(ii).i01 , inf(ii).x01 , inf(ii).y01 , inf(ii).cls  ] = ...
    seg_info(ii);

    %-------------------------------------------------------------------
    % Read XY-Data of Points

    [dat(ii).xyd,inf(ii).ins] = read_dat( dat(ii) , inf(ii).cls , ii );

    inf(ii).xyd = cell(dat(ii).snr,1);
    inf(ii).hst = cell(dat(ii).snr,1);

    if ~( dat(ii).snr == 0 )

        %----------------------------------------------------
        % Remove bad segments which are closed and not inside

        ok = ~( inf(ii).cls & ~inf(ii).ins );

        if ~all(ok)

            dat(ii).snr = sum(ok);

             ok = find(ok);

            dat(ii).pid = dat(ii).pid(ok);
            dat(ii).pnr = dat(ii).pnr(ok);
            dat(ii).xyd = dat(ii).xyd(ok);
            inf(ii).lev = inf(ii).lev(ok);
            inf(ii).skm = inf(ii).skm(ok);
            inf(ii).i01 = inf(ii).i01(:,ok);
            inf(ii).x01 = inf(ii).x01(:,ok);
            inf(ii).y01 = inf(ii).y01(:,ok);
            inf(ii).cls = inf(ii).cls(:,ok);
            inf(ii).ins = inf(ii).ins(:,ok);
            inf(ii).xyd = inf(ii).xyd(ok);
            inf(ii).hst = inf(ii).hst(ok);

        end

        %----------------------------------------------------
        % Set Data of border crossing segments inside Area to "inf"

        if any( ~inf(ii).cls )

           ok = find(~inf(ii).cls);

           inf(ii).xyd(ok) = dat(ii).xyd(ok);

        end

    end

    np(ii) = sum(dat(ii).pnr);

    % Checked by other segment in [ Entry ; Exit ]
    %  BinIndex + i * SegmentIndex

    inf(ii).chk = NaN*zeros(2,dat(ii).snr);

end

%***************************************************************

if sum(np) == 0 
   clear global cfg cnf bin inf hst
   return
end

%***************************************************************
% Search for connected segments

ww = warnstat;

if isequal(ww,'backtrace')
   warning('on');
end

rl0 = get( 0 , 'RecursionLimit' );
                     
set( 0 , 'RecursionLimit' , max(rl0,cnf.rcl+10) ); 

ns = zeros(n,1);

for ii = 1 : n

    if np(ii)

       for jj = find( inf(ii).cls == 0 );

           if ~isnan(inf(ii).cls(jj))

               xyd = inf(ii).xyd{jj};
               hst = ii + i*jj;
              mode = isnan(inf(ii).chk(2,jj));  

               ib = ii;
               is = jj;

              %***************************************************************************
              % RecursionLoop to take care before RecursionLimit

if Show
   f1 = figure; a1 = gca; hold on, 
   z = 0;
   line('parent',a1,'xdata',xyd(1,:),'ydata',xyd(2,:),'color','k');
   line('parent',a1,'xdata',xyd(1,1),'ydata',xyd(2,1),'marker','o');
   line('parent',a1,'xdata',xyd(1,end),'ydata',xyd(2,end),'marker','x');
end
              while 1
             
                 [xy,ok] = rec_seg( ib , is , zeros(2,0) , mode , 0 );

if Show
   z = z+1;
   patch('parent',a1,'xdata',[xy(1,:) NaN],'ydata',[xy(2,:) NaN], ...
      'cdata',[z+0*xy(1,:) NaN], ...
      'edgecolor','flat','facecolor','none', ...
      'marker','.','markerfacecolor','flat', ...
      'cdatamapping','direct');
   patch('parent',a1,'xdata',[xy(1,1) NaN],'ydata',[xy(2,1) NaN], ...
      'cdata',[z+0*xy(1,1) NaN], ...
      'edgecolor','flat','facecolor','none', ...
      'marker','o','markerfacecolor','none','markeredgecolor','flat', ...
      'cdatamapping','direct');
   patch('parent',a1,'xdata',[xy(1,end) NaN],'ydata',[xy(2,end) NaN], ...
      'cdata',[z+0*xy(1,end) NaN], ...
      'edgecolor','flat','facecolor','none', ...
      'marker','x','markerfacecolor','none','markeredgecolor','flat', ...
      'cdatamapping','direct');

   patch('parent',a1,'xdata',[xy(1,[1 end]) NaN],'ydata',[xy(2,[1 end]) NaN], ...
      'cdata',[z+0*xy(1,[1 end]) NaN], 'linestyle','--',...
      'edgecolor','flat','facecolor','none', ...
      'marker','x','markerfacecolor','none','markeredgecolor','flat', ...
      'cdatamapping','direct');
end

                 [xy,ins] = in_area(xy,cnf.xl,cnf.yl,cfg.bsi/cfg.sc);

if Show
   patch('parent',a1,'xdata',[xy(1,:) NaN],'ydata',[xy(2,:) NaN], ...
      'cdata',[z+0*xy(1,:) NaN], ...
      'edgecolor','flat','facecolor','none', ...
      'marker','.','markerfacecolor','flat', ...
      'cdatamapping','direct');
end

                 if mode
                    xyd = cat( 2 , xyd , xy );
                 else
                    xyd = cat( 2 , xy , xyd );
                 end

                 % No Segment found or Break by World-Surrounding Segment
                 brk = any( ok == [ -1  0 ] );

                 if mode & brk & isnan(inf(ii).chk(1,jj))
                    mode = 0;
                 elseif ~( ok == 9 )
                    break
                 end

                 ibs = hst( 1 + ( end - 1 ) * mode );
                 
                 ib = real(ibs);
                 is = imag(ibs);

              end
              %***************************************************************************

if Show
   cm = hsv(ceil(1.1*z));
   set(f1,'colormap',cm(1:z,:));
   % if ii == 1 & jj == 1, keyboard,end
end

               inf(ii).hst{jj} = hst;

               %--------------------------------------------------
               % Check for inside Area

               [xyd,inf(ii).ins(jj)] = in_area(xyd,cnf.xl,cnf.yl,cfg.bsi/cfg.sc);

               if ~inf(ii).ins(jj)

                   xyd = zeros(2,0);

               elseif ok == -1

                   % World-Surrounding Line: close polward

                   xyd  = xyd(:,[1 1:end end]);
                   xyd(2,[1 end])  = 90 * sign(mean(xyd(2,:)));

               end

               dat(ii).xyd{jj} = xyd;

           end

       end
       % jj

       ns(ii) = sum( ~isnan(inf(ii).cls) & inf(ii).ins );

    end                

end

set( 0 , 'RecursionLimit' , rl0 );

warning(ww);

%***************************************************************
% Close NetCDF-File

status = ncmex('close',fid);

%***************************************************************
% Extract good Data

nn = sum(ns);   
i0 = cat(1,0,cumsum(ns,1));  % StartIndex, ZERO-Based

xyd =  cell(nn,1);
lev = zeros(nn,1);
skm = zeros(nn,1);

ok  = zeros(nn,1);

for ii = 1 : n

    if ns(ii)

        jj = find( ~isnan(inf(ii).cls) & inf(ii).ins );

        ns(ii)   = min( ns(ii) , size(jj,2) );

        ind = i0(ii) + ( 1 : ns(ii) );

        xyd(ind) = dat(ii).xyd(jj);
        lev(ind) = inf(ii).lev(jj);
        skm(ind) = inf(ii).skm(jj);
 
         ok(ind) = 1;

    end

end

if ~all(ok)
    ok = find(ok);
    xyd = xyd(ok);
    lev = lev(ok);
    skm = skm(ok);
end

clear global cfg cnf bin inf hst

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,area,level,mskm] = checkin(file,vin);

msg = cell(0,1);

area  = [];
level = [];
mskm  =  NaN * zeros(1,2);

%********************************************************************
% Check File

if ~chkstr(file,1)
    msg = cat( 1 , msg , { 'FileName must be a nonempty String.' } );
end

ok = ( exist(file,'file') == 2 );

nf = size(file,2);

if ~ok & ( nf <= 2 )

    form = 'binned_shore_%s.cdf';

    % Get Resolution
      res = lower(file(1));

    if isequal(res,'s') & ( nf > 1 )
       if isletter(file(2))
          res = lower(file(2));
          file = file(2:end);
          nf = nf - 1;
       end
    end

    % Get Levels
    if nf == 2
       lev = lower(file(2));
       if ( double('0') <= double(lev) ) & ( double(lev) <= double('9') )
          level = eval(lev);
          level = level + 10 * ( level == 0 );
       end
    end

    file1 = sprintf(form,res);

    ok = ( exist(file1,'file') == 2 );
    if ok
       file = file1;
    else
       msg = cat( 1 , msg , { [ 'File not exists: ' file1 ] } );
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

function [bnr,x0,y0,xl,yl,nx,ny] = area2bin(area,bsi)

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

bx = [ floor(xl(1)/bsi)  ceil(xl(2)/bsi) ];

bx(2) = bx(2) + nx * ( bx(2) <= bx(1) );

bx = ( bx(1) : bx(2)-1 ) + 1;

bx = bx + nx * ( ( bx <= 0 ) - ( bx > nx ) );

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

x0 = ( bx - 1 ) * bsi;          % [   0   .. px   )

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

y0 = py/2 - by * bsi;

x0 = x0(ones(sb(1),1),:);
y0 = y0(:,ones(1,sb(2)));

%*********************************************
% Set to XBase

x0 = x0 + ( xb - xl(1) );
xl = xl + ( xb - xl(1) );

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

function [snr,pid,pnr,lev,skm,i01,x01,y01,cls] = seg_info(bid,mode)

global cfg cnf bin

if nargin < 2
   mode = 0;
end

[sid,status] = ncmex('vargetg',cfg.fid,cfg.vnr.sid,bin(bid,1)-1,1,1);
[snr,status] = ncmex('vargetg',cfg.fid,cfg.vnr.snr,bin(bid,1)-1,1,1);

if ( snr == 0 )
   pid = zeros(1,0);
   pnr = zeros(1,0);
   lev = zeros(1,0);
   skm = zeros(1,0);
   i01 = zeros(2,0);
   x01 = zeros(2,0);
   y01 = zeros(2,0);
   cls = zeros(1,0);
   return
end

%****************************************************************
% Get SegmentInfo

%------------------------------------------------------------------
% First ID of points and PointNumber of segment

[pid,status] = ncmex('vargetg',cfg.fid,cfg.vnr.pid,sid,snr,1);

if sid+snr < cfg.nsf
   [next_id,status] = ncmex('vargetg',cfg.fid,cfg.vnr.pid,sid+snr,1,1);
else
    next_id = cfg.npf;
end

pnr = diff(cat(2,pid,next_id),1,2);


%------------------------------------------------------------------
% GMT_SHORE.C
%
%   seg[s].level = (seg_info[s] >> 6) & 7;
%   seg[s].n = (seg_info[s] >> 9);
%   seg[s].entry = (seg_info[s] >> 3) & 7;
%   seg[s].exit  =  seg_info[s] & 7;
%
%  ">>" == bitshift right (divide)
%  "&"  == bitand 

[emb,status] = ncmex('vargetg',cfg.fid,cfg.vnr.emb,sid,snr,1);

%------------------------------------------------------------------
% [ Entry ; Exit ] and Level of segment

i01 = cat( 1 , bitand( bitshift( emb , -3 ) , 7 ) , ...
               bitand( emb , 7 )  );

lev = bitand( bitshift( emb , -6 ) , 7 );

%------------------------------------------------------------------
% Area of segment

[skm,status] = ncmex('vargetg',cfg.fid,cfg.vnr.skm,sid,snr,1);

%------------------------------------------------------------------
% Check segment with level and minimum area

if isempty(cnf.level)
   sok =  ones(1,snr);
else
   sok = zeros(1,snr);
   for ll = cnf.level
       sok = ( sok | ( lev == ll ) );
   end
end

if ~isnan(cnf.mskm(1))
    sok = ( sok & ( skm >= cnf.mskm(1) ) );
end

if ~isnan(cnf.mskm(2))
    sok = ( sok & ( skm <= cnf.mskm(2) ) );
end


cls = ( sum(i01==4,1) == 2 );   % True for not crosing border of bin

sok = ( sok & ( ~cls | isequal(mode,0) ) );

%------------------------------------------------------------------
% Remove bad segments

if ~all(sok)

   ok = find(sok);

   snr = sum(sok);

   pnr = pnr(ok);
   pid = pid(ok);
   lev = lev(ok);
   skm = skm(ok);
   cls = cls(ok);
   i01 = i01(:,ok);

end

%------------------------------------------------------------------
% Start and End-Points of segments, crossing border

x01 = NaN * zeros(snr,2);
y01 = NaN * zeros(snr,2);

for ii = find( ~cls )

    [x01(ii,:),status] =  ncmex('vargetg',cfg.fid,cfg.vnr.lon,pid(ii),2,pnr(ii)-1);
    [y01(ii,:),status] =  ncmex('vargetg',cfg.fid,cfg.vnr.lat,pid(ii),2,pnr(ii)-1);

end

x01 = ( x01 + (cfg.sc+1) * ( x01 < 0 ) ) / cfg.sc * cfg.bsi + bin(bid,2);       
y01 = ( y01 + (cfg.sc+1) * ( y01 < 0 ) ) / cfg.sc * cfg.bsi + bin(bid,3);       


x01 = permute(x01,[2 1]);
y01 = permute(y01,[2 1]);


%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xyd,ins] = read_dat(dat,cls,bid);

global cfg cnf bin

xyd = cell(dat.snr,1);

ins = zeros(1,dat.snr);

if dat.snr == 0
   return
end

for ii = 1 : dat.snr

    xyd{ii} = zeros(2,dat.pnr(ii));

    [x,status] = ncmex('vargetg',cfg.fid,cfg.vnr.lon,dat.pid(ii),dat.pnr(ii),1);
    [y,status] = ncmex('vargetg',cfg.fid,cfg.vnr.lat,dat.pid(ii),dat.pnr(ii),1);

    x = ( x + (cfg.sc+1) * ( x < 0 ) ) / cfg.sc * cfg.bsi + bin(bid,2);       
    y = ( y + (cfg.sc+1) * ( y < 0 ) ) / cfg.sc * cfg.bsi + bin(bid,3);       

    %--------------------------------------------------
    % Check for inside Area

    ins(ii) = any( ( cnf.xl(1) <= x ) & ( x <= cnf.xl(2) ) & ...
                   ( cnf.yl(1) <= y ) & ( y <= cnf.yl(2) )       );

    %--------------------------------------------------
    % Check for inside Area

    if ins(ii) 

       xyd{ii} = cat( 1 , x , y );

    elseif ~cls(ii)

       xyd{ii} = cat( 1 , x([1 dat.pnr(ii)]) , y([1 dat.pnr(ii)]) );

    end

end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xyd,ok] = rec_seg( ii , jj , xyd , mode , cnt );

global cfg cnf bin inf hst

%----------------------------------------------------

ok = 0;

px = 360;  % Period for X

%----------------------------------------------------
% Correspondings for [ Entry ; Exit ] - Flag
%
%   DIR    Flag  CorrFlag  BinOffs  dX0   dY0
%     
%   South   0      2        +nx      0    -bsi
%   East    1      3        +1      +bsi   0
%   North   2      0        -nx      0    +bsi
%   West    3      1        -1      -bsi   0
%
% mod([0 1 2 3],2) == [ 0  1  0  1 ]
%

%----------------------------------------------------
% Start at Exit

i0 = 1 + mode;
i1 = 2 - mode;

if ~isnan(inf(ii).chk(i0,jj))  % Allready checked
    ok = 2;
    return
end

f0 = inf(ii).i01(i0,jj);            % Exit Flag of segment

si = ( f0 > 1 );
md = mod(f0,2);

f1 = f0 + 2 - 4*si;                     % corresponding Entry  Flag

db = 1 + ( 1 - md ) * ( cfg.nx - 1 );   % BinDeviation

bn = bin(ii,1) + ( 1 - 2*si ) * db;     % Corresponding bin

bmd = mod( bin(ii,1) , cfg.nx ); 

bn = bn + md * cfg.nx * ( ( ( bmd == 1 ) & ( f0 == 3 ) ) - ...
                          ( ( bmd == 0 ) & ( f0 == 1 ) )          );

bn = bn + cfg.nbf * ( ( bn < 1 ) - ( bn > cfg.nbf ) );

dx =  ( 1 - 2*si ) *       md   * cfg.bsi;
dy = -( 1 - 2*si ) * ( 1 - md ) * cfg.bsi;

x00 = bin(ii,2) + dx;
y00 = bin(ii,3) + dy;


% Look for existing BinNumber with same XOffset !!!

sb = size(bin,1);

ok0 = zeros(sb,3);

ok0(:,1) = ( bin(:,1) == bn );
ok0(:,2) = ( bin(:,2) == x00 );

ok0(:,3) = ( sum(ok0(:,[1 2]),2) == 2 );

if any(ok0(:,3))

   ib = find(ok0(:,3));

else

   % Take care for World-Surrounding Lines (Antarctica), give a break
   if ( x00 < cnf.xl(1)-px-cfg.bsi ) | ( cnf.xl(2)+px+cfg.bsi < x00 ) 
      ok = -1;
      return
   end

   ib  = sb + 1;

   bin = cat( 1 , bin , [ bn x00 y00 ] );

   %------------------------------------------------------------------
   % BinNumber exists, change XOffset
   if any(ok0(:,1)) 
   %------------------------------------------------------------------
 
      bb = find(ok0(:,1));
      bb = bb(1);

      inf = cat( 1 , inf , inf(bb) );

      inf(ib).x01 = inf(ib).x01 + x00 - bin(bb,2);
      inf(ib).y01 = inf(ib).y01 + y00 - bin(bb,3);

      snr = size(inf(ib).lev,2);

      for kk = 1 : snr
          if ~isempty(inf(ib).xyd{kk})
              inf(ib).xyd{kk}(1,:) = inf(ib).xyd{kk}(1,:) + x00 - bin(bb,2);
              inf(ib).xyd{kk}(2,:) = inf(ib).xyd{kk}(2,:) + y00 - bin(bb,3);
          end
      end

      inf(ib).cls(find(isnan(inf(ib).cls))) = 0;

      inf(ib).chk = NaN*zeros(2,snr);

   %------------------------------------------------------------------
   % Extract new Bin
   else
   %------------------------------------------------------------------

      inf = cat( 1 , inf , inf(1) );

        [ snr , pid , pnr , inf(ib).lev , inf(ib).skm , ...
        inf(ib).i01 , inf(ib).x01 , inf(ib).y01 , inf(ib).cls ] = ...
       seg_info( ib , 1 );

       inf(ib).chk = NaN*zeros(2,snr);
       inf(ib).ins =     zeros(1,snr);
       inf(ib).xyd =  cell(snr,1);

    end

   %------------------------------------------------------------------

end

%----------------------------------------------------

jj0 = floor( (bin(ii,1)-1) / cfg.nx ) + 1;  % Column
ii0 = bin(ii,1) - (jj0-1) * cfg.nx;         % Row

jj1 = floor( (bn-1) / cfg.nx ) + 1;  % Column
ii1 = bn - (jj1-1) * cfg.nx;         % Row

form = 'Bin %.0f [ %.0f %.0f ] --%.0f-- [ %.4f %.4f ] --%.0f--> %.0f [ %.0f %.0f ]';

txt = sprintf(form,bin(ii,1),jj0,ii0,f0,inf(ii).x01(2,jj),inf(ii).y01(2,jj),f1,bn,jj1,ii1);

%----------------------------------------------------

if isempty(inf(ib).lev)
   msg = sprintf('%s\nNo Segments found in corresponding Bin.\n',txt);
   warning(msg);
   return
end

%-----------------------------------------------------------
% Corresponding Segments in Level and SquareKiloMeters

ok1 = ( ( inf(ib).lev == inf(ii).lev(jj) )  & ...
        ( inf(ib).skm == inf(ii).skm(jj) )        );

if ~any(ok1)
    msg = sprintf('%s\nNo corresponding Segments in Level and Area found.\n',txt);
    warning(msg);
    return
end
 
ok1 = find(ok1);

%-----------------------------------------------------------
% Corresponding Segments in Points with ExitPoint of Segment

x1 = inf(ii).x01(i0,jj);
x2 = inf(ib).x01(:,ok1);

x1 = x1 - px * floor((x1+px/2)/px);  % [ -px/2 .. px/2 )
x2 = x2 - px * floor((x2+px/2)/px);  % [ -px/2 .. px/2 )

n2 = prod(size(ok1));
o2 = ones(1,n2);

dx = abs( x2 - x1([1 1],o2) );
dy = abs( inf(ib).y01(:,ok1) - inf(ii).y01([i0 i0],jj*o2) );

acc = 1e-10;

for kk = [ 2  1 ] + [ -1  1 ]*mode  % [ Entry ; Exit ]
    ok2 = ( ( inf(ib).i01(kk,ok1) == f1 ) & ...
            ( dx(kk,:) <= acc ) & ( dy(kk,:) <= acc ) );
    if any(ok2)
       break
    end
end

if ~any(ok2)

    ok2 = ( ( inf(ib).i01(:,ok1) == f1 ) & isnan(inf(ib).chk(:,ok1)) & ...
            ( ( ( dx <= acc ) & md ) | ( ( dy <= acc ) & ~md ) ) );

    if ~any(ok2(:))
        msg = sprintf('%s\nNo corresponding Segments in Entry/Exit found.\n',txt);
    else
        msg = sprintf('%s\nNo corresponding Segments in Position found.\n',txt);
    end

    if 1

        warning(msg);
        return

    else

        ok2 = find(ok2);

        if md
        % East | West 
            [hh,kk] = min(dy(ok2));  [ md x1  hh ]
        else
        % South | North
            [hh,kk] = min(dx(ok2));  [md inf(ii).y01(i0,jj) hh ]
        end

        ok2 = ok2(kk);

        kk = mod(ok2-1,2) + 1;
        ok2 = floor( (ok2-1) / 2 ) + 1;

    end

else

    ok2 = find(ok2);

end

is = ok1(ok2);

if ~( isnan(inf(ib).chk(kk,is)) | isequal(inf(ib).chk(kk,is),ii+i*jj) )
    msg = sprintf('%s\nCorresponding Segment allready connected.\n',txt);
    warning(msg);
    return
end

% Found in Exit ==> Flip
if kk == i0

   inf(ib).i01([1 2],is) = inf(ib).i01([2 1],is);
   inf(ib).x01([1 2],is) = inf(ib).x01([2 1],is);
   inf(ib).y01([1 2],is) = inf(ib).y01([2 1],is);
   inf(ib).chk([1 2],is) = inf(ib).chk([2 1],is);

   if ~isempty(inf(ib).xyd{is})
       inf(ib).xyd{is} = inf(ib).xyd{is}(:,end:-1:1);
   end

   msg = sprintf('%s\nFlip Data of corresponding Segment.\n',txt);
   warning(msg);

end

%-----------------------------------------------------------
% Set CheckValue

ibs = ib + i * is;

inf(ib).chk(i1,is) = ii + i * jj;

inf(ii).chk(i0,jj) = ibs;

%-----------------------------------------------------------
% Check, if corresponding segment allready in History

if any( hst == ibs )
   ok = 1;
   return
end

%-----------------------------------------------------------
% Append Data

inf(ib).cls(is) = NaN;   % Data used from this segment

if inf(ib).ins(is)
   xy = inf(ib).xyd{is};
%%   ind = ( 1+mode : size(xy,2)-(~mode) );
%%   xy  = xy(:,ind);
else
   xy = cat(1,inf(ib).x01(i1,is),inf(ib).y01(i1,is));
end

if mode
   xyd = cat( 2 , xyd , xy  );
   hst = cat( 1 , hst , ibs );
else
   xyd = cat( 2 , xy  , xyd );
   hst = cat( 1 , ibs , hst );
end

cnt = cnt + 1;

if cnt  > cnf.rcl
   ok = 9;
   msg = sprintf('RecursionLimit of %.0f reached.',cnf.rcl);
   warning(msg);
   return
end

%************************************************************
% Recurse

[xyd,ok] = rec_seg( ib , is , xyd , mode , cnt );


%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xy,ins] = in_area(xy,xl,yl,off);

if isempty(xy)
   xy  = zeros(2,0);
   ins = 0;
   return
end

in = ( ( xl(1) <= xy(1,:) ) & ( xy(1,:) <= xl(2) ) & ...
       ( yl(1) <= xy(2,:) ) & ( xy(2,:) <= yl(2) )       );

ins = any(in);

%%% if ~ins
%%%     xy = xy(:,[1 end]);  % Area Crossing possible
%%% end

% Get Data inside Limits

if all(in)
   return
end

ll = [ xl(1) yl(1) xl(2) yl(2) ] + [ -1  -1  1  1 ]*off;
%      >=    >=    <=    <=

for ii = 1 : size(ll,2)

    jj = mod(ii-1,2) + 1;

    op = 1 - 2 * ( ii > 2 );

    if ins

       in = ( op*xy(jj,:) >= op*ll(ii) );

       in( find( diff(in) ==  1 ) + 0 ) = 1;
       in( find( diff(in) == -1 ) + 1 ) = 1;

       in([1 end]) = 1;

       if ~all(in)
           in = find(in);
           xy = xy(:,in);
       end

    end

    xy(jj,:) = op*max(op*xy(jj,:),op*ll(ii));

end

dd = ( abs(diff(xy,1,2)) < 1e-10 );

dd = ( sum(dd,1) == 2 );

if ~any(dd)
    return
end

xy( : , find(dd) + 1 ) = [];

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

