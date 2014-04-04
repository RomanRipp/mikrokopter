function [x,y,z,msg] = get_topo(cnf,xl,yl,stride,smt,step);

% GET_TOPO    Extract Topography-Data from NetCDF-File
%              using LOAD_CDF
%
%
% [ X , Y , Z , Msg ] = GET_TOPO( Config , XLim , YLim ); GridMode
%
% [ X , Y , Z , Msg ] = GET_TOPO( Config , XVec , YVec ); GridMode
%
% [ XYZD , Msg ] = GET_TOPO( Config , Pos1 , Pos2 );      SectionMode
%
% [ XYZD , Msg ] = GET_TOPO( Config , XVec , YVec );      SectionMode
%
%----------------------------------------------------------
% Input 
%
% Config      Name of MAT-File and optional VariableNames for X Y and Z
% 
% Config = { MAT_File [XVar] [YVar] [ZVar] }
%  
%            defaults:   'x'    'y'    'z'
%
%             X == Longitude[deg], Y == Latitude[deg]
% ------
%
% Config      Name of ConfigFile for NetCDF-File,
%              and optional DataVariable in CNF-File 
%              and a different NetCDF-FileName, to using LOAD_CDF
% 
% Config = { ConfigFile [DataVar] [FileName] }
%
%             A ConfigFile contains following KeyWords:
%
%               # FileName   % NetCDF-FileName
%               # XDim       % Name of Dimension for X
%               # YDim       % Name of Dimension for Y
%               # XVar       % Name of Variable  for X
%               # YVar       % Name of Variable  for Y
%               # <DataVar>  % Name of Variable  for TopographyData
%
%             default: DataVar = 'Topo'
%
%              use CNF_CDF( CDF_FileName )  to create this ConfigFile
%
% XLim        Range for X (Longitude)  [ XMin  XMax ] or Vector for X (Lon)
% YLim        Range for Y (Latitude)   [ YMin  YMax ] or Vector for Y (Lat)
%   
%  or use a single Input:  GET_TOPO( Config , [ XMin  XMax YMin YMax ] )
%
%
% Pos1        StartPosition of Section [ Lon1  Lat1 ]
% Pos2          EndPosition of Section [ Lon2  Lat2 ]
% 
%  or use a single Input:  GET_TOPO( Config , [ Lon1 Lat1 Lon2 Lat2 ] )
%
%
%----------------------------------------------------------
% Output:
%
%  X          Vector for X (Longitude)  [ Nx by 1 ]
%  Y          Vector for Y (Latitude)   [ Ny by 1 ]
%  Z          Matrice for Topography    [ Ny by Nx ]
%
%
%  XYZD       SectionData:  [ Lon Lat Topo Distance ]  [ N by 4 ]
%             Distance by SODANO
%
%----------------------------------------------------------
% additional Inputs:
%
% [ X , Y , Z ] = GET_TOPO( ... , Stride , SmoothWindow , Step );
%  
% Stride       Stride to extract Data form the NetCDF-File
%                 default:   Stride = 1
%                 example:   Stride = 2  ==> extract every second Value
%                                             from NetCDF-File (half size of Data)
%
% SmoothWindow Window to smooth the extracted TopographyData, using MEANIND2
%
%              Smooth | [ SmoothLand  SmoothWater ]
%
% Step         Step to (low)sample the extracted, (smoothed) Data
% 
%                 default:   Step = 1
%                 example:   Step = 2  ==> takes every second Value
%                                             from extracted (smoothed) Data
%
%
%----------------------------------------------------------
%
% see also:  LOAD_CDF, CNF_CDF, GEODIST
%

Nin  = nargin;
Nout = nargout;

x = [];
y = [];
z = [];

msg = '';

if Nin < 4, stride = []; end
if Nin < 5, smt    = []; end
if Nin < 6, step   = []; end


if isempty(stride), stride = 1; end
if isempty(smt),    smt    = 1; end
if isempty(step),   step   = 1; end

err = ~( mod(Nout,2) == 0 );  % True for Error
sct =  ( Nout <= 2 );         % True for SectionMode !!!

[msg,cnf,opt] = chkcnf(cnf);

if ~isempty(msg)
    msg = sprintf('Invalid Configuration.\n%s',msg);
    if err, error(msg), else, if sct, y = msg; end, return, end
end

if Nin < 3
   if prod(size(xl)) == 4
      yl = xl([3 4]);
      xl = xl([1 2]);
   else
      msg = 'XLim,YLim or Area required.';
      if err, error(msg), else, if sct, y = msg; end, return, end
   end
end

%----------------------------------------------
% Check for SectionMode

sct = ( Nout <= 2 );  % True for SectionMode !!!

if sct

   px = prod(size(xl));
   py = prod(size(yl));

   isp = ( ( px == 2 ) & ( py == 2 ) );  % True if 2 Point defined

   if isp

      % Get Section between 2 Points

      p1 = xl;
      p2 = yl;
      xl = [ p1(1)  p2(1) ];
      yl = [ p1(2)  p2(2) ];

      [msg,x,y] = loadcnf(cnf,opt,xl,yl,stride);

      if ~isempty(msg)
          if err, error(msg), else, y = msg; return, end
      end

      sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
      sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

      if ~isequal(sx,sy)
         if     vx & vy
            x = ones(py,1) * x(:)';
            y = y(:) * ones(1,px);
         elseif vx & ( px == sy(2) )
            x = ones(sy(1),1) * x(:)';
         elseif vy & ( py == sx(1) )
            y = y(:) * ones(1,sx(2));
         else
            msg = 'Size of X and Y must be agree.';
            if err, error(msg), else, if sct, y = msg; end, return, end
         end
      end

      if size(x,2) == 1
         dx = NaN;
      else
          x = abs(diff( x , 1 , 2 ));
          x = x .* cos( (y(:,1:(end-1))+y(:,2:end))/2 * pi/180 );
         dx = [ min(x(:))  mean(x(:)) ];
         dx = dx([1 2 2]) .* [ 1  1  1/2 ];
         dx = dx( 1 + 2 * ( dx(1) < 1/5 * dx(2) ) );
      end

      if size(y,1) == 1
         dy = NaN;
      else
          y = abs(diff(y,1,1));
         dy = [ min(y(:))  mean(y(:)) ];
         dy = dy([1 2 2]) .* [ 1  1  1/2 ];
         dy = dy( 1 + 2 * ( dy(1) < 1/5 * dy(2) ) );
      end

      d = min( dx , dy  );

      if isnan(d)

         n = 100;

      else

         d = d * 60 * 1852;

         dst = geodist( p1([2 1]) , p2([2 1]) );

         n = ceil( dst / (d+(d==0)) );
         n = n + ( n == 2 );   % Min 2 Points

      end

      d = NaN * zeros(4,n);

      [d(4,:),d(2,:),d(1,:)] = geodist( p1([2 1]) , p2([2 1]) , n );

      xl = [ d(1,:) xl ];
      yl = [ d(2,:) yl ];

      xl = [ min(xl) max(xl) ];
      yl = [ min(yl) max(yl) ];

   else

      if ~isequal(px,py)
          error('Size of Vectors X and Y must be agree in SectionMode.'),
      end
  
      d = cat(2,xl(:),yl(:),zeros(px,2));

      ii = ( 2 : px );

%sum(geodist(d(ii-1,[2 1]),d(ii,[2 1])),1)
%dd = d(ii,[1 2])-d(ii-1,[1 2]);
%dd(:,1) = dd(:,1) .* cos((d(ii,2)+d(ii-1,2))/2*pi/180);
%sum(sqrt(sum(dd.^2,2)),1)*60*1852
 
      d(ii,4) = cumsum(geodist(d(ii-1,[2 1]),d(ii,[2 1])),1);
     
      xl = [ min(d(:,1))  max(d(:,1)) ];
      yl = [ min(d(:,2))  max(d(:,2)) ];

      d = permute(d,[2 1]);

   end

end

%----------------------------------------------
% Load Data

fprintf(1,'\nGET_TOPO: Read Data from %s',cnf);

[msg,x,y,z] = loadcnf(cnf,opt,xl,yl,stride);

if ~isempty(msg)
    fprintf(1,' error\n%s',msg);
    if err, error(msg), else, if sct, y = msg; end, return, end
end

fprintf(1,' ok\n\n');

%----------------------------------------------
% Smooth

if any( smt > 1 )

  fprintf(1,'Smooth Data ...');

  if prod(size(smt)) >= 2

     z = min(-1,meanind2(z,smt(2),'cos')) .* ( z <  0 ) + ...
         max( 1,meanind2(z,smt(1),'cos')) .* ( z >= 0 );

  else

     is_water = ( z < 0 );

     z = meanind2(z,smt,1);

     z(find( is_water)) = min( z(find( is_water)) , -1 );
     z(find(~is_water)) = max( z(find(~is_water)) ,  1 );

  end

  fprintf(1,' done\n\n');

end

%----------------------------------------------
% LowSample

if step > 1

  sx = size(x); sy = size(y); sz = size(z);

  ind = { ':'  ':' };
  for ii = 1 : 2
      ind{ii} = ( 1 : step : sz(ii) );
  end

  x = x( ind{ 2-isequal(sx,sz) : 2 } );

  y = y( ind{ 1 : 1+isequal(sy,sz) } );

  z = z(ind{:});

end

%----------------------------------------------

if ~sct
    return
end

%----------------------------------------------
% SectionMode

if     prod(size(x)) == 1
    d(3,:) = interp1(y,z,d(2,:));
elseif prod(size(y)) == 1
    d(3,:) = interp1(x,z,d(1,:));
else
    d(3,:) = interp2(x,y,z,d(1,:),d(2,:));
end

x = permute(d,[2 1]);

y = msg;

%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,x,y,z] = loadcnf(cnf,opt,xl,yl,stride);


Nout = nargout;

msg = '';
x = [];
y = [];
z = [];

noz = ( Nout < 4 );

%**********************************************************
if prod(size(opt)) == 3       % MAT_File
%**********************************************************

  if noz
     opt = opt([1 2]);
  end

  try
     z = load(cnf,opt{:});
  catch
     msg = lasterr;
  end

  if ~isempty(msg)
      msg = sprintf('Error load Data using LOAD( %s )\n%s',cnf,msg);
      return
  end


  x = getfield(z,opt{1});
  y = getfield(z,opt{2});

  if noz
     return
  end

  z = getfield(z,opt{3});

  [msg,x,y,z] = chkgrid(x,y,z);

  if ~isempty(msg)
      msg = sprintf('Invalid Data in MAT-File"%s".\n%s',cnf,msg);
      return
  end

  if isempty(z) | ( stride == 1 )
     return
  end

  sx = size(x); sy = size(y); sz = size(z);

  ind = { ':'  ':' };
  for ii = 1 : 2
      ind{ii} = ( 1 : stride : sz(ii) );
  end

  x = x( ind{ 2-isequal(sx,sz) : 2 } );

  y = y( ind{ 1 : 1+isequal(sy,sz) } );

  z = z(ind{:});

%**********************************************************
else                   % NetCDF-Config
%**********************************************************

  if noz
     mode = {};
  else
     mode = { 'Data'  opt(1) };
  end

  if noz | ( isequal(size(xl),[1 2]) & isequal(size(yl),[1 2]) )
     mode = cat(2,mode,{'orig'});
  end

  [msg,z] = load_cdf( cnf , ...
                     'xvar' , xl , ...
                     'xdim' , stride  , ... 
                     'yvar' , yl , ...
                     'ydim' , stride , 'file' , opt{2} , mode{:} );

   if ~isempty(msg)
       z = [];
       return
   end

   x = z{1}';
   y = z{2};

   if noz
      return
   end

   z = z{5};

  [msg,x,y,z] = chkgrid(x,y,z);

  if ~isempty(msg)
      msg = sprintf('Invalid Data in NetCDF-File"%s".\n%s',cnf,msg);
      return
  end


%**********************************************************
end
%**********************************************************

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,x,y,z] = chkgrid(x,y,z);

% Checks / Sort Grid
%

msg = '';

%------------------------------------------------------------

sx = size(x); px = prod(sx); vx = ( px == max(sx) );
sy = size(y); py = prod(sy); vy = ( py == max(sy) );
sz = size(z); pz = prod(sz);

if all( [ px  py  pz ] == 0 )
   return
end

%------------------------------------------------------------

msg = cell(0,1);

%------------------------------------------------------------
% Check Dimensions

if ~all( [ size(sx,2) size(sy,2) size(sz,2) ] == 2 )

    msg = cat( 1 , msg , {'XYZ-Variables must have max. 2 Dimensions.'} );

elseif ~isequal(sx,sy,sz)

    ex = isequal(sx,sz);
    ey = isequal(sy,sz);

    vx = ( ~ex & vx & ( px == sz(2) ) );
    vy = ( ~ey & vy & ( py == sz(1) ) );
    
    if ~( ex | vx ) & ( px == 2 ) & ( sz(2) > 2 )
        x  = linspace(x(1),x(2),sz(2));
        sx = [ 1  sz(2) ];
        px = sz(2);
        vx = 1;
    end

    if ~( ey | vy ) & ( py == 2 ) & ( sz(1) > 2 )
        y  = linspace(y(1),y(2),sz(1));
        sy = [ sz(1)  1 ];
        py = sz(1);
        vy = 1;
    end

    if ~( ( ex | vx ) & ( ey | vy ) )

        msg = cat( 1 , msg , {'MatrixDimensions of XYZ must be agree.'} );

    elseif any( sz([1 2]) == 1 )

        msg = cat( 1 , msg , {'MatrixDimensions of Z must be nonsingular.'} );

    else

        %----------------------------------------------
        % Sort Grid

        %----------------------------------------------

        if vx
           x = x(:)';
           if px > 1
              sgn = sign(diff(x));
              if ( sgn(1) < 0 ) & all( sgn == sgn(1) )
                 ix = ( px : -1 : 1 );
                  x = x(ix);
                  z = z(:,ix);
                  if isequal(sy,sz)
                     y = y(:,ix);
                  end
              end
           end
        end

        %----------------------------------------------

        if vy
           y = y(:);
           if py > 1
              sgn = sign(diff(y));
              if ( sgn(1) < 0 ) & all( sgn == sgn(1) )
                 iy = ( py : -1 : 1 );
                  y = y(iy);
                  z = z(iy,:);
                  if isequal(sx,sz)
                     x = x(iy,:);
                  end
              end
           end
        end

        %----------------------------------------------
        
    end
   
end

%------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
else
    msg = '';
end


%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf,opt] = chkcnf(cnf);

% Check for NetCDF- or MAT-Configuration
%
% { MAT_File [ XVar  YVar  ZVar ] }
% { CNF_File [ DataVar  NC_File ] }


opt = {};

ok = chkstr(cnf);
if ok
   cnf = {cnf};
else
   if ischar(cnf)
      cnf = cellstr(cnf);
   elseif iscell(cnf);
      pc = prod(size(cnf));
      ok = chkcstr(cnf(1:min(3,pc)),0);    % Check for Config
   end
end

if ~ok
    msg = 'Configuration must be a String or CellArray of Strings.';
    return
end

pc = prod(size(cnf));

if pc >= 2
   opt = cnf(2:pc);
   opt = opt(:)';
end

cnf = cnf{1};

%------------------------------------
% Check for MAT_File

typ = { 'NetCDF'  'MAT' };

try
   v = whos('-file',cnf);
catch
   v = [];
end

mat = ~isempty(v);

typ = typ{ 1 + mat };

msg = cell(0,1);

%**********************************************************
if mat                 % MAT_File
%**********************************************************

 if ~( isempty(opt) | ( ( prod(size(opt)) >= 3 ) & chkcstr(opt,0) ) )

    msg = cat( 1 , msg , {'MAT-Config must be { MAT_File [ XVar  YVar  ZVar ] }.'} );

 else

    if isempty(opt)
       vr = { 'x' 'y' 'z' };
    else
       vr = opt(1:3);
    end

    cfg = { v.name  };
    cls = { v.class };
    siz = { v.size  };

    %-----------------------------------------------------------
    % Check for XYZ-Variables

    nv = size(vr,2);

    sv = NaN * zeros(nv,3);

    ok = zeros(nv,1);

    for ii = 1 : nv

        jj = strcmp(cfg,vr{ii});

        if any(jj)

           jj = find(jj);

           i3 = ( ii >= 3 );

           s2 = size(siz{jj},2);
           if ~any( s2 == 2+[0 i3] )  % DimensionNumber
               sv(ii,:) = -1;
           else
               sv(ii,1:s2) = siz{jj};
           end

           ck = ( strcmp(cls{jj},{'double' 'uint8'}) & [ 1  i3 ] );

           ok(ii) = 1 - 2 * ~any(ck);

        end

    end

    chk = ( 1 : nv );                        % X Y Z requested


       cc = { 'X' 'Y' 'Z' };

       mm = cell(1,3); mm(:) = {''};

       if  any( ok(chk) == 0 )
           jj = chk( find( ok(chk) == 0 ) );
            m = cat( 1 , cc(jj) , vr(jj) );
            m = sprintf(' %s:"%s" ',m{:});
            n = sprintf(' "%s" ',cfg{:});
            m = sprintf('Missing Variables in ##: %s',m);
            mm{1} = sprintf('%s\n   Available: %s',m,n);
       end

       if  any( ok(chk) == -1 )
           jj = chk( find( ok(chk) == -1 ) );
            m = cat( 1 , cc(jj) , vr(jj) );
            m = sprintf(' %s:"%s" ',m{:});
            mm{2} = sprintf('Variables in ## must be numeric: %s',m);
       end

       if  any( sv(chk,1) == -1 )
           jj = chk( find( sv(chk,1) == -1 ) );
            m = cat( 1 , cc(jj) , vr(jj) );
            m = sprintf(' %s:"%s" ',m{:});
            mm{3} = sprintf('Variables in ## must have 2 Dimensions: %s',m);
       end
 
       jj = ~strcmp(mm,'');
       if any(jj)
          jj = find(jj);
          mm = sprintf('%s\n',mm{jj});
       else
          mm = '';
       end


    %----------------------------------------------------------------------
    if ~isempty(mm)
    %----------------------------------------------------------------------

        msg = cat( 1 , msg , {mm} );

        ok(:) = 0;

    %----------------------------------------------------------------------
    else
    %----------------------------------------------------------------------

        ok = ( ok == 1 );

        xy = ( ok(1) & ok(2) );  % True for X AND Y
        ok([1 2]) = xy;


        jj = isnan(sv(:,3));
        if any(jj)
           jj = find(jj);
           sv(jj,3) = 1;
        end

        pv = prod(sv,2);
        vc = ( pv == max(sv,[],2) );

        %-----------------------------------------------------------------------------------
        % Check Size of Z with XY

        if ok(3) & ( pv(3) == 1 )
           m = sprintf('Z-Variable" %s" in ## must be nonsingular.',vr{3});
           msg = cat( 1 , msg , {m} );
        end

        if xy & ok(3) & ~( pv(3) == 0 )

           zk = ( ( isequal(sv(1,[1 2]),sv(3,[1 2])) | ( vc(1) & ( pv(1) == sv(3,2) ) ) ) & ...
                  ( isequal(sv(2,[1 2]),sv(3,[1 2])) | ( vc(2) & ( pv(2) == sv(3,1) ) ) )         );
           if ~zk
               if fc
                  ok(1:3) = 0;
               else
                  m = 'MatrixDimensions of XYZ-Variables in ## must be agree.';
                  msg = cat( 1 , msg , {m} );
               end
           end

        end

    %----------------------------------------------------------------------
    end % mm
    %----------------------------------------------------------------------

    opt = vr;

    jj = ( ok == 0 );
    if any(jj)
       jj = find(jj);
       opt(jj) = {''};    % Set Variables Empty !!!
    end


  end

%**********************************************************
else                   % NetCDF-Config
%**********************************************************

 if isempty(cnf)

    msg = cat( 1 , msg , {'Empty ConfigFile in NetCDF-Config.'} );

 else

   if isempty(opt)
      opt = { 'Topo'  '' };
   else
      opt = cat( 2 , opt(:)' , {''} );
   end

   opt = opt(1:2);

   try
       [m,cfg] = cnf_cdf(cnf,'check',opt{2});
   catch
        m = sprintf('Error call CNF_CDF to check ##.\n%s',lasterr);
   end

   if ~isempty(m)

       m = sprintf( 'Invalid ##.\n%s' , m );
       msg = cat( 1 , msg , {m} );

   else

       cfg  = cfg{2,1}(:,1);   % DataVariables

       if isempty(cfg)
          m = 'No DataVariables defined in ##.';
          msg = cat( 1 , msg , {m} );
       else
          if isempty(opt{1})
             opt{1} = cfg{1};   % Use First DataVariable !!!
          elseif ~any(strcmp(cfg,opt{1}))
             n = sprintf(' "%s" ',cfg{:});
             m = sprintf('Can''t find DataVariable "%s" in ##.',opt{1});
             m = sprintf('%s\n   Available: %s',m,n);
             msg = cat( 1 , msg , {m} );
          end
       end


   end

 end

%**********************************************************
end      % MAT | NetCDF
%**********************************************************


if ~isempty(msg)

    tp = '';
    if strcmp(typ,'NetCDF')
       tp = '-Config';
    end

    m = sprintf('%s%s-File "%s"',typ,tp,cnf);

    msg = strrep(msg,'##',m);

    msg = sprintf('%s\n',msg{:});

else

    msg = '';

end

