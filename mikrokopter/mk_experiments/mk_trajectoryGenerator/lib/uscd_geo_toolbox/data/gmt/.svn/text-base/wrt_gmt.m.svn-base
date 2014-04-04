function msg = wrt_gmt(xy,file,bsi,att,varargin);

% WRT_GMT  Write Lines to binned NetCDF-Files (GMT-Format)
%
% Msg = WRT_GMT( XY , File , BinSize , GlobalAttributes )
%
%------------------------------------------------------------------
%
% XY:      Matrice with [ X  Y ] Coordinates (2 Columns)
%           or optional [ X Y Level ]
%          Segments are separated by NaN-Rows in [ X  Y ]
% 
% File:    FileName for NetCDF-File
%          If the FileName has no Extension, ".cdf" will append
%
% BinSize: Size of Bins in degree
%          if no Value is given, BinSize will detect automaticly
%
% Attributes: { Name Type Value } Global Attributes for NetCDF-File
%            
%------------------------------------------------------------------
%
% WRT_GMT creates 2 temporary Files: "<File>.bin" and "<File>.dat"
%
% The temporary Files will removed, if the NetCDF-File is written
%  successfull. 
%
% If the Value of BinSize is NEGATIVE, the temporary File will not 
%  removed and no NetCDF-File will created.
%
% Create a temporary File only:
%
%    WRT_GMT( XY , File , -BinSize )  
%
% To append Data to the temporary Files use: 
%
%    WRT_GMT( XY , File , 'append' )  
%
% To append Data to the temporary Files and 
%    write the NetCDF-File use: 
%
%   WRT_GMT( XY , File , 'write' )  
% 
%------------------------------------------------------------------
% Example to create a NetCDF-DataSet for Europe
%
%  area = [-12 45 35 72];          % Europe
%  file = 'europe.cdf';
%  att = { 'title' 'char' ...
%          'CoastLines, Rivers, Political Boundaries of Europe' };
%
%  ini = { 'sc'   'rc'    'bc'     % Coast River Bound
%            1   [1 2 9]   1       % Level to read
%           'k'    'b'    'r' };   % Colors
%
%  wrt_gmt([],file,-20);           % Create the TempFiles
%
%  for ii = 1 : size(ini,2)
%      xy = read_gmt(ini{1,ii},area,ini{2,ii});
%      wrt_gmt(cat(2,xy,ii+0*xy(:,1)),file,'append');
%  end
%
%  wrt_gmt([],file,'write',att);       % Write NetCDF-File
%
%  figure, axis(area), hold on
%  
%  for ii = 1 : size(ini,2)
%      xy = read_gmt(file,area,ii);
%      plot(xy(:,1),xy(:,2),ini{3,ii});
%  end
%
%  mercator
%
%------------------------------------------------------------------
%
% See also: READ_GMT, SHORE_GMT, WRITE_CDF, CREATE_CDF, LOOK_CDF
%
%------------------------------------------------------------------
% Binned NetCDF-File-Format (GMT):
%
%   4 Dimensions:
%  --------------
%  1)            Dimension_of_scalar = 1       
%  2)        Dimension_of_bin_arrays = <162> = 360/BinSize * 180/BinSize    
%  3)    Dimension_of_segment_arrays = <547>       
%  4)      Dimension_of_point_arrays = <2123>      
%
%  12 Variables:
%  -------------
%  1)    Bin_size_in_minutes( Dimension_of_scalar ),   long
%  2)    Bin_size_in_degree( Dimension_of_scalar ),   long
%  3)    N_bins_in_360_longitude_range( Dimension_of_scalar ),   long
%  4)    N_bins_in_180_latitude_range( Dimension_of_scalar ),   long
%  5)    N_bins_in_file( Dimension_of_scalar ),   long
%  6)    N_segments_in_file( Dimension_of_scalar ),   long
%  7)    N_points_in_file( Dimension_of_scalar ),   long
%  8)    Id_of_first_segment_in_a_bin( Dimension_of_bin_arrays ),   long
%  9)    N_segments_in_a_bin( Dimension_of_bin_arrays ),   short
% 10)    Hierarchial_level_of_a_segment( Dimension_of_segment_arrays ),   short
% 11)    Id_of_first_point_in_a_segment( Dimension_of_segment_arrays ),   long
% 12)    Relative_longitude_from_SW_corner_of_bin( Dimension_of_point_arrays ),   short
%         units: "1/65535 of <Bin_size> degrees relative to SW corner of bin"
% 13)    Relative_latitude_from_SW_corner_of_bin( Dimension_of_point_arrays ),   short
%         units: "1/65535 of <Bin_size> degrees relative to SW corner of bin"
%

msg = '';
nl  = char(10);

Nin  = nargin;
Nout = nargout;

if Nin < 2
   msg = 'Not enough InputArguments.';
   return
end

lid = 1;  % STDOUT

%--------------------------------------------------------------
% Check XY

if ~isempty(xy)

    s  = size(xy);
    ok = ( isnumeric(xy) & ...
           ( prod(s) == s(1)*s(2) ) & any( s(2) == [ 2  3 ] ) );

    if ok
       ok = sum( isnan(xy(:,[1 2])) , 2 );
       ok = ( ~any( ok == 1 )  &  ~all( ok == 2 ) );
    end

    if ~ok
        msg = [ 'XY must be a 2- or 3-Column-Matrice [X Y] or [X Y L],'  nl ...
                ' with finite Values in X and Y, separated by NaN-Rows.' ];
    else
       xy = check_dat(xy);
    end

end

%--------------------------------------------------------------
% Check File

if ~chkstr(file,1)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'File must be a String' ];
else
   [p,n,e] = fileparts(file);

   if isempty(p)
      p = pwd;
   elseif ~( exist(p,'dir') == 7 )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              sprintf('Directory "%s" doesn''t exist.',p) ];
   end
  
   if isempty(n)
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'FileName not specified.' ];
   else
      if isempty(e)
         e = '.cdf';
      end
      file = cat(2,fullfile(p,n),e);
   end
end

%--------------------------------------------------------------
% Check BinSize

mode = 'new';

if Nin < 3
   bsi = [];
elseif chkstr(bsi,1)
   mode = bsi;
   bsi  = [];
elseif ~isempty(bsi)
   ok = ( prod(size(bsi)) == 1 );
   if ok
      ok = ( isfinite(bsi) & ( 0 < abs(bsi) ) & ( abs(bsi) <= 180 ) & ...
             ( mod( 180 , bsi ) == 0 ) );
   end
   if ~ok
       msg = [ msg nl(1:(end*(~isempty(msg)))) ...
               'Binsize must be a single numeric in range ( 0 180 ],' nl ...
               '180 must be a multiple of BinSize' ];
   end
end

%--------------------------------------------------------------
% Check Attributes

if Nin < 4
   att = [];
elseif ~isempty(att)
   ok = ( iscell(att) & ( size(att,2) == 3 ) & ( ndims(att) == 2 ) );
   if ok
      ok = chkcstr(att(:,[1 2]));
   end
   if ~ok
       msg = [ msg nl(1:(end*(~isempty(msg)))) ...
        'Attributes must be a 3-Column CellArray: { Name Type Value }.' ];
   end
end

if Nin > 4
   [m,ds,vs,as,f,dm,vr,at] = create_cdf('',varargin{:});
   if ~isempty(m)
       msg = [ msg nl(1:(end*(~isempty(msg)))) ...
               'Invalid following Inputs DIM/VAR/ATT.' ];
   elseif isempty(vr)
       dm = cell(0,2);
       vr = cell(0,5);
       at = cell(0,1);
   else
       dm = dm(:,1:2);
       v1 = size(vr,1);
       v2 = size(vr,2);
       at = at(1:min(end,v1),:);
       v2 = size(vr,2);
       if v2 > 5
          vr = vr(:,1:5);
       elseif v2 < 5
          a2 = cell(v1,5-v2);
          a2(:) = {[]};
          vr = cat(2,vr,a2);
       end
       for ii = 1 : v1
           if isnumeric(vr{ii,4})
              vr{ii,4} = dm(vr{ii,4}+1,1);
           end
       end
   end
else
   dm = cell(0,2);
   vr = cell(0,5);
   at = cell(0,1);
end


%--------------------------------------------------------------

if ~isempty(msg)
    if Nout == 0
       fprintf(lid,'\nInvalid Inputs.\n%s\n\n',msg);
       clear msg
    end
    return
end

fprintf(lid,'\n');

%***************************************************************
% Precision for Data

prc = 16;

%***************************************************************
% Check for Append/Write

BIN_File = cat(2,file,'.bin');
DAT_File = cat(2,file,'.dat');

no_write = 0;

mode = lower(mode(1));

if any(strcmp(mode,{'a' 'w'}))

   [msg,bsi] = check_file(BIN_File,DAT_File,prc,lid);

   if ~isempty(msg)
      if Nout == 0
         fprintf(lid,'%s\n\n',msg);
         clear msg
      end
      return
   end

   md = 'a';

   if isempty(xy) &  strcmp(mode,'a')
      return
   end

else

   md = 'w+';

   if ~isempty(bsi)
       no_write = ( bsi < 0 );
       bsi = abs(bsi);
   end

end

%***************************************************************
% AutoDetecting of BinSize

% GMT:  bsi = [ 60  120  300  600  1200 ]
%       
%       bsi/sqrt(median(dxy)) == 60*30
%

if isempty(bsi)

   if isempty(xy)
      msg = 'Nothing to do.';
      if Nout == 0
         fprintf(lid,'%s\n\n',msg);
         clear msg
      end
      return
   end

    n = size(xy,1);

   dd = sqrt( ( ( xy(2:n,1) - xy(1:n-1,1) ) .* ...
                cos( ( xy(2:n,2) + xy(1:n-1,2) ) / 2 * pi/180 ) ).^2 + ...
               ( xy(2:n,2) - xy(1:n-1,2) ) .^2 );

   dd = median(dd(find(~isnan(dd))));

   bsi = ceil( 30 * sqrt(dd) );

   while ~( mod(180,bsi) == 0 )
          bsi = bsi - 1;
   end

   fprintf(lid,'Autodetected BinSize: %.0f\n\n',bsi);

end

%***************************************************************
% BinResolution

px = 360;
nx = px / bsi;

py = 180;
ny = py / bsi;

%***************************************************************
% Write Data

if ~isempty(xy) | strcmp(md,'w+')

    [msg,ns,nb,np] = wrt_dat(BIN_File,DAT_File,md,prc,xy,bsi,px,nx,py,ny,lid);

    if ~isempty(msg)
        if Nout == 0
           fprintf(lid,'%s\n\n',msg);
           clear msg
        end
        return
    end

end

if strcmp(mode,'a') | no_write
   fprintf(lid,'\n');
   return
end

%***************************************************************
% Read Info

[msg,bsi,nx,ny,bin,did] = check_file(BIN_File,DAT_File,prc,lid);

if ~isempty(msg)
    if Nout == 0
       fprintf(lid,'%s\n\n',msg);
       clear msg
    end
    return
end

%***************************************************************
% Write NetCDF-File

if isempty(bin)
   msg = 'Nothing to do.';
else
   msg = write_file(file,did,bsi,nx,ny,bin,prc,att,dm,vr,at,lid);
   if ~isempty(msg)
       if Nout == 0
          fprintf(lid,'%s\n\n',msg);
          clear msg
       end
       return
   end
end

%***************************************************************
% Remove TempFiles


fprintf(lid,'Delete BIN-File: %s',BIN_File);
try
   delete(BIN_File); 
catch
   fprintf(lid,'   %s','error');
end
fprintf(lid,'\n');


fprintf(lid,'Delete DAT-File: %s',DAT_File);
try
   delete(DAT_File); 
catch
   fprintf(lid,'   %s','error');
end
fprintf(lid,'\n\n');

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = write_file(file,did,bsi,nx,ny,bin,prc,glb,dm,vr,at,lid);

msg = '';

nl = char(10);

[uprc,sprc,sc] = get_prec(prc);

[bin,sgm,dat,dim,var,att] = file_def(bsi,nx,ny,bin,sc);

if ~isempty(vr)
    dim = cat(1,dim,dm);
    var = cat(1,var,vr);
    att = cat(1,att,at);
end

if ~isempty(glb)
    att = cat(1,att,{glb});
end

fprintf(lid,'\nWrite NetCDF-File: %s',file);

[msg1,msg2,vm,ds,vs,as] = write_cdf(file,dim,var,att);

if ~( isempty(msg1) & isempty(msg2) )
   fprintf(lid,'   %s\n','error');
   msg = [ msg1 nl(1:(end*(~isempty(msg1)))) msg2 ];
   msg = sprintf('Error writing NetCDF-File "%s".\n%s',file,msg);
   return
end

fprintf(lid,'\n');

[fid,stat] = ncmex('open',file,'write');

if ~( ( fid > 0 ) & ( stat == 0 ) )
   msg = sprintf('%sError open NetCDF-File: %s',msg0,file);
   return
end
  
%------------------------------------------------------------------------
% Get [ Start Length ] of Segments      

form = '\r   %10.0f Points of %8.0f Segments in %6.0f Bins';

nb = size(sgm,1);

fprintf(lid,'\nWrite %.0f Segments (%.0f Points) in %.0f Bins ... \n',sum(sgm(:,2)),sum(bin(:,3)),nb);

np = 0;
ns = 0;

xid = strmatch( 'Relative_longitude' , var(:,1) ) - 1;
yid = strmatch( 'Relative_latitude'  , var(:,1) ) - 1;

for ii = 1 : nb

    ind = sgm(ii,1) - 1 + ( 1 : sgm(ii,2) );

    for jj = 1 : sgm(ii,2)

        kk = ind(jj);
 
        ip = dat(kk,1);
        lp = dat(kk,2);

        [x,y] = read_dat(did,prc,ip,lp);

        if ~isempty(x)
            ncmex('varput',fid,xid,bin(kk,4)-1,bin(kk,3),x);
            ncmex('varput',fid,yid,bin(kk,4)-1,bin(kk,3),y);
            ns = ns + 1;
            np = np + bin(kk,3);
            fprintf(lid,form,np,ns,ii);
        end

    end

end

ncmex('close',fid);

fclose(did);

fprintf(lid,'\n\n');

if ns < sum(sgm(:,2))
   msg = sprintf('Error read Data of %.0f Segments.',sum(sgm(:,2))-ns);
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,ns,nb,np] = wrt_dat(BIN_File,DAT_File,mode,prc,xy,bsi,px,nx,py,ny,lid)

msg = '';

ns = 0;
nb = 0;
np = 0;

%***************************************************************
% Open TMP-Files for writing

%---------------------------------------------------------------
% BIN_FILE: [ bsi  nx   ny  
%             bnr1 lev1 siz1 
%             bnr2 lev2 siz2 ... ]
 
fprintf(lid,'Open BIN_File for writing: %s',BIN_File);

bid = fopen(BIN_File,mode);

if bid == -1
   fprintf(lid,'   %s\n','error');
   msg = sprintf('Can''t open BIN-File "%s" for writing.',BIN_File);
   return
end

fprintf(lid,'\n');

%---------------------------------------------------------------
% DAT_FILE: [ bsi  nx   ny  
%             x1   y1   x2  y2 ... ]
%

fprintf(lid,'Open DAT_File for writing: %s',DAT_File);

did = fopen(DAT_File,mode);

if bid == -1

   fprintf(lid,' %s\n','error');

   if strcmp(mode,'w+')
      if ~( bid == -1 )
          fclose(bid);
      end
      fprintf(lid,'Delete BIN-File: %s',BIN_File);
      try
         delete(BIN_File); 
      catch
         fprintf(lid,'   %s','error');
      end
      fprintf(lid,'\n');
   end

   msg = sprintf('Can''t open DAT-File "%s" for writing.',DAT_File);

   return

end

fprintf(lid,'\n');

%***************************************************************
% Write Header

if strcmp(mode,'w+')
   wrt_tmp(bid,did,prc,bsi,nx,ny);
end

%***************************************************************
% Write Data

if isempty(xy)
   fclose(bid);
   fclose(did);
   return
end

is_lev = ( size(xy,2) == 3 );

form = '\r   %10.0f Points of %8.0f Segments in %6.0f BinSegments';

%***************************************************************
% Get StartIndex of Segments

n = size(xy,1);

i0 = cat( 1 , 1 , find(isnan(xy(:,1)))+1 );

if ~isnan(xy(n,1))
    i0 = cat( 1 , i0 , n+2 );
end

ns = size(i0,1) - 1;

l0 = diff(i0,1,1) - 1;

%***************************************************************
% Write Segments

% GMT-Area: ( 0 : bsi : 360 ), ( 90 : -bsi : -90 )

fprintf(lid,'\nWrite %.0f Segments (%.0f Points) ... \n',ns,sum(l0));

for ii = 1 : ns

    ind = i0(ii) - 1 + ( 1 : l0(ii) );

    if is_lev
       lv = xy(ind(1),3);
       if ~all( xy(ind,3)-lv == 0 )
           lv = 0;
       end
    else
       lv = 0;
    end

    [bx,by,x,y] = xy2bin(xy(ind,1),xy(ind,2),bsi,nx,ny,px,py);

    bn = bx + ( by - 1 ) * nx;

    %*******************************************************************
    if all( bn-bn(1) == 0 )
    %*******************************************************************

       [x0,y0] = binoff(bsi,py,bx(1),by(1));

       wrt_tmp(bid,did,prc,bsi,x,y,x0,y0,bn(1),lv)

       nb = nb + 1;
       np = np + size(x,1);

       fprintf(lid,form,np,ii,nb);

    %*******************************************************************
    else
    %*******************************************************************

       %-----------------------------------------
       % Make continous in X

       add = cat( 1 , 0 , diff(x,1,1) );
       add = -sign(add) .* ( abs(add) > px/2 );
       add = cumsum(add,1);

        x =  x + px * add;

       %-----------------------------------------
       % Get Values at Border
          
       c = interpg([x y+90]/bsi) * bsi;

       [bx,by,x,y] = xy2bin(c(:,1),c(:,2)-90,bsi,nx,ny,px,py);

       add = cat( 1 , 0 , diff(x,1,1) );
       add = -sign(add) .* ( abs(add) > px/2 );
       add = cumsum(add,1);

         x = x + px * add;

        bn = bx + ( by - 1 ) * nx;

       %-----------------------------------------

       bdd = cat( 1 , 0 , diff(bx,1,1) );
       bdd = -sign(bdd) .* ( abs(bdd) > nx/2 );
       bdd = cumsum(bdd,1);

        bx = bx + nx * bdd;

       %-----------------------------------------
       % StartIndex and Length of Segments

       n1 = size(x,1);

       i1 = cat( 1 , 1 , find(~( diff(bn) == 0 ))+1 );
       l1 = diff(cat(1,i1,n1+1));

       %-----------------------------------------
       % Points to insert at Border

       n2 = size(i1,1) - 1;

       i2 = i1( ( 1 : n2 ) + 1 );

       dx = bx(i2) - bx(i2-1);
       dy = by(i2) - by(i2-1);

       l2 = abs(dx) + abs(dy);

       %-----------------------------------------
       % New Index for Insertion

       % [i1             l1] = [Start Length] of Segment
       % [i1((1:n2)+1)-1 l2] = [Index Length] to Insert behind Index

       i1 = i1 + cumsum(cat(1,0,l2),1);

       ind = grp2ind(i1,l1);            % Expanded Origin

       bn(ind) = bn;
       bx(ind) = bx;
       by(ind) = by;
        x(ind) =  x;
        y(ind) =  y;

       ind = grp2ind(i1(1:n2)+l1(1:n2),l2);  % Index to Insert

       bn(ind) = NaN;
       bx(ind) = NaN;
       by(ind) = NaN;
        x(ind) = NaN;
        y(ind) = NaN;

       %-----------------------------------------
       % Insert BorderValues

       for jj = 1 : n2

           ik = i1(jj) + l1(jj) - 1 + ( 1 : l2(jj) );  % Index to insert
           jk = ik([1 l2(jj)]) + [ -1  1 ];            % Index of border

           %*******************************************************************
           % Cross in X and Y
           if     l2(jj) == 2

                  [xb,yb,xi,yi] = xy2bin(mean(x(jk)),mean(y(jk)),bsi,nx,ny,px,py);

                   bnr = xb + (yb-1) * nx;

                   dd = xi -  x(jk(1));
                   db = xb - bx(jk(1));

                   xi = xi - px * sign(dd) .* ( abs(dd) > px/2 );
                   xb = xb - nx * sign(db) .* ( abs(db) > nx/2 );

                  ins = cat( 1 , [ bn(jk(1)) bx(jk(1)) by(jk(1)) x(jk(1)) y(jk(1)) ] , ...
                                 [ bnr       xb        yb        xi       yi       ] , ...
                                 [ bn(jk(2)) bx(jk(2)) by(jk(2)) x(jk(2)) y(jk(2)) ] );

                   db = diff(ins(:,[2 3]),1,1);

                   for ib = [ 1  2 ]

                       dd = ~( db(ib,1) == 0 ) * db(ib,1) + ...
                            ~( db(ib,2) == 0 ) * db(ib,2);

                       kk = ib + [ 0  1 ];
                       kk = kk([1 2]+[1 -1]*( dd < 0 ));

                       bn(ik(ib)) = ins(kk(1),1);
                       bx(ik(ib)) = ins(kk(1),2);
                       by(ik(ib)) = ins(kk(1),3);
                        x(ik(ib)) = ins(kk(2),4);
                        y(ik(ib)) = ins(kk(2),5);

                   end

           %*******************************************************************
           % Cross in X xor Y
           elseif l2(jj) == 1

                  dd = ~( dx(jj) == 0 ) * dx(jj) + ...
                       ~( dy(jj) == 0 ) * dy(jj);

                  jk = jk([1 2]+[1 -1]*( dd < 0 ));

                  bn(ik) = bn(jk(1));
                  bx(ik) = bx(jk(1));
                  by(ik) = by(jk(1));
                   x(ik) =  x(jk(2));
                   y(ik) =  y(jk(2));

           end
           %*******************************************************************

       end

       %-----------------------------------------
       % Adjust bx

         bx = bx - nx * floor( bx / nx );
         bx = bx + nx * ( bx == 0 );

       %-----------------------------------------
       % StartIndex and Length of Segments

       n1 = size(x,1);

       i1 = cat( 1 , 1 , find(~( diff(bn) == 0 ))+1 );
       l1 = diff(cat(1,i1,n1+1));

       n1 = size(i1,1);

       for jj = ( 1+(l1(1)==1) : n1-(l1(n1)==1) )

           kk  = i1(jj);

           ind = kk - 1 + ( 1 : l1(jj) );

           [x0,y0] = binoff(bsi,py,bx(kk),by(kk));

           xi = x(ind) - px*floor(x(ind)/px);
           xi = xi + px * ( xi-x0 < 0 );

           wrt_tmp(bid,did,prc,bsi,xi,y(ind),x0,y0,bn(kk),lv)
          
           nb = nb + 1;
           np = np + l1(jj);

           fprintf(lid,form,np,ii,nb);
 
       end

    %*******************************************************************
    end
    %*******************************************************************

end

fprintf(lid,'\n\n');

fclose(bid);
fclose(did);

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bx,by,x,y] = xy2bin(x,y,bsi,nx,ny,px,py)

%*********************************************
% Transform Longitude

x = x - px * floor((x+px/2)/px);  % [ -px/2 .. px/2 )

x = x + px * ( x < 0 );           % [ 0 .. px )

bx = floor(x/bsi) + 1;

bx = bx + nx * ( ( bx <= 0 ) - ( bx > nx ) );

%*********************************************
% Transform Latitude

y = y + py/2;                          % [ -py/2 .. py/2 ] --> [ 0 .. py ]

y = y - 2*py * floor( y / (2*py) );    %  [ 0 .. 2*py ]

y = y + ( 2*py - 2*y ) .* ( y > py );  %  /\: [ 0 .. py ]

y = y - py/2;                          % [ -py/2  py/2 ]

by = -y + py/2;

by = floor(by/bsi) + 1;

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x0,y0] = binoff(bsi,py,bx,by);

x0 =  ( bx - 1 ) * bsi;
y0 = ( py/2 - by * bsi );

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [uprc,sprc,sc] = get_prec(prc);


uprc = sprintf('uint%.0f',2*prc);
sprc = sprintf( 'int%.0f',1*prc);

sc = 2^prc - 1;

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wrt_tmp(bid,did,prc,bsi,x,y,x0,y0,bn,lv)

[uprc,sprc,sc] = get_prec(prc);

if nargin == 6
% wrt_tmp(bid,did,prec,bsi,nx,ny)
  fwrite(did,[60*bsi x y],uprc);
  fwrite(bid,[60*bsi x y],uprc);
  return
end

% Remove Offset, Scale with BinSize and PrecisionScale

x = round( sc * ( x - x0 ) / bsi );
y = round( sc * ( y - y0 ) / bsi );

x = x - ( sc + 1 ) * ( x >= (sc+1)/2 );
y = y - ( sc + 1 ) * ( y >= (sc+1)/2 );

fwrite(did,x,sprc);
fwrite(did,y,sprc);

fwrite(bid,[bn lv size(x,1)],uprc);

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y] = read_dat(fid,prc,ip,lp);

[uprc,sprc] = get_prec(prc);

if ischar(fid)
   x = zeros(0,3);
   y = fopen(fid,'r');
   if y == -1
      x = -1;
      return
   end
   d = fread(y,3,uprc); % [ bsi*60 nx ny ]
   if size(d,1) < 3
      return
   end
   x = permute(d,[2 1]);
   x(1) = x(1)/60;
   if nargout < 2
      fclose(y);
   end
   return
end

fseek(fid,prc/8*(ip-1)+3*2*prc/8,-1);

d = fread(fid,lp,sprc);

if size(d,1) == lp

   n = lp/2;

   x = d(1:n);
   y = d(n+(1:n));

end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bsi,nx,ny,bin] = read_bin(file,prc);

bsi = [];
nx  = 0;
ny  = 0;
bin = zeros(0,3);

fid = fopen(file,'r');

if fid == -1
   bsi = NaN;
   return
end

prc = get_prec(prc);

%-----------------------------------
% Get Info: [ bsi nx ny ]

d = fread(fid,3,prc);

if size(d,1) < 3
   return
end

bsi = d(1) / 60;
 nx = d(2);
 ny = d(3);

%-----------------------------------
% Get BinSegments: [ bnr lev siz ]

if nargout == 4

   d = fread(fid,prc);

   if ~isempty(d)
       n = size(d,1);
       if ( mod(n,3) == 0 )
           bin = permute(reshape(d,3,n/3),[2 1]);
       else
           bin = n;
       end
   end

end

fclose(fid);

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,bsi,nx,ny,bin,did] = check_file(BIN_File,DAT_File,prc,lid)

Nout = nargout;

msg = '';
nl  = char(10);

bsi = [];
nx  = 0;
ny  = 0;
bin = zeros(0,3);
did = [];

%***************************************************************

if ~( exist(BIN_File,'file') == 2 )
    msg = sprintf('BIN-File "%s" doesn''t exist.',BIN_File);
end

if ~( exist(DAT_File,'file') == 2 )
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            sprintf('DAT-File "%s" doesn''t exist.',DAT_File) ];
end

if ~isempty(msg)
    return
end

%***************************************************************
% Read Info

%----------------------------------------------------------------------
% DatFile

fprintf(lid,'Read BinInfo from DAT-File: %s',DAT_File);

if Nout == 6
   [inf,did] = read_dat(DAT_File,prc);
else
    inf = read_dat(DAT_File,prc);
end

if isequal(inf,-1)
   msg = sprintf('Can''t open DAT-File "%s" for reading.',DAT_File);
elseif isempty(inf)
   msg = sprintf('Not enough HeaderData in DAT-FILE: %s',DAT_File);
end
   
if ~isempty(msg)
    fprintf(lid,'   %s\n','error');
    return
end

fprintf(lid,'\n');

%----------------------------------------------------------------------
% BinFile

fprintf(lid,'Read BinInfo from BIN-File: %s',BIN_File);

[bsi,nx,ny,bin] = read_bin(BIN_File,prc);

if     isempty(bsi)
   msg = sprintf('Not enough HeaderData in BIN-FILE: %s',BIN_File);
elseif isnan(bsi)
   msg = sprintf('Can''t open BIN-File "%s" for reading.',BIN_File);
elseif ~( size(bin,2) == 3 )
   msg = sprintf('No valid BinData in BIN-FILE: %s',BIN_File);
elseif ~isequal(inf,[bsi nx ny])
   msg = 'HeaderData of BIN-FILE and DAT-File must be equal.';
end   

if ~isempty(msg)
    fprintf(lid,'   %s\n','error');
    return
end

fprintf(lid,'\n');

%----------------------------------------------------------------------
% Check Size of DAT-File

s = prc/8 * ( 3*2 + 2*sum(bin(:,3)) );

d = dir(DAT_File);

if ~( s == d.bytes )
    msg = sprintf('Invalid Size of DAT-File: ',DAT_File);
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function xy = check_dat(xy);

%***************************************************************
% Remove following NaN-Values

ii = isnan(xy(:,1));

if all(ii)
   xy = zeros(0,size(xy,2));
   return
end

if any(ii)

   ii = find(ii);
   ii = ii( find(diff(ii)==1) + 1 );

   xy(ii,:) = [];

   if isnan(xy(1,1))
      xy = xy(2:end,:);
   end

end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bin,sgm,dat,dim,var,att] = file_def(bsi,nx,ny,bin,sc);

%-----------------------------------------------------------
% Add [ Start Length ] of Segment

ns = size(bin,1);

dat      = 2 * bin(:,[3 3]);
dat(:,1) = cumsum(cat(1,1,dat(1:ns-1,2))); % StartIndex in Length in DAT-File

[h,si] = sort(bin(:,1));

bin = bin(si,:);
dat = dat(si,:);

sgm = cat( 1 , 1 , find(~( diff(bin(:,1)) == 0 ))+1 );
sgm = cat( 2 , sgm , diff(cat(1,sgm,ns+1)) );

np = sum(bin(:,3));

nb = nx * ny;

bid = bin(sgm(:,1),1);
sid = zeros(nb,2);
sid(bid,1) = 1;
sid(bid,2) = sgm(:,2);

sid(:,1) = cumsum(sid(:,1));
    jj    = find( sid(:,1) > 0 );
sid(jj,1) = sgm(sid(jj,1),1) - 1;

bin = cat( 2 , bin , cumsum(cat(1,1,bin(1:ns-1,3))) ); % StartIndex

%-----------------------------------------------------------
% NetCDF-Definitions

dim = { ...
        'Dimension of scalar'              1
        'Dimension of bin arrays'         nb
        'Dimension of segment arrays'     ns
        'Dimension of point arrays'       np   };


var = { ...

'Bin size in minutes'                       'long'     [1]    [0]    bsi*60
'Bin size in degree'                        'long'     [1]    [0]    bsi
'N bins in 360 deg longitude range'         'long'     [1]    [0]     nx
'N bins in 180 deg latitude range'          'long'     [1]    [0]     ny
'N bins in file'                            'long'     [1]    [0]     nb
'N segments in file'                        'long'     [1]    [0]     ns
'N points in file'                          'long'     [1]    [0]     np
'Id of first segment in a bin'              'long'     [1]    [1]    sid(:,1)
'N segments in a bin'                       'short'    [1]    [1]    sid(:,2)
'Hierarchial level of a segment'            'short'    [1]    [2]    bin(:,2)
'Id of first point in a segment'            'long'     [1]    [2]    bin(:,4)-1
'Relative longitude from SW corner of bin'  'short'    [1]    [3]     []
'Relative latitude from SW corner of bin'   'short'    [1]    [3]     []  };



scl = sprintf('%.0f',sc);
uni = sprintf('1/%s of %.0f degrees relative to SW corner of bin',scl,bsi);

a0 = cell(0,3);

att = { ...

{ 'units'       'char'   'minutes'                 % 'Bin size in minutes'
  'short_name'  'char'   'bsi'     }

{ 'units'       'char'   'degree'                  % 'Bin size in degree'
  'short_name'  'char'   'dsi'
  'definition'  'char'   'dsi = bsi / 60' }

{ 'short_name'  'char'   'nx'                      % 'N bins in longitude'
  'definition'  'char'   'nx = 360 / dsi' 
  'bin_edges'   'float'  [ 0  360  bsi ]
  'bin_number'  'char'   'bx = ( 1 : nx )'
  'bin_offset'  'char'   'xb = ( bx - 1 ) * dsi' }

{ 'short_name'  'char'   'ny'                      % 'N bins in latitude'
  'definition'  'char'   'ny = 180 / dsi'
  'bin_edges'   'float'  [ 90 -90  -bsi ]
  'bin_number'  'char'   'by = ( 1 : ny )'
  'bin_offset'  'char'   'yb = 90 - by * dsi' }

{ 'short_name'  'char'   'nb'                      % 'N bins in file'
  'definition'  'char'   'nb = nx * ny'
  'bin_number'  'char'   'bn = bx + ( by - 1 ) * nx' }

a0   % 'N segments in file'
a0   % 'N points in file' 
a0   % 'Id of first segment in a bin' 
a0   % 'N segments in a bin'
a0   % 'Hierarchial level of a segment'
a0   % 'Id of first point in a segment'

{ 'units'       'char'  uni                        % 'Relative longitude'
  'short_name'  'char'  'x'
  'scale'       'long'  sc
  'conversion'  'char'  'lon = xb + dsi * ( x + (scale+1) * ( x < 0 ) ) / scale'
  'formula'     'char'  'lon = xb + dsi * ( x + (scale+1) * sign(x).*(sign(x)-1)/2 ) / scale' }

{ 'units'       'char'  uni                        % 'Relative latitude'
  'short_name'  'char'  'y'
  'scale'       'long'  sc
  'conversion'  'char'  'lat = yb + dsi * ( y + (scale+1) * ( y < 0 ) ) / scale'
  'formula'     'char'  'lat = yb + dsi * ( y + (scale+1) * sign(y).*(sign(y)-1)/2 ) / scale' }

};

dim(:,1) = strrep(dim(:,1),' ','_');
var(:,1) = strrep(var(:,1),' ','_');

%-----------------------------------------------------------
% Check DataType SHORT in VariableDefinition for valid Range
%  Convert to LONG if Range exeeded
%
% Limits for Short: [ -32768  32767 ]
%

lim = [ -2  2 ] .^ (16-1) + [ 0  -1 ];

ii = find(strcmp(var(:,2),'short'));

if ~isempty(ii)
    for jj = ii(:)'
        v = var{jj,5};
        if ~isempty(v)
            if ( min(v) < lim(1) ) | ( lim(2) < max(v) ) 
               var{jj,2} = 'long';
            end
        end
    end
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [z,v] = interpg(x,y)

% INTERPG  Interpolate coordinates of polygon on edges of grid
%
% Z = INTERPG( X )
%
% X = [ N by M ]  N Points with M Coordinates in units [1] of Grid
%
% Z = [ K by M ]  with ( K >= N ) if no duplicate Points of X exist
%
% Z contains the unique Coordinates off P, with added Points
%   where the lines of the polygon P crossing edges of the grid.
%   at the added Points:  any( Z == round(Z) )
%

if nargin < 1
   error('Not enough InputArguments.')
end

if ~isnumeric(x)
    error('X must be numeric.')
end

sx = size(x);
ps = prod(sx);
ns = size(sx,2);

if ps <= 1
   z = x;
   return
end

prm = ( ps == max(sx) );
if prm
   dim = sum(cumprod(double(sx==1))) + 1;
   prm = ( dim > 1 );
   if prm
      perm = cat( 2 , ( dim : ns ) , ( 1 : dim-1) );
         x = permute( x , perm );
        sx = sx(perm);
   end
end

nx = sx(1);
mx = sx(2);

if ~( ps == nx*mx )
    error('X must be a 2-dimensional numeric Matrice');
end

n1 = nx - 1;

i1 = ( 1 : n1 );
i2 = i1 + 1;

dx = x(i2,:) - x(i1,:);

%-----------------------------------------
% Remove duplicate Points

jj = ( sum( abs(dx) < 1e2*eps , 2 ) == mx );

if any(jj)

   if nx-sum(jj) < 2
      z = x;
      if prm
         z = permute(z,perm);
      end
      return
   end

   jj = find(jj);

   x(jj+1,:) = [];

   nx = size(x,1);

   n1 = nx - 1;

   i1 = ( 1 : n1 );
   i2 = i1 + 1;

   dx = x(i2,:) - x(i1,:);

end
   
%-----------------------------------------
% Check for BorderCrossings

sg = sign(dx);

x0 =  ceil(sg.*x(i1,:));
x1 = floor(sg.*x(i2,:));

x0 = x0 + ( x0 == sg.*x(i1,:) );
x1 = x1 - ( x1 == sg.*x(i2,:) );

nn = x1 - x0 + 1;
nn = max(nn,0);

if ~any(nn(:))
    z = x;
    if prm
       z = ipermute(z,perm);
    end
    return
end

x0 = sg.*x0;
x1 = sg.*x1;

ns = sum(nn,2) + 1;

i0 = cumsum(cat(1,1,ns));

nz = sum(ns) + 1;

%-----------------------------------------
% Add Distance

x  = cat( 2 , x , cumsum(cat(1,0,sqrt(sum(dx.^2,2)))) );

%-----------------------------------------

z = zeros(nz,mx+1);

z(i0(i1),:) = x(i1,:);
z(nz,:) = x(nx,:);

for ii = i1

    if ns(ii) > 1

         i1 = cumsum(cat(2,1,nn(ii,:))) + i0(ii) - 1;

         for jj = find( ~( nn(ii,:) == 0 ) )

             ix = ( x0(ii,jj) : sg(ii,jj) : x1(ii,jj) )' - x(ii,jj);
      
             z(i1(jj)+(1:nn(ii,jj)),:) = ix * ( x(ii+1,:)  - x(ii,:)  ) / ...
                                              ( x(ii+1,jj) - x(ii,jj) ) + ...
                                          x( ii*ones(1,size(ix,1)) , : );

         end

         ind = i0(ii) + ( 1 : ns(ii)-1 );

         [h,si] = sort(z(ind,mx+1));

         z(ind,:) = z(ind(si),:);

    end

end


%-----------------------------------------
% Remove duplicate Points

z = z(:,1:mx);

dx = diff(z,1,1);

jj = ( sum( abs(dx) < 1e2*eps , 2 ) == mx );

if any(jj)

   jj = find(jj);

   z(jj+1,:) = [];

end

%-----------------------------------------
   
if prm
   z = ipermute(z,perm);
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l,s);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%
% Index = GRP2IND( StartIndex , GroupLength , LowSampleStep  )
%

if isempty(i0);
   ii = [];
   return
end

if nargin < 3
   s = 1;
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

l = ceil( l / s );

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
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

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );

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

