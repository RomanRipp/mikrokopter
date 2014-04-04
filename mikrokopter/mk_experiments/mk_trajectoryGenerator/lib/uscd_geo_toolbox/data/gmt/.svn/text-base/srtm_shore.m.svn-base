function srtm_shore(area,pfad,res,gmt)

% SRTM_SHORE  Create LandMask-Image(s) for SRTM-HGT-Tiles
%
% SRTM_SHORE( Area , Directory )
%
% AREA = [ LonMin LonMax LatMin LatMax ]
%
% Creates Images of LandCovering in [ 1° x 1° ] - Tiles
%  with a Resolution of [ 1201 x 1201 ] (3arcsec). 
% 
% The Images are Gray-Scaled-PNG with the ColorValues: 
%   Land == 000, Ocean == 255, 000 < Lakes < 255;
%
% The FileName of the Image refers to the SW-Corner
%  of the [ 1° x 1° ] - Tile: "Directory/#LAT*LON.png", 
%  where "#" is "N" or "S", "*" is "E" or "W".
%
% The LandCovering is read from binned GMT-ShoreLines,
%  using SHORE_GMT, origin: WVS and CIA WDB-II by Wessel and Smith.
%
% The GMT-Data, extracted for the Area, will saved 
%  into a MAT-File:  "#LatMin#LatMax*LonMin*LonMax.mat",
% to use it for following calls of SRTM_SHORE.
%
%-----------------------------------------------------------------
%
% SRTM_SHORE( HGT_Directory , Directory )
%
% Creates the Images for all HGT-Tile ("*.hgt*"),
%  found in HGT-Directory
%
%-----------------------------------------------------------------
%
% see also: READ_HGT, SHORE_GMT (required)
%

Nin = nargin;

if Nin < 2
   pfad = '';
end

if Nin < 3
   res = [];
end

if Nin < 4
   gmt = '';
end

msg = cell(0,1);

if ~( ( isnumeric(area) & isequal(size(area),[1 4]) ) | chkstr(area,1) )
    msg = cat(1,msg,{'Area must be a DirectoryName or 4-Element Vector.'});
end

ok = chkstr(pfad,0);
if ~ok
    ok = ( isnumeric(pfad) & ( Nin <= 3 ) );
    if ok
       res = pfad;
       gmt = res;
       pfad = '';
    else
       msg = cat(1,msg,{'Invalid Input for Directory (or Resolution).'});
    end
end

if ~( isnumeric(res) & ( prod(size(res)) <= 1 ) )
    msg = cat(1,msg,{'Resolution must be empty or [ WDT + SIZ*i ]'});
end

ok = chkstr(gmt,0);
if ok
   ok = ( prod(size(gmt)) <= 1 );
end

if ~ok
    msg = cat(1,msg,{'GMT_ID must be empty or a single Character'});
elseif ~isempty(gmt)
    gmt = lower(gmt);
    if ~any( gmt == 'clihf' )
        msg = cat(1,msg,{'Invalid GMT_ID for ShoreLines.'});
    end
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%--------------------------------------------------------------
% Get Width and PixelSize

if isempty(res)
   res = 0;
end

wdt = abs(real(res));
if ~isfinite(wdt)
    wdt = 0;
end

wdt = wdt + 1 * ( wdt == 0 );

siz = abs(imag(res));
if ~isfinite(siz)
    siz = 0;
end

siz = siz + 1201 * ( siz == 0 );

if isempty(gmt)
   gmt = 'f';
end

%*****************************************************************
% Check for Directory, look for HGT-Files

if chkstr(area,1)

   cmd = sprintf('find "%s" -name "*.hgt*"',area);

   [s,w] = unix(cmd);

   if ~( s == 0 )
       fprintf(1,'Error call UNIX: %s\n%s',cmd,w);
       return
   elseif isempty(w)
       fprintf(1,'No HGT-Files found in %s',area);
       return
   end

   w = sepname(w,NaN,char(10));

   area = NaN * ones(1,4);

   for ii = 1 : prod(size(w))
       [p,w{ii}] = fileparts(w{ii});
       [m,v] = str2vec(w{ii});
       if isempty(m) & isequal(size(v),[1 2])
          v(1) = ( 1 - 2 * any(w{ii}=='S') ) * v(1);
          v(2) = ( 1 - 2 * any(w{ii}=='W') ) * v(2);
          area(1) = min(area(1),v(2));
          area(2) = max(area(2),v(2));
          area(3) = min(area(3),v(1));
          area(4) = max(area(4),v(1));
       end
   end

   area([1 3]) = floor(area([1 3]));
   area([2 4]) = floor(area([2 4])) + 1;

else

   w = {};

   area([1 3]) = floor(area([1 3]));
   area([2 4]) =  ceil(area([2 4]));

end

%*****************************************************************

ew = 'EW';
ns = 'NS';

%*****************************************************************
% Check for MAT-File

n1 = ns(1+(area(3)<0));
n2 = ns(1+(area(4)<0));
e1 = ew(1+(area(1)<0));
e2 = ew(1+(area(2)<0));

file = sprintf( '%s%2.2d%s%2.2d%s%3.3d%s%3.3d.mat'     , ...
                 n1 , abs(area(3)) , n2 , abs(area(4)) , ...
                 e1 , abs(area(1)) , e2 , abs(area(2))       );

ok = ( exist(file,'file') == 2 );
if ok
   try
      d = load(file);
      d = d.d;
   catch
     ok = 0;
   end
end
if ok
   ok = ( isfield(d,'s') & isfield(d,'l') );
end

if ok

    fprintf(1,'Use %s\n',file);

else

   a = area + [ -1  1 -1  1 ] * 10/60;

   d = struct( 's' , { shore_gmt([gmt '1'],a) } , ...      % Shore
               'l' , { shore_gmt([gmt '2'],a) } , ...      % Lake
               'i' , { shore_gmt([gmt '3'],a) }       );   % Islands in Lake

  fprintf(1,'Save %s',file);

  try
     save(file,'-mat','d');
     fprintf(1,' ok\n');
  catch
     fprintf(1,' error\n');
  end

end

%*****************************************************************
% Sort Shore by Size

nsh = size(d.s,1);

s2 = zeros(nsh,1);

for ii = 1 : nsh
    s2(ii) = size(d.s{ii},2);
end

ib = find( s2 > 10000 );

%*****************************************************************

fmt = 'png';

x0 = area(1);
y0 = area(3);

xl = area([1 2]) - x0;
yl = area([3 4]) - y0;

dg = 3600/1200;  % Intervall [sec]

scl = 3600/dg;

x = ( xl(1) : xl(2)-1 ) * scl;
y = ( yl(1) : yl(2)-1 ) * scl;

off =  20;  % PixelOffset
ext = scl/2 + 1;  % single ImageSize, [ 2 x 2 ] == [ 1201 x 1201 ] !!!

pix = ext + 2 * off;

fig = figure('position',[200 50 pix([1 1])], ...
       'toolbar','none', ...
       'menubar','none');

axe = axes('position',[0 0 1 1], ...
     'visible','off', ...
     'nextplot','add'    );

x01 = [0 pix] - 1/2 - off;
y01 = [0 pix] - 1/2 - off;

i1 = ( 2 : ext );

%----------------------------------------------------------------------

nsh = size(d.s,1);
nlk = size(d.l,1);
nil = size(d.i,1);

for ii = 1 : nsh
    d.s{ii}(1,:) = 0.5 * round( (d.s{ii}(1,:)-x0) * scl / 0.5 );
    d.s{ii}(2,:) = 0.5 * round( (d.s{ii}(2,:)-y0) * scl / 0.5 );
end

for ii = 1 : nlk
    d.l{ii}(1,:) = 0.5 * round( (d.l{ii}(1,:)-x0) * scl / 0.5 );
    d.l{ii}(2,:) = 0.5 * round( (d.l{ii}(2,:)-y0) * scl / 0.5 );
end

for ii = 1 : nil
    d.i{ii}(1,:) = 0.5 * round( (d.i{ii}(1,:)-x0) * scl / 0.5 );
    d.i{ii}(2,:) = 0.5 * round( (d.i{ii}(2,:)-y0) * scl / 0.5 );
end

%**********************************************************************

for iy = 1 : size(y,2)

    for ix = 1 : size(x,2)

%----------------------------------------------------------------------
% Check for File

        xx = round( x(ix)/scl + x0 );
        yy = round( y(iy)/scl + y0 );

        n = ns(1+(yy<0));
        e = ew(1+(xx<0));

fok = ~isempty(w);
if fok

        f = sprintf('%s%2.2d%s%3.3d',n,abs(yy),e,abs(xx));

        fok = ~isempty(strmatch(f,w));

end

if fok
        file = sprintf('%s%2.2d%s%3.3d.%s',n,abs(yy),e,abs(xx),fmt);
        fok = ~( exist(fullfile(pfad,file)) == 2 );
end

%----------------------------------------------------------------------
if fok
%----------------------------------------------------------------------


        delete(findobj(axe,'type','patch','tag','ISLAND'));
        delete(findobj(axe,'type','patch','tag','LAKE'));
        delete(findobj(axe,'type','patch','tag','ISL'));

        xlm = x01 + x(ix) + [ 0  0.5 ] * scl;
        ylm = y01 + y(iy) + [ 0  0.5 ] * scl;

        %---------------------------------------------------------------------
        % Shore

if ~isempty(d.s)

        ok = zeros(nsh,1);
        for ii = 1 : nsh
            ok(ii) = any( ( xlm(1) <= d.s{ii}(1,:) ) & ( d.s{ii}(1,:) <= xlm(2) ) & ...
                          ( ylm(1) <= d.s{ii}(2,:) ) & ( d.s{ii}(2,:) <= ylm(2) )       );
        end

        ok(ib) = 1;
 
        nn = sum(ok);

        ok = find(ok);
        for ii = 1 : nn
            xx = min(max(d.s{ok(ii)}(1,:),xlm(1)),xlm(2));
            yy = min(max(d.s{ok(ii)}(2,:),ylm(1)),ylm(2));
            jj = find( ( diff(xx,1,2) == 0 ) & ( diff(yy,1,2) == 0 ) );
            xx(jj+1) = [];
            yy(jj+1) = [];
            patch(xx,yy,'k','edgecolor','k','tag','ISLAND','parent',axe);
        end
end

        %---------------------------------------------------------------------
        % Lakes

if ~isempty(d.l)

        ok = zeros(nlk,1);
        for ii = 1 : nlk
            ok(ii) = any( ( xlm(1) <= d.l{ii}(1,:) ) & ( d.l{ii}(1,:) <= xlm(2) ) & ...
                          ( ylm(1) <= d.l{ii}(2,:) ) & ( d.l{ii}(2,:) <= ylm(2) )       );
        end
 
        nn = sum(ok);

        if nn > 253, keyboard, end

        cc = linspace(0,1,nn+2);
        cc = round(255*cc(2:end-1))/255;

        ok = find(ok);

        for ii = 1 : nn
            xx = min(max(d.l{ok(ii)}(1,:),xlm(1)),xlm(2));
            yy = min(max(d.l{ok(ii)}(2,:),ylm(1)),ylm(2));
            jj = find( ( diff(xx,1,2) == 0 ) & ( diff(yy,1,2) == 0 ) );
            xx(jj+1) = [];
            yy(jj+1) = [];
            patch(xx,yy,cc(ii)*[1 1 1],'edgecolor','none','tag','LAKE','parent',axe);
        end
end
        
        %---------------------------------------------------------------------
        % Islands

if ~isempty(d.i)

        ok = zeros(nil,1);
        for ii = 1 : nil
            ok(ii) = any( ( xlm(1) <= d.i{ii}(1,:) ) & ( d.i{ii}(1,:) <= xlm(2) ) & ...
                          ( ylm(1) <= d.i{ii}(2,:) ) & ( d.i{ii}(2,:) <= ylm(2) )       );
        end
 
        nn = sum(ok);

        ok = find(ok);

        for ii = 1 : nn
            xx = min(max(d.i{ok(ii)}(1,:),xlm(1)),xlm(2));
            yy = min(max(d.i{ok(ii)}(2,:),ylm(1)),ylm(2));
            jj = find( ( diff(xx,1,2) == 0 ) & ( diff(yy,1,2) == 0 ) );
            xx(jj+1) = [];
            yy(jj+1) = [];
            patch(xx,yy,'k','edgecolor','k','tag','ISL','parent',axe);
        end
end
        
        %---------------------------------------------------------------------

        set(axe,'xlim',x01+x(ix),'ylim',y01+y(iy)+0.5*scl);
        drawnow
        figure(fig)
        c1 = getframe(fig);
        c1 = c1.cdata(off+(1:ext),off+(1:ext),1);

        set(axe,'xlim',x01+x(ix)+0.5*scl,'ylim',y01+y(iy)+0.5*scl);
        drawnow
        figure(fig)
        c2 = getframe(fig);
        c2 = c2.cdata(off+(1:ext),off+(1:ext),1);

        set(axe,'xlim',x01+x(ix),'ylim',y01+y(iy));
        drawnow
        figure(fig)
        c3 = getframe(fig);
        c3 = c3.cdata(off+(1:ext),off+(1:ext),1);

        set(axe,'xlim',x01+x(ix)+0.5*scl,'ylim',y01+y(iy));
        drawnow
        figure(fig)
        c4 = getframe(fig);
        c4 = c4.cdata(off+(1:ext),off+(1:ext),1);

        ok = [ isequal(c1(:,end),c2(:,1)) isequal(c3(:,end),c4(:,1)) 
               isequal(c1(end,:),c3(1,:)) isequal(c2(end,:),c4(1,:))  ];

        c1(:,end) = max(c1(:,end),c2(:,1)); c1 = cat(2,c1,c2(:,i1));
        c3(:,end) = max(c3(:,end),c4(:,1)); c3 = cat(2,c3,c4(:,i1));

        c1(end,:) = max(c1(end,:),c3(1,:));
        
        c = cat( 1 , c1 , c3(i1,:) );

        if 0 % ~all(ok(:))
            file,ok, keyboard
        end

        fprintf(1,'%s',file);

        imwrite(c,fullfile(pfad,file),fmt);

        fprintf(1,'\n');

%----------------------------------------------------------------------
end  % FileCheck
%----------------------------------------------------------------------
         
    end

end