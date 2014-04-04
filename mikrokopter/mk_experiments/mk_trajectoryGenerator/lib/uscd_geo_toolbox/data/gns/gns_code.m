function  [cc,nn,msg] = gns_code(file,area);

% GNS_CODE  Returns FeatureCodes from GNS-NameFiles
%
% [ Codes , Number , Msg ] = GNS_CODE( File , Area )
%
%  Area   = [ LonMin  LonMax  LatMin  LatMax ]
%
%  Codes = ['F DSG']  CharacterArray
%
% F   = Feature Classification Code (single Character)
% DSG = Feature Designation Code
%
%****************************************************************
% TAB-separated Columns of GNS-NameFiles:
%
%  1. RC         Region Code (World)
%  2. UFI        Unique Feature Identifier
%  3. UNI        Unique Name Identifier
%  4. LAT        Latitude  [deg]
%  5. LONG       Longitude [deg]
%  6. DMS_LAT    Latitude  [ddmmss]
%  7. DMS_LONG   Longitude [ddmmss]
%  8. UTM        Universal Transverse Mercator coordinate grid reference.
%  9. JOG        Joint Operations Graphic reference.
% 10. FC         Feature Classification:
%                  A = Administrative region;
%                  P = Populated place;
%                  V = Vegetation;
%                  L = Locality or area;
%                  U = Undersea;
%                  R = Streets, highways, roads, or railroad;
%                  T = Hypsographic;
%                  H = Hydrographic;
%                  S = Spot feature.
% 11. DSG        Feature Designation Code
% 12. PC         Populated Place Classification (Scale)
% 13. CC1        Primary Country Code
% 14. ADM1       First-order administrative division
% 15. ADM2       Second-order administrative division
% 16. DIM        Dimension, elevation or population data.
% 17. CC2        Secondary Country Code
% 18. NT         Name Type:
%                   C = Conventional;
%                   D = Not verified;
%                   N = Native;
%                   V = Variant or alternate.
% 19. LC         Language Code
%
% 20. SHORT_FORM A specific part of the name that could substitute 
%                  for the full name.
%
% 21. GENERIC    The descriptive part of the full name
%
% 22. SORT_NAME  A form of the full name which allows for easy sorting 
%                 of the name into alpha-numeric sequence.
%
% 23. FULL_NAME  The full name is a complete name which identifies 
%                 the named feature.
%
% 24. FULL_NAME_ND  Same as the full name but the diacritics and 
%                   special characters are substituted with Roman characters.
%
% 25. MOD_DATE   (YYYY-MM-DD)
%
% 

if nargin < 1
   file = '';
end

if nargin < 2
   area = [];
end

msg = '';
cc  = char(ones(0,1));
nn  = ones(0,1);

if isempty(file)
   return
end

%*****************************************

q = char(39);  % "'"
b = char(32);  % Blank
n = char(10);  % NL
t = char(9);   % TAB

nt = 25;    % Number of TAB-spaced Columns

mm = 1e6;   % Number of Characters to read

%*****************************************

fc = 10;  % Feature Classification
fd = 11;  % Feature Designation Code

lat = 4;
lon = 5;

%*****************************************
% Check File

if isempty(file)
   return
end

[msg,file,fid,siz] = chkfile(file);

if ~isempty(msg)
    return
end

%*******************************************************************
% Start with Loop

act = 0;

nz = ceil(siz/mm);
sc = [ 50  50 ];

txt = sprintf('Read %.0f kB Blocks (%.1f MB)',mm./(2.^[10 20]));

loopdot(sc,nz,txt);

%*****************************************
for zz = 1 : nz
%*****************************************

nr = min(mm,siz-act);

v = fread(fid,nr,'char');

v = char(v(:)');

nr = max(find(v==n));

v = v(1:nr);

act = act + nr;

fseek(fid,act,'bof');

%*****************************************

v(find(v==13)) = [];

is = find( ( v == t ) | ( v == n ) );  % Seperators 25 per Row

ns = size(is,2);

if ~( mod(ns,nt) == 0 )
    msg = 'Invalid Number of Seperators.';
    fclose(fid);
    return
end

nr = ns/nt;

if ~all( v( is( nt * ( 1 : ns/nt ) ) ) == n )
    msg = 'Invalid Number of TAB-spaced Columns.';
    fclose(fid);
    return
end

%*****************************************

irw = nt * ( 0 : (nr-1) );   % RowIndex

%*****************************************
% Check for Number of RegionCode 
%  as First Character of Line:    0 .. 9

jj = cat(2,0,is(irw(2:nr))) + 1;

jj = find( ( '0' <= v(jj) ) & ( v(jj) <= '9' ) );

irw = irw(jj);

%*****************************************
% Check with Area

if ~isempty(area) & ~isempty(irw)

    irw = chkarea(v,is,irw,area,lat,lon,t);

end

%*****************************************

%=========================================
if ~isempty(irw)
%=========================================

nr = size(irw,2);

bb = char(ones(nr,0));

for ff = [ fc  fd ]

    im = ff-1 + irw;  % Index in IS of Marker

    ig = is(im);               % Start of Marker
    lg = is(im+1) - is(im);    % Length of FDM, incl. Start-TAB

    ig = grp2ind(ig,lg);

    vv = v(ig);

    vv = strrep(vv,q,[q q]);

    if all( lg == lg(1) )
       vv = cat(2,'[''',vv(2:end),''']');
       vv = strrep(vv,t,''';''');
    else
       vv = cat(2,'[{''',vv(2:end),'''}]');
       vv = strrep(vv,t,'''};{''');
    end

    vv = eval(vv);

    if iscellstr(vv)
       vv = char(vv);
    end

    bb = cat(2,bb,b(ones(1,nr),~(size(bb,2)==0)),vv);

end

bb = sortrows(bb);

ig = find( sum( diff(bb,1,1) == 0 , 2 ) < size(bb,2) );

ig = cat( 1 , 1 , ig+1 , size(bb,1)+1 );
lg = diff(ig,1,1);

bb = bb(ig(1:(end-1)),:);

s0 = size(cc,2);
s1 = size(bb,2);
if     s0 > s1
       bb = cat(2,bb,b(ones(1,size(bb,1)),ones(1,s0-s1)));
elseif s1 > s0
       cc = cat(2,cc,b(ones(1,size(cc,1)),ones(1,s1-s0)));
end

cc = cat(1,cc,bb);
nn = cat(1,nn,lg);

%=========================================
end   % ~isempty(irw)
%=========================================

loopdot(sc,nz,zz,1);

%*****************************************
end
%*****************************************


fclose(fid);


%*****************************************
% Check

if isempty(cc)
   return
end

[cc,si] = sortrows(cc);

nn = cumsum(nn(si),1);

ig = find( sum( diff(cc,1,1) == 0 , 2 ) < size(cc,2) );

ig = cat( 1 , 1 , ig+1 , size(cc,1)+1 );

ind = ( 1 : (size(ig,1)-1) );

cc = cc(ig(ind+0)-0,:);
nn = nn(ig(ind+1)-1);

nn = nn - cat(1,0,nn(1:(end-1)));


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function irw = chkarea(v,is,irw,area,lat,lon,t)

%*****************************************
% Check with Area

lb = { 'LON'  'LAT' };
il = [  lon    lat  ];

na = size(area,1);
nr = size(irw,2);

ok = zeros(2,nr);

for ii = [ 1  2 ]

    jj = il(ii) - 1 + irw;

    ig = is(jj);             % Start of Marker
    lg = is(jj+1) - is(jj);  % Length of Marker, without EndTab

    jj = find( lg > 0 );

    if ~isempty(jj)
        % Check 1. Marker, may be Header

        kk = ig(jj(1)) + ( 1 : ( lg(jj(1)) - 1 ) );
        kk = v(kk);

        try
           h = eval(kk);
        catch
           jj(1) = [];
        end

    end

    if ~isempty(jj)
 
        ig = ig(jj);
        lg = lg(jj);

        ig = grp2ind(ig,lg);

        vv = v(ig(2:end));  % Excl. StartTab

        vv = strrep(vv,t,',');

        kk = [ 1  2 ] + 2* ( ii - 1 );

        try
           vv        = eval(['[' vv ']']);
        catch
           warning(sprintf('Error evaluate %s.\n%s',lab{ii},lasterr));
           vv= [];
        end

        if isequal(size(vv),size(jj))
           for ll = 1 : na
                ok(ii,jj) = ( ok(ii,jj) | ...
                              ( ( area(ll,kk(1)) <= vv ) & ...
                                ( vv <= area(ll,kk(2)) )        ) );
           end
        end

    end

end 
 
ok = find( ok(1,:) & ok(2,:) );

irw = irw(ok);

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,fid,siz] = chkfile(file);

msg = '';
fid = [];
siz = 0;

%*****************************************
% Check File

file = which(file);

ok = ~isempty(file);
if ok
   ok = ( exist(file,'file') == 2 );
   if ok
      d = dir(file);
      ok = ( prod(size(d)) == 1 );
      if ok
         ok = ( ~d.isdir & ( d.bytes > 0 ) );
      end
   end
end

if ~ok
    msg = 'File not found.';
    return
end

siz = d.bytes;

%*****************************************

fid = fopen(file,'r');

if fid == -1
   msg = 'Can''t open File.';
   return
end
