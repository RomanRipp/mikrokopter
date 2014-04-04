function  [cc,msg] = gns_name(file,cnf,area,name);

% GNS_NAME  Returns Features from GIS-NameFiles
%
% [ Feature , Msg ] = GNS_NAME( File , Config , Area , Names )
%
% Config = { FC  {DSG} }
% Area   = [ LonMin  LonMax  LatMin  LatMax ]
% Names  = { Names } , WildCard: '*', 
%            extract Features, matching {Names}
%
% FC  = Feature Classification Code
% DSG = Feature Designation Code
%
% Use the FC  'x'   or '_' to match all Classes of Features
%
% Use the DSG 'xxx' or '_' to match all Features of Class.
%
% Use the WildCard '_' at end of DSG, to extract all Features, 
%  beginning with DSG(1:end-1)
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
   file = [];
end

if nargin < 2
   cnf = { 'x'  'xxx' };
end

if nargin < 3
   area = [];
end

if nargin < 4
   name = [];
else
   [ok,name] = chkcstr(name);
   if ~ok
       name = [];
   end
end

% area = [ 13  14  54  55 ]; %%% Ruegen


msg = '';

mm = 1e6;   % Number of Characters to read

% [msg,file,cnf,area,name] = checkin(file,varargin{:});

%*****************************************

q = char(39);  % "'"
b = char(32);  % Blank
n = char(10);  % NL
t = char(9);   % TAB

wc = '_';   % WildCard at End of Feature Designation Code

wfc = 'x';    % WildCard for Feature Classification Code
wfd = 'xxx';  % WildCard for Feature Designation Code

%*****************************************
% Initialisation of Fields, Indize

[ini,cc,fc,fd,lat,lon,nt] = get_ini(cnf,wc,wfc,wfd);

c0 = cc;

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

%**********************************************************
for zz = 1 : nz
%**********************************************************

nn = min(mm,siz-act);

v = fread(fid,nn,'char');

v = char(v(:)');

nn = max(find(v==n));

v = v(1:nn);

act = act + nn;

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

if ~all( v( is( nt * ( 1 : nr ) ) ) == n )
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

nrw = size(irw,2);

ifc = fc-1 + irw;  % Index in IS of Feature-Classification-Marker
ifd = fd-1 + irw;  % Index in IS of Feature-Designation-Marker

cf = fieldnames(cc);

for f1 = cf(:)'
  
    if strcmp(f1,wfc)
       ii = ( 1 : nrw );
    else
       ii = find( v(is(ifc)+1) == f1{1} );
    end

  %-------------------------------------------
  if ~isempty(ii)
  %-------------------------------------------

    d0 = getfield(c0,f1{1});
    dd = getfield(cc,f1{1});

    df = fieldnames(dd);

    for f2 = df(:)'

        ig = is(ifd(ii));                      % Start of Feature-Designation-Marker
        lg = is(ifd(ii)+1) - is(ifd(ii)) + 1;  % Length of FDM, incl. surr. TAB's

        jj = find( lg > 0 );

        if ~isempty(jj)

            ig = ig(jj);
            lg = lg(jj);

            ng = size(ig,2);

            if strcmp(f2,wfd)
               jj = cumsum( cat( 2 , 1 , lg(1:(ng-1)) ) , 2 );
            else
               ig = grp2ind(ig,lg);
               if f2{1}(end) == wc
                  jj = findstr(v(ig),cat(2,t,f2{1}(1:(end-1))));
               else
                  jj = sort(cat(2,findstr(v(ig),cat(2,t,f2{1},t)), ...
                                  findstr(v(ig),cat(2,t,f2{1},b))));
               end
            end

            if ~isempty(jj)

                ig = zeros(1,sum(lg));
                lg = cat( 2 , 1 , lg );
                ig( cumsum( lg(1:ng) ) ) = ( 1 : ng );

                ig = ig(jj);

                jj = ii(ig);

            end

        end

      %-------------------------------------------------------
      if ~isempty(jj)
      %-------------------------------------------------------
        % Extract

        tmp = getfield(d0,f2{1});

        tf  = fieldnames(tmp);

        % Search 'name' first
          kk = find(strcmp(lower(tf),'name'));
          if ~( kk == 1 )
              tf = tf([ kk (1:(kk-1)) ((kk+1):prod(size(tf))) ]);
          end

        ir = irw(jj);
 
        for f3 = tf(:)'

            kk = find(strcmp(ini(:,2),f3{1}));

            if ~isnan(ini{kk,1})

                [vv,lg] = get_value(v,is,ini(kk,:),ir,t,q);

                if strcmp(f3,'name') 
                   if ~isempty(name)
                       ok = strwcmp(cellstr(vv),name);
                       if ~any(ok)
                           break
                       end
                       ok = find(ok);
                       vv = vv(ok,:);
                       lg = lg(ok);
                       ir = ir(ok);
                   end
                   tmp.len = lg(:);  % NAME ==> Set Length
                end

                tmp = setfield(tmp,f3{1},vv);

            end  % ~isnan(ini{kk,1})

        end
        % f3

        if ~isempty(tmp.fid)

            %----------------------------------------------------
            % PreSort

            if size(tmp.fid,1) > 1

               ok = chkname(tmp.fid,tmp.typ,tmp.len);

               for f3 = tf(:)'
                   vv  = getfield(tmp,f3{1});
                   tmp = setfield(tmp,f3{1},vv(ok,:));
               end

            end

            %----------------------------------------------------
            % Append

            val = getfield(dd,f2{1});

            for f3 = tf(:)'
                v0 = getfield(val,f3{1});
                v1 = getfield(tmp,f3{1});
                if ischar(v0)
                   sp = char(32);
                else
                   sp = NaN;
                end
                s0 = size(v0,2);
                s1 = size(v1,2);
                if     s0 > s1
                   v1 = cat(2,v1,sp(ones(1,size(v1,1)),ones(1,s0-s1)));
                elseif s1 > s0
                   v0 = cat(2,v0,sp(ones(1,size(v0,1)),ones(1,s1-s0)));
                end
                val = setfield( val , f3{1} , cat(1,v0,v1) );
            end
            % kk

            dd = setfield( dd , f2{1} , val  );

        end

      %-------------------------------------------------------
      end % ~isempty(jj)
      %-------------------------------------------------------

    end
    % f2

    cc = setfield(cc,f1{1},dd);

  %-------------------------------------------
  end %  ~isempty(ii)
  %-------------------------------------------

end
% f1

%=========================================
end %  ~isempty(is)
%=========================================

loopdot(sc,nz,zz,1);

%**********************************************************
end
%**********************************************************


fclose(fid);


%*****************************************
% Check

cf = fieldnames(cc);

for f1 = cf(:)'
  
    dd = getfield(cc,f1{1});

    df = fieldnames(dd);

    for f2 = df(:)'

        tmp = getfield(dd,f2{1});

        if size(tmp.fid,1) > 1;

           ok = chkname(tmp.fid,tmp.typ,tmp.len);

           for ff = fieldnames(tmp)'
               vv  = getfield(tmp,ff{1});
               tmp = setfield(tmp,ff{1},vv(ok,:));
           end

           dd = setfield(dd,f2{1},tmp);

        end

    end

    cc = setfield(cc,f1{1},dd);

end
        
%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [vv,lg] = get_value(v,is,ini,irw,t,q)

% Get Value from Columns of V, specified by IS, INI and IRW

ic = ini{1} - 1 + irw;   % Index in IS of Column

ig = is(ic);                % Start of Column
lg = is(ic+1) - is(ic);     % Length of Column, incl. Start TAB's

i1 = ( lg == 1 );   % Start-TAB only
lg = lg + i1;       % Use End-TAB too             

ig = grp2ind(ig,lg);

vv = v(ig);

if any(i1)
   ig = cumsum(cat(2,1,lg),2);  % StartTAB in VV
   ig = ig(find(i1)) + 1;       %   EndTAB in VV
   if     strcmp(ini{3},'n')
          vv(ig) = char(0);
          vv     = strrep(vv,char(0),'NaN');
   elseif strcmp(ini{3},'c')
          vv(ig) = ' ';
   else
          vv(ig) = [];

   end
end

lg = lg - 1 - i1;   % Length of Values

if     strcmp(ini{3},'n')
       vv = cat(2,'[',vv(2:end),']');
       vv = strrep(vv,t,';');
elseif strcmp(ini{3},'c')
       vv = cat(2,'[''',vv(2:end),''']');
       vv = strrep(vv,t,''';''');
else
       vv = strrep(vv,q,[q q]);
       vv = cat(2,'[{''',vv(2:end),'''}]');
       vv = strrep(vv,t,'''};{''');
end

vv = eval(vv);

if iscellstr(vv)
   vv = char(vv);
end

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkname(id,typ,len);

% NT  Name Type: C = Conventional;
%                D = Not verified;
%                N = Native;
%                V = Variant or alternate.

typ(find(typ == 'N')) = 0;
typ(find(typ == 'C')) = 1;
typ(find(typ == 'V')) = 2;
typ(find(typ == 'D')) = 3;


[h,si] = sort(len);

[h,sj] = sort(typ(si));
si = si(sj);

[h,sj] = sort(id(si));
si = si(sj);

jj = find( diff(id(si),1,1) == 0 ) + 1;

ok = ( 1 : size(id,1) );

ok(si(jj)) = [];

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

f = which(file);

if ~isempty(f)
    file = f;
end

ok = ( exist(file,'file') == 2 );
if ok
   d = dir(file);
   ok = ( prod(size(d)) == 1 );
   if ok
      ok = ( ~d.isdir & ( d.bytes > 0 ) );
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

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,cnf,area,name] = checkin(file,varargin);

msg  = '';
cnf  = cell(0,2);
area = ones(0,4);
name = cell(0,1);

% area: 4 Element | 4 Columns , numeric
%  cnf: 2 Columns, 1. Column: CellString, ONE Character
%                  2. Column: CellString or CellArray of Character/CellString
% name: CharacterArray | CellString
%
 
%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,cc,fc,fd,lat,lon,nt] = get_ini(cc,wc,wfc,wfd);

%*****************************************
% Initialisation of Fields, Indize

nt  = 25;    % Number of TAB-spaced Columns

lat = 4;
lon = 5;

fc  = 10;     % Feature Classification
fd  = 11;     % Feature Designation Code

pc  = 12;     % PC Populated Place Classification (Scale), P:PPL only

fdc = 'fdc';  % FeatureDesignationCode-Field, *_ only
fcc = 'fcc';  % FeatureClassificationCode-Filed,  x only
pcf = 'scl';  % PC-Field, PPL* only 

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!! NAME is first Field !!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ini = {  24    'name'  ' ' % FULL_NAME_ND 
        NaN    'len'   'n' % Length of Name
         18    'typ'   'c' % NT         Name Type (search for N)
        lat    'lat'   'n' % LAT        Latitude  [deg]
        lon    'lon'   'n' % LONG       Longitude [deg]
         pc     pcf    'c' % PC         Populated Place Classification
          2    'fid'   'n' % UFI        Unique Feature Identifier
         fd     fdc    ' ' % DSG        Feature Designation Code
         fc     fcc    'c' % FC         Feature Classification
};

%---------------------------------------
% Structure

var = permute(ini(:,[2 1]),[2 1]);
var(2,:) = { {char(ones(0,0))} };

kk = find(strcmp(ini(:,3),'n'));
var(2,kk) = { {ones(0,1)} };

kk = find(strcmp(ini(:,3),'c'));
var(2,kk) = { {char(ones(0,1))} };

nc = size(cc,1);

for ii = 1 : nc
    if strcmp(cc(ii,1),wc)
       cc(ii,1) = {wfc};
    end  
    dd = cc{ii,2};
    if ischar(dd)
       dd = cellstr(dd);
    end
    dd = dd(:)';
    dd = dd([1 1],:);
    for jj = 1 : size(dd,2)
        if strcmp(dd(1,jj),wc)
           dd(1,jj) = {wfd};
        end
        vv = var;
        if ~( strcmp(dd(1,jj),wfd) |  ...
              ( strcmp(cc(ii,1),'P') & strncmp(dd(1,jj),'PPL',3) ) )
            vv(:,find(strcmp(vv(1,:),pcf))) = [];
        end
        if ~( strcmp(dd(1,jj),wfd) | strcmp(dd{1,jj}(end),wc) )
            vv(:,find(strcmp(vv(1,:),fdc))) = [];
        end
        if ~strcmp(cc(ii,1),wfc)
            vv(:,find(strcmp(vv(1,:),fcc))) = [];
        end
        dd(2,jj) = { struct(vv{:}) };
    end
    cc{ii,2} = { struct(dd{:}) };
end

cc = permute(cc,[2 1]);

cc = struct(cc{:});
