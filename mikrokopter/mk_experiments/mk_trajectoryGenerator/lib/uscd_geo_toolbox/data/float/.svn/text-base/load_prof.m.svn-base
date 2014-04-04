function [c,varargout] = load_prof(varargin);

% LOAD_PROF  Load Profiles from *_PROF.NC
%
% C = LOAD_PROF( Path , Area , Recurse )
%
%   Path: PathName(s) to search for ProfilesFiles
%
%   Area: [ LonMin LonMax LatMin LatMax ]
%         Get Profiles inside Area
%         empty: No check by Area
%
%   Recurse: '-r' to search recursivly in Path(s)
%            '-R' to follow DirectoryLinks (UNIX only)
%            see also WHICHFILE
%
% Read from ProfileFile the Variables: 
%
%   JULD, LONGITUDE, LATITUDE, INST_REFERENCE
%   PRES, TEMP, PSAL
%
%   JULD is based on 1950-01-01 00:00:00
%
% Note: No QualityCheck is done !!!
%    
% Returns Structure C, one Element per File:
%
%   name   FileName without '_TRAJ.NC'
%   pfad   Directory of File
%   inst   INST_REFERENCE
%   len    Number of Records
%   day    Matlab's DateNum from JULD
%   lon    LONGITUDE
%   lat    LATITUDE
%   dpt    Depth by P2Z80(PRES,LATITUDE)
%   pr     PRES  at Z
%   tp     TEMP  at Z
%   sl     PSAL  at Z
%
%------------------------------------------------------------
%
% [ C , Lon , Lat , Day , PRES , TEMP , PSAL ] = LOAD_PROF( Z , ... )
%
% Interpolates the Profile to the Levels of Z
%
% returns the Vectors for Day, Lon, Lat
%    and the Matrices for PRES, TEMP, PSAL
%
% LOAD_PROF( Z , Tolerance , ... )
%
% use a Tolerance-Intervall for Z-Levels.
%
%   Tolerance:  single numeric
%
%              use a imaginary Part for a linear Fit trough Z==0
%                dZ + i*Z0, Tolerance(z) = dZ * z/Z0
%
%              empty use default: 5 + 2000*i
%
% Use GRIDDATA to create a regular XY-Grid for TEMP and PSAL
%  per Z-Level. Take care on NaN-Values !!!
%
%------------------------------------------------------------
%
% See also: WHICHFILE, ASSIGN_CDF, LOOK_CDF, P2Z80
%

%********************************************************************

suf = '_prof.nc';

c = struct( 'name' , { '' } , ...
            'pfad' , { '' } , ...
            'inst' , { '' } , ...  % INST_REFERENCE
            'len'  , { 0  } , ...
            'day'  , { [] } , ...
            'lon'  , { [] } , ...
            'lat'  , { [] } , ...
            'dpt'  , { [] } , ...
            'pr'   , { [] } , ...
            'tp'   , { [] } , ...
            'sl'   , { [] }       );

vv  = { 'JULD' 'LONGITUDE' 'LATITUDE' };  % [ NProf ]        --> [ 1       NProf ]
ww  = { 'PRES' 'TEMP' 'PSAL' };           % [ NProf NLevels] --> [ NLevels NProf ]

ref = 'INST_REFERENCE';

bd = datenum(1950,01,01);  % BaseDay for JULD

out = { 'lon' 'lat' 'day' 'pr' 'tp' 'sl' };  % OutPutOrder

def = 5 + 2000*i;    % Default for Tolerance

tlim = [ -4  40 ];   % TempLimit for Profiles !!!

%********************************************************************

Nin  = nargin;
Nout = nargout - 1;

varargout = cell(1,Nout);

v = varargin(:);
if Nin < 5
   v = cat( 1 , v , cell(5-Nin,1) );
end

if chkcstr(v{1})
   zlev = [];
   ztol = [];
   pfad = v{1};
   area = v{2};
   rec  = v{3};
elseif chkcstr(v{2})
   zlev = v{1};
   ztol = [];
   pfad = v{2};
   area = v{3};
   rec  = v{4};
else
   zlev = v{1};
   ztol = v{2};
   pfad = v{3};
   area = v{4};
   rec  = v{5};
end

if ~isempty(zlev)
    if ~isnumeric(zlev)
        error('ZLevel must be numeric');
    else
        zlev = zlev(:);
    end
end

if isempty(ztol)
   ztol = def;
elseif ~( isnumeric(ztol) & ( prod(size(ztol)) == 1 ) )
   error('Tolerance must be a single Numeric.');
end

if isempty(pfad)
   pfad = {'.'};
else
   [ok,pfad] = chkcstr(pfad);
   if ~ok
       error('PathList must be Strings.');
   end
end

if ~( isempty(area) | isequal(size(area),[1 4]) )
   error('Area must have 4 Elements: [ LonMin LonMax LatMin LatMax ].');
end

if isempty(rec)
   rec = {};
else
   [ok,rec] = chkcstr(rec);
   if ~( ok & ( prod(size(rec)) == 1 ) ) 
       error('Option must be a String.');
   end
end

%********************************************************************
% Get Files

[f,p] = whichfile(pfad,suf,rec{:});

if isempty(f)
   fprintf(1,'No Targets found for "*%s".\n',suf);
   c = c(ones(0,1));
   return
end


n  = size(f,1);

c  = c(ones(n,1));

ok = zeros(n,1);

frm = size(char(p),2) + size(char(f),2);
frm = sprintf('%s-%.0fs ... ','%',frm);

nv = prod(size(vv));  % Number of Vectors to extract
sv = cell(nv,1);

nw = prod(size(ww));  % Number of Variables to extract
sw = cell(nw,1);

vv = cat( 2 , vv , ww , {ref} );

for ii = 1 : n

    file = [ p{ii}  f{ii} ];

    fprintf(1,frm,file);

    [m,v] = assign_cdf(file,vv,NaN);

    if ~isempty(m)

        fprintf(1,'error assign');

    else

        for jj = 1 : nv
           sv{jj} = size(getfield(v,vv{jj}));
        end
 
        for jj = 1 : nw
           sw{jj} = size(getfield(v,vv{jj+nv}));
        end
 
        if ~( isequal(sv{:}) & isequal(sw{:}) )

            fprintf(1,'invalid size');

        elseif ( prod(sv{1}) == 0 ) | ( prod(sw{1}) == 0 )

            fprintf(1,'empty variables');

        elseif ~all( prod(sv{1}) == [ max(sv{1}) sw{1}(2) ] )

            fprintf(1,'invalid size of vectors');

        elseif ~( prod(sw{1}) == sw{1}(1)*sw{1}(2) )

            fprintf(1,'no 2D-Variables');

        else

            xx = v.LONGITUDE;
            yy = v.LATITUDE;
            dd = v.JULD;

            ok(ii) = isempty(area);

            if ~ok(ii)
                jj = ( ( area(1) <= xx ) & ( xx <= area(2) ) & ...
                       ( area(3) <= yy ) & ( yy <= area(4) )       );
                ok(ii) = any(jj);
                if ok(ii) & ~all(jj)
                   jj = find(jj);
                   for kk = 1 : nw
                       w = getfield(v,ww{kk}); 
                       w = w(:,jj);
                       v = setfield(v,ww{kk},w);
                   end
                   xx = xx(jj);
                   yy = yy(jj);
                   dd = dd(jj);
                end               
            end

            if ~ok(ii)

                fprintf(1,'not in area');

            else

                pr = v.PRES;
                tp = v.TEMP;
                sl = v.PSAL;

                jj = ( (  tlim(1) <= min(tp,[],1) ) & ( max(tp,[],1) <= tlim(2) ) );

                ok(ii) = any(jj);

                if ~ok(ii)

                   fprintf(1,'no profiles in T-Range');

                else

                   if ~all(jj)
                       jj = find(jj);
                       pr = pr(:,jj);
                       tp = tp(:,jj);
                       sl = sl(:,jj);
                       xx = xx(jj);
                       yy = yy(jj);
                       dd = dd(jj);
                   end

                   dpt = p2z80(pr,yy(ones(1,size(pr,1)),:));

                   if ~isempty(zlev)
                       [xx,yy,dd,dpt,pr,tp,sl] = extract(xx,yy,dd,dpt,pr,tp,sl,zlev,ztol);
                   end

                   ok(ii) = ~isempty(xx);

                   if ~ok(ii)

                       fprintf(1,'not in Z-Range');

                   else

                       c(ii).lon = xx(:)';
                       c(ii).lat = yy(:)';
                       c(ii).day = dd(:)' + bd;

                       c(ii).dpt = dpt;
                       c(ii).pr  = pr;
                       c(ii).tp  = tp;
                       c(ii).sl  = sl;

                       c(ii).name = strrep(f{ii},suf,'');
                       c(ii).pfad = p{ii};

                       c(ii).len  = size(pr,2);

                       c(ii).inst = getfield( v , ref );

                       fprintf(1,'%4.0f records',c(ii).len);

                   end

                end

            end

        end

    end

    fprintf(1,'\n');

end

if ~any(ok)
    c = c(ones(0,1));
elseif ~all(ok)
    ok = find(ok);
    c  = c(ok);
end

if isempty(c) | ( Nout == 0 ) | isempty(zlev)
   return
end

for ii = 1 : min(Nout,size(out,2))

    v = eval(['{c.' out{ii} '}']);

    varargout{ii} = cat( 2 , v{:} );

end


%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xx,yy,dd,zp,pr,tp,sl] = extract(xx,yy,dd,z,p,t,s,zz,dz);

zp = [];
pr = [];
tp = [];
sl = [];

nz = size(zz,1);

if imag(dz) == 0
   dz = dz(ones(nz,1));
else
   dz = real(dz) * zz / imag(dz);
end

zl = cat( 1 , min(z,[],1) , max(z,[],1) );

np = size(xx,2);

op = ones(1,np);
oz = ones(1,nz);

ok = ( ( zl(1*oz,:)-dz(:,op) <= zz(:,op) ) & ...
       ( zl(2*oz,:)+dz(:,op) >= zz(:,op) )        );

if ~any(ok(:))
    xx = [];
    yy = [];
    dd = [];
    return
end

ip = ( sum(ok,1) > 0 );  % True  for good Profiles, matching ZZ
np =  sum(ip);           % Number of good Profiles
ip = find(ip);


xx = xx(ip);
yy = yy(ip);
dd = dd(ip);

ok = ok(:,ip);
zl = zl(:,ip);
z  =  z(:,ip);
p  =  p(:,ip);
t  =  t(:,ip);
s  =  s(:,ip);

zp = NaN * zeros(nz,np);
pr = zp;
tp = zp;
sl = zp;

mz = size(z,1);

qz = ones(mz,1);

for iz = 1 : nz

    nk =  sum(ok(iz,:),2);

    if nk > 1

       ik = find(ok(iz,:));

       zk = min(max(zl(1,ik),zz(iz)),zl(2,ik));

       i0 = ( z(:,ik) < zk(qz,:) );
       i0 = sum(i0,1);
       i1 = i0 + ( i0 < mz );
       i0 = i0 + ( i0 == 0 );

       i0 = i0 + mz * ( ik - 1 );
       i1 = i1 + mz * ( ik - 1 );

       mm = ( zk - z(i0) ) ./ ( z(i1) - z(i0) + ( i1 == i0 ) );

       zp(iz,ik) = zk;
       pr(iz,ik) = ( p(i1) - p(i0) ) .* mm + p(i0);
       tp(iz,ik) = ( t(i1) - t(i0) ) .* mm + t(i0);
       sl(iz,ik) = ( s(i1) - s(i0) ) .* mm + s(i0);

       if all(isnan(pr(iz,ik)))
          fprintf(1,'\nWarning: Invalid interpolation at %.0f. ',zz(iz));
       end

    end

end
