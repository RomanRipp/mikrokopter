function [c,x,y,d] = load_traj(pfad,area,rec);

% LOAD_TRAJ  Load Trajectories from *_TRAJ.NC
%
% C = LOAD_TRAJ( Path , Area , Recurse )
%
%   Path: PathName(s) to search for TrajectoryFiles
%
%   Area: [ LonMin LonMax LatMin LatMax ]
%         Get Trajectories which crossing Area
%         empty: No check by Area
%
%   Recurse: '-r' to search recursivly in Path(s)
%            '-R' to follow DirectoryLinks (UNIX only)
%            see also WHICHFILE
%
% Read from TrajectoryFile the Variables: 
%
%   JULD, LONGITUDE, LATITUDE, INST_REFERENCE
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
%
% [ C , Lon , Lat , Day ] = LOAD_TRAJ( ... )
%
% returns the Vectors for Lon, Lat, Day separated by NaN's
%
%------------------------------------------------------------
%
% See also: WHICHFILE, ASSIGN_CDF, LOOK_CDF
%

%********************************************************************

suf = '_traj.nc';

c = struct( 'name' , { '' } , ...
            'pfad' , { '' } , ...
            'inst' , { '' } , ...  % INST_REFERENCE
            'len'  , { 0  } , ...
            'day'  , { [] } , ...
            'lon'  , { [] } , ...
            'lat'  , { [] }       );

vv  = { 'JULD' 'LONGITUDE' 'LATITUDE' };

ref = 'INST_REFERENCE';

bd = datenum(1950,01,01);  % BaseDay for JULD

d = [];
x = [];
y = [];

%********************************************************************

Nin  = nargin;
Nout = nargout;

if Nin < 1
   pfad = {'.'};
else
   [ok,pfad] = chkcstr(pfad);
   if ~ok
       error('PathList must be Strings.');
   end
end

if Nin < 2
   area = [];
elseif ~( isempty(area) | isequal(size(area),[1 4]) )
   error('Area must have 4 Elements: [ LonMin LonMax LatMin LatMax ].');
end

if Nin < 3
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

nv = prod(size(vv));  % Number of Variables to extract
sv = cell(nv,1);

vv = cat( 2 , vv , {ref} );

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
 
        if ~isequal(sv{:})

            fprintf(1,'invalid size');

        elseif prod(sv{1}) == 0 

            fprintf(1,'empty variables');

        elseif ~( prod(sv{1}) == max(sv{1}) )

            fprintf(1,'no vectors');

        else


            xx = v.LONGITUDE;
            yy = v.LATITUDE;

            ok(ii) = isempty(area);

            if ~ok(ii)
                ok(ii) = any( ( area(1) <= xx ) & ( xx <= area(2) ) & ...
                              ( area(3) <= yy ) & ( yy <= area(4) )       );
            end

            if ~ok(ii)

                fprintf(1,'not in area');

            else

                c(ii).lon = xx(:);
                c(ii).lat = yy(:);
                c(ii).day = v.JULD(:) + bd;
                
                c(ii).name = strrep(f{ii},suf,'');
                c(ii).pfad = p{ii};

                c(ii).len  = prod(sv{1});

                c(ii).inst = getfield( v , ref );

                fprintf(1,'%4.0f records',c(ii).len);

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

if isempty(c) | ( Nout < 2 )
   return
end

n = size(c,1);

d = cell(2,n); d(1,:) = { c.day }; d(2,:) = {NaN};
d = cat(1,d{:});

x = cell(2,n); x(1,:) = { c.lon }; x(2,:) = {NaN};
x = cat(1,x{:});

y = cell(2,n); y(1,:) = { c.lat }; y(2,:) = {NaN};
y = cat(1,y{:});
