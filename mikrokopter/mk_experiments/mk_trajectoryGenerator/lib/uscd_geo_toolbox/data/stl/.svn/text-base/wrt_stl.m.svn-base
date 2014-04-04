function [files,init] = wrt_stl(ini,file,x,y,z,varargin)

% WRT_STL   Writes Data in Tiles to Files in STL-Format
%
% [ Files , Init ] = WRT_STL( INI , File , X , Y , Z , Format , Mode )
%
% INI = [ NRow NCol NRowZip ]
% 
% Writes NrowxNCol STL-Files, zip to ceil(NRow/NRowZip) ZIP-Files
%
% Init = [ RowNr ColNr  RowStart  RowEnd  ColStart  ColEnd ]
%
% see also: WRITE_STL
%


Nin = nargin;

if Nin < 5
   error('Not enough Input arguments.')
end

msg = {};

ok = ( isnumeric(ini) & ~isempty(ini) );
if ok
   ok = ( all( mod(ini(:),1) == 0 ) & all( ini(:) >= 0 ) );
end

if ~ok
    msg = cat( 1 , msg , {'Init must be positive Integers.'} );
end

vx = 0;
vy = 0;

try
   [vx,vy] = write_stl(file,x,y,z,varargin{:},'check');
catch
   msg = cat( 1 , msg , {lasterr} );
end


if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg);
end


%---------------------------------------------------------

files = {};

[pfd,name,ext] = fileparts(file);

file = fullfile(pfd,name);

if isempty(ext)
   ext = '.stl';
end

zxt = '.zip';

ini = ini(:)';

if size(ini,2) == 1
   ini = cat( 2 , ini , 1 );
end

if size(ini,2) == 2
   ini = cat( 2 , ini , 0 );
end

ini([1 2]) = max( ini([1 2]) , 1 );

zip = ( ini(3) > 0 );

ny = size(z,1);
nx = size(z,2);

ni = floor( [ ny  nx ] ./ ini([1 2]) ) - 1;

y01 = [ 0  ni(1) ];
x01 = [ 0  ni(2) ];

sc = ini([2 2]);
nn = prod(ini([1 2]));
 
%------------------------------------------------------------
% Write STL-Files

str = sprintf('Write Faces for "%s"',name)';
if zip
   str = sprintf('%s (ZIP each %.0f. row)',str,ini(3));
end

cl = loopdot(sc,nn,str);

init = zeros(nn,6);

for jy = 1 : ini(1)

    iy = (jy-1)*ni(1) + y01 + 1;

    if jy == ini(1)
       iy(2) = ny;
    end

    iy = ( iy(1) : iy(2) );

    for jx = 1 : ini(2)

        ix = (jx-1)*ni(2) + x01 + 1;
        if jx == ini(2)
           ix(2) = nx;
        end

        ix = ( ix(1) : ix(2) );

        jj = (jy-1) * ini(2) + jx;

        stl_file = sprintf('%s_%2.2d%s',file,jj,ext);

        init(jj,:) = [ jy jx iy([1 end]) ix([1 end]) ];

        if 1 % one_stl
           fprintf(1,'\r%2.2d %2.2d %2.2d %4.0f %4.0f %4.0f %4.0f\n',jj,init(jj,:));
        end

        %----------------------------------------------------------- 
        % STL

        if vx, xx = x(ix); else, xx = x(iy,ix); end
        if vy, yy = y(iy); else, yy = y(iy,ix); end

        write_stl( stl_file , xx , yy , z(iy,ix) , varargin{:} );

        %----------------------------------------------------------- 
        % ZIP
        %----------------------------------------------------------- 
        if zip
        %----------------------------------------------------------- 

           zn = ceil( jy / ini(3) );

           zip_file = sprintf('%s_%2.2d%s',file,zn,zxt);

           new = ( ( ( mod(jy,ini(3)) == 1 ) | ( ini(3) == 1 ) ) & ( jx == 1 ) );

           opt = '';
           if ( exist(zip_file,'file') == 2 )
              if new
                 delete(zip_file)
                 files = cat( 1 , files , zip_file );
              else
                 opt = '-u';
              end
           end

           cmd = sprintf('zip %s %s %s',opt,zip_file,stl_file);

           [s,w] = unix(cmd);

           if s == 0 
              delete(stl_file);
           end

        %----------------------------------------------------------- 
        else
        %----------------------------------------------------------- 

           files = cat( 1 , files , stl_file );
  
        %----------------------------------------------------------- 
        end
        %----------------------------------------------------------- 

        loopdot(sc,nn,jj,1,cl);

    end

end
