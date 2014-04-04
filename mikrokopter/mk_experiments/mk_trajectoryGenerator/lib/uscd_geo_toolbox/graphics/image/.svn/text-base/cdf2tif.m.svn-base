function cdf2tif(cnf,tif,varargin)

% CDF2TIF  converts NetCDF-ImageData to TIFF-File, using UNIX: PPM2TIFF
%           
% CDF2TIF( CNF , TIF , [DataVar] , [Options] , [Brighten] , [Range] )
%
% CNF       ConfigFile for NetCDF-File for using LOAD_CDF
% TIF       TIFF-File to create
%
% DataVar   DataVariable in CNF-File, which specifies the VarialeName
%            of the ImageData in the NetCDF-File.
%           DimLength in Z must be 3 !!!
%
% Options   Options for PPM2TIFF, String, see below for more info
%            The first Character must be  "-"   !!!
%                
% Brighten  Potenz to bright Colors
%
% Range     Range of Dimensions for X and Y to extract Data for conversion
%
%           2-Column CellArray, first Column: 'X' or 'Y'  DimensionCharacter
%                              second Column: Stride
%                                            [ Start Count ]  (Stride == 1 )
%                                            [ Start Count Stride ]
%           Note: NetCDF starts at ZERO !!!
%
% see also: CNF_CDF, LOAD_CDF, WRTPPM
% 
%------------------------------------------------------------
% PPM2TIFF
%------------------------------------------------------------
%
% usage: ppm2tiff [options] input.ppm output.tif
% where options are:
%  -r #           make each strip have no more than # rows
%  -R #           set x&y resolution (dpi)
% 
%  -c jpeg[:opts]  compress output with JPEG encoding
%  -c lzw[:opts]  compress output with Lempel-Ziv & Welch encoding
%                (no longer supported by default due to Unisys patent enforcement)
%  -c zip[:opts]  compress output with deflate encoding
%  -c packbits    compress output with packbits encoding
%  -c none        use no compression algorithm on output
%
% JPEG options:
%  #              set compression quality level (0-100, default 75)
%  r              output color image as RGB rather than YCbCr
% LZW and deflate options:
%  #              set predictor value
%
%------------------------------------------------------------
%

dat = '';
opt = '';
scl = 1;

Nin = nargin;

if ~isunix
    error('UNIX-Enviroment required.');
end

if Nin < 2
   error('Not enough InputArguments.');
end

if ~( chkstr(cnf,1) & chkstr(tif,2) )
    error('First to Inputs must be nonempty Strings.');
end

opt = cell(0,1);
dim = cell(0,2);

for ii = 1 : Nin-2
    v = varargin{ii};
    if chkstr(v,1)
       if v(1) == '-'
          opt = cat(1,opt,{v});
       else
          dat = v;
       end
    elseif isnumeric(v) & ~isempty(v)
       scl = v(1);
    elseif iscell(v)
       if ~( size(v,2) == 2 )
           error(sprintf('CellArray of %.0f. Input must have 2 Columns.',ii+2));
       elseif ~chkcstr(v(:,1))
           error(sprintf('1. Column of CellArray of %.0f. Input must be Strings.',ii+2));
       else
           v(:,1) = lower(v(:,1));
           for jj = 1 : size(v,1)
                w = v{jj,2};
                p = prod(size(w));
               ok = ( any(v{jj,1}=='xy') & any(p==[1 2 3]) & isnumeric(w) );
               if ~ok
                   error(sprintf('Invalid Elements in CellArray of %.0f. Input.',ii+2));
               else
                   if ~isempty(dim)
                       if any(strcmp(dim(:,1),v{jj,1}))
                          error(sprintf('Multiple Definition for Dimension %.0f. Input.',ii+2));
                       end
                   end
                   if p == 1
                      w = [ 0 NaN w ];  % Stride
                   elseif p == 2
                      w = [ w(:)' 1 ];  % [ Start Count ]
                   else
                      w =   w(:)';      % [ Start Count Stride ]
                   end
                   v{jj,2} = w;
                   v{jj,1} = [ v{jj,1} 'dim' ];
               end
           end
           dim = cat( 1 , dim , v );
       end
    else
       error(sprintf('Invalid %.0f. Input.',ii+2));
    end
end

opt = sprintf('%s ',opt{:});
   
%*****************************************************
% Check CNF-File

%----------------------------------------------------
% Check Variables

fprintf(1,'\nCheck CNF-File: %s ... ',cnf);

[m,v] = cnf_cdf(cnf,'check');

if isempty(m)
   m = cell(0,1);
   if ~( v{3,3} == 3 )
      m = cat(1,m,{'Z must have DimensionLength 3 (RGB).'});
   end
   if isempty(v{2,1})
       m = cat(1,m,{'No DataVariable defined.'});
   end

   if isempty(dat)
      if size(v{2,1},1) > 1 
          m = cat(1,m,{'Multiple DataVariable defined.'});
      end
      dat = v{2,1}{1};
   elseif ~any(strcmp(dat,v{2,1}(:,1)))
       m = cat(1,m,{'DataVariable not found.'});
   end
   if ~isempty(m)
       m = sprintf('%s\n',m{:});
   end
end

if ~isempty(m)
    fprintf(1,'error\n');
    error(m);
end

fprintf(1,'ok\n');

%----------------------------------------------------
% Check Requested Dimensions

w = v{1,3};
h = v{2,3};

if isempty(dim)
   w = v{1,3};
   h = v{2,3};
   dim = { 'xdim'  [ 0  w  1 ]
           'ydim'  [ 0  h  1 ] };
   dim = permute(dim,[2 1]);
else

   if size(dim,1) == 1
      if strcmp(dim(1),'xdim')
         dim = cat( 1 , dim , { 'ydim'  [ 0  h  1 ] } );
      else
         dim = cat( 1 , { 'xdim'  [ 0  w  1 ] } , dim );
      end
   end

   if isnan(dim{1,2}(2)), dim{1,2}(2) = length(1:dim{1,2}(3):w); end
   if isnan(dim{2,2}(2)), dim{2,2}(2) = length(1:dim{2,2}(3):h); end

   dim = permute(dim,[2 1]);

   fprintf(1,'\nCheck Dimensions: %s ... ',cnf);

   [m,d] = load_cdf(cnf,dim{:});

   if ~isempty(m)
       fprintf(1,'error\n');
       error(m);
   end

   fprintf(1,'ok\n');

   w = prod(size(d{1}));
   h = prod(size(d{2}));
   
end

fprintf(1,'ImageSize: %.0f x %.0f\n',w,h);

if any( [w h] == 0 ) 
   return
end

%*****************************************************
% Write PPM

tmp = [tif '.tmp.ppm'];

fprintf(1,'Open PPM-File: %s ... ',tmp);

fid = fopen(tmp,'w');

if fid == -1
   fprintf(1,'error\n');
   error(sprintf('Error open File "%s" for writing.',tmp));
end

fprintf(1,'ok\n');

fprintf(fid,'P6\n');
fprintf(fid,'%4d %4d\n',w,h);
fprintf(fid,'%4d\n',255);

dpr = 50;
cpr = h/10;                    % Max. 10 Rows
cpr = dpr * ceil( cpr / dpr );

sc = [ -dpr cpr ];

cl = loopdot(sc,h,'Write Data');

yd = dim{2,2};

for ii = 1 : h

    dim{2,2} = [ yd(1)+(ii-1)*yd(3)  1  yd(3) ];

    [m,d] = load_cdf(cnf,dim{:},'data',dat);

    if isempty(m)
       d = d{5};
       if ~( ndims(d) == 3 )
           m = 'Data must have 3 Dimensions.';
       end
    end

    if ~isempty(m)
        fclose(tmp); 
        delete(tmp);
        fprintf(1,'\n');
        error(m);
    end

    if ~( scl == 1 )
         d = rgb2hsv(d);
         d(:,:,3) = d(:,:,3) .^ scl;
         d = 255 * hsv2rgb(d);
    else
         d = double(d);
    end

    d = permute(d,[3 2 1]);

    fwrite(fid,d,'uchar');

    loopdot(sc,h,ii,0,cl);

end

fclose(fid);


%*****************************************************
% PPM2TIF

cmd = sprintf('ppm2tiff %s %s %s',opt,tmp,tif);

fprintf(1,'Call UNIX: %s ... ',cmd);

[s,w] = unix(cmd);

delete(tmp);

if ~( s == 0 )
    fprintf(1,'error\n');
    error(w);
else
    fprintf(1,'ok\n');
end
