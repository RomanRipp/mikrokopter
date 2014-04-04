function [file,fid] = cdf_file(file,ext);

% CHecks for NetCDF-File


fid = [];

if Nin == 0
   file = '';
   return
end

if Nin < 2
   ext = { ''  '.cdf'  '.nc'  '.cnf' };
else
   [ok,ext] = chkcstr(ext);
   if ~ok
       error('Following Input must be a String or CellArray of Strings.');
   end
end

for e = ext

    ff = cat(2,file,e{1});

    [ff,fid] = chkcdf(file);

    ok = ~isempty(fid);         % File exist
    if ok
       ok = ~isnan(fid);        % File readable
       if ok
          ok = ~( fid == -1 );  % NetCDF-File
          if ~ok

           % Try to load KeyWord "FileName" from ConfigFile, 
           %   but NOT if called by CNF_CDF
              ff = '';
              try   % NetCDF-ConfigFile
                 [m,v,ff] = load_asc(file,'FileName',{'#' ':'},'%');
              end
              ok = chkstr(ff,1);
              if ok
                 [ff,fid] = cdf_file(file);
              end

          end
       end
    end

    if ok
       file = ff;
       break
    end

end
 
%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [file,fid] = chkcdf(file)

% [ FileName , FID ] = CHKCDF( File )
%
%  FID    empty if File not exist
%         NaN   if File cann't be open
%         -1    if not a NetCDF-File, using NCMEX
%         ID    otherwise

fid = [];

if ~( exist(file,'file') == 2 );
    return
end

%---------------------------------------
% Get Full FileName

ff = which(file);
if ~isempty(ff)
    file = ff;
end

%---------------------------------------
% Try to Open File

fid = fopen(file,'r');

if fid == -1
   fid = NaN;
   return
end

fclose(fid);

%---------------------------------------
% Try to Open NetCDF-File

[fid,stat] = ncmex('open',file,'nowrite');

if ~( ( fid > 0 ) & ( stat == 0 ) )
    fid = -1;
end

