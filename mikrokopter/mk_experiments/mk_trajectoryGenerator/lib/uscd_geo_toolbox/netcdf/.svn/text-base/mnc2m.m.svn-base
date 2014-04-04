function  mnc2m(ModelFile,M_File,lid);

% MNC2M   Writes Matlab-Script-File to create NetCDF-File from ModelFile, using WRT_CDF
%
%  MNC2M( Model_NetCDF_File , M_File , LogID )
%
%  default: M_File = <ModeFileName>.m
%            LogID = 1   LogOutput to Matlab-CommandWindow
%
% see also: LOOK_CDF
%


Nin = nargin;

%---------------------------------------------------
% Check Inputs

if Nin < 1
   error('Input Model_NetCDF_File is missing.')
end

if Nin < 2
   M_File = '';
end

if isempty(M_File)
   [p,n,e] = fileparts(ModelFile);
   M_File = cat( 2 , fullfile(p,n) , '.m' );
end

if ~( ( ischar(ModelFile) & ~isempty(ModelFile) & ...
        ( prod(size(ModelFile)) == size(ModelFile,2) )  )  & ...
      ( ischar(M_File) & ...
        ( prod(size(M_File)) == size(M_File,2) ) )  )

   error('FileNames must be Strings.')

end
     
if Nin < 3
   lid = 1;
end

try
  fprintf(lid,'');
catch
  lid = 1;
end

%---------------------------------------------------
% Check M_File

if isempty(M_File)
   M_File = ModelFile;
end

[pfad,name,ext] = fileparts(M_File);

[Msg,name,M_File] = chk_name(name,pfad,'.m','','new');

if isempty(Msg) & isempty(M_File)
   Msg = 'Invalid M_File.';
end

if ~isempty(Msg)
   error(Msg)
end

%---------------------------------------------------
% Read Informations about Dimensions and Variables


fprintf(lid,'\n');
fprintf(lid,'Read NetCDF-ModelFile: %s ... ',strrep(ModelFile,'\','\\'));


[msg,dim,var] = look_cdf(ModelFile);

if ~isempty(msg)
    fprintf(lid,'error\n\n')
    error(msg)
end

fprintf(lid,'ok\n\n')
    
% dim = { DimName  DimLength  UnLimited    }
% var = { VarName VarType Ndim [dim] Nattr }

%---------------------------------------------------

fprintf(lid,'Write M-File: %s ...',strrep(M_File,'\','\\'));
 
fid = fopen(M_File,'wt');

if fid == -1
   fprintf(lid,'error\n\n')
   error(['Can''t open M_File: ' outfile ])
end

%--------------------------------------------------------------------------
% Help

[p,n,e] = fileparts(ModelFile);
  mname = upper(cat(2,n,e));

[p,n,e] = fileparts(M_File);
  fname = upper(n);

     bl = char(32*ones(1,size(fname,2)));

fprintf(fid,'%% %s   Creates NetCDF-File from ModelFile "%s", using WRT_CDF\n',fname,mname);
fprintf(fid,'%%\n');
fprintf(fid,'%% %s   Modify the Value  for "OutFile", the NetCDF-File to create,\n',bl);
fprintf(fid,'%% %s          the Length of the Dimensions in second column of "dim",\n',bl);
fprintf(fid,'%% %s          the Value  of the Variables  in second column of "var".\n',bl);
fprintf(fid,'%%\n\n');
fprintf(fid,'%%----------------------------------------------------------------------\n');

%--------------------------------------------------------------------------
% Write Files to use

fprintf(fid,'%% Modify the Value for "OutFile", the NetCDF-File to create\n\n');
fprintf(fid,'  OutFile = ''%s'';\n\n\n' , 'test.nc' );

fprintf(fid,'%% Never modify the Value for "ModelFile", the NetCDF-DefaultFile to use\n\n');
fprintf(fid,'ModelFile = ''%s'';\n\n' , strrep(ModelFile,'\','\\') );
 
%--------------------------------------------------------------------------
% Write Dimensions

fprintf(fid,'%%-------------------------------------------------------------\n');
fprintf(fid,'%% Definition of Dimensions, modify the Length in second Column\n\n');

fprintf(fid,'%% dim = {  Name   Length  } \n\n');
fprintf(fid,' dim = { ... \n');


dl = size(char(dim(:,1)),2);   % Maximum Length of DimensionNames

bl1 = char(32*ones(1,8));      % First Blanks

for ii = 1 : size(dim,1);

    bl2 =  char( 32 * ones( 1 , 4 + dl-size(dim{ii,1},2) ) );

    fprintf( fid , '%s%s%s%s%s%8.0f\n' , bl1 , '''' , dim{ii,1} , '''' , ...
                                         bl2 ,        dim{ii,2} );

end

fprintf( fid , '%s};\n\n' , bl1 );

%--------------------------------------------------------------------------
% Write Variables

fprintf(fid,'%%-----------------------------------------------------------\n');
fprintf(fid,'%% Definition of Variables, modify the Value in second Column\n\n');
fprintf(fid,'%% Note: The Order of Dimensions are reverse\n');
fprintf(fid,'%%       Variables can be missing or EMPTY\n');
fprintf(fid,'%%       The Length of a Variable in a Dimension can be less\n');
fprintf(fid,'%%        than the defined DimensionLength\n\n');

fprintf(fid,'%% var = {  Name   Value  }    %% Requested DimensionOrder\n\n');
fprintf(fid,' var = { ... \n');

vl = size(char(var(:,1)),2);   % Maximum Length of VariableNames

bl1 = char(32*ones(1,8));      % First Blanks
bl3 = char(32*ones(1,4));      % Third Blanks

for ii = 1 : size(var,1);

    bl2 =  char( 32 * ones( 1 , 4 + vl-size(var{ii,1},2) ) );

    % Dimensions for Variable
    dv = '<SCALAR>';
    if var{ii,3} > 0

       id = var{ii,4}(var{ii,3}:-1:1) + 1;   % Reverse Order

       dv = strhcat( dim(id,1) , ', ' );
   
    end

    fprintf( fid , '%s%s%s%s%s%s%s%% %s\n'        , ...
                    bl1 , '''' , var{ii,1} , '''' , ...
                    bl2 , '[]' , bl3 , dv  );

end

fprintf( fid , '%s};\n\n' , bl1 );

%--------------------------------------------------------------------------
% Create NC_File

fprintf(fid,'%%------------------------------------------------------\n');
fprintf(fid,'%% Write NetCDF-File using WRT_CDF\n\n');

fprintf(fid,'Msg = wrt_cdf( OutFile , ModelFile  , dim , var );\n\n');

fprintf(fid,'if ~isempty(Msg)\n');
fprintf(fid,'   error(Msg);\n');
fprintf(fid,'end\n');

%--------------------------------------------------------------------------

fclose(fid);

fprintf(lid,'ok\n\n')
