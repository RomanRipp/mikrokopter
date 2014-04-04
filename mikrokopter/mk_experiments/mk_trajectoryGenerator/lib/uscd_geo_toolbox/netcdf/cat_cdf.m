function [ Msg, ReadMsg, WriteMsg ] = cat_cdf(filename,LoopDim,ncfiles,LoopFiles,NoVar)


% CAT_CDF  Concatenates NetCDF-Files
%
% [ Msg, ReadMsg, WriteMsg ] = CAT_CDF( FileName, LoopDim, NC_Files, LoopFiles, NoVar );
%
%  FileName is a existing NetCDF-File, to cat Data from NC_Files into.
%
%  NC_Files are the FileNames of the NetCDF-Files to concatenate.
%
% LoopDim is a 2-Column-CellArray with 
%          the DimensionNames (1. Column)  and 
%          the    StartValues (2. Column)
%         along this Dimensions the Data from NC_Files will be filled 
%          into filename, starting at the StartValue.
%         In NC_Files this Dimensions must have the Length 1 !!!!!!
% 
% LoopFiles is a 2-Column-CellArray where the RowNumber corresponds to
%            NC_Files, with the Name of the LoopDimension in the 1. Column and
%            the Value of this LoopDimensions in the 2. Column
%              
%
% NoVar   Names of Variables from NC_Files, 
%           which data should NOT to be cat !!!
%
%
% Note:  The NetCDF-File, specified by FileName, must have the full Length
%          of LoopDimensions, to fill the Data in!!!
%
%
%       Msg is the ErrorMessage
%   ReadMsg is the ErrorMessage using  READ_CDF( NC_Files )
%  WriteMsg is the ErrorMessage of write the Data into FileName
%
% See also:  LOOK_CDF, READ_CDF, WRITE_CDF , ADD_CDF , WRT_CDF
% 

Nin = nargin;

      Msg = '';
  ReadMsg = cell(0,2);
 WriteMsg = cell(0,2);

nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);
nl2 = char([ 10 32*ones(1,nm0+2) ]);
nl4 = char([ 10 32*ones(1,nm0+4) ]);


%*****************************************************************
% Check Inputs

if Nin < 3
 Msg = [ Msg0 'Not enough InputArguments.' ];
 return
end


%-----------------------------------------------------------------

if ~ischar(filename)
  Msg = [ 'FILENAME must be a String.' ];
end


%-----------------------------------------------------------------

if ischar(LoopDim)
  LoopDim = cellstr(LoopDim);
end

if iscellstr(LoopDim)
  LoopDim = LoopDim(:);
end

if ~iscellstr(LoopDim(:,1))

  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          '1. Column of LoopDim must contain Strings.'  ];

else

  if size(LoopDim,2) < 2
    LoopDim      = LoopDim(:,[1 1]);
    LoopDim(:,2) = { 1 };            % StartValues
  else

    try
      dl = cat(1,LoopDim{:,2});
      ok = ( isnumeric(dl)  &  ( size(dl,1) == size(LoopDim,1) ) );
    catch 
      ok = 0;
    end

    if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
              '2. Column of LoopDim must contain single numerics.'  ];
    end

  end
end

%-----------------------------------------------------------------

if ischar(ncfiles)
  ncfiles = cellstr(ncfiles);
end

if ~iscellstr(ncfiles)
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'ncfiles must be a CharArray or CellStringArray.'  ];
end

%-----------------------------------------------------------------

if ischar(LoopFiles)
  LoopFiles = cellstr(LoopFiles);
end

if iscellstr(LoopFiles)
  LoopFiles = LoopFiles(:);
end

if ~iscellstr(LoopFiles(:,1))

  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          '1. Column of LoopFiles must contain Strings.'  ];

else

  if size(LoopFiles,2) < 2
    LoopFiles      = LoopFiles(:,[1 1]);
    LoopFiles(:,2) = { 1 };            % StartValues
  else

    try
      dl = cat(1,LoopFiles{:,2});
      ok = ( isnumeric(dl)  &  ( size(dl,1) == size(LoopFiles,1) ) );
    catch 
      ok = 0;
    end

    if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
              '2. Column of LoopFiles must contain single numerics.'  ];
    end

  end
end

%-----------------------------------------------------------------

if Nin < 5
 NoVar = {''};
end

if isempty(NoVar)
 NoVar = {''};
end

if ischar(NoVar)
 NoVar = cellstr(NoVar);
end

if ~iscellstr(NoVar)
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'NoVar must be a CharArray or CellStringArray.'  ];
end


%-----------------------------------------------------------------

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


%*****************************************************************

LoopStart = cat(1,LoopDim{:,2});

  LoopDim = LoopDim(:,1);

    NLoop = size(LoopDim,1);

%-----------------------------------------------------------------
% Check LoppFiles with LoopDim

nlf = size(LoopFiles,1);

ok = zeros( nlf , 1 );

for ll = 1 : NLoop

  ok = ( ok  |  strcmp( LoopDim{ll} , LoopFiles(:,1) ) );

end

if any(~ok)
  Msg = [ Msg0  'DimensionNames in 1. Column of LoopFiles ' ...
                'must correspond with LoopDim.' ];
  return
end


%*****************************************************************
% Read Data from filename

[Msg, DIM, VAR, ATT] = look_cdf( filename );
 
if ~isempty(Msg)
 Msg = [Msg0 'Error using LOOK_CDF( '  filename ' )'  nl  Msg ];
 return
end


%*****************************************************************
% Search for LoopDimension 


      ok = zeros(NLoop,1);
LoopDimN = zeros(NLoop,1);  % Actual LoopDimLength

for ll = 1 : NLoop
  ok(ll) = any( strcmp(LoopDim{ll},DIM(:,1)) );
  if ok(ll)
      dd = find( strcmp(LoopDim{ll},DIM(:,1)) );
      dd = dd(1);
     LoopDimN(ll)   = DIM{dd,2};
          DIM{dd,2} = DIM{dd,2} + ll*i ; 
     % Imaginary Part is Index in LoopDim !!!
  end
end

if any( ~ok )

   jj = find( ~ok );
      
   Msg = [ Msg0  'Invalid NetCDF-File: '   filename ...
                nl0 'Didn''t found LoopDimensions:'   ...
                nl4  strhcat(LoopDim(jj),', ',4,nl4) ];
   return

end

LoopStart = ones(size(VAR,1),1) * LoopStart(:)';



%*****************************************************
% Store Actual Length of VariableDimensions in CNT

CNT = cell(size(VAR,1),1);

for ii = 1 : size(VAR,1)

  CNT(ii) = {  cat(2,DIM{VAR{ii,4}+1,2})  };

end


%*****************************************************
% Open filename to Write

[FID,STAT] = ncmex('open',filename,'write');
 
if ~( ( FID > 0 ) & ( STAT == 0 ) )
  Msg = [ Msg0  'Error open NetCDF-File: ' filename ];
  return
end



  ncfiles = ncfiles(:);

NFiles = size(ncfiles,1);

 ReadMsg = cell(NFiles,2);
WriteMsg = cell(NFiles,2);

 ReadMsg(:,1) = { 1};
 ReadMsg(:,2) = {''};
WriteMsg(:,1) = { 1};
WriteMsg(:,2) = { cell(0,4) };


 
for ii = 1 : NFiles

 [msg,dim,var,att]=read_cdf(ncfiles{ii},'dim',LoopFiles(ii,[1 2]),'miss',[],'fill',[]);

 % dim = { DimName  DimLength }
 % var = { VarName VarType Ndim [dim] Nattr  []   Data  VarGetStatus}
  
 ReadMsg{ii,1} = isempty(msg);

 if ~ReadMsg{ii,1}

   ReadMsg{ii,2} = msg;

 else

   % Check if LoopDim exist

   if any( strcmp( LoopFiles{ii,1} , dim(:,1) ) )

          kk    = find( strcmp( LoopFiles{ii,1} , dim(:,1) ) );
      dim{kk,2} = 1;

      Nvar = size(var,1);

     WriteMsg{ii,2}      = cell(Nvar,4);
     WriteMsg{ii,2}(:,1) = { 1};
     WriteMsg{ii,2}(:,2) = {''};
     WriteMsg{ii,2}(:,3) = var(:,1);
     WriteMsg{ii,2}(:,4) = ncfiles(ii);
          
     for jj = 1 : Nvar

       %------------------------------------------------
       % Check VariableName

       kk = find(strcmp(var{jj,1},VAR(:,1)));
        
       
       if ~isempty(kk);

         if  any( strcmp( LoopFiles{ii,1} , dim(var{jj,4}+1,1) ) ) & ...
            ~any( strcmp(var{jj,1},NoVar) )  ;
         % NOT in NoVar

           %---------------------------------------
           % Check DimensionNames

           WriteMsg{ii,2}{jj,1} = isequal( dim(var{jj,4}+1,1) , ...
                                           DIM(VAR{kk,4}+1,1)       );
           if ~WriteMsg{ii,2}{jj,1}

               WriteMsg{ii,2}(jj,2) = { 'Dimensions must be agree.' };

           else

             %------------------------------------------- 
             % Check Length of Dimensions

             ndvar  = var{jj,3};
        
                cnt = cat( 2 , dim{var{jj,4}+1,2} );

              start = ones(1,ndvar);

              LoopNr = imag(CNT{kk}(1:ndvar));   
             is_loop = find( LoopNr  );
 
             if ~isempty(is_loop)

                start(is_loop) = LoopStart( kk , LoopNr(is_loop) );
     
               WriteMsg{ii,2}{jj,1} = all( start-1+cnt <= real(CNT{kk}(1:ndvar)) );
               if ~WriteMsg{ii,2}{jj,1}

                  WriteMsg{ii,2}(jj,2) = { 'Length of Dimensions exceeded.' };

               else

                 LoopStart( kk , LoopNr(is_loop) ) = ...
                 LoopStart( kk , LoopNr(is_loop) ) + 1;
        
  
                 status = 0;
                   msg1 = '';
                 try
                   status = ncmex( 'varput' , FID , kk-1 , start-1 , cnt , var{jj,7} );
                 catch
                   msg1 = lasterr;
                 end
 
                 WriteMsg{ii,2}{jj,1} = ( isempty(msg1)  & ...
                                          ( status ~= -1 ) );
 
                 if  ~WriteMsg{ii,2}{jj,1}
                  if status == -1
                    msg1 = [ 'NCMEX:   STATUS  -1 ' ];
                  end  
                  WriteMsg{ii,2}(jj,2) = { msg1 };
                 end

               end
               % ~isempty(is_loop)
             end
             % Check for Length of Dimensions
           end
           % Check DimensionNames
         end
         % NOT in NoVar
       end
       % Check VariableName

     end
     % jj
    
     if Nvar
        WriteMsg{ii,1} = all(cat(1,WriteMsg{ii,2}{:,1}));
     end

   end
   % Check for Length( DIM at LoopDim )  == 1 !!!!!!
   %  &  found_loop

 end
 % Read Ok

end
% ii = 1 : NFiles

status = ncmex('close',FID);

if status == -1
 Msg = [ 'Error close NetCDF-File: ' filename ];
end


% Cat ReadMsg and WriteMsg into Msg

  jj = find( [ cat(1,ReadMsg{:,1}) == 0 ]);
  if ~isempty(jj)
    Msg = [ Msg nl0(1:(end*(~isempty(Msg))))  ...
            'Invalid NetCDF-Files to concatenate:'           ...
             nl4 strhcat(ReadMsg(jj,2),nl4)       ];
  end

  jj = find( ( cat(1,WriteMsg{:,1}) == 0 ) );
  if ~isempty(jj) 

    % Messages for single Variables  
    %    {  0|1  Messg  VarName  FileName }
    msg  = cat(1,WriteMsg{jj,2});

    jj = find( cat(1,msg{:,1}) == 0 );
    if ~isempty(jj)
      Msg = [ Msg nl0(1:(end*(~isempty(Msg))))  ...
             'Invalid Data in NetCDF-Files to concatenate:' ... 
              nl4 strhcat( msg(jj,[4 3 2])' , ': ' , 3 , nl4 ) ];
   end
  end


if ~isempty(Msg)

  Msg = [ Msg0  Msg ];

end

 

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ~( ischar(str)  |  iscellstr(str) )
   error('StringArray must be a CharArray or CellStringArray.');
end

if iscellstr(str)
  str = char(str);
end

str = double(str);

    jj    = find( sum( ( str == 32 ) , 2 ) == size(str,2) );
str(jj,:) = [];
 
str = cellstr(char(str));

str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


