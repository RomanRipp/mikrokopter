function Msg = add_cdf(filename,dim,var);

% ADD_CDF  Writes Data to NetCDF-File
%
% Msg = ADD_CDF( FileName , DIM , VAR )
%
%  DIM : New Definition for N Dimensions
%  DIM = { DimName  DimLength }       [ N by 2 ] CellArray
%
%  VAR : New Definition for M Variables
%  VAR = { VarName  VarValue  }       [ M by 2 ] CellArray 
%
%
% Msg contains ErrorMessages.
%
%
% Reads Data from NetCDF-File, specified by FileName,
%  resets the Length of Dimensions and Data of Variables
%  and write a New NetCDF-File, using WRITE_CDF
%
% See also: LOOK_CDF, READ_CDF, WRITE_CDF, CREATE_CDF, CAT_CDF
%
 
Nin = nargin;

Msg  = '';
nl   = char(10);

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



%-----------------------------------------------------------------
% Check Inputs

if Nin < 3
 Msg = [ Msg0 'Not enough InputArguments.' ];
 return
end


if ~ischar(filename)
  Msg = [ 'FILENAME must be a String.' ];
end


if ~iscell(dim) | ~iscell(var) 
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'Inputs DIM and VAR must be CellArrays.'];
end

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


if ( size(dim,2) < 2 )  |  ( ndims(dim) ~= 2 )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
           'CellArray DIM must contains 2 Columns.'  ...
            nl4 'DIM = { DimName  DimLength }'           ];
else
  if ~iscellstr(dim(:,1))
    Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
            '1. Column of CellArray DIM must contain Strings.' ];
  end
end

if ( size(var,2) < 2 )  |  ( ndims(var) ~= 2 )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'CellArray VAR must contains 2 Columns.'  ...
          nl4 'VAR = { VarName Value  } '   ];
else
  if ~iscellstr(var(:,1))
    Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
            '1. Column of CellArray VAR must contain Strings.'  ];
  end
end


if ~isempty(Msg)
 Msg = [ Msg0 Msg ];
 return
end


%*****************************************************************
% Read Data from CDF_File

[Msg, DIM, VAR, ATT] = read_cdf( filename );
 
if ~isempty(Msg)
 Msg = [Msg0 'Error using READ_CDF( '  filename ' )'  nl  Msg ];
 return
end

VAR = VAR(:,[1 2 3 4 7]);   % Data in 5. Column now


%*****************************************************************
% Check    dim  and  var  
% Fill up  DIM  and  VAR


dim_ok = zeros(size(dim,1),1);
for ii = 1 : size(dim,1)
  dim_ok(ii) = any(strcmp(dim{ii,1},DIM(:,1)));
  if dim_ok(ii)
     jj = find(strcmp(dim{ii,1},DIM(:,1)));
     DIM{jj,2} = dim{ii,2};
  end
end

if any(~dim_ok)
  jj = find(~dim_ok);
  Msg = [ nl2 'Didn''t found following Dimensions: ' ...
          nl4 strhcat(dim(jj,1),', ',4,nl4) ];
end


var_ok = zeros(size(var,1),1);
for ii = 1 : size(var,1)
  var_ok(ii) = any(strcmp(var{ii,1},VAR(:,1)));
  if var_ok(ii)
     jj = find(strcmp(var{ii,1},VAR(:,1)));
     VAR{jj,5} = var{ii,2};
  end
end

if any(~var_ok)
  jj = find(~var_ok);
  Msg = [ Msg  nl2(1:(end*(~isempty(Msg)))) ...
          'Didn''t found following Variables: ' ...
          nl4 strhcat(var(jj,1),', ',4,nl4) ];
end

if ~isempty(Msg)
  Msg = [ Msg0  'Invalid Data in NetCDF-File: ' filename ...
          Msg ];
  return
end




%*****************************************************************
% Create New File and store 'Append'-Data
 
 [ CreateMsg , WriteMsg , VarMsg , dstat, vstat, astat ] =  ...
     write_cdf(filename,DIM,VAR,ATT);

  Msg = [ CreateMsg  nl(1:(end*(~isempty(CreateMsg)))) ...
          WriteMsg ];

  if ~isempty(Msg)
    Msg = [ Msg0 'Error write NetCDF-File: '  filename nl ...
           Msg ];
    return
  end
  

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


