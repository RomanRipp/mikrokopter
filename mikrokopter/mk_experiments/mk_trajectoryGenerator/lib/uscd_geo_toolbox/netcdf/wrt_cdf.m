function  Msg = wrt_cdf( filename , defaultfile , dims , vars , atts );

% WRT_CDF  Writes Data into a new NetCDF-File, using a DefaultFile
%
% Msg = WRT_CDF( FileName , DefaultFileName  , DIM , VAR , ATT );
%
%  DefaultFileName  is the Default-NetCDF-File
%
%  A new NetCDF-File with the Name "FileName" is created, using 
%    the Dimensions, Variables and Attributes of the DefaultFile.
%
%  DIM : New Definition for N Dimensions
%  DIM = { DimName  DimLength }       [ N by 2 ] CellArray
%
%  VAR : New Definition for M Variables
%  VAR = { VarName  VarValue  }       [ M by 2 ] CellArray 
%
%  ATT : New Definition for M+1 Attributes (single or last are global Att)
%  ATT = { { AttName  AttValue } } [ [M][+1] by 1 ] CellArray 
%
%  Msg  contains ErrorMessages, empty if all was succesfull.
%
%  If new Dimensions are defined by the Input DIM, 
%     WRT_CDF doesn't copies the Values of the Variables of the DefaultFile.
%
%  WRT_CDF( FileName , DefaultFileName ) copies the entire DefaultFile.
%
%  WRT_CDF( FileName , DefaultFileName , VAR , ATT ) copies the entire 
%   DefaultFile, but replace the Variables and Attributes like defined 
%   by the Inputs VAR and ATT
%
% see also:  WRITE_CDF  LOOK_CDF  CREATE_CDF  ADD_CDF  CAT_CDF
%
                     

Nin = nargin;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

nl = char(10);

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);
nl2 = char([ 10 32*ones(1,nm0+2) ]);
nl4 = char([ 10 32*ones(1,nm0+4) ]);


Msg = '';

%-----------------------------------------------------------------
% Check Inputs

if Nin < 2
 Msg = [ Msg0 'Not enough InputArguments.' ];
 return
end


if ~ischar(filename)
  Msg = [ 'FILENAME must be a String.' ];
end

if ~ischar(defaultfile)
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
         'DefaultFileName must be a String.' ];
end

if Nin == 2
   dims = {};
   vars = {};
end

if Nin < 5
   atts = {};
end

if     ( Nin == 3 )
       vars = dims;
       dims = {};
elseif ( Nin == 4 ) & ( size(vars,2) == 1 )
       atts = vars;
       vars = dims;
       dims = {};
end

if isempty(dims)
   dims = cell(0,2);
end

if isempty(vars)
   vars = cell(0,2);
end


if ~iscell(dims) | ~iscell(vars) 
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'Inputs DIM and VAR must be CellArrays.'];
end

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


if ( size(dims,2) < 2 )  |  ( ndims(dims) ~= 2 )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
           'CellArray DIM must contains 2 Columns.'  ...
            nl4 'DIM = { DimName  DimLength }'           ];
elseif ~isempty(dims)
  if ~iscellstr(dims(:,1))
    Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
            '1. Column of CellArray DIM must contain Strings.' ];
  end
end

if ( size(vars,2) < 2 )  |  ( ndims(vars) ~= 2 )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
          'CellArray VAR must contains 2 Columns.'  ...
          nl4 'VAR = { VarName Value  } '   ];
elseif ~isempty(vars)
  if ~iscellstr(vars(:,1))
    Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
            '1. Column of CellArray VAR must contain Strings.'  ];
  end
end


ok = isempty(atts);

if  ok
    atts = { {} };
else
    s1 = size(atts,1);
    s2 = size(atts,2);
    ok = ( any( s1 == [ 1 size(vars,1) ] ) & ...
           any( s2 == [ 1  2 ] )  &  ( ndims(atts) == 2 ) );
    if ok 
       if ( s2 == 1 )
          try
             at = cat(1,atts{:});
          catch
             ok = 0;
          end
          if ok & ~isempty(at)
             ok = ( ( size(at,2) == 2 )  &  ( ndims(at) == 2 ) );
          end
       else
          at = atts;
          for ii = ( 1 : s1 )
              atts{ii,1} = atts(ii,:);
          end
          atts = atts(:,1);  
       end
    end

    if ~ok
        Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
               'Elements CellArray ATT must contains 2 Columns ATT = { { AttName Value } }.'  nl0 ...
        'Number of Elements in CellArray ATT must correspond to Number of Variables and/or Global' ];
    else
      ok = iscellstr(at(:,1));
      if ~ok
          Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
                 '1. Column of Elements in CellArray ATT must contain Strings.'  ];
      end
      for ii = 1 : size(at,1)
          kk = chkstr(at{ii,2});
          if ~kk
              kk = isnumeric(at{ii,2});
          end
          if ~kk
              ok = 0;
              Msg = [ Msg nl0(1:(end*(~isempty(Msg))))          ...
                      '2. Column of Elements in CellArray ATT must contain Strings or Numerics.'  ];
              break
          end
      end
    end

end

na = size(atts,1);
nv = size(vars,1);

if ok & ( na < nv+1 )

   atts = cat( 1 ,  atts , cell(nv+1-na,1) );
   atts( (na:nv) + 1 ) = { {} };

end


if ~isempty(Msg)
    Msg = [ Msg0 Msg ];
    return
end

%-------------------------------------------------------
% Read defaultfile
 

[msg,dim,var,att]=look_cdf(defaultfile);
 
if ~isempty(msg)
    Msg = [ Msg0 'Error Read default NetCDF-File: '  defaultfile nl msg ];
    return
end


 % Check for Unlimited Dimension
 ul = find(strcmp(lower(dim(:,3)),'unlimited'));
 if ~isempty(ul)
     dim(ul,2) = { 'unlimited' };
 end

 dim = dim(:,1:2);
 var = var(:,1:5);

 var(:,5) = {[]};


%--------------------------------------------------------
% Check for Dimensions in  dim

msg = '';

if ~isempty(dims)

    dim_ok = zeros(size(dims,1),1);

    for vv = 1:size(dims,1)

        jj = find(strcmp(dims{vv,1},dim(:,1)) );

        dim_ok(vv) = ~isempty(jj);
     
        if dim_ok(vv)
           dim{jj,2} = dims{vv,2};
        end

    end
 
    jj = find(~dim_ok);
    if ~isempty(jj)

        msg = [ msg nl2(1:(end*(~isempty(msg))))       ...
                'Didn''t found following Dimensions:'  ...
                nl4 strhcat(dims(jj,1),', ',4,nl4)          ];
 
    end

end

%------------------------------------------------------------------
% Check for Variables in var and Attributes

mv = size(var,1);

% True if Variable to read from DefaultFile
rv = isempty(dims) * ones(mv,1);  

if ~( isempty(vars) & isempty(atts) )

    nv = size(vars,1);

    var_ok = zeros(nv+1,1);    % One more for global Atts

    var_ok(nv+1) = 1;          % Global is true by default

    for vv = 1 : size(var_ok,1)

        if ~var_ok(vv)
            jj = find(strcmp(vars{vv,1},var(:,1)) );
            var_ok(vv) = ~isempty(jj);
            if var_ok(vv)
               var{jj,5} = vars{vv,2};
                rv(jj)   = 0;  % Do NOT read from DefaultFile
            end
        else
            jj = size(var,1) + 1;
        end

        if var_ok(vv)
           nat = atts{vv};           % New Attribute
           vat =  att{jj};           % Existing Attributes
           if ~isempty(nat)
               if isempty(vat)
                  nat = nat(:,[1 2 2]);
                  for ii = 1 : size(nat,1)
                      [nat{ii,2},nat{ii,3}] = attype(nat{ii,2});
                  end
                  att{jj} = nat;
               else
                  for ii = ( 1 : size(nat,1) )
                      kk = find(strcmp(nat{ii,1},vat(:,1)));
                      if isempty(kk)
                         [typ,val] = attype(nat{ii,2});
                         vat = cat( 1 , vat , { nat{ii,1} typ val '' } );
                      else
                         vat(kk,3) = nat(ii,2);
                      end
                  end
                  att{jj} = vat;
              end
           end
        end

    end

    jj = find(~var_ok);

    if ~isempty(jj)
  
        msg = [ msg nl2(1:(end*(~isempty(msg))))       ...
                'Didn''t found following Variables:'    ...
                nl4  strhcat(vars(jj,1),', ',4,nl4)        ];
 
    end

end

if ~isempty(msg)
    Msg = [ Msg0 ' Invalid defaultfile: ' ...
             defaultfile  nl2  msg ];
    return
end


if any(rv)
   ii = find(rv);
   [m,d,v] = read_cdf(defaultfile,'var',var(ii,1));
   var(ii,5) = v(:,7);
end


 [ CreateMsg , WriteMsg , VarMsg , dstat, vstat, astat ] =  ...
     write_cdf(filename,dim,var,att);

  Msg = [ CreateMsg  nl(1:(end*(~isempty(CreateMsg)))) ...
          WriteMsg ];

  if ~isempty(Msg)
    Msg = [Msg0 'Error write NetCDF-File: '  filename nl ...
           Msg ];
  end
  


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [typ,val] = attype(val)

%   DataTypes:  1. 'byte'      [        -128        127 ]      'int8'
%               2. 'char'      [           0        255 ]     'uint8'
%               3. 'short'     [      -32768      32767 ]     'int16'
%               4. 'long'      [ -2147483648 2147483647 ]     'int32'
%               5. 'float'       32 bit floating point
%               6. 'double'      64 bit floating point
  

ok = chkstr(val);

if ok
   typ = 'char';
elseif classchk(val,'int8')
   typ = 'byte';
elseif classchk(val,'int16')
   typ = 'short';
elseif classchk(val,'int32')
   typ = 'long';
else
   typ = 'float';
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


