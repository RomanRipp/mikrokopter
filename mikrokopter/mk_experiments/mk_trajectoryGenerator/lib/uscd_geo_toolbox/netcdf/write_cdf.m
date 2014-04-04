function  [msg,WrtMsg,VarMsg,ds,vs,as] = write_cdf(varargin);

% WRITE_CDF  Creates a NetCDF-File, using CREATE_CDF, and writes Data into
%
% [ CreateMsg , WriteMsg ] = WRITE_CDF( FILENAME, DIM, VAR, ATT );
%
%------------------------------------------------
% Inputs:
%
%  DIM : Definition for N Dimensions
%  DIM = { DimName  DimLength  }               [ N   by 2 ] CellArray
%
%  VAR : Definition for M Variables
%  VAR = { VarName VarType Ndim [dim] Value }  [ M   by 5 ] CellArray 
%
%  ATT : Definition for Variable- and Global-Attributes
%  ATT = { VarAtt ; GlbAtt }                   [ M+1 by 1 ] CellArray
%
%  Att = { AttName AttType AttValue }          [ K   by 3 ] CellArray
%
%------------------------------------------------
%
%  DataTypes:  1. 'byte'      [        -128        127 ]      'int8'
%              2. 'char'      [           0        255 ]     'uint8'
%              3. 'short'     [      -32768      32767 ]     'int16'
%              4. 'long'      [ -2147483648 2147483647 ]     'int32'
%              5. 'float'       32 bit floating point
%              6. 'double'      64 bit floating point
%
%------------------------------------------------
% OutPuts: 
%
%  CreateMsg   Errors detected creating, 
%                open or close the new NetCDF-File 
%
%  WriteMsg    Errors detected writing the VariableValues
%
%
% [...,VarMsg,DimStatus,VarStatus,AttStatus] = WRITE_CDF( ... )
%
%  VarMsg      Errors detected writing the VariableValues
%                for each Variable (CellArray)
%
%  DimStatus   Status of define the Dimensions
%  VarStatus   Status of define the Variables
%  AttStatus   Status of define the Attributes
%
%------------------------------------------------
%
%  see also: CREATE_CDF, WRT_CDF, LOOK_CDF, READ_CDF
%

Nin = nargin;

msg = '';

VarMsg = cell(0,2);
WrtMsg = '';

ds  = [];
vs  = [];
as  = [];


nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));

nlm = sprintf('%s%s',nl,char(32*ones(size(msg0))));

%***********************************************************
% Create NetCDF-File

file = '';

try
   [msg,ds,vs,as,file,dim,var,att] = create_cdf(varargin{:});
catch
   msg = lasterr;
end

if isempty(file) & isempty(msg)
   msg = cat(2,'Empty FileName.',nl(1:(end*(~isempty(msg)))),msg);
end

if ~isempty(msg)
   msg = sprintf('%sError call CREATE_CDF.%s%s',msg0,nlm,msg);
   return
end

%------------------------------------------------------------

nv = size(var,1);

VarMsg      = cell(nv,2);
VarMsg(:,1) = {  0 };
VarMsg(:,2) = { '' };

%------------------------------------------------------------
% Check Status for Msg

msg = cell(0,1);

jj = ( ds == -1 );
if any(jj)
   jj = find(jj);
   msg = cat(1,msg,{'Error define Dimensions:'} , ...
                   { strhcat(dim(jj,1),', ',4,nlm) }          );
end

jj = ( vs == -1 );
if any(jj)
   jj = find(jj);
   msg = cat(1,msg,{'Error define Variables:'} , ...
                   { strhcat(var(jj,1),', ',3,nlm) }          );
end

if any( cat(1,as{:}) == -1 )
   jj = find(jj);
   msg = cat(1,msg,{'Error define Attributes.'});
end

if ~isempty(msg)
    msg = sprintf('%sError create NetCDF-File.%s%s',msg0,nlm,strhcat(msg,nlm));
    return
end

%------------------------------------------------------------
% No Data ==> return

VarMsg(:,1) = { 1 };

msg = '';

if size(var,2) < 5
   return
end

%***********************************************************
% Open NetCDF-File to write

[fid,stat] = ncmex('open',file,'write');
 
if ~( ( fid > 0 ) & ( stat == 0 ) )
   msg = sprintf('%sError open NetCDF-File: %s',msg0,file);
   return
end

%***********************************************************
% Check for Unlimited Dimensions

for ii = 1 : size(dim,1)
    if ischar(dim{ii,2})
       dim(ii,2) = { inf };
    end
end

%***********************************************************
% Check Values

is_char = strcmp( var(:,2) , 'char' );

ok = zeros(nv,1);

ww = warnstat;
     warning('off');

for ii = 1 : nv

    val = var{ii,5};

    ok(ii) = ( isnumeric(val) | ( is_char(ii) & ischar(val) ) );

    if ok(ii) & is_char(ii) & ~ischar(val) & ~isempty(val)
       lastwarn('')
       var{ii,5} = char(double(val));
       ok(ii) = isempty(lastwarn);
    end

    if ~ok(ii)
        VarMsg(ii,[1 2]) = { 0   'Invalid Type of Variable' };
    else
        ok(ii) = ~isempty(val);
    end

end

warning(ww)

%***********************************************************
% Write Variables

ok = find(ok);

for ii = ok(:)'

    val = var{ii,5};

    if ~is_char(ii) & ~strcmp(class(val),'double')
        val = double(val);
    end

    %------------------------------------------------------------
    % Check with Attributes

    if ~isempty(att{ii}) &  ~is_char(ii)

        mv = find( strcmp( att{ii}(:,1) , 'missing_value' ) );
        fv = find( strcmp( att{ii}(:,1) , '_FillValue'    ) );
        mn = find( strcmp( att{ii}(:,1) , 'valid_min'     ) );
        mx = find( strcmp( att{ii}(:,1) , 'valid_max'     ) );
        rg = find( strcmp( att{ii}(:,1) , 'valid_range'   ) );
        ad = find( strcmp( att{ii}(:,1) , 'add_offset'    ) );
        sc = find( strcmp( att{ii}(:,1) , 'scale_factor'  ) );

        %------------------------------------------------
        % Check for Valid_Min, Valid_Max ofr ValidRange,
        %   Set to FillValue or MissingValue

        is_out = [];         % Values outside valid Range

        if ~isempty(mn)
            at = att{ii}{mn(1),3};
            if prod(size(at)) == 1
               is_out = cat( 1 , is_out , find( val < at ) );
            end
        end

        if ~isempty(mx)
            at = att{ii}{mx(1),3};
            if prod(size(at)) == 1
               is_out = cat( 1 , is_out , find( val > at ) );
            end
        end

        if ~isempty(rg)
            at = att{ii}{rg(1),3};
            if prod(size(at)) == 2
               is_out = cat( 1 , is_out , ...
                        find( ( val < min(at) ) | ( max(at) < val ) ) );
            end
        end

        %------------------------------------------------
        % Check for Scale and Offset

        if ~isempty(ad)
            at = att{ii}{ad(1),3};
            if ( prod(size(at)) == 1 )  & ~isequal(at,0)
               val = val - at;
            end
        end

        jj = strcmp( att{ii}(:,1) , 'scale_factor' );
        if ~isempty(sc)
            at = att{ii}{sc(1),3};
            if ( prod(size(at)) == 1 )  & ~isequal(at,1)
               val = val / at;
            end
        end

        %------------------------------------------------
        % Check for Miss and Fill

        is_nan = find(isnan(val));  % Missing Values

        % Check for MissingValue
        if ~isempty(mv)
            at = att{ii}{mv(1),3};
            if ( prod(size(at)) == 1 )
                mv = at
            else
                mv = []; 
            end
        end

        % Check for FillValue
        if ~isempty(fv)
            at = att{ii}{fv(1),3};
            if ( prod(size(at)) == 1 )
               fv = at; 
            else
               fv = [];
            end
        end

        if     isempty(mv), mv = fv;
        elseif isempty(fv), fv = mv;
        end

        if ~isempty(is_nan) & ~isempty(mv)
            val(is_nan) = mv;
        end

        if ~isempty(is_out) & ~isempty(fv)
            val(is_out) = fv;
        end

    end
    % Attributes

    %-----------------------------------------------------------------
    % Write Value

    si = size(val);
    
    ndval = size(si,2);
    ndvar = var{ii,3};

    is_single = ( ( ndvar == 1 )  &  ( sum(si>1) <= 1 )  );
    is_scalar = ( ( ndvar == 0 )  &  ( sum(si>1) == 0 )  );

    VarMsg{ii,1} = ( ( ndvar >= ndval ) | is_single | is_scalar  );

    if ~VarMsg{ii,1}
        VarMsg{ii,2} = 'Number of Dimensions of Variable must be agree';
    else
        if ~is_scalar
            dim_si = cat( 2 , dim{var{ii,4}+1,2} ); % !!!!!!!!!!!
            var_si = ones(1,ndvar);
            var_si(1:ndval) = size(val); 
            if is_single
               VarMsg{ii,1} = ( dim_si >= max(var_si) );
               dim_si = max(var_si);  
            else 
               nd = size(dim_si,2);
               VarMsg{ii,1} = all( dim_si(nd:-1:nd-ndvar+1) >=  var_si(1:ndvar) ); 
               dim_si(nd:-1:nd-ndvar+1) = var_si(1:ndvar);
            end
        end

        if ~VarMsg{ii,1}
            VarMsg{ii,2} = 'Length of Dimensions must be agree';
        else

 % status = ncmex('VARPUT1' , cdfid, varid , count, value )
 % status = ncmex('VARPUT' , cdfid, varid , start, count, value, autoscale)
          
            if is_scalar
               action = 'VARPUT1';
               status = ncmex('varput1',fid,ii-1,length(val),val);
            else
               action = 'VARPUT';
               status = ncmex('varput',fid,ii-1,zeros(1,ndvar),dim_si,val);
            end
if min(val(:)) == 7, keyboard, end
            VarMsg{ii,1} = ~( status == -1 );
            if ~VarMsg{ii,1}
                VarMsg{ii,2} = sprintf('Error using NCMEX(''%s'',''%s''), STATUS = %.0f ', ...
                                        action , var{ii,1} , status );
            end

        end % Dimensions ok

    end  % ndims ok

end % ii 

%***********************************************************

status = ncmex('close',fid);

if status == -1
   msg = sprintf('%sError close NetCDF-File: %s',file);;
end

%***********************************************************
% Cat VarMsg to WriteMsg

ok = cat( 1 , VarMsg{:,1} );

if all(ok)
   return
end

jj = find(~ok);

WrtMsg = sprintf( '%sError write Data.%s%s' , msg0 , nlm , ...
         strhcat( permute(cat(2,var(jj,1),VarMsg(jj,2)),[2 1]) , ': ' , 2 , nlm ) );


%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );



