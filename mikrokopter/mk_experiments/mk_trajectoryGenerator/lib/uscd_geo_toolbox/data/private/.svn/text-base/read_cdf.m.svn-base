function [msg,dim,var,att,txt] = read_cdf(varargin)

% READ_CDF   reads Data from a NetCDF-File, using LOOK_CDF
%
% [ Msg, DIM, VAR, ATT, TXT ] = READ_CDF(FILE)
%
%  DIM is an CellArray with Informations about the Dimensions:
%  DIM = { DimName  DimLength  UNLIMITED  [ Start Count Stride ]  }
%
%  VAR is an 8-Column CellArray with 
%       Informations and Data of the Variables:
%  VAR = { VarName VarType Ndim [dim] Nattr  []   Data  VarGetStatus}
%   if Ndim == 1, Var{:,6} = [First Last]Value, if  VarType~=2
%                 Var{:,6} =   String         , if  VarType==2
%
%  ATT is an CellArrays with Information about the Attributes:
%  ATT contains one CellArray for each Variable 
%        and one more for Global Attributes 
%  ATT = { { AttName AttType AttValue AttString } }
%       
%  TXT is an CellArray of Strings,
%       containing the Informations as ASCII-Text.
%
%  MSG contains the ErrorMessages.
%
%----------------------------------------------------------------
%
% additional Options:
%
% [ ... ] = READ_CDF( CDF_FILENAME , Property1 , Property1Value , ... )
%
% Valid Properties:
%
%  'var' ,  VarNames   %  Read's only the Variables and their Attributes,
%                      %   specified by the CellStringArray VarNames.
%                      %  If the VarName wasn't found in the CDF-File
%                      %   the Value for the Type of this Variable in the 
%                      %   2. Column of VAR is set to 'none'.
%
%  'dim' , { DimName DimVal } % Select the Variables at the specified
%                             %  Dimensions only at the specified Values
%                             %   DimName   ... String, 
%                             %   DimVal    ... (complex) IntegerVector
%
%      DimVal = [ i*Start ]                % Extract at single Value
%      DimVal = [ Stride + i*Start      ]
%      DimVal = [ Start  Stride         ]
%      DimVal = [ Start  Count   Stride ]
%
%   Note: The first Value refers to Start == 1
%         A Negative or ZERO Start will measured from End
%         A Negative Value of Stride will extract to Begin
%
%
%  'miss' , MissingValue  % MissingValues will set to this Value, default: NaN  
%  'fill' ,    FillValue  %    FillValues will set to this Value, default: EMPTY  
%
%       Empty Values ==> no replacing of this Values
%
%----------------------------------------------------------------
%
%  see also: ASSIGN_CDF, LOOK_CDF, WRITE_CDF, CREATE_CDF 
%


Nout = nargout;

msg = '';
dim = [];
var = [];
att = [];
txt = [];

nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));




nlm = sprintf('%s%s',nl,char(32*ones(size(msg0))));

%***********************************************************

[msg,file,vin,din,mval,fval] = checkin(varargin{:});

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    return
end

is_var = ~isempty(vin);
is_dim = ~isempty(din);

%***********************************************************
% Read NetCDF-File

if Nout < 5
   [msg,dim,var,att] = look_cdf(file);
else
   [msg,dim,var,att,txt] = look_cdf(file);
end

if ~isempty(msg)
    msg = sprintf('%sError using LOOK_CDF(''%s'').%s%s',msg0,file,nlm,msg);
    return
end

[fid,stat] = ncmex('open',file,'nowrite');

if ~( ( fid > 0 ) & ( stat == 0 ) )
   name = which(file);
   [fid,stat]  = ncmex('open',name,'nowrite');
   if ~( ( fid > 0 ) & ( stat == 0 ) )
      msg = sprintf('%sError open NetCDF-File "%s".',msg0,file);
      return
   end
end

%***********************************************************
% Check for VAR

if is_var

   nv = size(vin,1);

   v = cell(nv,size(var,2));
   v(:,1) = vin;
   v(:,2) = {'none'};

   a = cell(nv,1);
   a(:) = { cell(0,3) };

   for ii = 1 : nv
       jj = strcmp(var(:,1),vin{ii});
       if any(jj)
          jj = find(jj);
          v(ii,:) = var(jj,:);
          a(ii)   = att(jj);
       end
   end

   var = v;
   att = a;

end

nv2 = size(var,2);
var = var(:,cat(2,(1:nv2),nv2,nv2));
var(:,nv2+1) = {  [] };
var(:,nv2+2) = { NaN };

iv = nv2+1;

%***********************************************************
% Check for DIM

nd2 = size(dim,2);
dim = dim(:,cat(2,(1:nd2),nd2));
dim(:,nd2+1) = { [ 0  NaN  1 ] };  %  [ Start Count Stride ] for each Dimension

id = nd2+1;

ok = ones(size(din,1),1);

for ii = 1 : size(dim,1)

  len = dim{ii,2};     % Counts

  dim{ii,id}(2) = len;  % Counts

  if is_dim

     jj = strcmp(din(:,1),dim{ii,1});

     if any(jj)

        jj = find(jj);
        jj = jj(end);

        dd = din{jj,2};
        nd = prod(size(dd));

        scs = NaN * zeros(1,3);

        %----------------------------------------
        if nd == 1
        %  Stride | Start*i
           if dd == 0
              scs = [ 0  1  1 ];
           else
              scs(3) = real(dd);   % Stride
              scs(1) = imag(dd);   % Start
              if scs(3) == 0
                 scs(2) = 1;       % Single Count
              end
           end
        %----------------------------------------
        elseif nd == 2
        % Start Stride
           scs([1 3]) = dd;
        %----------------------------------------
        else
        % Start Count Stride
           scs = dd;
        end
        %----------------------------------------

        scs(1) = scs(1) + 1;  % !!! Start with 1

        scs(3) = scs(3) + ( scs(3) == 0 );
        scs(1) = scs(1) - len * floor( scs(1) / len );
        scs(1) = scs(1) + len * ( scs(1) == 0 );

        if isnan(scs(2))
           if scs(3) > 0
              scs(2) = floor( ( len - scs(1) ) / scs(3) ) + 1;
           else
              scs(3) = abs(scs(3));
              scs(2) = floor( ( scs(1) - 1 ) / scs(3) ) + 1;
              scs(1) = scs(1) - scs(3) * ( scs(2) - 1 );
           end
        else
           if scs(3) < 0
              scs(3) = abs(scs(3));
              scs(1) = scs(1) - scs(3) * ( scs(2) - 1 );
           end
        end

        ok(jj) = ( ( 1 <= scs(1) ) &  ...
                   ( scs(1)+scs(3)*( scs(2) - 1 ) <= len ) );

        if ok(jj)
           dim{ii,id} = scs + [ -1  0  0 ]; %  !!! Start with 0
        end

     end

  end

end

if is_dim
   if any(~ok)
      jj = find(~ok);
      msg = sprintf( '%sOut of Range Values for Dimensions to read:%s%s' , ...
                     msg0 , nlm , strhcat(din(jj,1),', ',4,nlm) );
      return
   end
end

scs = cat(1,dim{:,id})';   % [ Start ; Count ; Stride ]

%***********************************************************

% Values to replace

%         AttName         Value    CharValue  Index
repl = { 'missing_value'  mval       ' '        []  
         '_FillValue'     fval       ' '        []  };

nr = size(repl,1);

is_char = strcmp(var(:,2),'char');
is_byte = strcmp(var(:,2),'byte');
is_int  = 0 * is_char;               % True for UINT8


ind = find(~strcmp(var(:,2),'none'));

for ii = ind(:)'

    %--------------------------------------------------------------
    % Get Data

    if ~var{ii,3}
    % Variable without any Dimension

      [ val , stat ] = ncmex( 'varget1' , fid , var{ii,1} , 1  );

    else
    
      did    = var{ii,4}+1;
      start  = scs(1,did);
      count  = scs(2,did);
      stride = scs(3,did);
  
      if all(count)
        [ val , stat ] = ncmex( 'vargetg' , fid , var{ii,1} , start , count , stride );
      else
        stat = 0;
        val  = zeros(count(var{ii,3}:-1:1));
      end

    end

    %--------------------------------------------------------------
    % Check Attributes

    if ~isempty(att{ii})  &  ~isempty(val)  &  ~( stat == -1 )

        %--------------------------------------------------------------
        % Search for Miss and Fill

        repl(:,4) = { [] };   % Reset Index to Replace

        for aa = 1 : nr
            jj = strcmp( att{ii}(:,1) , repl{aa,1} );
            if any(jj) & ( ~isempty(repl{aa,2})  |  is_char(ii) )
               jj = find(jj);
               at = att{ii}{jj,3};
               if  strcmp(att{ii}(jj,2),'char')  &  isequal(at,'\0') 
                   at = char(0);
               end
               if prod(size(at)) == 1
                  repl{aa,4} = find( val == at );
               end
            end
        end

        %--------------------------------------------------------------
        % Scale and Add Offset, Transform CHAR to BYTE(UINT8)

        scl = 1;
        jj = strcmp( att{ii}(:,1) , 'scale_factor' );
        if any(jj)
           is_int(ii)  = ( is_int(ii) | is_char(ii) );
           is_char(ii) = 0;
           jj = find(jj);
           scl = att{ii}{jj,3};
           if ( prod(size(scl)) == 1 ) & ~isequal(scl,1)
              val = val * scl;
              if is_char(ii)
                 val = char(val);
              end
           end
        end

        off = 0; 
        jj = strcmp( att{ii}(:,1) , 'add_offset' );
        if any(jj)
           is_int(ii) = ( is_int(ii) | is_char(ii) );
           is_char(ii) = 0;
           jj = find(jj);
           off = att{ii}{jj,3};
           if ( prod(size(off)) == 1 ) & ~isequal(off,0)
              val = val + off;
              if is_char(ii)
                 val = char(val);
              end
           end
        end

        is_int(ii) = ( is_int(ii) | ...
                     ( is_byte(ii) & isequal(scl,1) & isequal(off,128) ) );

        %--------------------------------------------------------------
        % Replace Miss and Fill

        rr = 2 + is_char(ii);
        for aa = 1 : nr
            if ~isempty(repl{aa,4})
                val(repl{aa,4}) = repl{aa,rr};
            end                
        end

    end

    %--------------------------------------------------------------
    % Check for ZERO-Characters

    if ~isempty(val) & ~( stat == -1 ) & is_char(ii)
            jj  = find( double(val) == 0 );
        val(jj) = char(32);
    end

    %--------------------------------------------------------------
    if ~( stat == -1 )
       if is_int(ii)
          val = uint8(val);
       end
       var{ii,iv} = val; 
    end

    var{ii,iv+1} = stat;
   
end

%***********************************************************

msg = '';

status = ncmex('close',fid);

if status == -1
   msg = sprintf('%sError close NetCDF-File "%s".',msg0,file);
end

  
%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,var,dim,mv,fv] = checkin(varargin)

msg = cell(0,1);

file = '';

var = cell(0,1);
dim = cell(0,2);

mv = NaN;   % MissingValue
fv =  [];   %    FillValue

Nin = nargin;

if Nin == 0
   msg = {'Input FileName is missing.'};
   return
end

%---------------------------------------------------------
% Check for File

val = varargin{1};
if ~chkstr(val,1)
    msg = cat(1,msg,{'FileName must be a String.'});
elseif ~( exist(val,'file') == 2 )
    msg = cat(1,msg,{sprintf('File "%s" not found.',val)});
else
    file = val;
end
 
if Nin == 1
   return
end

%-----------------------------------------------------------
% Basic

nv = Nin - 1;

if ~( mod(nv,2) == 0 ) | ~chkcstr(varargin(2:2:nv))
   msg = cat(1,msg,{'Additional Inputs must Property-Value-Pairs.'}, ...
                   {' Properties must be Strings.'});
   msg = strhcat(msg,'',1);
   return
end

prop = varargin(2:2:nv);

%-----------------------------------------------------------
% Check for Variables

jj = strcmp(prop,'var');
if any(jj)
   jj = 2 * find(jj) + 1;
   [ok,val] = chkcstr(varargin{jj(end)});
   if ~ok
       msg = cat(1,msg,{'VariableNames must be a CharArray or CellStringArray.'});
   else
       var = val(:);
   end
end

%-----------------------------------------------------------
% Check for Dimensions

jj = strcmp(prop,'dim');
if any(jj)
   jj  = 2 * find(jj) + 1;
   val = varargin{jj(end)};
   s = size(val); p = prod(s);
   ok = ( iscell(val) & ( p == s(1)*s(2) ) & ~( p == 0 ) & ( s(2) == 2 ) );
   if ok
      ok = chkcstr(val(:,1));
   end
   if ~ok
       msg = cat(1,msg,{'Dimension must be a 2-Column-CellArray: { DimName  DimNumber }.'});
   else
      for ii = 1 : s(1)
          dd  = val{ii,2};    % [ Start + Stride*i ];
          ok = ( ok & isnumeric(dd) & ~isempty(dd) );
          if ok
             dd = dd(:)';
             nd = min( 3 , size(dd,2) );
             dd = dd(1:nd);
             ok = all(   isfinite(dd)           & ...
                       ( mod(real(dd),1) == 0 ) & ...
                       ( mod(imag(dd),1) == 0 )        );
             ok = ( ok & ( ~( dd(nd) == 0 ) | ( nd == 1 ) ) );
          end
          if ~ok
              break
          end
      end
      if ~ok
          msg = cat(1,msg,{'DimNumbers must be a single (complex) Integer with NonZero last Value.'});
      end
   end
   if ok
      dim = val;
   end
end

%-----------------------------------------------------------
% Check for Miss and Fill

jj = strcmp(prop,'miss');
if any(jj)
   jj  = 2 * find(jj) + 1;
   val = varargin{jj(end)};
   ok = ( isnumeric(val) & ( prod(size(val)) == 1 ) );
   if ~ok
       msg = cat(1,msg,{'Missing_Value must be a single numeric.'});
   else
       mv = val;
   end
end

jj = strcmp(prop,'fill');
if any(jj)
   jj  = 2 * find(jj) + 1;
   val = varargin{jj(end)};
   ok = ( isnumeric(val) & ( prod(size(val)) == 1 ) );
   if ~ok
       msg = cat(1,msg,{'Fill_Value must be a single numeric.'});
   else
       fv = val;
   end
end


%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


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


