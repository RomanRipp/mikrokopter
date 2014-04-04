function [msg,dim,var,att,txt] = read_hdf(varargin)

% READ_HDF   reads Data from a HDF-File, using LOOK_HDF
%
% [ Msg, DIM, VAR, GlobalATT, TXT ] = READ_HDF(FILE)
%
%  VAR : Information for Variables,  [ Nvar by 15 ] CellArray
%
%  VAR =   1. Name
%          2. Rank
%          3. Type
%          4. Ndim
%          5. Size
%          6. Natt 
%          7. Range       : valid_range = [ Min Max ]
%          8. Fill        : _FillValue
%          9. DIM         : { Name  Length  Type  Scale  ATT    }
%         10. ATT         : { Name  Count   Type  Value  String }
%         11. Compressed  : 0 | 1
%         12. Chunked     : 0 | 1
%         13. ChunkLenght
%         14. Value
%         15. ReadStatus
%
%  DIM : Information for Dimensions, [ Ndim by  6 ] CellArray
%
%  DIM = { Name  Length  Type  Scale  ATT  [Start Count Stride] }
%
%  ATT : Information for Attributes, [ Natt by  5 ] CellArray
%
%  ATT = { Name  Count   Type  Value  String }
%
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
% [ ... ] = READ_HDF( HDF_FILENAME , Property1 , Property1Value , ... )
%
% Valid Properties:
%
%  'var' ,  VarNames   %  Read's only the Variables and their Attributes,
%                      %   specified by the CellStringArray VarNames.
%                      %  If the VarName wasn't found in the HDF-File,
%                      %   the Typ of this Variable in the 2. Column of VAR
%                      %   is set to 'none'.
%
%  'dim' , { DimName DimVal } % Select the Variables at the specified
%                             %  Dimensions only at the single DimNumber
%                             %   DimName   ... String, 
%                             %   DimVal    ... (complex) Integer
%                             % DimVal = DimNumber +      0*i
%                             % DimVal = DimStart  + Stride*i
%                             % DimVal = DimEnd    - Stride*i
%                             % Note: the first Value is equal to 1
%                             % ( real(DimVal) <= 0 )  ==> measured from End
%      
%  'miss' , MissingValue  % MissingValues will set to this Value, default: NaN  
%  'fill' ,    FillValue  %    FillValues will set to this Value, default: EMPTY  
%
%       Empty Value, no replacing of this Values
%
%----------------------------------------------------------------
%
%  see also: LOOK_HDF, ASSIGN_HDF, HDFSD
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
% Read HDF-File

if Nout < 5
   [msg,dim,var,att] = look_hdf(file);
else
   [msg,dim,var,att,txt] = look_hdf(file);
end

% var = { 1.Name 2.Rank 3.Type 4.Ndim 5.[Size] 6.Natt 7.Fill 8.[Range] 9.{dim} 10.{att} ...
%         11.Compressed 12. Chunked  13.ChunkLenght }
% dim = { Name  Length  Type  Scale  ATT   }
% att = { Name  Count   Type  Value  String }

if ~isempty(msg)
    msg = sprintf('%sError using LOOK_HDF(''%s'').%s%s',msg0,file,nlm,msg);
    return
end

fid = hdfsd('start',file,'read');

if fid == -1
   name = which(file);
  fid = hdfsd('start',name,'read');
   if fid == -1
      msg = sprintf('%sError open HDF-File "%s".',msg0,file);
      return
   end
end

%***********************************************************
% Check for VAR

vnr = ( 1 : size(var,1) ) - 1;

if is_var

   nv = size(vin,1);

   vn = zeros(1,nv);

   v = cell(nv,size(var,2));
   v(:,1) = vin;
   v(:,3) = {''};

   v(:,[9 10]) = {cell(0,5)};

   for ii = 1 : nv
       jj = strcmp(var(:,1),vin{ii});
       if any(jj)
          jj = find(jj);
          v(ii,:) = var(jj,:);
          vn(ii)  = vnr(jj);
       end
   end

   var = v;
   vnr = vn;

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

  len = dim{ii,2};      % Counts

  dim{ii,id}(2) = len;  % Counts

  if is_dim

     jj = strcmp(din(:,1),dim{ii,1});
     if any(jj)

        jj = find(jj);
        jj = jj(end);

        dd = din{jj,2};
        d0 = real(dd);    % Start
        d1 = imag(dd);    % Stride

        d0 = d0 + len * ( d0 <= 0 );   % Check for Measuring from End

        ok(jj) = ( ( 1 <= d0 ) & ( d0 <= len ) );

        if ok(jj)

            if     ( d1 == 0 )
               dim{ii,id}(1) = d0 - 1;                                % Start
               dim{ii,4}(2) = 1;                                      % Count
            elseif ( d1 > 0 )
               dim{ii,id}(1) = d0 - 1;                                % Start
               dim{ii,id}(2) = floor( ( len - d0 ) / d1 ) + 1;        % Count
               dim{ii,id}(3) = d1;                                    % Stride
            elseif ( d1 < 0 )
                          d1 = abs(d1);
               dim{ii,id}(3) = d1;                                    % Stride
               dim{ii,id}(2) = floor( ( d0 - 1 ) / d1 ) + 1;          % Count
               dim{ii,id}(1) = d0 - d1 * ( dim{ii,id}(2) - 1 ) - 1;   % Start
            end

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

msg    = cell(size(var,1),1);
msg(:) = {''};

ind = find(~strcmp(var(:,3),''));

for ii = ind(:)'

% var = { 1.Name 2.Rank 3.Type 4.Ndim 5.[Size] 6.Natt 7.Fill 8.[Range] 9.{dim} 10.{att} ...
%         11.Compressed 12. Chunked  13.ChunkLenght }

    is_char = type_cmp(var{ii,3},'char');

    %--------------------------------------------------------------
    % Get Data

    vid = hdfsd('select',fid,vnr(ii));

    nd  = var{ii,4};
    start  = zeros(1,nd);
    count  = zeros(1,nd);
    stride = zeros(1,nd);

    for jj = 1 : nd
        kk = find(strcmp(dim(:,1),var{ii,9}{jj,1}));
        start(jj)  = scs(1,kk);
        count(jj)  = scs(2,kk);
        stride(jj) = scs(3,kk);
    end

      if all(count)
        [ val , stat ] = hdfsd( 'readdata' , vid , start , stride , count );
      else
        stat = 0;
        val  = zeros(count(var{ii,3}:-1:1));
      end

    %---------------------------------------------------

    status = hdfsd('endaccess',vid);

    %--------------------------------------------------------------
    % Check Attributes
    % att = { Name  Count   Type  Value  String }

    att = var{ii,10};

    if ~isempty(att)  &  ~isempty(val)  &  ~( stat == -1 )

        %--------------------------------------------------------------
        % Search for Miss and Fill

        repl(:,4) = { [] };   % Reset Index to Replace

        for aa = 1 : nr
            jj = strcmp( att(:,1) , repl{aa,1} );
            if any(jj) & ( ~isempty(repl{aa,2})  |  is_char )
               jj = find(jj);
               at = att{jj,4};
               if  type_cmp(att(jj,3),'char')  &  isequal(at,'\0') 
                   at = char(0);
               end
               if prod(size(at)) == 1
                  repl{aa,4} = find( val == at );
               end
            end
        end

        %--------------------------------------------------------------
        % Scale and Add Offset

        jj = strcmp( att(:,1) , 'scale_factor' );
        if any(jj)

           jj = find(jj);
           at = att{jj,4};

           if ~strcmp(class(at),class(val))
              if type_cmp(class(at),'char')
                 cl = 'char';
              elseif type_cmp(class(at),'int','byte','short','long')
                 switch att{jj,3}
                    case 'byte'
                        cl = 'int8';
                    case 'short'
                        cl = 'int16';
                    case 'long'
                        cl = 'int32';
                    otherwise
                        cl = att{jj,3};
                 end
              else
                 cl = 'double';
              end
              try
                 val   = feval(cl,val);
              catch
                 msg{ii} = sprintf('Error convert "%s" from %s to %s.', ...
                                 var{ii,1},upper(class(val)),upper(cl));
              end
           end

           if ( prod(size(at)) == 1 ) & ~isequal(at,1) & isempty(msg{ii})
              val = val * at;
              if is_char
                 val = char(val);
              end
           end

        end

        jj = strcmp( att(:,1) , 'add_offset' );
        if any(jj) & isempty(msg{ii})

           jj = find(jj);
           at = att{jj,4};

           if ~strcmp(class(at),class(val))
              if type_cmp(class(at),'char')
                 cl = 'char';
              elseif type_cmp(class(at),'int','byte','short','long')
                 switch att{jj,3}
                    case 'byte'
                        cl = 'int8';
                    case 'short'
                        cl = 'int16';
                    case 'long'
                        cl = 'int32';
                    otherwise
                        cl = att{jj,3};
                 end
              else
                 cl = 'double';
              end
              try
                 val = feval(cl,val);
              catch
                 msg{ii} = sprintf('Error convert "%s" from %s to %s.', ...
                                 var{ii,1},upper(class(val)),upper(cl));
              end
           end

           if ( prod(size(at)) == 1 ) & ~isequal(at,0) & isempty(msg{ii})
              val = val + at;
              if is_char
                 val = char(val);
              end
           end

        end

        %--------------------------------------------------------------
        % Replace Miss and Fill

        rr = 2 + is_char;
        for aa = 1 : nr
            if ~isempty(repl{aa,4})
                val(repl{aa,4}) = repl{aa,rr};
            end                
        end

    end

    %--------------------------------------------------------------
    % Check for ZERO-Characters

    if ~isempty(val) & ~( stat == -1 ) & is_char
            jj  = find( double(val) == 0 );
        val(jj) = char(32);
    end

    %--------------------------------------------------------------
    if ~( stat == -1 )
       var{ii,iv} = val; 
    end

    var{ii,iv+1} = stat;
   
end

%***********************************************************

msg = msg(find(~strcmp(msg,'')));

status = hdfsd('end',fid);

if status == -1
   msg = cat(1,msg,{sprintf('Error close HDF-File "%s".',file)});
end

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
else
    msg = '';
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = type_cmp(typ,varargin);

ok = 0;
if nargin < 2
   return
end

[vok,cmp] = chkcstr(varargin(:));

if ~vok
    return
end

for ii = 1 : size(cmp,1)
    if size(typ,2) >= size(cmp{ii},2)
       ok = ~isempty(findstr(typ,cmp{ii}));
       if ok
          break
       end
    end
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
          ok = ( ok & isnumeric(dd)  &  ( prod(size(dd)) == 1 ) );
          if ok
             dd = [ real(dd)  imag(dd) ];
             ok = ( ok & all( isfinite(dd)  & ( mod(dd,1) == 0 ) ) );
          end
          if ~ok
              break
          end
      end
      if ~ok
          msg = cat(1,msg,{'DimNumbers must be a single (complex) Integer.'});
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


