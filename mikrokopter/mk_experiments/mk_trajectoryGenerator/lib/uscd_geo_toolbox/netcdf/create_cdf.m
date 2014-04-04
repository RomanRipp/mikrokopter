function [msg,ds,vs,as,file,dim,var,att] = create_cdf(file,dim,var,att)

% CREATE_CDF  Creates NetCDF-File
%
% MSG = CREATE_CDF( FILENAME , DIM , VAR , ATT )
%
%------------------------------------------------
% Inputs:
%
%  DIM : Definition for N Dimensions
%  DIM = { DimName  DimLength  }         [ N   by 2 ] CellArray
%
%  VAR : Definition for M Variables
%  VAR = { VarName VarType Ndim [dim] }  [ M   by 4 ] CellArray 
%
%  ATT : Definition for Variable- and Global-Attributes
%  ATT = { VarAtt ; GlbAtt }             [ M+1 by 1 ] CellArray
%
%  Att = { AttName AttType AttValue }    [ K   by 3 ] CellArray
%
%------------------------------------------------
%
%  DataTypes:  1. 'byte'      [        -128        127 ]      'int8'
%              2. 'char'      [           0        255 ]     'uint8'
%              3. 'short'     [      -32768      32767 ]     'int16'
%              4. 'long'      [ -2147483648 2147483647 ]     'int32'
%              5. 'float'
%              6. 'double'
%
%------------------------------------------------
% OutPuts: 
%
% [ MSG, DimStatus, VarStatus, AttStatus ] = CREATE_CDF( ... )
%
%  returns also the NCMEX-Status 
%    from the Defining of Dimensions, Variables, Attributes
%    "-1" means, that the actions wasn't successfull.
%
%    DimStatus = [ N by 1 ];
%    VarStatus = [ M by 1 ];
%    AttStatus = { [ Natt by 1 ] }; [M+1 by 1] CellArray
%
%
% The following Outputs returns the Inputs:
%
%    [ ... , FILENAME , DIM , VAR , ATT ] = CREATE_CDF( ... )
%
%
% In case of empty FILENAME the Inputs will checked only.
%
%------------------------------------------------
%
%  see also: WRITE_CDF, WRT_CDF, LOOK_CDF, READ_CDF
%


Nin = nargin;

msg = '';
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

if Nin < 3
   dim = cell(0,2);
   var = cell(0,4);
   msg = cat(2,msg0,'Not enough InputArguments');
   return
end

if Nin < 4
   att = cell(0,1);
end

msg = cell(0,1);

if isempty(file)
   file = '';
elseif ~chkstr(file)
   msg = cat(1,msg,{'FileName must be a String.'});
end

if ~iscell(dim) | ~iscell(var) | ~iscell(att)
   msg = cat(1,msg,{'Inputs DIM VAR and ATT must be CellArrays.'});
end

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    return
end

%---------------------------------------------------------------------
% Check DIM

sd = size(dim); pd = prod(sd);

if ~( ( sd(2) >= 2 ) & ( pd == sd(1)*sd(2) ) & ~( pd == 0 ) )
    msg = cat(1,msg,{'CellArray DIM must contains 2 Columns: { DimName  DimLength }'});
else
   if ~chkcstr(dim(:,1))
       msg = cat(1,msg,{'1. Column of CellArray DIM must contain Strings.'});
   end
end

%---------------------------------------------------------------------
% Check VAR

sv = size(var); pv = prod(sv);

if ~( ( sv(2) >= 4 ) & ( pv == sv(1)*sv(2) ) & ~( pv == 0 ) )

   msg = cat(1,msg,{'CellArray VAR must contains 4 Columns: { VarName VarType Ndim [dim]  }'});

else

   if ~chkcstr(var(:,1))
       msg = cat(1,msg,{'1. Column of CellArray VAR must contain Strings.'});
   end

   for ii = 1 : sv(1)

       %----------------------------------------------------
       % Check VarType

       vt = var{ii,2};
       ok = ( isnumeric(vt) & ( prod(size(vt)) == 1 ) );
       if ok
          ok = ( ( mod(vt,1) == 0 ) & ( 1 <= vt ) & ( vt <= 6 ) );
       else
          ok = chkstr(vt,1);
          if ok
             vt = lower(vt);
             if strcmp(vt,'int')
                vt = 'long';
             end
             var{ii,2} = vt;
          end
       end
       if ~ok
           str = 'single Integers between 1 and 6 or nonempty Strings';
           str = sprintf('2. Column of CellArray VAR must contain %s.',str);
           msg = cat(1,msg,{str});
       end

       %----------------------------------------------------
       % Check VarDim

       vd = var{ii,4};
       ok = isnumeric(vd);
       if ok
          ok = all( ( mod(vd,1) == 0 ) & ( 0 <= vd ) & ( vd <= sd(1)-1 ) );
          if ok
             vd = vd(:)'; 
             var{ii,4} = vd;
             var{ii,3} = size(vd,2);
          end
       else
          [ok,vd] = chkcstr(vd);
          if ok
             vd = vd(:)';
             nd = size(vd,2);
             for jj = 1 : nd
                 kk = strcmp(lower(dim(:,1)),lower(vd{jj}));
                 if ~any(kk)
                     str = sprintf('%.0f. Dimension "%s" of %.0f. Variable',jj,vd{jj},ii);
                     str = sprintf('Can''t find %s.',str);
                     msg = cat(1,msg,{str});
                     vd{jj} = NaN;
                 else
                     kk = find(kk);
                     vd{jj} = kk(1) - 1;
                 end
              end
              var{ii,4} = cat(2,vd{:});
              var{ii,3} = nd;
           end
       end
       if ~ok
           str = 'Integers in Range of DIM or Strings with DimNames';
           str = sprintf('4. Column of CellArray VAR must contain %s.',str);
           msg = cat(1,msg,{str});
       end

   end

end

%---------------------------------------------------------------------
% Check ATT

na = sv(1) + 1;

if Nin < 4

   att = cell(na,1);
   att(:) = { cell(0,3) };

else

   att = att(:);
   sa  = size(att,1);

   if     sa > na
          att = att(1:na);
   elseif sa < na
          att = cat(1,att,cell(na-sa,1));
   end

   ok = 1;
   for ii = 1 : min(sa,na)
       ok = isempty(att{ii});
       if ~ok
           s = size(att{ii}); p = prod(s);
           ok = ( ( s(2) >= 3 ) & ( p == s(1)*s(2) ) & ~( p == 0 ) );
       end
       if ~ok
           break
       end
   end

   if ~ok
       msg = cat(1,msg,{'Each Cell in ATT must contain a CellArray with 3 Columns: { AttName AttType AttValue }'});
   end

end

%-----------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    return
end

%***********************************************************

nd = sd(1);
nv = sv(1);

%***********************************************************
% Check for UNLIMITED Dimensions

ul = zeros(nd,2);   %  [ IsChar  IsUnlimited ]
for ii = 1 : nd
    ul(ii,1) = ischar(dim{ii,2});
    if ul(ii,1)
       ul(ii,2) = strcmp(lower(dim{ii,2}),'unlimited');
    else
       ok = isnumeric(dim{ii,2});
       if ok
          ok = ( prod(size(dim{ii,2})) == 1 );
          if ok
             ok = ( ( mod(dim{ii,2},1) == 0 )  &  ( dim{ii,2} >= 1 ) );
          end
       end
       ul(ii,1) = ~ok;
    end   
end

%-----------------------------------------------------------
% 1. Check for   IsChar  &  ~IsUnlimited

  BadDim = ( ( ul(:,1) == 1 )  &  ( ul(:,2) == 0 ) );

  if any(BadDim)
     jj = find( BadDim );
     msg = cat(1,msg,{'Invalid Values for Dimensions:'} , ...
                     { strhcat(dim(jj,1),', ') }        , ...
                     {'Values could be a positive, nonzero Integer or ''unlimited''.'} );
  end


%-----------------------------------------------------------
% 2. Check for   Duplicate Unlimited

  if sum(ul(:,2),1) > 1 
     jj = find(ul(:,2));
     msg = cat(1,msg,{'Only ONE Unlimited Dimension allowed, found for:'}, ...
                     { strhcat(dim(jj,1),', ')} );
  end


if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    return
end


%-----------------------------------------------------------
% 3. Check for Unlimited Dimension in VariableDimension

  if any(ul(:,2));

     dul = find(ul(:,2)) - 1;  % ID of Unlimited Dimension
 
     vok = ones(nv,1);
     for ii = 1 : nv
         if ~isempty(var{ii,4})
             if any( var{ii,4} == dul ) 
                vok(ii) = ( var{ii,4}(1) == dul );
             end
         end
     end

     if any(~vok)
        jj = find(~vok);
        msg = sprintf('%sUnlimited Dimension must be the 1. in Variables:%s%s', ...
                       msg0 , nlm , strhcat(var(jj,1),', ',3,nlm) );
        return
     end

  end
  % dul


%***********************************************************

if isempty(file)
   msg = '';
   return
end
    
%***********************************************************
% CREATE File

fid = ncmex('create',file,'clobber');

if fid <= 0
   msg = sprintf('%sError create NetCDF-File "%s".',msg0,file);
   return
end


ds = NaN*ones(nd,1);  % DimensionStatus
vs = NaN*ones(nv,1);  % VariableStatus

as    = cell(na,1);  % AttributeStatus
as(:) = { NaN*zeros(0,1) };


%***********************************************************
% Store Global Attributes

if ~isempty(att{na})

    s = size(att{na},1);

    as{na} = NaN*ones(s,1);

    for jj = 1 : s

        name =  att{na}{jj,1};
        typ  =  att{na}{jj,2};
        val  =  att{na}{jj,3};

        if strcmp(typ,'int')
           typ = 'long';
        end

        if ~strcmp(typ,'char') & ~strcmp(class(val),'double')
           val = double(val);
        end

        if ~isempty(name)
            as{na}(jj) = ncmex('attput',fid,'global',name,typ,length(val),val);
        end

    end 

end

%***********************************************************
% Define Dimensions
%
% status = ncmex('DIMDEF',fid, 'name'  , length)
%

for ii = 1 : nd
    ds(ii) = ncmex('dimdef',fid,dim{ii,1},dim{ii,2});
end


%***********************************************************
% Define Variables
%
% status = ncmex('VARDEF',fid,'name'   , datatype, ndims   , [dim])
%

for ii = 1 : nv

  if strcmp(var{ii,2},'int')
     var{ii,2} = 'long';
  end

  vs(ii) = ncmex('vardef',fid,var{ii,1},var{ii,2},var{ii,3},var{ii,4});
 

  if ~isempty(att{ii})

      s = size(att{ii},1);

      as{ii} = NaN*ones(s,1);

      for jj = 1 : s

          name =  att{ii}{jj,1};
          typ  =  att{ii}{jj,2};
          val  =  att{ii}{jj,3};

          if strcmp(typ,'int')
             typ = 'long';
          end

          if ~strcmp(typ,'char') & ~strcmp(class(val),'double')
              val = double(val);
          end

          if ~isempty(name)

              if      strcmp( typ , 'char' )  &  isequal(val,'\0')  &  ...
                  any(strcmp( name , {'missing_value' '_FillValue'} ))
                  val = char(0);
              end

              as{ii}(jj) = ncmex('attput',fid,(ii-1),name,typ,length(val),val);

          end

      end 

  end

end

%***********************************************************
% Thats all

 msg = '';

 status = ncmex('close',fid);

 if status == -1
    msg = sprintf('%sError close NetCDF-File "%s".',msg0,file);
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

