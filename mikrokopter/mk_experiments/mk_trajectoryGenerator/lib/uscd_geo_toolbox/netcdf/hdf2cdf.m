function [msg,dim,var,att,cdf,ds,as,vs,ws] = hdf2cdf(varargin);

% HDF2CDF  Converts HDF-Files into NetCDF-Files
%
% Msg = HDF2CDF( HDF_File , CDF_File , BlockSize , '@<LOG_File>' )
%
%  HDF_File    Name of HDF-File to convert
%
%  CDF_File    Name of CDF-File to create
%               optional, default: <HDF-FileName>.cdf
%
%  BlockSize   defines the maximum Size of DataCounts
%               to read and write at same Time, 
%               usefull to reduce the consumption of Memory
%               optional, default: 1,000,000
%
%  LOG_File    specifies a File for LogOutput, 
%               empty (just an '@') ==> no LogOutput
%               optional, default: LogOutput to MatlabCommandWindow
%
%  Msg contains ErrorMessages, empty if all was successfull
%
%-------------------------------------------------------------------
%
% [ Msg , DIM , VAR , ATT ] = HDF2CDF( ... )
%
%  returns the Definitions for the created NetCDF-File.
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
%
% NOTE: The Values of DimName, VarName, AttName will checked 
%        for valid Characters: 'a' .. 'z'
%                              'A' .. 'Z'
%                              '0' .. '9'
%                              '_'
%
%       Invalid Characters will replaced by '_' !!!
%
%-------------------------------------------------------------------
%
% [ ... , CDF_File ] = HDF2CDF( ... )
% 
% returns the Name of the created NetCDF-File, usefull in case of
%  CDF_File is not given in Inputs or empty.
%
%-------------------------------------------------------------------
%
% [ ... , DimStat , VarStat , AttStat , WrtStat ] = HDF2CDF( ... )
%
%  returns also the NCMEX-Status 
%    from the Defining of Dimensions, Variables, Attributes
%    "-1" means, that the actions wasn't successfull.
%
%    DimStat = [ N by 1 ];
%    VarStat = [ M by 1 ];
%    AttStat = { [ Natt by 1 ] }; [M+1 by 1] CellArray
%
%    WrtStat = [ M by 2 ];  [ Read_HDF_Status  Write_CDF_Status ]
%
%-------------------------------------------------------------------
%
% see also: LOOK_HDF, READ_HDF, ASSIGN_HDF, HDFSD, 
%           LOOK_CDF, READ_CDF, ASSIGN_CDF, NCMEX, READ_CDF, WRITE_CDF 
%
%


msg = '';
dim = [];
var = [];
att = [];
ds  = [];
vs  = [];
as  = [];
ws  = [];
 
no_out = ( nargout == 0 );


msg = '';

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));

nl   = char(10);

nlm = sprintf('%s%s',nl,char(32*ones(size(msg0))));

%************************************************************

[msg,hdf,cdf,bsi,lid] = checkin(varargin{:});

if ~isempty(msg)
    msg = sprintf('%sInvalid Inputs.%s%s',msg0,nlm,strhcat(msg,nlm));
    if no_out
       error(msg);
    end
end

if isempty(hdf)
   return
end

%***********************************************************

log_file = ( lid > 1 );

if log_file
   setmark(lid,'start',msg0);
else
   fprintf(lid,'\n%s',msg0);
end

fprintf(lid,'Convert  %s  --->  %s\n',hdf,cdf);
 
%***********************************************************
% Read HDF-File

fprintf(lid,'\nCall LOOK_HDF(''%s'') ... ',hdf);

[msg,dim,var,att] = look_hdf(hdf);

% var = { 1.Name 2.Rank 3.Type 4.Ndim 5.[Size] 6.Natt 7.Fill 8.[Range] 9.{dim} 10.{att} ...
%         11.Compressed 12. Chunked  13.ChunkLenght }
% dim = { Name  Length  Type  Scale  ATT   }
% att = { Name  Count   Type  Value  String }

if ~isempty(msg)
    fprintf(lid,'error\n\n');
    msg = strrep(msg,nl,nlm);
    msg = sprintf('%sError using LOOK_HDF(''%s'').%s%s',msg0,hdf,nlm,msg);
    if no_out | log_file
       fprintf(lid,'%s\n\n',msg);
       if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
       if no_out, clear msg, end
    end
    return
end

fprintf(lid,'ok\n');

dim = dim(:,[1 2]);      % { Name Length }

att = cat(1,var(:,10),{att});

var = var(:,[1 3 4 9]);  % { Name Type Ndim [dim] }


nvar = size(var,1);
natt = size(att,1);

%***********************************************************
% Check for Names

fprintf(lid,'\nCheck Parameter ... ');

msg = cell(0,1);

for ii = 1 : size(dim,1)

    dim{ii,1} = check_name(dim{ii,1});

end

hdf_org = var(:,[1 2]);  % Store original Name and Type

for ii = 1 : size(var,1)

    var{ii,1} = check_name(var{ii,1});

    var{ii,3} = size(var{ii,4},1);

    for jj = 1 : var{ii,3}
        var{ii,4}{jj,1} = check_name(var{ii,4}{jj,1});
    end

end

for ii = 1 : natt

    for jj = 1 : size(att{ii},1)
        att{ii}{jj,1} = check_name(att{ii}{jj,1});
    end

    d = char(att{ii}(:,1));
    d = ( diff(d,1,1) == 0 );
    d = ( sum(d,2) == size(d,2) );
        
   if any(d)
      d = find(d);
      if size(d,1) > 1
         d(find(diff(d,1,1)==1)+1) = [];
      end
      if ii < natt
         str = sprintf('of Variable "%s"',var{ii,1});
      else
         str = 'global Attributes';
      end
      str = sprintf( 'Duplicate AttributeNames of %s:%s%s' , ...
                       str , nlm , strhcat(var{ii,10}(d,1),', ',4,nlm) );
      msg = cat(1,msg,{str});
   end

end

for ii = { dim  var }

   if size(ii{1},1) > 1

      d = char(ii{1}(:,1));
      d = ( diff(d,1,1) == 0 );
      d = ( sum(d,2) == size(d,2) );

      if any(d)
         d = find(d);
         if size(d,1) > 1
            d(find(diff(d,1,1)==1)+1) = [];
         end
         if size(ii{1},2) == 2
            str = 'Dimension';
         else
            str = 'Variable';
         end
         str = sprintf('Duplicate %sNames:%s%s',str,nlm,strhcat(ii{1}(d,1),nlm));
         msg = cat(1,msg,{str});
      end

   end

end

%***********************************************************
% VarDim and Types

for ii = 1 : size(var,1)

    typ = var{ii,2};

    var{ii,2} = check_type(typ);

    if strcmp(var{ii,2},'none')
       msg = cat( 1 , msg , ...
            {sprintf('Can''t translate Type "%s" for Variable "%s".',typ,var{ii,1})} );
    end

    dn = var{ii,4}(:,1);

    var{ii,4} = zeros(1,var{ii,3});

    for jj = 1 : var{ii,3}

        kk = find(strcmp(dim(:,1),dn{jj}));

        if isempty(kk)
            str = sprintf( 'Can''t find Dimension "%s" for Variable "%s".', ...
                           dn{jj} , var{ii,1} );
            msg = cat( 1 , msg , {str} );
        else
           var{ii,4}(jj) = kk - 1;
        end

    end

end

for ii = 1 : natt

    att{ii} = att{ii}(:,[1 3 4]);  % { Name Type Value }

    for jj = 1 : size(att{ii},1)

        typ = att{ii}{jj,2};

        att{ii}{jj,2} = check_type(typ);

        if strcmp(att{ii}{jj,2},'none')

           if ii < natt
              str = sprintf('of Variable "%s"',var{ii,1});
           else
              str = 'global Attributes';
           end
           str = sprintf('Can''t translate Type "%s" for Attribute "%s" of %s.', ...
                         typ , att{ii}{jj,1} , str );
           msg = cat( 1 , msg , {str} );

        elseif ~type_cmp(att{ii}{jj,2},'char') & ...
               ~strcmp(class(att{ii}{jj,3}),'double')

            att{ii}{jj,3} = double(att{ii}{jj,3});

        end

    end

end

if ~isempty(msg)
    fprintf(lid,'error\n\n');
    msg = strhcat(msg,nlm);
    msg = sprintf('%sError convert Parameter.%s%s',msg0,nlm,msg);
    if no_out | log_file
       fprintf(lid,'%s\n\n',msg);
       if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
       if no_out, clear msg, end
    end
    return
end

fprintf(lid,'ok\n');

%***********************************************************
% Create NetCDF-File

%-------------------------------------------------------------
% Quiet Mode

ncmex('setopts',0);

%-------------------------------------------------------------

fprintf(lid,'\nCall CREATE_CDF(''%s'') ... ',cdf);

try
   [msg,ds,vs,as] = create_cdf(cdf,dim,var,att);
catch
   msg = lasterr;
end

if ~isempty(msg)
    fprintf(lid,'error\n\n');
    msg = strrep(msg,nl,nlm);
    msg = sprintf('%sError call CREATE_CDF.%s%s',msg0,nlm,msg);
    if no_out | log_file
       fprintf(lid,'%s\n\n',msg);
       if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
       if no_out, clear msg, end
    end
    return
end

fprintf(lid,'ok\n\n');

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
    if no_out | log_file
       fprintf(lid,'%s\n\n',msg);
       if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
       if no_out, clear msg, end
    end
    return
end

%***********************************************************
% Open Files to read/write

%-------------------------------------------------------------

fprintf(lid,'Open    HDF-File to  read: "%s" ... ',hdf);

fid = hdfsd('start',hdf,'read');

if fid == -1
  name = which(hdf);
  fid = hdfsd('start',name,'read');
   if fid == -1
      fprintf(lid,'error\n\n');
      msg = sprintf('%sError open HDF-File "%s".',msg0,hdf);
      if no_out | log_file
         fprintf(lid,'%s\n\n',msg);
         if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
         if no_out, clear msg, end
      end
      return
   end
end

fprintf(lid,'ok\n');

%-------------------------------------------------------------

fprintf(lid,'Open NetCDF-File to write: "%s" ... ',cdf);

[nid,stat] = ncmex('open',cdf,'write');
 
if ~( ( nid > 0 ) & ( stat == 0 ) )
   fprintf(lid,'error\n\n');
   msg = sprintf('%sError open NetCDF-File: %s',msg0,cdf);
   hdfsd('end',fid);
   if no_out | log_file
      fprintf(lid,'%s\n\n',msg);
      if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
      if no_out, clear msg, end
   end
   return
end

fprintf(lid,'ok\n');

%***********************************************************
% Read/Write Data

fprintf(lid,'\nRead and Write %.0f Variables, Blocksize = %.0f\n',nvar,bsi);

ws = zeros(nvar,2);  % [ ReadStatus  WriteStatus ]

for ii = 1 : nvar

    vid = hdfsd('select',fid,ii-1);

    nd  = var{ii,3};

    siz = cat(2,dim{var{ii,4}+1,2});

    psz = prod(siz);

    str = sprintf(' %.0f ',siz);

    fprintf(lid,'\n%2.0f) %s : #%s [%s]',ii,var{ii,1},var{ii,2},str);

    if ~isequal(hdf_org(ii,:),var(ii,[1 2]))
        fprintf(lid,'   (%s : #%s)',hdf_org{ii,:});
    end

    fprintf(lid,'\n\n%11.0f Counts',psz);

    if any( siz == 0 )

       fprintf(lid,'\n\n');

    else

        is_char = strcmp(var{ii,2},'char');

        %---------------------------------------------------------------------
        % Counts per Dimension

        cnt = floor( siz ./ max(cumprod(siz)/bsi,1) );
        cnt = cnt + ( cnt == 0 );
        
        nc  = ceil( siz ./ cnt );  % Cycles per Dimension

        ncl  = prod(nc);
        csi  = prod(cnt);

        str = sprintf(' %.0f ',cnt);

        fprintf(lid,', Blocksize = %.0f = [%s], Cycles = %.0f',csi,str,ncl);

        %--------------------------------------------------------
        % Formats

        p10 = ( nc - 1*( nc > 1 ) ) .* cnt;
        p10 = floor( log(p10) / log(10) ) .* ( nc > 1 ) + 1;
        p10 = sprintf(' %%%.0f.0f ',p10);

        pcl = floor( log(ncl) / log(10) ) + 1;
        pcl = max(pcl,3);
        pcl = sprintf('%s%.0f.%.0fd','%',pcl,pcl);

        psc = floor( log(csi) / log(10) ) + 1;
        psc = sprintf('%s%.0f.0f','%',psc);

        psi = floor( log(psz) / log(10) ) + 1;
        psi = sprintf('%s%.0f.0f','%',psi);

        %--------------------------------------------------------

        if lid == 1
           fprintf(lid,'\n');
           ret = '\r';
        else
           ret = '\n';
        end

        form = sprintf('%s   #%s: %s %s Counts,  Start = [%s];  %s Counts total (%s%s) written ... ', ...
                       ret , pcl , '%s' , psc , p10 , psi , '%3.0f' , '%%' );

        %---------------------------------------------------------------------

        start = ones(1,nd);
        strd  = ones(1,nd);  % Stride

        nv = 0;
       
        for jj = 1 : ncl

            for kk = nd : -1 : 2
                if start(kk) > siz(kk)
                   start(kk)   = 1;
                   start(kk-1) = start(kk-1) + cnt(kk-1);
                end
            end

            count = min(cnt,siz-start+1);

            pc = prod(count);

            fprintf(lid,form,jj,'Read ',pc,start-1,nv,floor(100*nv/psz));

            [ val , ws(ii,1) ] = hdfsd( 'readdata' , vid , start-1 , strd , count );
            if ~( ws(ii,1) == -1 )
                if ~is_char & ~strcmp(class(val),'double')
                    val = double(val);
                end
                fprintf(lid,form,jj,'Write',pc,start-1,nv,floor(100*nv/psz));
                ws(ii,2) = ncmex('varput',nid,ii-1,start-1,count,val);
                if ~( ws(ii,2) == -1 )
                     nv = nv + pc;
                     fprintf(lid,form,jj,'Wrote',pc,start-1,nv,floor(100*nv/psz));
                end
            end

            if any( ws(ii,:) == -1 )
               break
            end

            start(nd) = start(nd) + cnt(nd);

        end
            
        if any( ws(ii,:) == -1 )
           fprintf(lid,'error');
        else
           fprintf(lid,'ok');
        end

        fprintf(lid,'\n');


    end

    %---------------------------------------------------

    status = hdfsd('endaccess',vid);

end

fprintf(lid,'\n');

%***********************************************************
% Close Files

status = hdfsd('end',fid);

if status == -1
   msg = cat(1,msg,{sprintf('Error close HDF-File "%s".',hdf)});
end

status = ncmex('close',nid);

if status == -1
   msg = cat(1,msg,{sprintf('Error close NetCDF-File: "%s"',cdf)});
end

%***********************************************************
% Check Status

str = { 'read'  'write' };
for ii = [ 1  2 ]
    jj = ( ws(:,ii) == -1 );
    if any(jj)
       jj = find(jj);
       msg = cat(1,msg,{sprintf('Error %s Variables:',str{1+(ii==2)})} , ...
                       { strhcat(var(jj,1),nlm) }          );
    end
end

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    if no_out | log_file
       fprintf(lid,'%s\n\n',msg);
       if log_file, setmark(lid,'error',msg0); else, fprintf(lid,char(7)); end
       if no_out, clear msg, end
    end
    return
end

msg = '';

if log_file
   setmark(lid,'end',msg0);
end

if no_out
   clear msg
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setmark(lid,mode,msg0)

 cl = clock;
str = datestr(datenum(cl(1),cl(2),cl(3),cl(4),cl(5),cl(6)),0);

ini = { '>>'  '<<'  '--' };

ind = [ 1  2 ] + [ 1 -1 ] * strcmp(mode,'end');
ind =   ind + ( 3 - ind ) * strcmp(mode,'error');

ini = ini(ind);

fprintf(lid,'\n%s %s%s %s %s\n\n',ini{1},msg0,upper(mode),str,ini{2});

if ~strcmp(mode,'start')
    try
       fclose(lid);
    end
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = check_name(str)

  fill = '_';

  df = double(fill);

  %---------------------------------------------------
  % Check for valid Characters
  % 0 .. 9  |  A .. Z   |  a .. z  |  _  |  FileSeparator

  str = double(str); 

  is_fill = ( str == df );

  ok = ( (  48 <= str  &  str <= 57  )  | ...
         (  65 <= str  &  str <= 90  )  | ...
         (  97 <= str  &  str <= 122 )  | ...
         is_fill );

  if all(ok)
     str = char(str);
     return
  end

  str(find(~ok)) = df;

  str(find(is_fill)) = 0;

  str = rmblank(char(str),2,fill);

  str(find(double(str)==0)) = fill;

  ok = ( double(str) == df );

  if any(ok)

     ok = find(ok);
     ok = ok(find(diff(ok)==1)+1);

     str(ok) = [];

  end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function typ = check_type(typ);

% CDF:  1. 'byte'      [        -128        127 ]      'int8'
%       2. 'char'      [           0        255 ]     'uint8'
%       3. 'short'     [      -32768      32767 ]     'int16'
%       4. 'long'      [ -2147483648 2147483647 ]     'int32'
%       5. 'float'
%       6. 'double'
%
% HDF:  char8
%       int*
%       uint*
%       float
%       double
%

if any(strcmp(typ,{'byte' 'char' 'short' 'long' 'float' 'double'}))
   return
end

if type_cmp(typ,'char');
   typ = 'char';
   return
end

if ~type_cmp(typ,'int')
    typ = 'none';
    return
end

s = ~type_cmp(typ,'uint');  % True for Signed Integer

v = '0123456789';
v = [ min(double(v))  max(double(v)) ];
v = ( ( v(1) <= double(typ) ) & ( double(typ) <= v(2) ) );

typ = typ(find(v));

try
   typ = eval(sprintf('[%s]',typ));
   ok  = ( prod(size(typ)) == 1 );
catch
   ok = 0;
end

if ~ok
    typ = 'none';
    return
end


% CDF:  1. 'byte'      [        -128        127 ]      'int8'
%       2. 'char'      [           0        255 ]     'uint8'
%       3. 'short'     [      -32768      32767 ]     'int16'
%       4. 'long'      [ -2147483648 2147483647 ]     'int32'

lim = [ 0  2^typ-1 ] - 2^(typ-1)*s;

ini = [ 8 ; 16 ; 32 ]; typ = { 'byte' 'short' 'long' };

ini = ( 2.^ini - 1 ) * [ 0  1 ] - 2.^(ini(:,[1 1])-1);

ok  = ( ( ini(:,1) <= lim(1) ) & ( lim(2) <= ini(:,2) ) );

if ~any(ok)
    typ = 'double';
else
    ok = find(ok);
    ok = ok(1);
    typ = typ{ok};
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
  
%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,hdf,cdf,bsi,lid] = checkin(varargin);

Nin = nargin;

msg = cell(0,1);

hdf = '';
cdf = '';
bsi = 1e6;
lid =  1;  % Matlab Command

%********************************************************
% Check HDF-File

if Nin < 1
   msg = cat(1,msg,{'Input HDF_File is missing.'});
   return
else
   v = varargin{1};
   if ~isempty(v)
       ok = chkstr(v);
       if ok
          ok = ( exist(v,'file') == 2 );
       end
       if ~ok
           msg = cat(1,msg,{'HDF_File must be a Name of an existing File.'});
       else
           hdf = v;
       end
    end
end


for ii = 2 : Nin

    m = '';

    v = varargin{ii};

    %----------------------------------------------
    if chkstr(v)

       if ~isempty(v)
           if strcmp(v(1),'@')
              lid = v(2:end);
           elseif isempty(cdf)
              cdf = v;
           else
              m = 'Duplicate Input for CDF_File.';
           end
       end

    %----------------------------------------------
    else

       ok = ( ( prod(size(v)) == 1 ) & isnumeric(v) );
       if ok
          ok = ( ( v > 0 ) & ( mod(v,1) == 0 ) & isfinite(v) );
       end

       if ok
          bsi = v;
       else
          m = 'Blocksize must be a positive Integer.';
       end

    end

    %----------------------------------------------
    if ~isempty(m)
        m = sprintf('Invalid %.0f. Input: %s',m);
        msg = cat(1,msg,{m});
    end

end

if ~isempty(msg) | isempty(hdf)
    return
end

%********************************************************

if isempty(cdf)
   [p,n,e] = fileparts(hdf);
   for ext = { '.cdf' '.nc' cat(2,e,'.cdf') }
       if ~strcmp(e,ext{1})
           cdf = cat(2,n,ext{1});
           break
       end
   end
end

%********************************************************

if isempty(lid)
   lid = 0;
elseif chkstr(lid)
   mode = { 'wt+' 'a' };
   exst = ( exist(lid,'file') == 2 );
   mode = mode{ 1 + exst };
   try
      lid = fopen(lid,mode);
   catch
      warning(sprintf('Can''t open LOG_File: %s',lid));
      lid = 1;
   end
end


%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%    A positive complex Value for DIM, to removes Blanks only from Start,
%    A negative complex Value for DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 160  32  13  10  9  0 ];  % [ NBSP Space CR LF TAB ZERO ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:);
    if ~all( ( abs(dim) == 1 ) |  ( abs(dim) == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must define 1. or 2. Dimension.' ];
    end
  end 
end

if Nin < 3
  cc = [ 160  32  13  10  9  0 ];  % [ NBSP  Space CR LF TAB ZERO ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = double(str);
  return
end


  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for ii = 1 : size(dim,1)

    d = dim(ii);

    s = sign(imag(d));  % Remove from wich Side:  1  0  -1 
 
    d = abs(d);

    jj = find( sum(blank,3-d) == si(3-d) );  % Columns with full Blanks

    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);

        %--------------------------------------------
        % Blank at Begin

        ind = ( 1 : nb );
        jj1 = find( ( ( jj == ind ) & ( s >= 0 ) ) );

        %--------------------------------------------
        % Blank at End

        ind = ind + si(d) - nb;
        jj2 = find( ( ( jj == ind ) & ( s <= 0 ) ) );

        %--------------------------------------------

        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);


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



%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
   str = cellstr(str);
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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

   
