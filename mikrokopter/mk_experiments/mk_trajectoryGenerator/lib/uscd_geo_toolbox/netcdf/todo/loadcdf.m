function [msg,cnf] = loadcdf(varargin);

% New, UNCOMPLETED version of LOAD_CDF

Nout = nargout;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));

nl = char(10);

nlm = sprintf('%s%s',nl,char(32*ones(size(msg0))));

%***********************************************************

[msg,file,cnf] = checkin(varargin{:});

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    if Nout == 0
       error(msg)
    end
    return
end

%***********************************************************

[msg,cnf] = checkfile(file,cnf);

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,strhcat(msg,nlm));
    if Nout == 0
       error(msg)
    end
    return
end

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,cnf] = checkin(varargin)

msg = cell(0,1);

file = '';

Nin = nargin;

%---------------------------------------------------------
% ConfigurationStructure

ord = { 'x' 'y' 'z' 't' };
prm = [  2   1   3   4  ];
        
d = struct(  'id' , {[]} , ...    % DimID
            'len' , {[]} , ...    % DimLength
            'scs' , {[]}       ); % StartCountStride

dim      = ord([1 1],:);
dim(2,:) = { {d} };
dim      = struct(dim{:});

v = struct(  'id' , {[]} , ...    % VarID
            'val' , {[]} , ...    % Value NetCDF
            'req' , {[]} , ...    % Value Requested by Input
            'itp' , { 0}       ); % TRUE for Interpolate

var      = ord([1 1],:);
var(2,:) = { {v} };
var      = struct(var{:});

dat = struct( 'name' , { '' } , ...    % KeyName
              'var'  , { '' } , ...    % VarName (NetCDF)
               'id'  , { [] } , ...    % VarID
              'did'  , { {} } , ...    % VarDimIDs
              'scs'  , { [] } , ...    % VarDimStartCountStride
              'prm'  , { [] } , ...    % Permutation
              'siz'  , { [] } , ...    % Size after Permutation
              'itp'  , { [] } , ...    % TRUE for Interpolate
              'scl'  , {  1 } , ...    % scale_factor
              'off'  , {  0 } , ...    % add_offset
              'miss' , { [] } , ...    % missing_value
              'fill' , { [] }       ); % _FillValue

cnf = struct( 'file' , {  '' } , ...
               'ord' , { ord } , ...
               'prm' , { prm } , ...
               'dim' , { dim } , ...
               'var' , { var } , ...
               'dat' , { dat(zeros(0,1)) } , ...
              'miss' , { NaN } , ...
              'fill' , {  [] } , ...
              'mode' , { 'interp' } );

%---------------------------------------------------------

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

if ~isempty(msg) | ( Nin < 2 )
    return
end

%---------------------------------------------------------

iv = 1;

% Requested Inputs
sets = { 'd'   'v'
         'dim' 'var'
         'scs' 'req' };

while iv < Nin

      iv = iv + 1;

      val = varargin{iv};

        m = '';

      if ~chkstr(val,1)
          m = 'String required.';
      else

          v1 = lower(val(1));
          if size(val,2) == 1
             v2 = '';
          else
             v2 = lower(val(2));
          end

          %------------------------------------------------------------------------
          if any(strcmp(v1,{'o' 's' 'i'}))  %  orig | sect | interp

            cnf.mode = v1;

          %------------------------------------------------------------------------
          elseif iv == Nin

            m = sprintf('Missing Value for Option "%s".',val);

          %------------------------------------------------------------------------
          else

            iv  = iv + 1;
              v = varargin{iv};

            switch v1

              %--------------------------------------------------------------------
              case {'m' 'f'}   %  miss | fill
                if ( isnumeric(v) & ( prod(size(v)) == 1 ) ) | isempty(v)
                   if strcmp(v1,'m')
                      cnf.miss = v;
                   else
                      cnf.fill = v;
                   end
                else
                   m = ' must be a single numeric';
                end

              %--------------------------------------------------------------------
              case 'd'         %  data

                if ~isempty(v);
                    [ok,v] = chkcstr(v,0);
                    if ok
                        v = v(:);
                        n = size(v,1);
                        d = dat(ones(n,1));
                        for ii = 1 : n
                            d(ii).name = v{ii};
                        end
                        cnf.dat = cat( 1 , cnf.dat , d );
                    else
                        m = ' must be a StringArray';
                    end
                end

              %--------------------------------------------------------------------
              case ord         %  x | y | z | t

                if ~any(strcmp(v2,sets(1,:)))

                    m = sprintf('Unkown Option "%s".',val);

                elseif ~isempty(v)

                    ok = isnumeric(v);
                    if ok
                       ok = all( isfinite(v) );
                    end

                    if ~ok
                        m = ' must be finite numeric';
                    elseif strcmp(v2,'d')
                        v = v(:)'
                        if ~( any( size(v,2) == [ 1  3 ] ) & ( round(v) == v ) )
                            m = ' must have 1 or 3 Integer';
                        elseif size(v,2) == 1
                            if abs(v) == 0
                               m = ' must have nonzero Stride';
                            end
                        elseif size(v,2) == 3
                            if abs(v(3)) == 0
                               m = ' must have nonzero Stride';
                            end
                            if ~( v(2) >= 0 )
                                if ~isempty(m)
                                    m = sprintf('%s and Counts >= 0',m);
                                else
                                   m = ' must have Counts >= 0';
                                end
                            end
                        end
                    elseif strcmp(v2,'v')
                        if ~chkmono(v)
                            m = ' must be a Vector';
                        else
                            v = v(:);
                        end
                    end

                    if isempty(m)
                       jj = find(strcmp(sets(1,:),v2));
                        c = getfield(cnf,sets{2,jj});
                        d = getfield(c,v1);
                        d = setfield(d,sets{3,jj},v);
                        c = setfield(c,v1,d);
                        cnf = setfield(cnf,sets{2,jj},c);
                    end

                end

              %--------------------------------------------------------------------
              otherwise

                m = sprintf('Unkown Option "%s".',val);

            end
            % v1

            if ~isempty(m)
                if strcmp(m(1),' ')
                   m = sprintf('Value for Option "%s" %s.',val,m);
                end
            end

          end
      end

      if ~isempty(m)
          msg = cat(1,msg,{sprintf('Invalid %.0f. Input: %s',iv,m)});
      end

end

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [msg,cnf] = checkfile(file,cnf);

%**********************************************************
% Check Configuration: CNF_CDF

[msg,c,dim,var,att] = cnf_cdf(file,'check');

if ~isempty(msg)
    return
end

cnf.file = c{1,1};

%**********************************************************
% Check Dim-Var-Configuration

ok = ~strcmp(c(:,[2 4]),'');
jj = xor(ok(:,1),ok(:,2));
if any(jj)
   jj = find(jj);
   msg = sprintf('Dimension AND Variable for %s must be defined' , ...
                  strhcat(upper(ord(jj)),' ') );
   return
end

%**********************************************************
% Open NetCDF-File

fid = ncmex('open',cnf.file,'nowrite');
if fid == -1
   msg = sprintf('Can''t open NetCDF-File "%s" for reading.',cnf.file);
   return
end

%**********************************************************
% Get Dimensions and Variables for Order

ind = find(ok(:,1));
ind = ind(:)';

n = size(ind,2);

m    = cell(2,n);  % Messages
m(:) = { '' };

for ii = ind

    ord = cnf.ord{ii};

    str = upper(ord);

    %**********************************************************
    % Dimension

    d = getfield(cnf.dim,ord);

    d.id = find(strcmp(dim(:,1),c{ii,2})) - 1;

    len   = c{ii,3};
    scs   = d.scs;

    auto_dim = ( size(scs,2) <= 1 );

    if isempty(scs)
       scs = [ 1  len  1 ];
    elseif auto_dim
       s = sign(scs);
       scs = abs(scs);
       scs = [ 1  ceil(len/scs) scs ];
       if s == -1
          scs(1) = len - (scs(2)-1)*scs(3);
       end
    else
       scs(1) = scs(1) + len * ( scs(1) <= 0 );
       if scs(3) < 0
          scs(3) = abs(scs(3));
          scs(1) = scs(1) - (scs(2)-1)*scs(3);
       end
    end

    scs(3) = scs(3) + ( 1 - scs(3) ) * ( scs(2) == 1 );

    if ~( ( 1 <= scs(1) )  &  ( scs(1) <= len ) & ...
          ( scs(1)+(scs(2)-1)*scs(3) <= len )  )
        m{1,ii} = sprintf('Requested Dimension %s out of range [ %.0f  %.0f ].',str,1,len);
    end

    scs(1) = scs(1) - 1;  % Start with ZERO

    %**********************************************************
    % Variable

    v = getfield(cnf.var,ord);

    v.id = find(strcmp(var(:,1),c{ii,2})) - 1;

    is_req = ~isempty(v.req);

    if isempty(m{1,ii})

       %-------------------------------------------------------
       % Check for Periodic

       pr = ( strcmp(ord,'x') & ( is_req | auto_dim ) & ( len > 1 ) );

       %-------------------------------------------------------
       % Get Value

       if pr
           [val,stat] = ncmex('vargetg',fid,v.id,0,len,1);
       else
           [val,stat] = ncmex('vargetg',fid,v.id,scs(1),scs(2),scs(3));
       end

       if stat == -1

          m{2,ii} = sprintf('Error read NetCDF-Variable %s.',str);

       else

          %----------------------------------------------------
          % Check Value

          [val,at] = check_att(val(:),att{v.id+1});

          [ok,sgn] = chkmono(val);
          if ~ok
              m{2,ii} = sprintf('NetCDF-Variable %s must be monotonic.',str);
          else
              v.val = val;
              if pr
                 p  = 360;
                 vr = abs( val(1) - val(len) );
                 dv = diff(val);
                 dv = max(abs(dv([1 len-1])));
                 pr = ( ( p-2*dv < vr ) & ( vr <= p ) );
                 ie = ( vr == p );
              end
          end

       end

    end

    %**********************************************************
    % Compare NetCDF-Variable with requested Variable

    if isempty(m{1,ii}) & isempty(m{2,ii}) & ( pr | ( is_req & ( len > 1 ) )

       req = v.req;

       %-------------------------------------------------------
       % Periodic
       %-------------------------------------------------------
       if pr

          if ~is_req

             scs = cat( 2 , scs , 0 );
             if auto_dim
                dn = scs(1) + [ -1  scs(2) ]*scs(3) + [ 1  -1 ] *  len + [ 0  ie ];
                scs = cat( 1 , [ dn(1)  1  1  -p ] , ...      % One before
                               [ scs             ] , ...
                               [ dn(2)  1  1   p ]        );  % One behind
             end

          else
%auto_dim ???

             % Sort Value
             ivl = sgn * ( 1 : len ) + ( len + 1 ) * ( sgn == -1 );

             vlb = val(ivl(1));  % ValBase == min(val)
             val = val - vlb;    % [ 0 .. p ]

             vrb = req(1);

             req = req - vlb;

             req = req - p * floor(req/p);  % [ 0 .. p )  

             if size(req,1) == 1

             else

                mvr = sign( req(2) - req(1) );

                 rb = cat( 1 , 0 , ~( sign(diff(req,1,1)) == mvr ) );

                req = req + cumsum( mvr*p*rb , 1 );

                req 
                i0  = find( ~( sign(diff(req,1,1)) == mvr ) );
                i0  = cat( 1 , 1 , i0+1 , size(req,1)+1 );
                n0  = size(i0,1) - 1;

                scs = scs(ones(n0,1),:);

                add = cat( 1 , 0 , ~( sign(diff(req,1,1)) == mvr ) );

                add = cumsum(add,1);

             end

          end

          %----------------------------------------------------
          % Extract Value

          val = zeros(0,1);
          for jj = 1 : size(scs,1)
              vn = ncmex('vargetg',fid,v.id,scs(jj,1),scs(jj,2),scs(jj,3));
              val = cat( 1 , val , check_att(vn(:),at) + scs(jj,4) );
          end

       %-------------------------------------------------------
       % NonPeriodic
       %-------------------------------------------------------
       elseif is_req

           jj = ( 1 : (len-1) );

           vr = [ min(req)  max(req) ];
           ok = ( ( vr(1) <= val ) & ( val <= vr(2) ) );

           if any(ok)
              dk = diff(ok);
              ok(jj+0) = ( ok(jj+0) | ( dk ==  1 ) );  % One before
              ok(jj+1) = ( ok(jj+1) | ( dk == -1 ) );  % One behind
           else
              vr = ( vr(1) + vr(2) ) / 2;
              dv = abs( val - vr );
              ok = ( dv == min(dv) );
              if ~( ( sum(ok) == 2 ) | any(ok([1 len])) )
                  dk = diff(ok);
                  ok(jj+0) = ( ok(jj+0) | ( ( dk ==  1 ) & ( sgn*val(jj+1) > sgn*vr ) ) );  % One before
                  ok(jj+1) = ( ok(jj+1) | ( ( dk == -1 ) & ( sgn*val(jj+0) < sgn*vr ) ) );  % One behind
              end
           end

           scs(2) = sum(ok);                          % Count
           ok     = find(ok);
           scs(1) = scs(1) + ( ok(1) - 1 ) * scs(3);  % Start
           val    = val(ok);

       end
       %-------------------------------------------------------

       v.val   = val;

       if is_req
          v.itp = ~isequal(val,req);
       end

    end

    cnf.var = setfield(cnf.var,ord,v);

    d.len   = len;
    d.scs   = scs;

    cnf.dim = setfield(cnf.dim,ord,d);

end

%**********************************************************

fid = ncmex('close',cnf.file);

m  = m(:);
jj = ~strcmp(m,'');
if any(jj)
   jj = find(jj);
   msg = m(jj);
   return
end

%**********************************************************

if isempty(cnf.dat)
   return
end

v = c{2,1};  % { KeyName VarName }

if isempty(v)
   msg = sprintf('No DataVariables defined in ConfigFile: %s',file);
   return
end

%**********************************************************
% Get DimensionIDs for Order

ord = cnf.ord(cnf.prm);

no = size(ord,2);

did = NaN*zeros(1,no);

for ii = 1 : no
    d = getfield(cnf.dim,ord{ii});
    if ~isempty(d.id)
        did(ii) = d.id;
    end
end

%**********************************************************
% Check Data

%----------------------------------------------------------
% Check for periodic X-Dimension

xid = cnf.dim.x.id;

if isempty(xid)
   nxd = 1;
   xid = NaN;
else
   nxd = size(cnf.dim.x.scs,1);
end

isx = ( nxd > 1 );

%----------------------------------------------------------

dat = cnf.dat;

nn = size(dat,1);

m    = cell(nn,1);  % Messages
m(:) = { '' };

for ii = 1 : nn

    d = dat(ii);

    str = sprintf('DataVariable "%s"',d.name);

    jj = strcmp(v(:,1),d.name);
    if ~any(jj)
       m{ii} = sprintf('%s not defined in ConfigFile: %s',str,file);
    else
       jj    = find(jj);
       d.var = v{jj(1),2};
       jj    = strcmp(var(:,1),d.var);
       if ~any(jj)
           m{ii} = sprintf('NetCDF-Variable "%s" for %s not defined in NetCDF-File: %s', ...
                            d.var , str , cnf.file );
       else
          jj    = find(jj);
          d.id  = jj - 1;
          d.did = var{jj,4};  % DimensionIDs
          if ~isempty(att{jj})
              [h,at] = check_att([],att{jj});
              d.scl  = at.scl;
              d.off  = at.off;
              d.miss = at.miss;
              d.fill = at.fill;
          end
       end
    end

    if isempty(m{ii})

       nd = size(d.did,2);

       d.scs =  cell(1,nd);
       d.siz = zeros(1,nd);
       d.prm = ( 1 : nd ) + no;
       d.itp = zeros(1,nd);

       exd = ( isx & any( d.did == xid ) );
       for jj = 1 : nd
           jp = nd + 1 - jj;
           kk = ( did == d.did(jj) );
           if any(kk)
              kk = find(kk);
              d.scs{jj} = getfield(getfield(cnf.dim,ord{kk}),'scs');
              d.prm(jp) = kk;
              d.itp(jp) = getfield(getfield(cnf.var,ord{kk}),'itp');
           else
              d.scs{jj} = [ 0  dim{d.did(jj)+1,2} 1 ];
           end
           d.siz(jp) = sum(d.scs{jj}(:,2));
           if ~( d.did(jj) == xid ) & exd
               d.scs{jj} = d.scs{jj}(ones(1,nxd),:);
           end
       end

       [h,d.prm] = sort(d.prm);

       d.siz = d.siz(d.prm);
       d.itp = d.itp(d.prm);

    end

    dat(ii) = d;

end

cnf.dat = dat;

jj = ~strcmp(m,'');
if any(jj)
   jj = find(jj);
   msg = m(jj);
end

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,cnf] = check_att(x,cnf,miss,fill)

if nargin < 3
   miss = NaN;
end

if nargin < 4
   fill = [];
end

if iscell(cnf)

   att = cnf;

   cnf = { 'miss'  {[]}  'missing_value'
           'fill'  {[]}  '_FillValue'
           'scl'   {1}   'scale_factor'
           'off'   {0}   'add_offset'    };

   for ii = 1 : size(cnf,1)
       jj = strcmp( att(:,1) , cnf{ii,3} );
       if any(jj)
          jj = find(jj);
          cnf{ii,2} = att(jj,3);
       end
   end

   cnf = permute(cnf(:,[1 2]),[2 1]);
   cnf = struct(cnf{:});

end

is_miss = ~( isempty(cnf.miss) | isempty(miss) );
is_fill = ~( isempty(cnf.fill) | isempty(fill) );

if ( ( cnf.scl == 1 ) & ( cnf.off == 0 ) & ...
     ~is_miss & ~is_fill ) | isempty(x)
   return
end

if is_miss
   jm = find( x == cnf.miss );
end

if is_fill
   jf = find( x == cnf.fill );
end

if ~( cnf.scl == 1 )
    x = x * cnf.scl;
end

if ~( cnf.off == 0 )
    x = x + cnf.off;
end

if is_fill
   x(jf) = fill;
end

if is_miss
   x(jm) = miss;
end

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,m] = chkmono(x)

m = 0;

s = size(x); p = prod(s);

ok = ( p <= 1 );
                 
if ~ok
    ok = ( p == max(s) );
    if ok
       s = sign(diff(x));
       ok = ( all( s == 1 )  |  all( s == -1 ) );
       if ok
          m = s(1);
       end
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


