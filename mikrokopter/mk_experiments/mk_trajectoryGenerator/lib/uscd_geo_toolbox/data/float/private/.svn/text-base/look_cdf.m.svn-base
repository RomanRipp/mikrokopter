function [msg,varargout] = look_cdf(varargin);

% LOOK_CDF  Reads Basic Informations from a CDF-File
%
% [ Msg , DIM , VAR , ATT , TXT ] = LOOK_CDF( FILENAME )
%
%  Reads Informations from NetCDF-File, specified by FILENAME.
%
%  DIM : Information for N Dimensions
%  DIM = { DimName  DimLength  '<UNLIMITED>' }    [ N   by 3 ] CellArray
%
%  VAR : Information for M Variables
%  VAR = { VarName VarType Ndim [dim] Natt Val }  [ M   by 6 ] CellArray 
%
%  Val is a  ShortValue: - a single Value in case of Scalar
%                        - the [ First Last ] - Values in case of Vector
%                        - a String in case of short String
%
%  ATT : Information for Variable- and Global-Attributes
%  ATT = { VarAtt ; GlbAtt }                      [ M+1 by 1 ] CellArray
%
%  Att = { AttName AttType AttValue AttString }   [ K   by 4 ] CellArray
%
%
%  TXT is an CellArray of Strings,
%        containing the Informations as ASCII-Text.
%
%  Msg contains the ErrorMessages.
%
%-----------------------------------------------------------------------
% DisplayMode
%
%  LOOK_CDF( ... , MODE ) , where Mode is a scalar, specifies the Display:
%
%  MODE == 0   displays the TEXT on the Terminal,
%  MODE ~= 0   displays the TEXT in a ListBox of a Figure.
%
%  In case of ( MODE ~= 0 ), the second Output is the Handle 
%   of the LOOK_CDF-Figure: [ Msg , FIG , ... ] = LOOK_CDF( ... )
%
%  Use the Value of FIG for MODE, to reset the LOOK_CDF-Figure with
%    the Informations about a new NetCDF-File.
%
%-----------------------------------------------------------------------
% FileMode
%
%  LOOK_CDF( ... , OUTFILENAME )  saves the TEXT into an ASCII-File,
%                                 specified by OUTFILENAME.
%
%  By default, a new File is created, an existing File with the same Name
%    will be overwritten. If the OUTFILENAME begins with an '@',
%    the TEXT will be append on an existing File.
%
%  The Inputs MODE and OUTFILENAME are optional 
%   and could be used together.
%
%-----------------------------------------------------------------------
%
%  See also: READ_CDF, ASSIGN_CDF, WRITE_CDF, CREATE_CDF, NCMEX
%
%-----------------------------------------------------------------------
%
%  DataTypes:  1. 'byte'      [        -128        127 ]      'int8'
%              2. 'char'      [           0        255 ]     'uint8'
%              3. 'short'     [      -32768      32767 ]     'int16'
%              4. 'long'      [ -2147483648 2147483647 ]     'int32'
%              5. 'float'       32 bit floating point
%              6. 'double'      64 bit floating point
%


Nin  = nargin;
Nout = nargout;

NoOut = ( Nout == 0 );

Nout = Nout - 1;

msg = '';

varargout = cell(1,Nout);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));



nl   = char(10);

%*************************************************************
% Check for CallBack

if Nin > 1
   fcn = varargin{1};
   if chkstr(fcn,1)
      if ( double(fcn(1)) == double('#') ) & ( size(fcn,2) > 1 )
         try
            feval(fcn(2:end),varargin{2:end});
         catch
            warning(sprintf('%s Error Callback: %s\n%s',msg0,fcn,lasterr));
         end
         return
      end
   end
end

%*************************************************************
% Check Inputs

[msg,file,mode,outfile] = checkin(varargin{:});

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,msg);
    if NoOut
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
       clear msg
    end
    return
end

if NoOut  &  isempty(outfile)  &  isempty(mode)
   mode = -1;    % Default: ListBoxWindow
end

%*************************************************************
% Get Info

is_mode = ~isempty(mode);
is_file = ~isempty(outfile);

is_text  = ( is_mode | is_file | ( Nout > 3 ) );

[msg,dim,var,att,file] = look( file , max(Nout,3*is_text) );

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,msg);
    if NoOut
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
       clear msg
    end
    return
end

%*************************************************************
% Get Text

if is_text
   [txt,vnr] = gettext(file,dim,var,att);
else
   txt = cell(0,1);
   vnr = zeros(0,2);
end

%*************************************************************
% Write File

if is_file
   msg = write(outfile,txt);
   if ~isempty(msg)
       msg = sprintf('%s%s',msg0,msg);
   end
end

%*************************************************************
% Display

fig = [];
if is_mode
   if mode == 0
   % Terminal
     fprintf(1,char(10));
     fprintf(1,strhcat(txt,char(10)));
   else
   % Figure with ListBox
     fig = list(file,dim,var,txt,vnr,mode);
   end
end

%*************************************************************
% OutPut

if NoOut
   if ~isempty(msg)
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
   end
   clear msg
end

if Nout <= 0
   return
end

if  ~isempty(fig)
   v = { fig dim var att txt vnr };
else
   v = { dim var att txt vnr };
end

n = min( Nout , prod(size(v)) );

varargout(1:n) = v(1:n);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [msg,file,mode,outfile] = checkin(varargin)

% CHECKIN  First Check of Inputs

Nin = nargin;

msg     = '';
file    = '';
mode    = [];
outfile = '';

%-------------------------------------------------------------
% Check FileName

if Nin < 1
   file = '';
else
   file = varargin{1};
end

if isempty(file)
   [f,p] = uigetfile('*.nc','Select a NetCDF-File ...');
   if isequal(f,0)
      return
   end
   file = fullfile(p,f);
end
       
if ~chkstr(file,0)
    msg = 'Input FileName must be a String.';
    return
end

for e = { ''  '.cdf'  '.nc'  '.cnf' };
    f  = cat( 2 , file , e{1} );
    ok = ( exist(f) == 2 );
    if ok
       break
    end
end

if ok
   file = f;
else
   msg = sprintf('File "%s" doesn''t exist.',file);
   return
end

%-------------------------------------------------------------
% Check Options

for ii = 2 : Nin
    if chkstr(varargin{ii},1)
       outfile = varargin{ii};
    elseif isnumeric(varargin{ii}) & ( prod(size(varargin{ii})) == 1 )    
          mode = varargin{ii};
    end
end


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,dim,var,att,file] = look(file,mode);

% LOOK  Get Information about DIM, VAR, ATT

msg = '';
dim = cell(0,3);
var = cell(0,6);
att = cell(0,1);

%-------------------------------------------------------------
% DataTypes

types = { 'byte'  'char'  'short'  'long'  'float'  'double' };

%-------------------------------------------------------------
% Quiet Mode

ncmex('setopts',0);

%-------------------------------------------------------------
% Open NetCDF file

[fid,stat] = ncmex('open',file,'nowrite');

ok = ( ( fid > 0 ) & ( stat == 0 ) );
if ~ok
    ff = which(file);
    [fid,stat] = ncmex('open',ff,'nowrite');
    ok = ( ( fid > 0 ) & ( stat == 0 ) );
end

if ~ok

    % Try to load KeyWord "FileName" from ConfigFile, 
    %   but NOT if called by CNF_CDF

    try
       c  = caller;
       ok = ~any(strcmp(c(:,1),'cnf_cdf'));  % !!!!!!
    end

    if ok
       ff = '';
       try   % NetCDF-ConfigFile
           [m,v,ff] = load_asc(file,'FileName',{'#' ':'},'%');
       end
       ok = chkstr(ff,1);
    end

    if ok
       file = ff;
       [fid,stat] = ncmex('open',file,'nowrite');
       ok = ( ( fid > 0 ) & ( stat == 0 ) );
       if ~ok
           ff  = which(file);
           [fid,stat] = ncmex('open',ff,'nowrite');
           ok = ( ( fid > 0 ) & ( stat == 0 ) );
       end
    end

end

if ~ok   
    msg = sprintf('Error open File "%s" as NetCDF-File.',file);
    return
end

if mode == 0
   status = ncmex('close',fid);
   return
end

%-------------------------------------------------------------
% NetCDF-Inquire

[ndim,nvar,natt,recdim,status] = ncmex('inquire',fid);

if status == -1
   status = ncmex('close',fid);
   msg = 'NCINQUIRE: status = -1 '; 
   return
end


%--------------------------------------------------------------
% get Information about Dimensions
%
dim      = cell(ndim,3);     % { Name  Size "UNLIMITED" }
dim(:,3) = { '' };

for ii = 1 : ndim
    [dim{ii,1},dim{ii,2},status] = ncmex('diminq',fid,ii-1);
end

if recdim >= 0
   dim{recdim+1,3} = 'UNLIMITED';
end

if mode == 1
   status = ncmex('close',fid);
   return
end

%--------------------------------------------------------------
% get Information about Variables
%

var   = cell(nvar,6); % { Name Type Ndim [dim] Natt [Range] }

att  = cell(nvar+1,1);  % Attribute, nvar+1 == global

for ii = 1 : nvar

 % Name      Type      Ndim      [dim]     Natt 
 [var{ii,1},var{ii,2},var{ii,3},var{ii,4},var{ii,5},status] = ncmex('varinq',fid,ii-1);

  var(ii,2) = types(var{ii,2});

  %-----------------------------------------------------------
  % Read Attributes
  att{ii} = cell(var{ii,5},4);   % { Name Type Value String }

  for jj = 1 : var{ii,5}

      [ att{ii}{jj,1} ,       status ] = ncmex('attname',fid,ii-1,jj-1); % Name
      [ att{ii}{jj,2} , len , status ] = ncmex('attinq' ,fid,ii-1,jj-1); % Type
      [ att{ii}{jj,3} ,       status ] = ncmex('attget' ,fid,ii-1,jj-1); % Value 

      att{ii}(jj,2) =   types( att{ii}{jj,2}   );
      att{ii}{jj,4} = att2str( att{ii}{jj,1:3} );

      if      strcmp( att{ii}(jj,2) , 'char' )  &  isequal(att{ii}{jj,3},'\0')  &  ...
          any(strcmp( att{ii}{jj,1} , {'missing_value' '_FillValue'} ))

          att{ii}{jj,3} = char(0);
    
      end

  end

  %-----------------------------------------------------------
  % Check for UINT8

  is_num = any( strcmp(att{ii}(:,1),'scale_factor') | ...
                strcmp(att{ii}(:,1),'add_offset')          );
 
  is_char = strcmp(var{ii,2},'char');

  is_byte = ( is_char  &   is_num );  % UINT8
  is_char = ( is_char  &  ~is_num );

  if is_byte
     jj = ( strcmp(att{ii}(:,2),'char') & ...
            ( strcmp(att{ii}(:,1),'valid_range')  | ...
              strcmp(att{ii}(:,1),'valid_min')    | ...
              strcmp(att{ii}(:,1),'valid_max')    | ...
              strcmp(att{ii}(:,1),'step_size')    | ...
              strcmp(att{ii}(:,1),'add_offset')   | ...
              strcmp(att{ii}(:,1),'scale_factor') | ...
              strcmp(att{ii}(:,1),'_FillValue')   | ...
              strcmp(att{ii}(:,1),'missing_value')      ) );
     if any(jj)
        jj = find(jj);
        for kk = jj(:)'
            att{ii}{kk,3} =   uint8(att{ii}{kk,3});
            att{ii}{kk,4} = att2str(att{ii}{kk,1:3});
        end
     end
  end

  %-----------------------------------------------------------
  % Variable with single Dimension !!!!!!!!!!!!!!!!

  val = [];

  if var{ii,3} == 0

     val = ncmex('varget1',fid,ii-1,1);

  elseif ( sum( ( cat(2,dim{var{ii,4}+1,2}) == 1 ) , 2 ) >= (var{ii,3}-1) )

     dl = cat(2,dim{var{ii,4}+1,2});

     if all(dl) & ~any(strcmp(dim(var{ii,4}+1,3),'UNLIMITED'))

        if is_char
           cnt = dl;
           str = ones(1,var{ii,3});
        else
           cnt = ones(1,var{ii,3}) + ( dl > 1 );
           str = dl - 1 * ( dl > 1 );  
        end

        val = ncmex('vargetg',fid,ii-1,...
                    zeros(1,var{ii,3}),cnt,str);
        val = val(:)';

     end

  end
 
  %-----------------------------------------------------------
  if ~isempty(val)

      %---------------------------------------------------
      if is_char & ( cnt <= 1024 )
      %---------------------------------------------------
      % Char  

        if ~isempty(att{ii})
            kk = find( strcmp( att{ii}(:,1) , 'missing_value' ) );
            if ~isempty(kk)
                if prod(size(att{ii}{kk(1),3})) == 1
                   val( find( ( val == att{ii}{kk(1),3} ) ) ) = 32;
                end
            end

            kk = find( strcmp( att{ii}(:,1) , '_FillValue' ) );
            if ~isempty(kk) 
                if prod(size(att{ii}{kk(1),3})) == 1
                   val( find( ( val == att{ii}{kk(1),3} ) ) ) = 32;
                end
            end
        end

        val = checkchar(val,1);   % incl. RMBLANK

        ok = isempty(val);
        if ~ok
            jj = find( double(val) == 10 );
            jj = cat( 2 , 0 , jj , size(val,2)+1 ); 
            ok = ( max(diff(jj)) <= 121 );
        end 

        if ok
           var{ii,6} = val;
        end

      %---------------------------------------------------
      elseif ~is_char
      %---------------------------------------------------

      %  [ Byte  |  UINT8  |  Short  |  Long  |  Float  |  Double ]
      %  Get [ First  Last ] Value

        if is_byte
           val = double(val);
        end

        fl = NaN; 
        if ~isempty(att{ii})
            kk = find( strcmp( att{ii}(:,1) , 'missing_value' ) );
            if ~isempty(kk)
                if prod(size(att{ii}{kk(1),3})) == 1
                   val( find( ( val == att{ii}{kk(1),3} ) ) ) = NaN;
                end
            end
            kk = find( strcmp( att{ii}(:,1) , '_FillValue' ) );
            if ~isempty(kk) 
                if prod(size(att{ii}{kk(1),3})) == 1
                   fl = att{ii}{kk(1),3};
                end
            end
        end

        if all( isnan(val) | ( val == fl ) ), val = []; end

        if ~isempty(val)
            ind = [ 1  size(val,2) ];
            val = val( ind( 1 : ( 1 + (ind(1)<ind(2)) ) ) );
            if is_byte
               val = uint8(val);
            end
            var{ii,6} = val;
        end

      %-----------------------------------------------------------
      end
      %-----------------------------------------------------------
  
  end

end


%--------------------------------------------------------------
% Global Attributes

if mode >= 3

        ii  = nvar+1;
    att{ii} = cell(natt,4);   % { Name Type Value String }

    for jj = 1:natt
        [ att{ii}{jj,1} ,       status ] = ncmex('attname',fid,'global',jj-1); 
        [ att{ii}{jj,2} , len , status ] = ncmex('attinq' ,fid,'global',jj-1);
        [ att{ii}{jj,3} ,       status ] = ncmex('attget' ,fid,'global',jj-1);
          att{ii}(jj,2) =   types(att{ii}{jj,2});
          att{ii}{jj,4} = att2str(att{ii}{jj,1:3});
    end

end

%--------------------------------------------------------------
% Close File

 status = ncmex('close',fid);

 
%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [txt,vnr] = gettext(file,dim,var,att);

% GETTEXT   Create Text

ndim = size(dim,1);
nvar = size(var,1);

vnr = zeros(nvar,2);  % [ ListNr  DescrNr ]

tab  = char(32*ones(1,4));

splt  = 60.2 + i*240;  % [ SplitLength MaxLength ] of AttributeString

app   = 's';

%----------------------------------------------------------------
% FileName

f = which(file);
if ~isempty(f)
    file = f;
end

siz = dir(file);
siz = siz.bytes;

scl = 10 * [ 0  1  2 ];
scl = 2 .^ scl;

ii = sum( siz >= scl );
ii = max(ii,1);

pre = { ''  'k'  'M' };
fmt = sprintf('%s.%.0ff','%',3*(ii>1));
fmt = sprintf(' %s %sBytes',fmt,pre{ii});

siz = sprintf(fmt,siz/scl(ii));

txt = { ... 
  ''
 sprintf( 'File: %s' , file )
 sprintf( 'Size: %s' , siz  )
 sprintf( '------%s' , char('-'+0*file) )
  '' };

%----------------------------------------------------------------
% Global Attributes

if ~isempty(att{nvar+1});

    natt = size(att{nvar+1},1);

    txt = cat( 1 , txt , ...
              { sprintf(' Global Attribute%s:',app(1:(end*(natt>1))))
                ' ------------------' 
                ''   } );

    s0 = size( char(att{nvar+1}(:,1)) , 2 );
    nt = sprintf( '\n %s%s' , tab , char(32*ones(1,s0+3)) );

    for jj = 1 : natt

        s2 = size(att{nvar+1}{jj,1},2);

        tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

        str = att{nvar+1}{jj,4};

        str = splitstr(str,splt,'; ');

        str = strrep( str , char(10) , ...
                    [ nt char(32*ones(1,strcmp(str(1),'"'))) ] );

        str = sprintf( '%s%s:%s  %s' , tb , att{nvar+1}{jj,1}    , ...
                                            att{nvar+1}{jj,2}(1) , str );

        [m,str] = char2cell(str);  
   
        txt = cat( 1 , txt , str );

    end

    txt = cat( 1 , txt , ...
              { '' 
                '' } );

end


%----------------------------------------------------------------
% Dimensions

txt = cat( 1 , txt , ...
{ sprintf(' %2.0f Dimension%s:',ndim,app(1:(end*(ndim>1))))
  ' --------------' 
  ''  } );

s0 = size( char(dim(:,1)) , 2 );

for ii = 1 : ndim

    s2 = size(dim{ii,1},2);

    tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

    str = sprintf('%3.0f)%s%s = %.0f%s   %s',ii,tb,dim{ii,1},dim{ii,2},tab,dim{ii,3});

    txt = cat( 1 , txt , {str} );

end



%----------------------------------------------------------------
% Variables

txt = cat( 1 , txt , ...
{ ''
  ''
  sprintf(' %2.0f Variable%s:',nvar,app(1:(end*(nvar>1))))
  ' -------------'
  '' } );

s0 = size( char(var(:,1)) , 2 );

for ii = 1 : nvar

    vnr(ii,1) = size(txt,1) + 1;
 
    s2 = size(var{ii,1},2);

    tb = sprintf( '%s%s' , tab(1:end-2) , char(32*ones(1,s0-s2)) );

    var_txt = sprintf('%3.0f %s)%s%s',ii,var{ii,2}(1),tb,var{ii,1});

    if var{ii,3}
       % Number of Dimensions not Zero
       dim_txt = { dim{ var{ii,4}+1 , 1 } };
       var_txt = sprintf( '%s( %s )' , var_txt , strhcat(dim_txt,' , ') );
    end
 
    txt = cat( 1 , txt , { var_txt } );

end

txt = cat( 1 , txt , ...
{ ''
  ' Description:'
  ' ------------'
  '' } );

   s0 = cat(1,att{1:nvar});
   if ~isempty(s0)
       s0 = size( char(s0(:,1)) , 2 );
   else
       s0 = 0;
   end

   s0 = max(s0,6);
   nt = sprintf( '\n %s%s%s' , tab , tab , char(32*ones(1,s0+3)) );

   s1 = size( char(var(:,2)) , 2 ) + size(tab,2);  % Length of TypeString

   l1 = 9;                                         % Length of 1. VarValue

sep = tab(1:end-1);

for ii = 1 : nvar

   vnr(ii,2) = size(txt,1) + 1;

   var_txt = sprintf('%3.0f)%s%s',ii,sep,var{ii,1});

   if var{ii,3}
      % Number of Dimensions not Zero
      dim_txt = { dim{ var{ii,4}+1 , 1 } };
      var_txt = sprintf( '%s( %s )' , var_txt , strhcat(dim_txt,' , ') );
      siz     = cat(2,dim{var{ii,4}+1,2});
   else
      siz     = 1;
   end
 
   txt = cat( 1 , txt , { var_txt ; '' } );

   typ = var{ii,2};

   if ~isempty(var{ii,6}) 

       %-----------------------------------------------------------
       if strcmp(typ,'char') & ischar(var{ii,6})
       %-----------------------------------------------------------

            str       = sprintf('"%s"',var{ii,6});
%%%%%%      str       = strrep( str , char(10) , [ nt ' ' ] );  %%% below for ATT
            att{ii}   = cat( 1 , { 'String' '_' [] str } , att{ii} );
            var{ii,5} = var{ii,5} + 1;

       %-----------------------------------------------------------
       else
       %-----------------------------------------------------------
       % ~CHAR ==> [ First  Last ] Value 
       % 5 zuverlaessige Ziffern

           if var{ii,3}
              len = max(cat(1,dim{var{ii,4}+1,2}));
           else
              len = 1;
           end

           if     any(strcmp(typ,{'float' 'double'}))
              prc = '.5g';
           elseif any(strcmp(typ,{'byte' 'char'}))
              prc = '.3d';
           else
              prc = '.0f';
           end

           if strcmp(typ,'char')
              typ = cat( 2 , typ , ' (UINT8)' );
           end

           bl = s1 - size(typ,2);          % NBlanks after Type

           f1 = sprintf('%s%.0f%s','%',l1+min(bl,0),prc);
           f2 = sprintf('%s%s','%',prc);

           if     len == 1
                  form = f1;
           elseif len == 2
                  form = sprintf('%s    %s',f1,f2);
           else
                  form = sprintf('%s .. %s',f1,f2);
           end

           str = sprintf(form,double(var{ii,6}),[]);

           typ = cat( 2 , typ , ',' , char(32*ones(1,max(bl,0))) , str );

       %-----------------------------------------------------------
       end
       %-----------------------------------------------------------

   end

   att{ii}   = cat( 1 , {'Size' '_' [] sprintf('%.0f ',siz) } , att{ii} );
   att{ii}   = cat( 1 , {'Type' '_' [] typ } , att{ii} );

   var{ii,5} = var{ii,5} + 2;

   for jj = 1 : var{ii,5}

      s2 = size(att{ii}{jj,1},2);

      tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

      str = att{ii}{jj,4};

      str = splitstr(str,splt);

      str = strrep( str , char(10) , ...
                    [ nt char(32*ones(1,strcmp(str(1),'"'))) ] );

      str = sprintf( '%s%s%s:%s  %s' , tab , tb , att{ii}{jj,1}    , ...
                                                  att{ii}{jj,2}(1) , str );

      [m,str] = char2cell(str);

      txt = cat( 1 , txt , str );

  end

  if var{ii,5}
     txt = cat( 1 , txt , {''} );
  end

  txt = cat( 1 , txt , {''} );

end

 txt = strrep(txt,char(0),'\0');


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = splitstr(str,splt,ins);

% SPLITSTR Split String
%
% SPLITSTR( String , Split , Insert )
%
% Split = RowLength + i*MaxLength
%
% RowLength  Maximum Length of Row,
% MaxLength  Maximum Length of String
%
% RowLength < 0  ==> call of RMBLANK
% 
% RowLength + OverSize/RowLength with  OverSize <= RowLength
%
% Insert     Strings, after each a NewLine will insertet 
%

Nin = nargin;

nl  = char(10);  % NewLine  (SplitCharacter)

%-----------------------------------------
% Check Inputs

if isempty(str)
   return
end

if Nin < 2
   splt = [];
end

if Nin < 3
   ins = [];
end

if ~chkstr(str)
    error('Input must be a String.');
end

if isempty(splt)
   splt = 60;
end

nmax = abs(imag(splt(1)));
splt =     real(splt(1));

if ~isempty(ins)
    [ok,ins] = chkcstr(ins);
    if ~ok
        ins = [];
    end
end

%-----------------------------------------
% Check for Quota

n = size(str,2);

is_quot = all( str([1 n]) == '"' );

if nmax == 0
   nmax =  n - 2*is_quot;
end

%-----------------------------------------
% RMBLANK if negative SPLIT

if splt <= 0
   if is_quot
      str = cat(2,str(1),rmblank(str(2:n-1),2),str(n));
   else
      str = rmblank(str,2);
   end
end

%-----------------------------------------
% Check with Length of String

nmax = nmax + sum( str == nl );

if all( n <= [ splt nmax ]+2*is_quot )
   return
end

app = '';

if n > nmax+2*is_quot

   app = ' ...';

   ind = cat( 2 , ( 1 : nmax+is_quot ) ,  n*ones(1,is_quot) );

   str = str(ind);

   n   = nmax + 2*is_quot;

   if splt <= 0
      if is_quot
         str = cat(2,str(1),rmblank(str(2:n-1),2-i),str(n));
      else
         str = rmblank(str,2-i);
      end
   end
   
end

%-----------------------------------------
% Check with Insert

if ~isempty(ins)

    n  = size(str,2);
    i1 = zeros(1,0);

    for cc = ins(:)'
        jj = findstr( str , cc{1} );
        if ~isempty(jj)
            i1 = cat( 2 , i1 , jj+size(cc{1},2)-1 );
        end
    end

    if ~isempty(i1)
        i1 = sort(i1);
       str = insert(str,i1,nl);
    end

end

%-----------------------------------------
% No SPLIT ==> return

splt = abs(splt);

ovs  = splt;
splt = floor(splt);
ovs  = ovs - splt;
ovs  = 1e-10 * round(ovs/1e-10);
ovs  = splt + ceil( ovs * splt );

if splt == 0
   if ~isempty(app)
       na  = size(app,2);
       ind = ( 1 : na ) + n-na-is_quot;
       str(ind) = app;
   end
   return
end

%-----------------------------------------
% Check Length of Segments between NewLine

n  = size(str,2);

i0 = find( str == nl );
lg = diff(cat(2,1,i0+1,n+1));

if ~any( lg > splt )
    if ~isempty(app)
        na  = size(app,2);
        ind = ( 1 : na ) + n-na-is_quot;
        str(ind) = app;
    end
    return
end

%-----------------------------------------
% Check for Characters to Split

cc = { '; '  ', '  '. '   ' ' };

i1 = cat( 2 , 0 , i0 , n+1 );
ok = 0 * i1;

ok(end) = NaN;

for ii = cc

    jj = findstr( str , ii{1} );

    if ~isempty(jj)

        len = size(ii{1},2);

        j0 =  1 + ( jj( 1) == 1 );

        j1 = size(jj,2);
        j1 = j1 - ( jj(j1) >= n-(len-1) );

        jj = jj( j0 : j1 );

        if ~isempty(jj)
     
            % NewLine   before                behind
            kk = ~( ( str(jj-1) == nl ) | ( str(jj+len) == nl ) );
 
            if any(kk)
               if ~all(kk)
                   kk = find(kk);
                   jj = jj(kk);
               end
               jj = jj + len - 1;        % End of Character to Split
               i1 = cat( 2 , i1 , jj );
               ok = cat( 2 , ok , ones(size(jj))+strcmp(ii{1},' ') );
            end

        end

    end

end

if all( ok == 0 )
   if ~isempty(app)
       na  = size(app,2);
       ind = ( 1 : na ) + n-na-is_quot;
       str(ind) = app;
   end
   return
end


[i1,si] = sort(i1);
 ok     = ok(si);

jj = ( diff(i1,1,2) == 0 );   % Same EndIndex !!!
if any(jj)
   jj = find(jj) + 1;
   i1(jj) = [];
   ok(jj) = [];
end

m  = size(i1,2);

l1 = diff(cat(2,i1(1),i1),1,2);    % Length between Characters to Split, incl.

l1 = cumsum(l1,2);

ind = ( 2 : m );

while 1

   % Remove Length-Offset

   n1 = zeros(1,m);
   i0 = find( ok <= 0 );
   n1(i0) = 1;
   n1 = cumsum(n1,2);
   n1 = l1(i0(n1));

   ls = ( l1 - n1 ) .* ~isnan(ok);

   ns = floor(ls/splt);

   jj = ( ( ns(ind-1) == 0 ) & ( ns(ind) > 0 ) );

   if ~any(jj)
       break
   end

   jj = find(jj) + 1;

   % NewLine in Range of OverSize or Split after Blank in Range of OverSize
   kk = ( ( ( ok(jj+1) == 0 ) & ( ls(jj+1) <= ovs ) ) | ...
          ( ( ok(jj+1) == 1 ) & ( ls(jj+1) <= ovs ) & ( ok(jj) == 2 ) ) );

   %%% ovs
   %%% lk = 0*ok; lk(jj) = 1-2*kk; [i1 ; l1 ; ok ; n1 ; ls ; ns ; lk] 

   ok(jj) = -1;

   if any(kk)
        m    = m - sum(kk);
        ind  = ind( 1 : (m-1) );
         kk  = jj(find(kk));
      i1(kk) = [];
      l1(kk) = [];
      ok(kk) = [];
   end

end

ok = ( ok == -1 );
if any(ok)
   ok = find(ok);
   i1 = i1(ok);
   str = insert(str,i1,nl);
end

if ~isempty(app)
    n   = size(str,2);
    na  = size(app,2);
    ind = ( 1 : na ) + n-na-is_quot;
    str(ind) = app;
end


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = write(file,txt)

% WRITE  Write Text to OutFile
 

msg = '';

 if file(1) == '@'
    file = file(2:end);
    mode = 'a';
 else
    mode = 'wt';
 end

 fid  = fopen(file,mode);

  if fid ~= -1
     fprintf( fid , strhcat(txt,char(10)) );
     fclose(fid);
  else
     msg = sprintf('Error open OutFile "%s" for writing.',file);
  end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  fig = list(file,dim,var,txt,vnr,fig)

% LIST  Display Info in ListBoxFigure

Nin = nargin;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

tag = upper(fcn);

new = ( Nin < 6 );
if ~new
    new = ~chkhndl(fig,'figure',tag);
end

if new
   fig = newfig(tag,fcn);
end

set( fig , 'name' , file );

hl = findobj(fig,'type','uicontrol','tag','LIST');

set( hl , 'listboxtop' , 1 , ...
               'value' , 1 , ...
              'string' , txt , ...
            'userdata' , vnr    );

ha = findobj(fig,'type','uimenu','tag','ASSIGN');


ch = get(ha,'children');
fl = strcmp( get(ch,'tag') , 'FLIP' );

hf = ch(find( fl));
ch = ch(find(~fl));

delete(ch);

set(hf,'checked','on');  % !!!

SelectCB = sprintf('%s(''#select'',gcbo);',fcn);
AssignCB = sprintf('%s(''#assign'',gcbo);',fcn);

set( ha , 'userdata' , file );

nn = 0;

sep = { 'off'  'on' };

for ii = 1 : size(var,1)

    h0 = uimenu( 'parent'   , ha , ...
                 'label'    , sprintf(' %s',var{ii,1}) , ...
                'separator' , sep{1+(ii==1)} , ...
                 'callback' , SelectCB , ...
                 'userdata' , var{ii,1}  );

    if     var{ii,3} == 0
       si = [ 1  1 ];
    elseif var{ii,3} == 1
       si = [ dim{ var{ii,4}+1 , 2 }  1  ];
    else
       si = cat( 2 , dim{ var{ii,4}+1 , 2 } );
    end

    si = si(end:-1:1);

    lab = sprintf( '%.0f x ' , si(end:-1:1)  );
    lab = lab(1:end-3);

    h1 = uimenu( 'parent'   , h0  , ...
                 'label'    , lab , ...
                 'userdata' , si  , ...
                 'callback' , AssignCB );

    nn = nn + prod(si);

end

    h0 = uimenu( 'parent'   , ha         , ...
                 'label'    , ' All ... ' , ...
                'separator' , 'on' , ...
                 'callback' , ''   , ...
                 'userdata' , var(:,1)   );

    lab = sprintf( '%.1f MBytes' , nn*8 / 2^20  );

    h1 = uimenu( 'parent'   , h0  , ...
                 'label'    , lab , ...
                 'userdata' , si  , ...
                 'callback' , AssignCB );

drawnow

refresh(fig);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = newfig(tag,fcn);

% NEWFIG   Create New ListBoxFigure


  %--------------------------------------
  scr_uni = get(0,'units');      set(0,'units','pixels')
  scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);
  
  ppi    = get(0,'ScreenPixelsPerInch');
        
  %--------------------------------------

   is_tall = -1 + ( scr_si(4) >=  480 ) + ...
                  ( scr_si(4) >=  600 ) + ...
                  ( scr_si(4) >= 1024 );

   is_win = strcmp( upper(computer) , 'PCWIN' );

   fontsize =  8 + 2 * is_tall - 2 * is_win;

 
   if is_win
      fontname = 'courier';
   else
      fontname = { 'arrial'  get(0,'fixedwidthfontname') };
      fontname = fontname{ 1 + ( scr_si(4) >= 1050 ) } ;
   end

  %--------------------------------------
  
  fs = floor( 25/18 * ppi/100 * fontsize );  % Points --> Pixels

  ww = 80;                   % Character
  hh = 50;
                          
  fig00 = ceil([ 0.55*ww  1.2*hh ] * fs);                         
  fig11 = floor([ 1/2  2/3 ].*scr_si(3:4));
  
  figpos = NaN*ones(1,4);
  
  figpos(3:4)= fig00 + ( fig11 - fig00 ) .* ( fig11 < fig00 );
  
  voffs = max( 60 , min( ceil(1/6*scr_si(4)) , 80 ) );

  figpos(1) = 50;
  figpos(2) = scr_si(4)-voffs-figpos(4);
  

 fig  = figure('position'   , figpos , ...
               'numbertitle', 'off'  , ...
             'integerhandle', 'off'  , ...
               'menubar'    , 'none' , ...
               'toolbar'    , 'none' , ...
               'name'       , tag    , ...
               'tag'        , tag    , ...
               'createfcn'  , ''     , ...
          'handlevisibility','callback' );


  hl = uicontrol( 'parent'    , fig        , ...
                  'style'     , 'listbox'  , ...
            'backgroundcolor' , [1 1 1]    , ...
            'foregroundcolor' , [0 0 0]    , ...
                  'units'     ,'normalized', ...
                  'position'  , [0 0 1 1] , ...
                  'min'       , 0         , ...
                  'max'       , 1         , ...
                  'fontunits' , 'points'  , ...
                  'fontsize'  , fontsize  , ...
                  'fontname'  , fontname  , ...
                  'string'    , ''        , ...
                  'tag'       , 'LIST'    , ...
        'horizontalalignment' , 'left'           );


    ParentCB = 'get(gcbo,''parent'')';

       NewCB = sprintf('%s(''#new'',%s);',fcn,ParentCB);
      FlipCB = sprintf('%s(''#flip'',gcbo);',fcn);
      DownCB = sprintf('%s(''#down'',gcbo);',fcn);
    SelectCB = sprintf('%s(''#select'',gcbo);',fcn);
   RefreshCB = sprintf('refresh(%s);',ParentCB);
     CloseCB = sprintf('delete(%s);',ParentCB);


   set( hl , 'callback' , DownCB );

   hn = uimenu( 'parent'      , fig   , ...
                'label'       , '&New' , ... 
                'tag'         , 'NEW' , ...
                'callback'    , NewCB       );

   ha = uimenu( 'parent'      , fig      , ...
                'label'       , '&Assign' , ...
                'tag'         , 'ASSIGN' , ...
                'callback'    ,  ''            );
 
   hf = uimenu( 'parent'      , ha      , ...
                'label'       , '&Flip Dimensions' , ...
                'tag'         , 'FLIP'  , ...
                'checked'     , 'on'    , ...      % !!!!!!
                'callback'    ,  FlipCB         );
 
   hr = uimenu( 'parent'      , fig      , ...
                'label'       , '&Refresh' , ...
                'tag'         , 'REFRESH' , ...
                'callback'    , RefreshCB     );
 
   hh = uimenu( 'parent'      , fig     , ...
                'label'       , '&Help' , ...
                'tag'         , 'HELP' , ...
                'callback'    , ''           );

        uimenu( 'parent' , hh , ...
                'label'  , '&DoubleClick on VariableNames to' , ...
                'tag'    , 'Help01' , ...
              'callback' , ''           );

        uimenu( 'parent' , hh , ...
                'label'  , ' &switch between List and Description' , ...
                'tag'    , 'Help01' , ...
              'callback' , ''           );

   hb = uimenu( 'parent'      , fig     , ...
                'label'       , '   '   , ...
                'tag'         , 'BLANK' , ...
                'callback'    , ''           );

   hc = uimenu( 'parent'      , fig     , ...
                'label'       , '&Close' , ...
                'tag'         , 'CLOSE' , ...
                'callback'    , CloseCB       );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function new(fig)

% CallBack for New-Menu

   [f,p] = uigetfile('*.nc','Select a NetCDF-File ...');
   if isequal(f,0)
      return
   end

   file = fullfile(p,f);

   [msg,dim,var,att,txt,vnr] = look_cdf(file);

   if ~isempty(msg)
       msg = cat( 2 , 'Error using LOOK_CDF' , char(10) , msg );
       warndlg(msg,'Error','warn');
       return
   end

   list(file,dim,var,txt,vnr,fig);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function select(h0)

% Select  CallBack for Select-Variable-Menu
%
%         Set Label of SubMenu by Size of Label 
%         Takes care by activated FLIP
%

   ha = get(h0,'parent');
   hf = findobj(ha,'type','uimenu','tag','FLIP');

   flip = ~isempty(hf);
   if flip
      flip = strcmp( get(hf,'checked') , 'on' );
   end

   h1 = get(h0,'children');
   si = get(h1,'userdata');

   if flip
      si = si(end:-1:1);
   end

    lab = sprintf( '%.0f x ' , si );
    lab = lab(1:end-3);

    set( h1 , 'label' , lab );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function assign(h1)

% ASSIGN  CallBack for Assign-Menu

   h0 = get( h1 , 'parent' );
   ha = get( h0 , 'parent' );

   file = get( ha , 'userdata' );
   name = get( h0 , 'userdata' );

   hf = findobj(ha,'type','uimenu','tag','FLIP');

   flip = ~isempty(hf);
   if flip
      flip = strcmp( get(hf,'checked') , 'on' );
   end

   [ok,name] = chkcstr(name,0);

   if ~ok
       msg = 'Invalid UserData of UIMenu.';
   else
       try
          [msg,d,v] = read_cdf(file,'var',name,'fill',NaN);
          if isempty(msg)
             for ii = 1 : size(v,1)
                 if flip
                    v{ii,7} = permute(v{ii,7},(ndims(v{ii,7}):-1:1));
                 end
                 assignin('base',strrep(v{ii,1},' ','_'),v{ii,7});
             end
          end
       catch
          msg = lasterr;
       end 
   end
  
   if ~isempty(msg)
       msg = cat( 2 , 'Error Assign Variable' , char(10) , msg );
       warndlg(msg,'Error','warn');
       return
   end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function flip(hf)

% FLIP  CallBack for Flip-Menu

sets = { 'on'  'off' };

sets = sets{ 1 + strcmp(get(hf,'checked'),'on') };

set( hf , 'checked' , sets );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function down(hl)

sel = get(get(hl,'parent'),'selectiontype');

if ~strcmp(sel,'open');
    return
end

vnr = get(hl,'userdata');

sv = size(vnr);
pv = prod(sv);

ok = ( isnumeric(vnr) & ~isempty(vnr) );
if ok
   sv = size(vnr);
   pv = prod(sv);
   ok = ( ( sv(2) == 2 ) & ( sv(1)*sv(2) == pv ) );
   if ok
       v = vnr(:);
      ok = all( ( v > 0 ) & ( mod(v,1) == 0 ) );
      if ok
         v  = max(v(:));
         ok = ( v <= size(get(hl,'string'),1) );
      end
   end
end

if ~ok
    return
end
 
cc = ( vnr == get(hl,'value') );

ok = any(cc,1);

if ~any(ok) | all(ok)
    return
end

jj = 1 + ok(2);

cc = find( cc(:,jj) );

val = vnr(cc,3-jj);

if jj == 1
   top = val - 1;
else
   top = vnr(1,1) - 3;
end

set( hl , 'value' , val , 'listboxtop' , top );

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = att2str(name,typ,val)

%  STRING = ATT2STR( DATATYPENUMBER , VALUE )
%
%  Converts AttributeValues to String
%
%
%  DataTypeNumber:    6     |    5    |   4    |    3    |   2    |    1
%  DataType      : 'double' | 'float' | 'long' | 'short' | 'char' | 'byte'
%

str = '';

switch typ

 %-----------------------------------------------------------
 case 'byte'

  ok = ( ~any( strcmp( name , {'missing_value' '_FillValue'} ) ) & ...
          ( size(val,1) == 1 )  &  ( size(val,2) > 3 ) );
  if ok
     ok = ( all( val >= 0 ) & any( val >= 28 ) );
     if ok
        ok = all( ( val ==  0 ) | ( val ==  9 ) | ( val == 10 ) | ...
                  ( val == 13 ) | ( val >= 28 ) );
     end
  end

  if ok
      val = checkchar(val,1);
      str = sprintf('"%s"',val);
  else
      form = '%3.3d  ';
      form = permute(form(ones(1,8),:),[2 1]);
      form = cat(2,permute(form(:),[2 1]),'\n');
      str = sprintf(form,val);
      str = str( 1 : ( end - ( mod(prod(size(val)),8) == 0 ) ) );
  end

 %-----------------------------------------------------------
 case  'char'

  if ~ischar(val)         % UINT8
      form = '%3.3d  ';
      form = permute(form(ones(1,8),:),[2 1]);
      form = cat(2,permute(form(:),[2 1]),'\n');
      str = sprintf(form,double(val));
      str = str( 1 : ( end - ( mod(prod(size(val)),8) == 0 ) ) );
  else
      
      if isequal(val,0)
         val = '\0';
      end

      val = checkchar(val,~any(strcmp(name,{'missing_value' '_FillValue'})));

      str = sprintf('"%s"',val);

  end

 %-----------------------------------------------------------
 case  { 'short'  'long' }

  str = int2str(val);

 %-----------------------------------------------------------
 case  { 'float'  'double' }

  str = num2str(val);

end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = checkchar(str,mode)

% CHECKCHAR  Checks for valid Characters

if nargin < 2
   mode = 1;
end

  str = double(str);

  bad = find( ~( ( str ==  9 ) |  ...
                 ( str == 10 ) |  ...
                 ( str == 13 ) |  ...
                 (  28 <= str  &   str <= 126 ) | ...
                 ( 160 <= str  &   str <= 255 )        ) );

str(bad) = 32;


if mode
  str = rmblank(str,2);
else
  str = char(str);
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,ind,ix] = insert(x,ind,y,dim)

% INSERT Inserts Values into Array
%
% [ Z , IY , IX ] = INSERT( X , Index , Y , DIM )
%
% inserts Y into X at the Index in the Dimension DIM.
%
% defaults: DIM = 1, for Vectors X the first non-singleton Dimension
%
% Index >  0  ==>  Y will inserted behind X(Index)
% Index <  0  ==>  Y will inserted before X(Index)
% Index == 0  ==>  Y will inserted before X(1)
%
% Make sure, that Index match the Size of X in DIM.
% 
% Y can be Empty, a Scalar, Vector or Matrice.
%
%  Make sure, that the Size of Y match the Size of X and 
%   the Length of Index in DIM.
%
% If Y is empty or not given, X(Index) will inserted (duplicated).
%
%--------------------------------------------------------------
% The Output-IndexVectors IY and IX refers to the inserted Y 
%   and the original X.
%
% IY = [ N by 2 ];  with ( N <= length(Index) )  
%
%   The 1. Column of IY is the Index in Z of the Inserted  Y
%   The 2. Column of IY is the Index in Y of the Inserted  Y
% 
% IX = [ size(X,DIM) by 1 ] is the Index in Z of the Original X
%
%--------------------------------------------------------------
% for Vectors X and Vectors or Scalars Y
%
%  Z( IY(:,1) ) == Y( IY(:,2) ) 
%
%  Z(IX)        == X
% 
%--------------------------------------------------------------
% for Matrice X and Y
%
%  Z( : , ... , IY(:,1) , ... , : ) == Y( : , ... , IY(:,2) , ... , : ) 
%
%  Z( : , ... , IX , ... , : )      == X
%


Nin = nargin;

Nout = nargout;

ix = [];

if Nin < 2
   error('Not enough Input Arguments.');
end

if Nin < 3
   y = [];
end

if Nin < 4
   dim = [];
end

%--------------------------------------------------------
% Check Class of X and Y

if ~strcmp(class(x),class(y))
    error('X and Y must be of same Class.');
end

%--------------------------------------------------------
% Check Size of X and Index

sx = size(x);

if isempty(dim)
   [ms,dim] = max(sx);
   if ~( ms == prod( sx + ( sx == 0 ) ) )  % Not a Vector
      dim = 1;
   end
else
   sx = cat( 2 , sx , ones(1,dim-size(sx,2)) );
end

if all( sx == 0 )
   sx      = ones( 1 , dim+(dim==1) );
   sx(dim) = 0;
   x       = zeros(sx);
end

ind = floor(ind(:));
 ni =  size(ind,1);

ind = abs(ind) - 1 * ( ind < 0 );  % Negative Inserts before !!!

if ~all( ( 0 <= ind ) & ( ind <= sx(dim) ) );
    error(sprintf('Index exceeds Matrix Dimension %.0f.',dim));
end

%--------------------------------------------------------
% Check Size of Y

sy = size(y);

not_y =  all( sy == 0 );   % ZERO-Size
one_y =  all( sy == 1 );   % Scalar

vec_y = ( max(sy) == prod( sy + ( sy == 0 ) ) ); % Scalar | Vector
vec_y = ( vec_y & ~one_y );                      % Vector

mat_y = ~( not_y | one_y );                      % Vector | Matrice

%--------------------------------------------------------
% Vector ==> reshape into Dimension DIM

if vec_y 

   ss = cat( 2 , sy , ones(1,size(sx,2)-size(sy,2)) );
   ss(dim) = sx(dim);

   if ~isequal(ss,sx)

       [ny,dy] = max(sy);

       if ~( dy == dim )

           p      = ( 1 : size(sy,2) );
           p(dim) = dy;
           p(dy ) = dim;

                y = permute( y , p );

               sy = size(y);

       end

    end

end

%--------------------------------------------------------
% Check Size of Y with Size of X and Length of Index at DIM

if mat_y

    sy = cat( 2 , sy , ones(1,dim-size(sy,2)) );

    s      = sx;
    s(dim) = ni;

    if ~isequal( s , sy )
       error(sprintf('Size of Y must match Size of X and Length of Index in Dimension %.0f.',dim));
    end

end


%--------------------------------------------------------
% Permute and Reshape X and Y to 1. Dimension at DIM

ns = size(sx,2);

perm = cat( 2 , dim , ( 1 : dim-1 ) , ( dim+1 : ns ) );

x = reshape( permute(x,perm) , sx(perm(1)) , prod(sx(perm(2:ns))) );

if mat_y
   y = reshape( permute(y,perm) , sy(perm(1)) , prod(sy(perm(2:ns))) );
end

%--------------------------------------------------------
% Remove Multiple Indize

[ind,si] = sort(ind,1);

    bad  = find( diff(ind,[],1) == 0 );

ind(bad) = [];

ni = size(ind,1);

if mat_y
     y        = y(si,:);
     y(bad,:) = [];
    si(bad)   = [];
else
    si        = ones(ni,1);
    if not_y
       si = NaN * si;
    end
end

%--------------------------------------------------------
% Build IndexVector to Insert

ii = ones( sx(dim)+ni , 1 );

ind = ind + ( 1 : ni )';

ii(ind) = 0;

if Nout > 2
   ix = find(ii);
end

if ~isempty(ii)

  ii    = cumsum(ii,1);

  ii(1) = ii(1) + ( ii(1) == 0 );   % Check for ZERO-Index

end

% Insert X (duplicate)
if ~isempty(x)
   x = x(ii,:);
end

% Insert Y
if ~not_y
   x(ind,:) = y;
end

if Nout > 1
   ind = cat( 2 , ind , si );
end

%--------------------------------------------------------
% Reshape and Permute back

sx(dim) = sx(dim) + ni;

x = reshape( x , sx(perm) );

perm = cat( 2 , ( 1 : dim-1 ) + 1 , 1 , ( dim+1 : ns ) );

x = permute( x , perm );

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
    if ~all( ( real(dim) == 1 ) |  ( real(dim) == 2 ) )
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
 
    d = real(d);

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,bb,ref] = char2cell(bb,mnl,ref01);

% CHAR2CELL Converts CharArray to CellStringArray
%
% [Msg,CellString] = CHAR2CELL( String )
%
% Converts the CharacterArray String into a CellStringArry.
%
%---------------------------------------------------------
% CHAR2CELL( String , MarkerNewLine )  
%
% Replace the String's defined by MarkerNewLine with NewLine,
%   default: EMPTY
%
%---------------------------------------------------------
% [Msg,CellString,Reference] = ...
%    CHAR2CELL( String , MarkerNewLine , ReferenceCell)
%  
% Returns References in String, using GET_REF, 
%   ReferenceCell = { S1 S2 C1 C2 }, type: >> help get_ref
%


 Msg = '';
 
 Msg0 = 'CHAR2CELL: ';

 ref = cell(0,2);

 Nin  = nargin;
 Nout = nargout;
 
 nl = char(10);

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  % Marker for NewLine
  if Nin < 2
    mnl = '';
  end

%*****************************************************
% Check Inputs

  if iscellstr(bb)
     bb = cellstr(char(bb));
     bb = strhcat(bb,'',1);
  end

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  ok = ( isnumeric(bb) | ischar(bb) );
  if ok

    if ischar(bb)
       bb = double(bb);
    end

    ok = all( ( mod(bb,1) == 0 )  & ( bb >= 0 ) & isfinite(bb)  );

  end

  if ~ok

      Msg = [ Msg0 ...
              'Input String must be a CharacterArray or ASCII-Codes.'];
 
      bb = cell(0,1);

      return

  end

  %---------------------------------------------------
  % Check MarkerNewLine

  if ~( isempty(mnl) |  ...
        ( ischar(mnl) &  ( prod(size(mnl)) == size(mnl,2) ) ) )

    Msg = [ Msg0  'Input MarkerNewLine must be a String.'];

    return

  end

  %---------------------------------------------------
  % Check Reference

  if ( Nin == 3 )  &  ( Nout == 3 )

     [Msg,ref] = get_ref('',ref01{:});

    if ~isempty(Msg)
       Msg = [ Msg0  'Invalid Input for Reference.' nl Msg ];
       return
    end

  end

%*****************************************************

  if ( size(bb,1) > 1 )  &  ( size(bb,2) > 1 )
     bb = cat( 2 , bb , 10*ones(size(bb,1),1) );
     bb = permute( bb , [ 2 1 ] );
     bb = bb(:);
  end

  if ( size(bb,1) > 1 ) 
     bb = permute( bb , [ 2 1 ] );
  end

  %---------------------------------------------------
  % Check Characters

  ok = all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

  if ~ok
    Msg = [Msg0 'Invalid Characters in String.' ];
    return
  end


%*****************************************************
 

  %---------------------------------------------------
  % Remove CR
  bb( find( bb == 13 ) ) = [];


  bb = char(bb);


  %---------------------------------------------------
  % TAB  --> 8 Blanks 
  bb = strrep( bb , char(9) , char(32*ones(1,8)) ); 

  %---------------------------------------------------
  % mnl --> NewLine   % !!!!!
  if ~isempty(mnl)
    bb = strrep( bb , mnl , char(10) ); 
  end

  %---------------------------------------------------
  % Reference
  if ( Nin == 3 )  &  ( Nout == 3 )

     [MsgR,ref] = get_ref(bb,ref01{:});

  end

  %---------------------------------------------------
  % Form CellString

  % 1. "'"     --> "''"
  bb = strrep( bb , char(39) , char([39  39]) ); 

  % 2. NL --> "';'"
  bb = strrep( bb , char(10) , char([39  59  39]) );
  

  bb = eval([  '{'''  bb  '''}' ]);

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,typ,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%


Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

is_typ = chkstr(typ,0);

if is_typ

  ok = isempty(typ);
  if ~ok
      ok = strcmp( get(h,'type') , typ  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = typ;
  
end


%-------------------------------------------------------------------------
% Check Tag

[ok,tag] = chkcstr(tag,0);

if ~ok
    msg = 'Input Tag must be a CharArray or CellStringArray.';
    return
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = strcmp(t,tag);


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

   
