function [ msg , v ,a ] = assign_hdf(file,varargin)

% ASSIGN_HDF   assign Data of a HDF-File in Workspace of Caller
%
% ASSIGN_HDF( FILENAME , VarNames , ... , { DimName DimVal } , ... )
%
% VarNames can be CharacterArrays or CellArrays of Strings 
%  for Variables to assign.
%
% The 2-Column CellArrays { DimName DimVal } defines special
%  Values for the Dimensions to extract.
%
%  DimName is a String for DimensionName, DimVal a (complex) Integer:
%
%   DimVal = DimNumber +      0*i
%   DimVal = DimStart  + Stride*i
%   DimVal = DimEnd    - Stride*i
%   
%    Note: the first Value is equal to 1
%    ( real(DimVal) <= 0 )  ==> measured from End
%   
% optional Outputs:
%
%   MSG    = ASSIGN_HDF( ... )  returns the ErrorMessages.
%
% [MSG,V,A] = ASSIGN_CDF( ... )  doesn''t assign the Variables 
%
%  in the Workspace of caller, the Values and the Attributes are
%  returned in the Structures V and A.
%
% [MSG,A,D] = ASSIGN_CDF( FILENAME , {} )  inquires for all Variables
%  
%  the Attributes and DimensionInfo are returned in the Structures A and D.
%
%
%  see also: READ_HDF, LOOK_HDF
%

Nin  = nargin;
Nout = nargout;

msg = '';
  v = [];
  a = [];

nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));


if Nin < 1
   msg = sprintf('%sNot enough InputArguments.',msg0);
   if Nout == 0
      fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
      clear msg
   end
   return
end

%------------------------------------------------------------
% Check Inputs

Nin = Nin - 1;

inq = ( Nin == 1 );
if inq
   inq = isempty(varargin{1});  % Inquire Only
end

ReadIn = cell(0,1);    % Inputs for READ_HDF

%------------------------------------------------------------
% Check for Inputs for READ_HDF

if ~inq & ~isempty(varargin)

    for ii = 1 : Nin

        val = varargin{ii};
    
        [ok(ii),str] = chkcstr(val,0);
        if ok(ii)
           varargin{ii} = str(:);
        elseif iscell(val) & ~isempty(val)
           s = size(val);
           n = prod(s);
           if ( s(2) == 2 ) & ( n == s(1)*s(2) )
               ok(ii) = 2 * chkcstr(val(:,1));
           end
        else
           ok(ii) = 3 * ( isnumeric(val) & ( prod(size(val)) == 1 ) );
        end
    
    end

    if any(~ok)
       msg = 'Inputs must be Character- or CellStringArrays for VariablesNames or DimensionDefinitions.';
       msg = sprintf('%sInvalid Inputs.\n%s',msg0,msg);
       if Nout == 0
          error(msg);
          clear msg
       end
       return
    end

    %------------------------------------------------------------
    % Inputs for READ_HDF

    ReadIn = cell(0,1);

    if any( ok == 2 )
       ReadIn = cat( 1 , ReadIn , { 'dim' ; cat( 1 , varargin{find(ok==2)} ) } );
    end

    if any( ok == 3 )
       ReadIn = cat( 1 , ReadIn , { 'fill' ; varargin{max(find(ok==3))} } );
    end

    if any( ok == 1 )

       var = cat( 1 , varargin{find(ok==1)} );
   
       %------------------------------------------------------------
       % Remove Duplicate VariableNames

       [str,si] = sort(var);

       str = double(char(str));
        jj = find( sum( diff(str,1,1) == 0 , 2 ) == size(str,2) );
   
       if ~isempty(jj)
           var(si(jj)) = []; 
       end
   
       ReadIn = cat( 1 , ReadIn , { 'var' ; var } );
   
    end

end

%------------------------------------------------------------
% Call READ_HDF or LOOK_HDF

if inq
   fcn = 'look_hdf';
else
   fcn = 'read_hdf';
end

[msg,d,v,a] = feval(fcn,file,ReadIn{:});

if ~isempty(msg)
    msg = sprintf('%sError call %s.\n%s',msg0,upper(fcn),msg);
end

%------------------------------------------------------------
% Check for Valid Output

if isempty(v) | ( size(v,2) < 14-4*inq )
   if isempty(msg)
      msg = sprintf('%sInvalid output from %s.',msg0,upper(fcn));
   end
   if Nout == 0
      error(msg);
      clear msg
   elseif Nout >= 2
      v = [];
      a = [];
   end
   return
end

%------------------------------------------------------------
% Check VariableNames for strange Characters

v(:,1) = chkname(v(:,1));

nv = size(v,1);

%------------------------------------------------------------
% Check, Prepend Attributes

a = cat(1,v(:,10),{a});      % { VarAtt ; GlbAtt }

if inq | ( Nout == 3 )

   % DimName --> DimIndex 
   for ii = 1 : nv
       vd = v{ii,9};
       nd = size(vd,1);
       di = zeros(1,nd);
       for jj = 1 : nd
           di(jj) = find(strcmp(vd{jj,1},d(:,1))) - 1;
       end
       v{ii,9} = di;
   end

   d(:,1) = chkname(d(:,1));  % Check DimNames first

   na = size(a,1);

   if     na < nv
      a = cat( 1 , a , cell(nv-na,1) );
   elseif na > nv
      na = nv + 1 - ( isempty(a{nv+1}) | ( Nout == 1 ) );
       a = a(1:na);
   end

   for ii = 1 : nv

       if isempty(a{ii})
          a{ii} = cell(0,2);
       else
          a{ii} = a{ii}(:,[1 4]);  %  { Name  Value }
       end

       sz = 0;
       dn = {};

       di = v{ii,9};
       if isempty(di) 
          if isequal(v{ii,4},0)
             sz = 1;            % Single without Dimension
          end
       else
          di = di(end:-1:1) + 1;
          if inq
             sz = cat(2,d{di,2});
          else
             sz = size(v{ii,14});
          end
          dn = d(di,1);
          dn = dn(:)';
       end

       a{ii} = cat( 1 , { 'Type' v{ii,3} } , ...
                        { 'Size' sz      } , ...
                        { 'Dims' {dn}    } , ...
                        a{ii}                    );

   end

   if na == nv+1
       a{na} = a{na}(:,[1 4]);  %  Global: { Name  Value }
   end
   
   % Check AttributeNames, build Structure

   for ii = 1 : na
       a{ii}(:,1) = chkname(a{ii}(:,1));
       a{ii} = permute(a{ii},[2 1]);
       a{ii} =  struct(a{ii}{:});
   end

   % Structure

   iv = min( ( 1 : na ) , nv );

   a = cat( 2 , v(iv,1) , a(1:na) );
   if ( na == nv+1 )
      a{na,1} = 'Global';
      a(:,1)  = chkname(a(:,1),2);
   end

   if inq
      v = a;
      a = d(:,[1 2]);
   end

end

%------------------------------------------------------------
% StructureOutput

if ~inq
    v = v(:,[1 14]);
end

if Nout > 1
   v = permute(v,[2 1]);
   v =  struct(v{:});
   if Nout > 2
      a = permute(a,[2 1]);
      a =  struct(a{:});
   end      
   return
end

%------------------------------------------------------------
% Assign Variables in Workspace of Caller

for ii = 1 : size(v,1)
    if ~isempty(v{ii,1})
        assignin( 'caller' , v{ii,1} , v{ii,2} );
    end
end

if Nout == 0
   clear msg
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function v = chkname(v,mode)

[ok,v] = chkcstr(v,0);

if ~ok
    return
end

n = prod(size(v));

if nargin < 2
   mode = 0;    % Both
end

%---------------------------------------------------------
if ~( mode == 2 )
%---------------------------------------------------------

% Replace invalid Characters by "_"

for ii = 1 : n
    w  = v{ii};
    jj = ~( ( ( '0' <= w ) & ( w <= '9' ) ) | ...
            ( ( 'A' <= w ) & ( w <= 'Z' ) ) | ...
            ( ( 'a' <= w ) & ( w <= 'z' ) )       );
    if any(jj)
         jj  = find(jj);
       w(jj) = ' ';
       w = rmblank(w,2);
       w =  strrep(w,' ','_');
    end
    v{ii} = w;
end

%---------------------------------------------------------
end
%---------------------------------------------------------

%---------------------------------------------------------
if ~( mode == 1 )
%---------------------------------------------------------

% Append duplicate Names with '_'

for ii = n : -1 : 1
    jj = cat(2,(1:(ii-1)),((ii+1):n));
    while any(strcmp(v{ii},v(jj)))
          v{ii} = cat(2,v{ii},'_');
    end
end

%---------------------------------------------------------
end
%---------------------------------------------------------

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
