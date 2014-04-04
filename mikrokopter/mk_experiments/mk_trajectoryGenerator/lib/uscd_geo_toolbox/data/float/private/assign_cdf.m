function [ msg , v , a ] = assign_cdf(file,varargin)

% ASSIGN_CDF   assign Data of a NetCDF-File in Workspace of Caller
%
% ASSIGN_CDF( FILENAME , VarNames , ... , { DimName DimVal } , ... )
%
% VarNames can be CharacterArrays or CellArrays of Strings 
%  for Variables to assign.
%
% The 2-Column CellArrays { DimName DimVal } defines special
%  Values for the Dimensions to extract.
%
%  DimName is a String for DimensionName, DimVal a (complex) IntegerVector:
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
% ASSIGN_CDF( FILENAME , ... , FillValue , ... )
%
%  defines the FillValue to use in READ_CDF
%
%
% optional Outputs:
%
% MSG = ASSIGN_CDF( ... )  
%
%   returns the ErrorMessages in MSG, the VariableValues are assigned.
%
%
% [MSG,V,A] = ASSIGN_CDF( ... )  doesn''t assign the Variables,
%
%  the Values and the Attributes or the requested Variables are 
%  returned in the Structures V and A.
%
% [MSG,A,D] = ASSIGN_CDF( FILENAME , {} )  inquires for all Variables
%  
%  the Attributes and DimensionInfo are returned in the Structures A and D.
%
%
%  see also: READ_CDF, LOOK_CDF
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
      error(msg);
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

ReadIn = cell(0,1);    % Inputs for READ_CDF

%------------------------------------------------------------
% Check for Inputs for READ_CDF

if ~inq & ~isempty(varargin)

    ok = double(zeros(Nin,1));

    for ii = 1 : Nin

        val = varargin{ii};
    
        [ok(ii),str] = chkcstr(val,0);
        if ok(ii)
           varargin{ii} = str(:);
        elseif iscell(val) & ~isempty(val)
           s = size(val);
           n = prod(s);
           if ( s(2) == 2 ) & ( n == s(1)*s(2) )
               ok(ii) = 2 * double(chkcstr(val(:,1)));
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
    % Inputs for READ_CDF

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
% Call READ_CDF or LOOK_CDF

if inq
   fcn = 'look_cdf';
else
   fcn = 'read_cdf';
end

[msg,d,v,a] = feval(fcn,file,ReadIn{:});

if ~isempty(msg)
    msg = sprintf('%sError call %s.\n%s',msg0,upper(fcn),msg);
end

%------------------------------------------------------------
% Check for Valid Output

if isempty(v) | ( size(v,2) < 7-inq )
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

if inq | ( Nout == 3 )

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
          a{ii} = a{ii}(:,[1 3]);  %  { Name  Value }
       end

       sz = 0;
       dn = {};

       di = v{ii,4};
       if isempty(di) 
          if isequal(v{ii,3},0)
             sz = 1;            % Single without Dimension
          end
       else
          di = di(end:-1:1) + 1;
          if inq
             sz = cat(2,d{di,2});
          else
             sz = size(v{ii,7});
          end
          dn = d(di,1);
          dn = dn(:)';
       end

       a{ii} = cat( 1 , { 'Type' v{ii,2} } , ...
                        { 'Size' sz      } , ...
                        { 'Dims' {dn}    } , ...
                        a{ii}                    );

   end

   if na == nv+1
       a{na} = a{na}(:,[1 3]);  %  Global: { Name  Value }
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
    v = v(:,[1 7]);
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
%%% size(v)

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
%    A positive complex Part of DIM, to removes Blanks only from Start,
%    A negative complex Part of DIM, to removes Blanks only from End.
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
