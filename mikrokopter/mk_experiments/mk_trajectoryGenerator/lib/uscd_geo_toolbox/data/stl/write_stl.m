function [z,kk] = write_stl(file,x,y,z,varargin);

% WRITE_STL   Writes gridded Data into STL-Format (triangulated)
%
% [XYZ,IND] = WRITE_STL( File , X , Y , Z , Format , Mode )
% [XYZ,IND] = WRITE_STL( File , V , IND , Format , Mode )
%
% Input for Grid:
%
% X, Y, Z  define a Grid, X and Y can be Vectors matching the size of Z
%
% V are Vertices and IND the Index for triangular Faces in V
%
% V   = [ NPoint by 3 ]  [ X Y Z ]
% IND = [ 3 by NFaces ]
%
% Output for triangulated Data:
%
% XYZ = [ 3 by NFaces by XYZ ]
% IND = [ 3 by NFaces ]        linear Index for triangular Faces in Z
%
% FID = WRITE_STL( File )                             Open   File
% XYZ = WRITE_STL( FID , X , Y , Z , Format , Mode )  Append to File
%       WRITE_STL( FID  )                             Close  File
%
% Format  String for NumberFormat, using FPRINTF
%         '%XYZ' | '%XY %Z' | '%X %Y %Z'
%
% Mode defines the Triangulation by GRID2TRI
%
% [VX,VY] = WRITE_STL( File , X , Y , Z , ... , 'check' )
%
% Checks Inputs only, doesn't writes STL-File,
% Returns VX and VY which are True of Vectors X and Y
%
% see also: GRID2TRI, READ_STL
%

Nout = nargout;

Nin  = nargin;

vin = varargin;

if Nin == 0
   error('Not enough InputArguments.');
end

fid = [];

nl = char([13 10]);

frm  = '%13.6e';
mode = 1;

%*******************************************************
% Check File

ok = ( ischar(file) & ~isempty(file) & ...
       ( prod(size(file)) == size(file,2) ) );
if ~ok
    ok = ( isnumeric(file) & ( prod(size(file)) == 1 ) );
    if ok
       ok = ~isempty(fopen(file));
       if ok
          fid = file;
       end
    end
end

if ~ok
    error('First Input must be a FileName or FileIdentifer.');
end

%*******************************************************
% Check Size of MatriceInputs

tri = 0;  % Check for already triangulated Input

if Nin > 1

   if Nin >= 3
      tri = ( isnumeric(x) & isnumeric(y) );
      if tri
         tri = ( ( size(x,2) == 3 ) & ( size(y,1) == 3 ) );
         if tri
            tri = all( ( mod(y(:),1) == 0 ) & ...
                       ( y(:) >= 1 ) & ( y(:) <= size(x,1) ) );
         end
      end
   end

end

if tri

   if Nin > 3
      vin = cat( 2 , vin , z );
   end

elseif ( Nin > 1 )
 
   if Nin < 4
      error('X, Y AND Z or valid V & IND are required.');
   end

   if ~( isnumeric(x) & isnumeric(y) & isnumeric(z) )
      error('X, Y AND Z must be numeric.');
   end

   sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
   sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );
   sz = size(z);  pz = prod(sz);

   if vx & ~( px == sx(2) )
      x = x(:)';
   end

   if vy & ~( py == sy(1) )
      y = y(:);
   end

   if vx & vy
      ok =   isequal( sz , [ py  px ] );
   elseif vy
      ok = ( isequal( sz , sx ) & ( sz(1) == py ) );
   elseif vx
      ok = ( isequal( sz , sy ) & ( sz(2) == px ) );
   else
      ok = ( isequal( sz , sx ) & isequal( sz , sy ) );
   end

   if ~ok
       error('Matrix Dimension must be agree');
   elseif ~( ndims(z) == 2 )
       error('Matrix must be 2-dimensional.');
   elseif ~all( sz >= 2 )
       error('Single Matrix Dimensions.');
   end

end

%*******************************************************
% Check for Format and Mode

chk = 0;

if Nin > 4-tri

   chk = isequal(vin{end},'check');

   if chk
      Nin = Nin - 1;
      vin = vin(1:(end-1));
   end

end

if Nin > 4-tri

   for v = vin(:)'
       if ~isempty(v{1})
           if ischar(v{1})
              frm  = v{1};
           elseif isnumeric(v{1})
              mode = v{1};
           else
              error('Invalid Input for Format or Mode.');
           end
       end
   end

end

fi = '%';

ok = ( ischar(frm) & ( size(frm,2) > 1 ) & ...
       ( prod(size(frm)) == size(frm,2) ) );

if ok
   np = sum( frm == fi );
   ok = any( np == [ 1 2 3 ] );
   if ok
      frm = separate(frm,fi);
      ok = ( ~any(strcmp(frm,'')) & ( size(frm,2) == np ) );
      if ok
         if     np == 1
            frm = frm([1 1 1]);
         elseif np == 2
            frm = frm([1 1 2]);
         end
         frm = frm([1 1],:);
         frm(1,:) = {fi};
         frm = sprintf(' %s%s',frm{:});
         d = [ -1  0  1 ];
         c = eval(['[' sprintf(frm,d) ']'],'0');
         ok = isequal(size(d),size(c));
         if ok
            ok = all( sign(d) == sign(c) );
         end
      end  

   end
         
end

if ~ok
    error('Format must be a String for using FPRINTF, returning a 3-Element RowVector.');
end

if ~( isnumeric(mode) & ( prod(size(mode)) == 1 ) )
    error('Mode must be a single Numeric.');
end

%*******************************************************

if chk
   z = vx; kk = vy; return
end

%*******************************************************

if ischar(file)

   fid = fopen(file,'wt');
   if fid == -1
      error('Cann''t open File.');
   end

   fprintf(fid,'solid%s',nl);

elseif Nin == 1

   fprintf(fid,'endsolid%s',nl);
   fclose(fid);

   fid = [];

end
 
if Nin == 1

   z  = fid;
   kk = [];

   if Nout == 0
      clear z
   end

   return

end

%*******************************************************
% Triangulation

if tri

    kk = y;

    z  = x(:,3);
    y  = x(:,2);
    x  = x(:,1);

    vx = 0;        % Don't forget, See below !!!
    vy = 0;        % Don't forget, See below !!!

elseif vx | vy

   [ii,jj,kk] = grid2tri(sz,mode);

else

   kk = grid2tri(sz,mode);

end

%*******************************************************
% Format for FACE

str = sprintf(' %.0f',[0 0 0]);

frm = sprintf( [     ...
'  facet normal%s%s' ...
'    outer loop%s'   ...
'      vertex%s%s'   ...
'      vertex%s%s'   ...
'      vertex%s%s'   ...
'    endloop%s'      ...
'  endfacet%s'           ] , str,nl,nl,frm,nl,frm,nl,frm,nl,nl,nl);

%*******************************************************
% Matrices of FaceVertices

if vx & vy
   z = cat(3,x(jj),y(ii),z(kk));
elseif vx
   z = cat(3,x(jj),y(kk),z(kk));
elseif vy
   z = cat(3,x(kk),y(ii),z(kk));
else
   z = cat(3,x(kk),y(kk),z(kk));
end

%*******************************************************
% Write Faces, check for NaN-Vertices
%
% [ 3 by NFaces by XYZ ] --> [ XYZ by 3 by NFaces ]
%

if any(isnan(z(:)))

   ok = find( sum(sum(isnan(z),1),3) == 0 );

   fprintf(fid,frm,permute(z(:,ok,:),[3 1 2]));
   
else

   fprintf(fid,frm,permute(z,[3 1 2]));

end

%*******************************************************

if ischar(file)
   fprintf(fid,'endsolid%s',nl);
   fclose(fid);
end

%*******************************************************

if Nout == 0

   clear z

end

%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ii,jj,kk] = grid2tri(si,mode)

% GRID2TRI  Triangulation of Grid
%
% [ I1 , I2 , IND ] = GRID2TRI( SIZ , Mode )
%
%             IND   = GRID2TRI( SIZ , Mode )
%
%  SIZ defines the Size of an 2-dimensional Grid
%  SIZ = [ NRow  NColumn ]
%
% A Grid defined by SIZ has PROD( SIZ - 1 ) 
%  rectangular Faces. Each of them is devided into 2
%   triangular Faces. The Indices of the Vertices
%  of the triangular Faces are returned.
%
%  The Outputs have 3 Rows for the Vertices of the Triangular
%   and 2 * PROD( SIZ - 1 ) Columns for the Faces.
%
%  I1 and I2 are Subscripts for RowIndex and ColumnIndex
%   of the Vertices.
%
%  IND is the linear Index of the Vertices.
%  IND = I1 + ( I2 - 1 ) * SIZ(1)
%
%------------------------------------------------------------
%
% Mode is a single (imaginary) Number to define the Order
%  and the Shape of the Triangles, following the Columns
%  and switching at each Row.
%
% Mode = Swap + Shape * i
%
%  Swap defines the Mode to change the Direction of the Columns:
%
%  Swap = 0  no switch
%         1  switch at the even     Rows
%        -1  swicth at the oddinary Rows
%
%         0: [1,1] .. [1,N] , [2,1] .. [2,N] , ... , [M,1] .. [2,N]
%         1: [1,1] .. [1,N] , [2,N] .. [2,1] , ...
%        -1: [1,N] .. [1,1] , [2,1] .. [2,N] , ...
% 
%        where [ M  N ] is the Size of the rectangular Faces:
%              [ M  N ] = SIZ - 1;  M = NRow - 1;  N = NColumn - 1;
%
%  Shape defined the Triangles
%
%  Shape = 0  "V" - Shape
%          1  "\" - Shape
%         -1  "/" - Shape
%
%  default: Mode = 1;  switch at the even Rows, "V"-Shape
%
%*******************************************************
%
% Example:
%
%%-----------------------------
%% Matlab Logo
%
%  zz = membrane(1,8);
%  si = size(zz);
%
%  x = linspace(-1,1,si(2));
%  y = linspace(-1,1,si(1));
%
%  xx = sign(x) .* (abs(x).^1.5);
%  yy = sign(y) .* (abs(y).^1.5);
%
%  [xx,yy] = meshgrid(xx,yy); % NonRegular Grid
%
%  % Triangulation
%  [i1,i2,ind] = grid2tri(si);
%
%%-----------------------------
%% One Patch per Face
%
%  figure, hold on, view([-37 30]);
%
%  patch( 'xdata' ,  x(i2)  , ...
%         'ydata' ,  y(i1)  , ...
%         'zdata' , zz(ind) , ...
%     'edgecolor' , 'k'    , ...
%     'facecolor' , 'r'          );
%
%%-----------------------------
%% Vertices and Faces
%
% figure, hold on, view([-37 30])
%
% patch( 'vertices'        , [ xx(:) yy(:) zz(:) ] , ...
%        'faces'           ,   ind'   , ...
%        'facevertexcdata' ,  zz(:)   , ...
%        'cdatamapping'    , 'scaled' , ...
%        'edgecolor'       , 'k'      , ...
%        'facecolor'       , 'interp'         );
%
%

Nout = nargout;

ok = ( isnumeric(si) & ( prod(size(si)) == 2 ) );
if ok
   ok = all( isfinite(si) & ( si > 1 ) & ( mod(si,1) == 0 ) ); 
end

if ~ok
    error('Size must be a 2-Element Vector with Integers larger 1.');
end

if nargin < 2
   mode = 1;
else
   if ~( isnumeric(mode) & ( prod(size(mode)) == 1 ) )
       error('Mode must be a single Numeric.');
   end
end

mode = [ real(mode)  imag(mode) ];
mode = sign(mode);

swap = mode(1);
mode = mode(2);

%----------------------------------------
% Indices, run over Columns

si = ( si(:)' - 1 ) .* [ 1  2 ];

ii = ( 1 : prod(si) );

jj = ii - si(2) * floor( ii / si(2) );
jj = jj + si(2) * ( jj == 0 );
jj = ceil( jj / 2 );

ii = ceil( ii / si(2) );

%----------------------------------------
% Swap JJ if even/odd II ans NonZero Mode(1)

if swap == 0
   im = ones(size(ii));     % Allways Odd II
else
   im = ii - 2 * floor( ii / 2 );           % True for Odd II
   jm = ( 1 + swap ) / 2 - swap * im;       % IM | ~IM
   jj = jj + ( 1 + si(2)/2 - 2*jj ) .* jm;
end

if mode == 0
   jm = jj - 2 * floor( jj / 2 );  % True for Odd JJ
else 
   % Mode:    1  ==>  JM ==  IM
   %         -1  ==>  JM == ~IM
   jm = ( 1 - mode ) / 2 + mode * im;
end

%****************************************
% Three Points

ii = ii([1 1 1],:);
jj = jj([1 1 1],:);

%----------------------------------------
% 1. Point

ii(1,:) = ii(1,:) + ( 1 - jm );
jj(1,:) = jj(1,:) + ( 1 - im );

%----------------------------------------
% 2. Point

kk = ( 1 : 2 : si(1)*si(2)-1 );

ii(2,kk) = ii(2,kk) + jm(kk);
jj(2,kk) = jj(2,kk) + ( 1 - im(kk) );

ii(2,kk+1) = ii(2,kk+1) + ( 1 - jm(kk+1) );
jj(2,kk+1) = jj(2,kk+1) + im(kk+1);

%----------------------------------------
% 3. Point

ii(3,:) = ii(3,:) + jm;
jj(3,:) = jj(3,:) + im;


%****************************************

if     Nout < 2
       ii = ii + ( jj - 1 ) * (si(1)+1);
elseif Nout > 2
       kk = ii + ( jj - 1 ) * (si(1)+1);
end


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = separate(name,sep);

% SEPARATE  separates Name
%
% String = SEPNAME( Name , Seperator  )
%
% String is a CellStringArray with removed surrounded Blanks
%

Nin = nargin;

if Nin < 1
   name = '';
end

if Nin < 2
   sep = '%';
end


%********************************************

str = cell(1,0);

if isempty(name)
   return
end

if ~chkstr(name)
    error('Name must be a String.')
end

if ~( chkstr(sep,1) & ( prod(size(sep)) == 1 ) )
   error('Seperator must be a single Character.')
end

n = size(name,2);

%---------------------------------------------
% Find Seperator in Name

is = ( double(name) == double(sep) );

if all(is)
   str    = cell(1,n-1);
   str(:) = { '' };
   return
end

%---------------------------------------------

i0 = ~is(1);
i1 = ~is(n);

is = cat( 2 , ones(1,i0) , is , ones(1,i1) );

is = find( is );

is = is(:);

ni = size(is,1) - 1;

if ni == 0 
   return
end
     
%---------------------------------------------
% [ Start  End ]

ind = ( 1 : ni ) ;

is  = cat( 2 , is(ind)+1 , is(ind+1)-1 ) - i0;

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------
   
   ind = [ 1  ni ];

   ind = ( 1 : ni );

   is = is(ind,:);

   str = cell(1,ni);

   for ii = 1 : ni

       str{ii} = rmblank( name( is(ii,1) : is(ii,2) ) , 2 );

   end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
