function [Msg,c,b] = readppm(file,area,out);

% READPPM  Reads ImageData from PPM-Format
%
% [ Msg , C ] = READPPM(FileName)
%
%    returns the 3-dimensional UINT8-ColorMatrice C, contains
%     the RGB-ColorValues ( 0 <= C <= 255 ).
%
%    Multiple frames will added to the 4. Dimension of C.
%
%
% [ Msg , C , ColorMap ] = READPPM( FileName , [Area] )
%
%    returns the 2-dimensional Indexed-ColorMatrice C, 
%     and the corresponding RGB-ColorMap with ColorValues ( 0 <= C <= 1 ).
%
%    Multiple frames will added to the 3. Dimension of C.
%
% Use the second Input Area to read only a selected Area from Image:
%
% Area = [ X1 Y1   X2  Y2  ]
%        [ X1 Y1  Width+i Height+i ]
%        [ +XStart  +YStart   ]  Read to end
%        [ -XEnd    -YEnd     ]  Read from Start
%        [ +Width+i +Height+i ]  Read from Start
%        [ -Width+i -Height+i ]  Read from End
%
% Negative or ZERO Coordinates are measured from End.
% Negative  Width or Hight looks backward.
%
% READPPM( FileName , Area , OutFile ) writes the selected Image directly
%  to OutFile in PPM-Format.
%
%  see also: WRTPPM, READXPM, WRTXPM, IMREAD, IMWRITE
%
%
% Description:  PPM format, a widely-used, simple, 
%      portable, but non-compressed format.  PPM images can be converted
%      to gif, jpg, tif, bmp, pict, and nearly every other image format
%      know to man (or nerd).   Look for the 'netpbm' distribution on 
%      the internet.
%
%

Nin  = nargin;
Nout = nargout;

Msg = '';
c   = [];
b   = [];

nl  = char(10);

lmax = 100;    % Max Number of HeaderLines
cm   = '#';    % CommentMark
hd   = { 'P6' 'P3' };   % First Line

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

if Nin < 1
   Msg = sprintf('%sNot enough InputArguments.',Msg0);
   return
end

Msg = cell(0,1);

if ~( ischar(file) & ~isempty(file) & ...
      ( prod(size(file)) == size(file,2) ) );
    Msg = cat(1,Msg,{'File must be a String.'});
end

if Nin < 2, area = []; end

if ~isempty(area)
    area = area(:)';
    ok = ( isnumeric(area) & any( size(area,2) == [ 2  4 ] ) );
    if ok
       ok = all( mod(real(area),1) == 0 );
    end
    if ~ok
        Msg = cat(1,Msg,{'Area must be a 2 or 4 element Vector of Integer.'});
    end
end

if Nin < 3, out = ''; end

if ~isempty(out)
   if ~( ischar(out) & ( prod(size(out)) == size(out,2) ) );
       Msg = cat(1,Msg,{'OutFile must be a String.'});
   elseif exist(out,'file') == 2
       Msg = cat(1,Msg,{'OutFile allready exist.'});
   end
end
  
%--------------------------------------------------

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%s Invalid Inputs.\n%s',Msg0,Msg);
    return
end

Msg = '';

%**************************************************
% Open File

fid = fopen(file,'r');

if isempty(fopen(fid))
   Msg = sprintf('%sError open File "%s" for reading.',Msg0,file);
   return
end

%*****************************************************************
% Get Header: hd  Width Height MaxColorValue

str0 = sprintf('%sEnd of File reached during read Header.',Msg0);
str1 = sprintf('%sTo much CommentLines.',Msg0);

cc = 0;

inf = zeros(1,0);

%************************************
% FirstLine

while ( cc < lmax )  &  ( size(inf,2) < 3 )

   m  = 3 - size(inf,2);

   cc = cc + 1;
   bb = fgetl(fid);
   if isequal(bb,-1)
      Msg = str0;
      fclose(fid);
      return
   end
   
   bb = rmblank(bb,2*i);

   %------------------------------------
   % Check FirstLine

   if cc == 1

      bb = rmblank(bb,2*i);
      ok = ( size(bb,2) >= 2 );
      if ok
         ok = any(strcmp(bb(1:2),hd));
      end
      if ~ok
          Msg = sprintf('%sFirst HeaderLine "%s" or "%s" expected in PPM-File "%s".',Msg0,hd{:},file);
          fclose(fid);
          return
      end

      hd = bb(1:2);

      bb = rmblank(bb(3:end),2);

   end

   if ~isempty(bb)
       if ~strcmp(bb(1),cm);
           try
              v = eval(sprintf('[%s]',bb));
           catch
              Msg = sprintf('%sError evaluate HeaderValues: "%s".',Msg0,bb);
              fclose(fid);
              return
           end
           if ~isempty(v)
               ok = isnumeric(v);
               if ok
                  v = v(:)';
                  v = v( 1 : min(size(v,2),m) );
                  ok = all( isfinite(v) & ( v > 0 ) & ...
                            ( round(v) == v ) );
               end
               if ~ok
                   Msg = sprintf('%sInvalid HeaderValues: "%s".',Msg0,bb);
                   fclose(fid);
                   return
               end
               inf = cat( 2 , inf , v );
           end
       end
   end

end

if ( size(inf,2) < 3 )
   Msg = str1;
   fclose(fid);
   return
end

%*****************************************************************

sz = inf([1 2]);   % [ Width Height ]
mc = inf(3);       % MaxColorValue

l = ( 1 + floor( log(mc)/log(2) / 8 ) );  % Length of Color in Bytes

b = 8 * l;                                % Length of Color in Bits

fmt = sprintf('ubit%.0f',b);              % Binary format to read ColorValues

l = 3 * l;                                % Length of Pixel in Bytes RGB

is_plain = strcmp(hd,'P3');

%*****************************************************************
% Get Area

% Area = [ X1 Y1   X2  Y2  ]
%        [ X1 Y1  Width+i Height+i ]
%        [ +XStart  +YStart   ]  Read to end
%        [ -XEnd    -YEnd     ]  Read from Start
%        [ +Width+i +Height+i ]  Read from Start
%        [ -Width+i -Height+i ]  Read from End

if ~isempty(area)

    img  = imag(area);
    area = real(area);

    sgn =    sign(area);
    neg =  ( sgn <= 0 );
    img = ~( img  == 0 );

    n = size(area,2);

    def = [ 1  1  sz ];
  
    ii = [ 1  2 ];
    jj = [ 3  4 ];

    if n == 2
       area = cat(2,area,NaN*zeros(1,2));
       kk = ( ( neg & img ) | ~( neg | img ) );
       area(jj) = def(ii+2*kk);
       area(ii) = area(ii) .* ( 1 - 2 * ( ~img & neg ) );
       area(ii) = real(area(ii)) + (area(jj)-sgn).*img;
     else
       area = area + def([3 4 3 4]) .* ( ~img & neg );  % From End
       area(jj) = area(jj) + (area(ii)-sgn(jj)).*img(jj);
     end

     for ii = 1 : 2
         jj = ii + [ 0  2 ];
         area(jj) = min(max(area(jj),def(jj(1))),def(jj(2)));
         area(jj) = area( jj + 2 * [ 1 -1 ] * ( area(jj(1)) > area(jj(2)) ) );
     end

else

     area = [ 1  1  sz ];

end

wh = area([3 4]) - area([1 2]) + 1;

pz = prod(sz);  % Number of Pixel per Frame

read_all = isequal(wh,sz);

%*****************************************************************
% Open Outfile

is_out = ~( isempty(out) | isempty(area) );

wrtout = ( is_out & ~( read_all | is_plain ) );  % Write directly

if wrtout
   oid = fopen(out,'w');
   if isempty(fopen(oid))
      Msg = sprintf('%sError open File "%s" for writing.',Msg0,out);
      return
   end
else
   oid = 0;
end

%*****************************************************************
% Read Image

reqout = ( Nout > 1 );   % OutputImage requested

if is_plain

   if ( reqout | is_out ), c = fscanf(fid,'%f'); end

else

   p = ftell(fid); fseek(fid,0,'eof');
   s = ftell(fid); fseek(fid,p,'bof'); 

   s = s - p;                          % Number of Bytes until EOF
   m = l * pz;                         % Number of Bytes per Frame

   if ~( mod(s,m) == 0 )

       wrn = sprintf('Size of %.0f-bit ImageData doesn''t correspond with HeaderSize: %.0f %.0f.', ...
                      b , sz );
       warning(wrn);

   end

   n = ceil(s/m);                      % Number of Frames

   if read_all

      if ( reqout | is_out ), c = fread(fid,fmt); end

   else

      c = uint8(zeros(0,1));
      w = l * sz(1);                   % Bytes per Row
      p = p + ( area(1) - 1 ) * l;     % Offset incl. XStart

      ok = ( wh(1) == sz(1) );         % Read full Rows
      nr = wh(2) + ( 1 - wh(2) ) * ok; % Number of Rows to Read
      nb = 3 * prod(wh(1:(1+ok)));     % Number of Bytes to Read once

      if wrtout
         dt = clock;
         dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
         dt = datestr(dt,0);
         fprintf(oid,'P6\n');
         fprintf(oid,'# %s\n',dt);
         fprintf(oid,'# %s: Matlab %s %s\n',fcn,version,computer);
         if n > 1
            fprintf(oid,'# %.0f Frames\n',n);
         end
         fprintf(oid,'%4d %4d\n',wh);
         fprintf(oid,'%4d\n',255);
      end

      for ii = 1 : n
          q  = p + ( ii - 1 ) * m;     % Offset incl. Frame
          for jj = 1 : nr
              f  = q + ( area(2) + jj - 2 ) * w; 
              fseek( fid , f , 'bof' );
              [d,cnt] = fread( fid , [ nb  1 ] , fmt );
              if ( cnt < nb )
                 d = cat( 1 , d , zeros(nb-cnt,1) );
              end
              if ~( mc == 255 )
                  d = 255 * d / mc;
              end
              if reqout
                 c = cat( 1 , c , uint8(d) );
              end
              if wrtout, fwrite(oid,d,'uchar'); end
          end
      end

      sz = wh;
      pz = prod(sz);  % Number of Pixel per Frame

   end

end

fclose(fid);

if wrtout, fclose(oid); end

%------------------------------------------------------------------------

wrtimg = ( is_out & ~wrtout ); % Write Image using WRTPPM

if ~( reqout | wrtimg )
    return
end

%------------------------------------------------------------------------
% Check, Transform Image

if ( isequal(c,-1) | isempty(c) )
   Msg = sprintf('%sEnd of File reached during read ImageData.',Msg0);
   return
end

m = 3 * pz;         % Number of ColorValues
s = size(c,1);

if ~( mod(s,m) == 0 )
    wrn = sprintf('Size of %.0f-bit ImageData doesn''t correspond with HeaderSize: %.0f %.0f.', ...
                    b , sz );
    warning(wrn);
    if s < m
       c = cat( 1 , c , zeros(m-s,1) );
    else
       c = cat( 1 , c , zeros(floor(s/m)-s,1) );       
    end
end

n = size(c,1) / m;  % Number of Frames

c = reshape(c,3,sz(1),sz(2),n);

c = permute(c,[3 2 1 4]);

if is_plain & ~isequal(wh,sz)
   c = c(area(2):area(4),area(1):area(3),:,:);
end

if read_all | is_plain
   if ~( mc == 255 )
       c = 255 * c / mc;
   end
   c = uint8(c);
end

if wrtimg
   Msg = wrtppm(c,out);
   if ~isempty(Msg)
       Msg = sprintf('%sError call WRTPPM.\n%s',Msg0,Msg);
   end
end

if Nout < 3
   return
end

%*****************************************************************
% Transform to Indexed

c = permute(c,[1 2 4 3]);
c = reshape(c,sz(2),sz(1)*n,3);

[c,b] = mat2ind(c);

c = reshape(c,sz(2),sz(1),n);

b = b / 255;

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,m] = mat2ind(c)

% MAT2IND  Converts Matrice to Indexed
%
% [ Y , Map ] = MAT2IND( X )
%
%   X   numeric Matrice with max. 3 Dimensions
%
%   Y   Indexed Matrice, refering to Values of Map
%
%  Map  Unique Values of X: [ N  by  size(X,3) ] 
%
% Convert back:
%
%   s1 = size(Y,1); s2 = size(Y,2); s3 = size(Map,2);
%   X  = feval(class(Map),zeros(s1,s2,s3));
%   for ii = 1 : s3
%       X(:,:,ii) = reshape(Map(Y,ii),s1,s2);
%   end
%

if ~( isnumeric(c) & ( ndims(c) <= 3 ) )
    error('Input must be numeric with max. 3 Dimensions.');
elseif isempty(c)
    c = zeros(size(c));
    m = zeros(0,size(c,3));
    return
end

cl = class(c);
if ~strcmp(cl,'double')
    c = double(c);
end

if ~all( isfinite(c(:)) | isnan(c(:)) )
    error('Input must have finite Values or NaN.');
end

s1 = size(c,1);
s2 = size(c,2);
s3 = size(c,3);

n = s1 * s2;

c = reshape(c,n,s3);

si = ( 1 : n )';

for ii = s3 : -1 : 1
    [m,jj] = sort(c(si,ii));
       si  = si(jj);
end

c = c(si,:);

is_nan = find(isnan(c));
if ~isempty(is_nan)
    nn = ceil(max(c(:))) + 1;
    c(is_nan) = nn+1;
    is_nan = 1;
end

c(2:n,:) = c(2:n,:) - c(1:(n-1),:);

if is_nan    
    n = ~( sum( abs(c) < 1e3*eps , 2 ) == s3 );
else
    n = ~( sum(     c  < 1e3*eps , 2 ) == s3 );
end

n(1) = 1;

m = cumsum(c(find(n),:),1);

if is_nan
   m( find( m > nn ) ) = NaN;
end

n = cumsum(n);

c = reshape(n,s1,s2);

c(si) = c;

if ~strcmp(cl,'double')
    m = feval(cl,m);
end

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
