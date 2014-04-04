function [msg,file,cnf,frm] = isave(file,x,prc)

% ISAVE   Save Variable using UBIT# Format
%
% ISAVE( [File] , Variable , [Precision+i*Mode] )
%
% The Values of Variable are rounded to multiples of Precision.
% Mode defines the Format to Write the Values:
%
%  Mode >= 0  : BIT#
%       <  0  : INT#
%
%  "#" is the Number of Bits to use, it will determined by
%  the Range of the Values.
%
% The nonzero Mode defines the multiple number  of Bytes (8 Bit)
%   to use for Format:  1=8, 2=16, 3=32, 4=64
%
% ISAVE without or empty File-Input calls UIPUTFILE.
%
% [Msg,File,Config,Format] = ISAVE( ... )
%
% returns ErrorMessages, FileName and Configuration
%
%   Config = [ Precision  NBit  NGroup ]
%
% ISAVE is recommended for large Matrices to save Diskspace.
% ISAVE writes only Groups of not-NaN-Values
%
%--------------------------------------------------
% Example:
%
%   a = 5 * rand(10,1);
%   isave('test.int',a,0.1);
%   b = iload('test.int')
%
%   [ a   b ]
%
%--------------------------------------------------
% Header of File                      bytes
%
%  1     UINT32    NHeaderBytes         4    NH
%  1     UINT8     NString              1    NS
%  NS    UCHAR     String (Name)     NS*1    
%  1     UINT8     NDim                 1    ND
%  ND    UINT32    DimSize           ND*4
%  1     UINT32    NGroup               4    NG
%  NG    UINT32    StartIndex        NG*4 
%  NG    UINT32    Length            NG*4    LG
%  1     UINT8     NFormatString        1    NF
%  NF    UCHAR     FormatString      NF*1
%  1     UINT8     NBit                 1    NB
%  1      INT32    Offset               4
%  1     DOUBLE    Precision            8
%
%  NN     BIT#     Values        ceil(NN*NB/8)
%  
%  NH = 24 + NS + 4*ND + 8*NG + NF;
%  NN = sum(LG)
%
%  MachineFormat: ieee-be
%
%--------------------------------------------------
%
% see also: ILOAD, SAVE, FWRITE
%

machineformat = 'ieee-be';  % MachineFormat

bm  = 24;    % BasicNumber of HeaderBytes

acc = 1e-3;  % Accuracy in case of unknown Precision

%****************************************************

Nin =   nargin;
out = ( nargout > 0 );

msg = cell(0,1);

%****************************************************
% Check Inputs

if Nin == 0
   msg = 'Not enough InputArguments.';
   if out
      return
   end
   error(msg);
end

if Nin < 2
   x = [];
end

if Nin < 3
   prc = [];
end

nb = [];

%---------------------------------------------------
% File

ix = 2;

ok = ( strcmp(class(file),'char') & ...
       ( prod(size(file)) == size(file,2) ) );

if ~ok
    if Nin == 3
        msg = cat( 1 , msg , {'File must be a String.'} );
    else
        x = file;  ix = 1;
        if Nin == 1
           prc = [];
        else
           prc = x;
        end
        file = '';
    end
end

nm = inputname(ix);  % Name

%----------------------------------------------------------------
% X

if ~isnumeric(x)

    msg = cat( 1 , msg , {'Variable must be numeric.'} );

else


    sx = size(x);
    nd = size(sx,2);
    nx = prod(sx);

    ix = 1;     % Start
    lx = nx;    % Length

    if nx > 0

       x  = x(:);

       if ~strcmp(class(x),'double')
           x = double(x);
       end

       lx = isnan(x);

       if ~all( isfinite(x) | lx )

           msg = cat( 1 , msg , {'Values must be finite or NaN.'} );

       else

           % Groups of not NaN

           if  all(lx)

             lx = 0;

           elseif any(lx)

             % IND2GRP

             ix = find(~lx);

             lx = cat( 1 , 1 , find( diff(ix,1,1) > 1 )+1 , size(ix,1)+1 );
             ix = ix(lx(1:end-1));
             lx = diff(lx,1,1);

           else

             lx = nx;

           end

       end

    end

end

%----------------------------------------------------------------
% Precision

rb = 0;  % RoundBit

if ~isempty(prc)
    ok = ( isnumeric(prc) & ( prod(size(prc)) == 1 ) );
    if ok
       prc = double(prc);
        rb = imag(prc);
       prc = real(prc);
       ok  = ( isfinite(prc) & ( prc >= 0 ) );
    end
    if ~ok
        msg = cat( 1 , msg , {'Precision must be a single finite positive  numeric.'} );
    end
else
    prc = 0;
end

%----------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    if out
       return
    end
    error(msg);
end

msg = '';

%*********************************************************
% Scale X

frm = { 'bit'  'int' };
frm = frm{ 1 + ( rb < 0 ) };

rb = abs(rb);
ib = ( rb > 0 );
rb = floor(rb);
rb = rb + ib * ( rb == 0 );

rb = 8 * ib * 2^(rb-1) + ( rb == 0 );

off = 0;
nb  = rb;  % NBit

nx = sum(lx);
ng = size(ix,1);

if nx > 0

   if ng == 1

      nb = ix + ( 1 : lx ) - 1;

   else

      % GRP2IND

      nb = ones(nx,1);
      jj = cumsum( cat(1,1,lx) , 1 );

      nb(jj(1:ng)) = ix;

      if ng > 1
         nb(jj(2:ng)) = nb(jj(2:ng))-(ix(1:ng-1,1)+lx(1:ng-1)-1);
      end

      nb = cumsum(nb,1);

   end

   if ( prc == 0 )

      if all( x(nb) == round(x(nb)) )
         prc = 1;
      else
         prc = max(x(nb)) - min(x(nb));  % Intervall
         if prc == 0
            prc = min(x(nb));
         end
         prc = 10 ^ floor( log(acc*prc) / log(10) + 2*eps );
      end

   end

   x = round( x / prc );

   off = min(x(nb));

   x   = x - off;

   nb  = max(x(nb));

   if nb == 0
      nb = 1;
   else
      nb = floor( log(nb) / log(2) ) + 1;
   end

   nb  = rb * ceil( nb / rb );

   nb = min( nb , 64 );
   if strcmp(frm,'int')
      nb = ceil( log(nb/8) / log(2) );
      nb = 8 * 2^nb; % [ 8 16 32 64 ]
   end

   x   =  x  - 2^(nb-1);

   off = off + 2^(nb-1);
   
elseif ( prc == 0 )

   prc = 1;

end

cnf = [ prc  nb  ng ];

frm = sprintf('%s%.0f',frm,nb);

%*********************************************************
% Check / Open File

if isempty(file)
   if ~isempty(nm)
       str = sprintf('Save "%s" as ...',nm);
   else
       str = 'Save as ...';
   end
   [file,pfad] = uiputfile('*.int',str);
   if isequal(file,0)
      file = '';
      if ~out
          clear msg
      end
      return
   end
   file = fullfile(pfad,file);
end

[fid,msg] = fopen(file,'w+',machineformat);

if fid == -1
   msg = sprintf('Error open File "%s" to write.\n%s',file,msg);
   if out
      return
   end
   error(msg);
end

%*********************************************************
% Write Data

%  1     UINT16    NHeaderBytes         2    nh
%  1     UINT8     NString              1    ns
%  NS    UCHAR     String (Name)     NS*1    nm
%  1     UINT8     NDim                 1    nd
%  ND    UINT32    DimSize           ND*4    sx
%  1     UINT32    NGroup               4    ng
%  NG    UINT32    StartIndex        NG*4    ix
%  NG    UINT32    Length            NG*4    lx
%  1     UINT8     NFormatString        1    NF
%  NF    UCHAR     FormatString      NF*1
%  1     UINT8     NBit                 1    nb
%  1      INT32    Offset               4    off
%  1     DOUBLE    Precision            8    prc
%
%  NN    UBIT#     Values       ceil(sum(LG)*NB/8)
%  
%  NH = 21 + NS + 4*ND + 8*NG;
%  NN = sum(LG)

ns = size(nm,2) * ~isempty(nm);
nf = size(frm,2);

nh = bm + ns + 4*nd + 8*ng + nf;

fwrite(fid,nh,'uint32');
fwrite(fid,ns,'uint8');
fwrite(fid,nm,'uchar');
fwrite(fid,nd,'uint8');
fwrite(fid,sx,'uint32');
fwrite(fid,ng,'uint32');
fwrite(fid,ix,'uint32');
fwrite(fid,lx,'uint32');
fwrite(fid,nf,'uint8');
fwrite(fid,frm,'uchar');
fwrite(fid,nb,'uint8');
fwrite(fid,off,'int32');
fwrite(fid,prc,'double');

if nx > 0
   for ii = 1 : ng
       fwrite(fid,x(ix(ii)+(1:lx(ii))-1),frm);
   end
end

fclose(fid);

if ~out
    clear msg
end