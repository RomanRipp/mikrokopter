function [x,msg,prc] = iload(file)

% ILOAD  Loads Variable from UBIT# formatted File
%
% X = ILOAD( File ) loads the Variable X from File, 
%                    written by ISAVE.
%
% ILOAD without InputArguments try to assign the
% Variable und their Name, written in the File,
% in the Workspace of caller.
%
% ILOAD without Inputs or empty Input calls UIGETFILE.
%
% [X,Msg,Precision] = ILOAD( File )
%
% returns ErrorMessages and the Precision of X
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
%  NN    UBIT#     Values        ceil(NN*NB/8)
%  
%  NH = 24 + NS + 4*ND + 8*NG + NF;
%  NN = sum(LG)
%
%  MachineFormat: ieee-be
%
%--------------------------------------------------
%
% see also: ISAVE, LOAD, FREAD
%

machineformat = 'ieee-be';  % MachineFormat

bm  = 24;    % BasicNumber of HeaderBytes

%****************************************************

x   = [];
msg = '';
prc = NaN;

Nin  = nargin;
Nout = nargout;

nox = ( Nout == 0 );
out = ( Nout == 2 );

%****************************************************
% Check / Open File

if Nin == 0

   file = '';

elseif ~( strcmp(class(file),'char') & ...
          ( prod(size(file)) == size(file,2) ) )

   msg = 'File must be a String.';

   if out
      return
   end

   error(msg);

end

%----------------------------------------------------

if isempty(file)

   str = 'Select a File to load ...';
   [file,pfad] = uigetfile('*.int',str);

   if isequal(file,0)
      file = '';
      if nox
         clear x
      end
      return
   end

   file = fullfile(pfad,file);

else

   if ~( exist(file,'file') == 2 )
       msg = 'File doesn''t exist.';
   else
       f = which(file);
       if ~isempty(f)
           file = f;
       end
   end

end

%----------------------------------------------------

if isempty(msg)

   bt = dir(file);
   bt = bt.bytes;

   if bt < bm
      msg = sprintf('File must have min. %.0f Bytes for Header.',bm);
   end

end

%----------------------------------------------------

if isempty(msg)

   [fid,msg] = fopen(file,'r',machineformat);
 
   if fid == -1
      msg = sprintf('Error open File "%s" to read.\n%s',file,msg);
   else
      msg = '';
   end

end

%----------------------------------------------------

if ~isempty(msg)

    if out
       return
    end

    error(msg);

end

%*********************************************************
% Read / Check Header
%
%  1     UINT16    NHeaderBytes         2    nh
%  1     UINT8     NString              1
%  NS    UCHAR     String (Name)     NS*1    name
%  1     UINT8     NDim                 1    ng
%  ND    UINT32    DimSize           ND*4    sx
%  1     UINT32    NGroup               4    ng
%  NG    UINT32    StartIndex        NG*4    ix
%  NG    UINT32    Length            NG*4    lx
%  1     UINT8     NFormatString        1    nf
%  NF    UCHAR     FormatString      NF*1    frm
%  1     UINT8     NBit                 1    nb
%  1      INT32    Offset               4    off
%  1     DOUBLE    Precision            8    prc
%
%  NN    UBIT#     Values    ceil(NN*NB/8)   nx
%  

nh  = fread(fid,1,'uint32');

hb = bm;

if nh < bm

   msg = sprintf('Number of HeaderBytes less then %.0f.',bm);

elseif bt < nh

   msg = sprintf(['Checksum Error in Header (Bytes).\n' ...
             'File must have min. %.0f Bytes for Header.'],nh);

else

   ns  = fread(fid,1,'uint8');
   hb  = hb + ns;

   if bt < hb
      msg = sprintf(['Checksum Error in Header (Name).\n' ...
               'File must have min. %.0f Bytes for Header.'],hb);
   else
   
      nm  = fread(fid,[1 ns],'uchar');
      nd  = fread(fid,1,'uint8');
      hb  = hb + 4*nd;

      if bt < hb
         msg = sprintf(['Checksum Error in Header (Size).\n' ...
                  'File must have min. %.0f Bytes for Header.'],hb);
      else

         sx  = fread(fid,[1 nd],'uint32');
         ng  = fread(fid,1,'uint32');
         hb  = hb + 8*ng;

         if bt < hb
            msg = sprintf(['Checksum Error in Header (Group).\n' ...
                 'File must have min. %.0f Bytes for Header.'],hb);
         else

            msg = cell(0,1);

            ix  = fread(fid,[ng 1],'uint32');
            lx  = fread(fid,[ng 1],'uint32');
            nf  = fread(fid,1,'uint8');
            hb  = hb + nf;

            if bt < hb
               msg = sprintf(['Checksum Error in Header (Format).\n' ...
                    'File must have min. %.0f Bytes for Header.'],hb);
            else

               msg = cell(0,1);

               frm = fread(fid,[1 nf],'uchar');
               nb  = fread(fid,1,'uint8');
               off = fread(fid,1,'int32');
               prc = fread(fid,1,'double');

               px = prod(sx);
               nx =  sum(lx);

               pos = ftell(fid);

               if ~( ( nh == pos ) & ( nh == hb ) )
                   msg = cat(1,msg,{'Checksum Error in HeaderBytes.'});
               end

               if px < ix(ng)+lx(ng)-1
                  msg = cat(1,msg,{'Invalid Size in Header.'});
               end

               sz = pos + ceil(nx*nb/8);  % Expected Size of File !!!

               if ~( bt == sz ) %%% bt < sz
                  msg = cat(1,msg,{sprintf(['File must have %.0f Bytes' ...
                                         ' for Header and Data.'],sz)});
               end

               if ~isempty(msg)
                   msg = sprintf('%s\n',msg{:});
               else
                   msg = '';
               end

            end

         end

      end

   end

end

%----------------------------------------------------

if ~isempty(msg)
    fclose(fid);
    msg = sprintf('Invalid File.\n%s',msg);
    if out
       return
    end
    error(msg);
end

%*********************************************************
% Read Variable

frm = char(frm); %%%  sprintf('bit%.0f',nb);

x = NaN * zeros(px,1);

if nx > 0
   for ii = 1 : ng
       x(ix(ii)+(1:lx(ii))-1) = fread(fid,lx(ii),frm);
   end
end

fclose(fid);

%----------------------------------------------------

if nx > 0
   x = x + off;
   x = x * prc;
end

x = reshape(x,sx);

%----------------------------------------------------

if nox & ~isempty(nm)

   nm = char(nm(:)');

   try
      assignin('caller',nm,x);
   catch
      msg = sprintf('Cann''t assign Variable "%s".',nm);
   end

   if isempty(msg)
      clear x
   else
      warning(msg);
   end

end
