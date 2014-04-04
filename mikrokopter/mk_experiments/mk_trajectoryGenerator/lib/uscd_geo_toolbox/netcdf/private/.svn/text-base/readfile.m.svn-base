function  [Msg,bb,ref] = readfile(file,mnl,ref01);

% READFILE  Read Text from File and Converts to CellStringArray
%
% [Msg,CellString] = READFILE( File )
%
%---------------------------------------------------------
% READFILE( File , MarkerNewLine )  
%
% Replace the String's defined by MarkerNewLine with NewLine,
%   default: EMPTY
%
%---------------------------------------------------------
% [Msg,CellString,Reference] = ...
%    READFILE( File , MarkerNewLine , ReferenceCell)
%  
% Returns References in String, using GET_REF, 
%   ReferenceCell = { S1 S2 C1 C2 }, type: >> help get_ref
%
%---------------------------------------------------------
%
% see also CHAR2CELL, GET_REF, LOADFILE
%

 MaxSize = 2*2^20; % 500000;  % MaximumFileSize [Bytes]

 Msg = '';
 
 Msg0 = 'READFILE: ';

 bb  = cell(0,1);
 ref = cell(0,2);

 nl  = char(10);

 Nin = nargin;

%------------------------------------------------------

 if Nin == 0
    Msg = [Msg0 'Input File is missing.'];
    return
 end

 if isempty(file)
    return
 end

 % Marker for NewLine, using in char2cell
 if Nin < 2
   mnl = '';
 end

 if Nin < 3
   ref01 = { '' '' };
 end

 %----------------------------------------------

 [Msg,bb] = loadfile(file,MaxSize,'char');

 if ~isempty(Msg)
     Msg = [ Msg0 'Error call LOADFILE.' nl Msg ];
     if isempty(bb)
        return
     end
     warning(Msg);
 end

 %----------------------------------------------
 % Char to CellString

  try

    [Msg,bb,ref] = char2cell(bb,mnl,ref01);

  catch

     Msg = [Msg0 'Error using CHAR2CELL.' ...
             nl lasterr ];  

  end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,bb] = loadfile(file,varargin);

% LOADFILE  Load binary Data from File, using FREAD
%
% [ Msg , V ] = LOADFILE( FileName , MaxSize , Precision )
%
%   MaxSize     MaximumFileSize [Bytes],   default: 2097152 == 2MBytes
%   Precision   Precision for using FREAD, default: 'char'
%
%  For more Informations  type: >> help fread
%
%  In case of Precision 'char', a CharacterString will returned,
%   if all Bytes are valid Characters for conversion, if some
%   Characters are invalid, Msg is not empty.
%

msg = '';
bb  = [];

msg0 = 'LOADFILE: ';

%------------------------------------------------
% Defaults

m = 2*2^20; % 500000;  % MaximumFileSize [Bytes]
p = 'char';

Nin = nargin;

%************************************************
% Check File

if Nin == 0
   msg = [msg0 'Input File is missing.'];
   return
end

if isempty(file)
   return
end

if ~chkstr(file,1)
   msg = [msg0 'Input File must be a String.'];
   return
end

%------------------------------------------------
% Get MaxSize and Precision

for ii = 1 : Nin-1

    v = varargin{ii};

    if chkstr(v,1)

       p = v;

    elseif ( isnumeric(v)  &  ( prod(size(v)) == 1 ) )

       if ( isfinite(v)  &  ( v >= 0 ) )
          m = v;
       end

    end

end

%************************************************
% Open and Read File

  fid = fopen(file,'r');

  if fid == -1  
     msg = [ msg0  'Can''t open File.' ];
     return
  end

 %----------------------------------------------
 % Check Size of File

  d = dir(file);

  if isempty(d)
   d = dir( which(file) );  % which(file) gives the full Name
                            %  [ PathName FileName ]
  end

  if d.bytes > m
    msg = [ msg0 'File too large, Limit = '  ...
            sprintf('%.0f Bytes',m) '.' ];
    fclose(fid);
    return
  end


 %----------------------------------------------

  try
     bb = fread(fid,p);
  catch
     msg = [ msg0 'Error call FREAD.' char(10) lasterr ];
  end

  fclose(fid);

 %----------------------------------------------
 % Precision: 'char'  ==>  Transform to String

  if isempty(msg) & isequal( p , 'char' )
    
    % Check Characters

    if any( ( 126 < bb ) & ( bb < 160 ) );

       % Old DOS !!!
       old = [ 132   148   129   142   153   154 ];
       new = [ 'ä'   'ö'   'ü'   'Ä'   'Ö'   'Ü' ];

       for ii = 1 : size(old,2)
           bb(find(bb==old(ii))) = double(new(ii));
       end 

    end

    bb(find(bb==12)) = [];   % FormFeed

     ok = ( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

     if ~all(ok)

         msg = [ msg0 'Invalid Characters in File.' ];

         if ~any(ok)
             bb = [];
         else
                ok  = find(~ok);
             bb(ok) = 32;
         end

     end

     bb = char(bb(:)');

  end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
