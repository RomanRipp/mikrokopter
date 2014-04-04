function [msg,bb] = loadfile(file,varargin);

% LOADFILE  Load binary Data from File, using FREAD
%
% [ Msg , V ] = LOADFILE( FileName , MaxSize , Precision , Extension )
%
%   MaxSize     MaximumFileSize [Bytes],   default: 2097152 == 2MBytes
%   Precision   Precision for using FREAD, default: 'char'
%   Extension   optional Extension to expand FileName,id neccessary
%               '.ext'
%
%  For more Informations to Precision  type: >> help fread
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

m   = 2*2^20; % 2 MBytes;  % MaximumFileSize [Bytes]
p   = 'char';
ext = '';

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

    if chkstr(v,1) & ~isempty(v)

       if strcmp(v(1),'.')
          ext = v;
       else
          p = v;
       end
 
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
     if ~isempty(ext)
        ie = 1 + strcmp( file(end) , ext(1) );
        ff = cat( 2 , file , ext(ie:end) );
        fid = fopen(ff,'r');
     end
     if fid == -1   
        msg = [ msg0  'Can''t open File.' ];
        return
     end
     file = ff;
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
       old = [ 132   148   129   142   153   154  225 ];
       new = [ 'ä'   'ö'   'ü'   'Ä'   'Ö'   'Ü'  'ß' ];

       n  = size(old,2);
       ok = zeros(1,n);
 
       for ii = 1 : n
           if ( ( ii < n ) | any(ok(1:(n-1*(ii>1)))) )
              jj = find(bb==old(ii));
              ok(ii) = ~isempty(jj);
              if ok(ii)
                 bb(jj) = double(new(ii));
              end
           end
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
                ok  = find(~ok); abs(bb(ok))
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
