function [msg,varargout] = load_asc(file,key,km,cm);

% LOAD_ASC  Load Variables from an ASCII-File
%
% [Msg,V] = LOAD_ASC( File , KEY, KeywordMarker, CommentMarker )
%
% Load Keywords and their Values from an ASCII-File.
% 
% The Keywords are marked in the File with the KeywordMarker.
%   A Keyword ends normally with NL, CommentMarker or a additional Marker,
%    specified as second Value for KeywordMarker in a 2-Element CellString:
%
%     KeywordMarker = { KeywordStart KeywordEnd }
%
% The Input KEY specifies the Keywords to load.
%  KEY can be a CharacterArray or CellArray of Strings.
%  Use an empty Value for KEY to load all KeyWords from File.
%  The WildCard '*' can used in the KeywordStrings,
%   i.e. do NOT use "*" as Marker!
%
% The Comment's are marked in the File with the CommentMarker.
%   A Comment is valid until End of Line (NL).
%
% If the first Comment, which follow a KeyWord, starts like:
%
%             % EVAL VariableName
%
%    the Value of the KeyWord will interpreted as MatlabExpression,
%    in which the Value of "VariableName" returns the Value for the KeyWord.
%   Example:
%             # ColorMap
%             %  EVAL cmap
%                cmap = jet(256);
%
% Use "" around a CommentMarker or KeywordMarker to mask them.
%
% Defaults:
%
%    KeywordMarker  = { '#'  ',' }
%    CommentMarker  =   '%'
%
%
%----------------------------------------------------------------------------
% OUTPUT
%
%  Msg contains the ErrorMessage. Msg is empty, if LOAD_ASC works successful.
%
%  V   returns an [ N by 5 ]-CellArray for N KeyWords:
%
%     { KeyWord  Value  Class  NDim  Size }
%
%     Class is a String with the class of Value, 
%      it's set to 'none', if the KeyWord can't evaluate successfull,
%      it's empty, if a requested KeyWord was not found in the File.
%
%     The first Column of V is equal to the nonempty Input KEY, 
%       but may differ, if the Wildcard "*" is used in the Strings for KEY.
%
% Multiple Output:
%
%  [ Msg , V , v1 , v2 , .. , v# ] = LOAD_ASC( File , ... )
%
%  returns the Values of the KeyWords in v#, 
%
%----------------------------------------------------------------------------
% Please note: 
%
%   To define the Value of an KeyWord (Variable), use a normal Matlab-Syntax.
%
%   KeyWords, defined in the File, can used as VariableNames in following
%    Statements in the File.
%
%   Do NOT use Names for KeyWords in Syntax "evl_kwd_???"
%    or MatlabCommands like: 
%       "for"   "if"      "try"     "catch"    "break"    "end"       "all"
%       "cat"   "eval"    "round"   "isempty"  "warning"  "lastwarn"  "class"
%       "char"  "ischar"  "strcmp"  "strrep"   "findstr"  "rmblank"   "size"
%
%----------------------------------------------------------------------------
%
% see also:  WRT_ASC
%
%----------------------------------------------------------------------------
% Example:
%  
%  cmap  = jet(12);            % 12 by 3 RGB-Matrice
%  cdat  = rand(5,6);
%  txt   = ['Lon';'Lat'];
%  keys  = {'ColorMap','CData %Random','Text'};
%  form  = ['%6.2f';'%8.3f'];
%  comm1 = ['Jet         ';'use COLORMAP'];
%  comm3 = 'Labels for Axis'; 
%  
%  % Write the Data into ASCII-File
%  msg = wrt_asc('test.asc',keys,{cmap;cdat;txt},form,{ comm1 ; '' ;comm3}) 
%
%
%  key = {'ColorMap','CData','Text'};       
%
%  % Load the Data from ASCII-File
%  [msg,V]  = load_asc('test.asc',key); 
%
%   % produce an 3 by 3 CellArray, contains KeyWord, Value, Class
%
%  [msg,V,v1,v2,v3] = load_asc('test.asc',key); 
%
%   % Returns the Values in v1 .. v3
%



Nin  = nargin;
Nout = nargout - 1;

varargout = cell(Nout,1);

msg = cell(0,1);

nl  = char(10);

ini = cell(0,3);  % { KeyWord Value Class }


if ~( Nout == 0 )
   varargout(:) = {  [] };
   varargout(1) = { ini };
end

quota = '"';

%*********************************************************************
% Get Inputs

if Nin < 1
   msg = catmsg('Input File is missing.','',i);
   return
end
if Nin < 2
   key = cell(0,1);
end
if Nin < 3
 km = { '#'  ',' };
end
if Nin < 4
 cm = '%';
end

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------------
% File

is_file = 0;

s1 = size(file,1);

if chkcstr(file,1)
   str = sprintf('%s\n',file{:});
elseif ~( ischar(file) & ( ndims(file) == 2 ) )
   msg = cat( 1 , msg , {'Input File must be a CharacterArray or CellArray of Strings.'} );
else
  if isempty(file)
     str = '';
  elseif ( s1 > 1 ) | any( file == 10 ) | any( file == 13 )
     str = file;
     if s1 > 1
        str = cat( 2 , str , char(10*ones(s1,1)) );
        str = permute(str,[2 1]);
        str = str(:);
        str = permute(str,[2 1]);
     end
  else
     is_file = 1;
  end
end

%---------------------------------------------------------------------
% KeyWords

ok = isempty(key);
if ok
   key = cell(0,1);
else
   [ok,key] = chkcstr(key);
   if ok
      key = key(:);
      ok  = ~any(strcmp(key,''));
   end
end

if ~ok

   msg = cat( 1 , msg , ...
              {'Input KeyWords must be a CharacterArray or CellArray of Strings.'} , ...
              {' An empty String is only valid for a single KeyWord.' } );

elseif ~( Nout == 0 )

  varargout{1} = cell(size(key,1),3);

  if ~isempty(key)

     varargout{1}(:,1) =  key;
     varargout{1}(:,2) = { [] };  % Value
     varargout{1}(:,3) = { '' };  % Class

  end

end

%---------------------------------------------------------------------
% KeyMarker

[ok,km] = chkcstr(km);
if ok
   ok = ~isempty(km{1});
   if ok
      km = km(:);
      km = cat( 1 , km , {''} );
      km = km([1 2]);
      k0 = rmblank(km{1},2);    % KeyStart
      k1 = rmblank(km{2},2);    % KeyEnd
      ok = ~isempty(k0);
   end
end

if ~ok
   msg = cat( 1 , msg , ...
              {'Input KeyMarker must be a CharacterArray or CellArray of Strings.'} , ...
              {' 2 KeyMarker defines { KeyWordStart KeyWordEnd }.' } , ...
              {' KeyWordStart must be nonempty.' }  );
end

%---------------------------------------------------------------------
% CommentMarker

if ~chkstr(cm)
   msg = cat( 1 , msg , {'Input CommentMarker must be a String.'} );
end

%---------------------------------------------------------------------

if ~isempty(msg)
   msg = catmsg( 'Invalid Inputs,' , msg , i );
   return
end

%*********************************************************************
% Read File

if is_file

   [msg,str] = loadfile(file,'char');

   if ~isempty(msg)
       msg = catmsg( 'Invalid File.' , msg , i );
       if isempty(str)
          return
       end
       warning(msg)
   end

end

msg = '';

if ( Nout == 0 ) | isempty(str)
   return
end

str = strrep(str,char([13 10]),nl);  % DOS: CRLF --> LF  
str = strrep(str,char(13),nl);       % MAC: CR   --> LF

%*********************************************************************
% Mask

msk = { k0  char(1)
        k1  char(2)
        cm  char(3)
      quota char(4) };

str = mask(str,msk,quota);

msk = msk(:,[2 1]);  % Flip to remask

%*********************************************************************
% Initialisation for KeyWords
%
% { KeyWord  Value  Class }
%

ini = get_ini(str,k0,k1,cm,nl,msk);

if isempty(ini) | isempty(key)
   if ~isempty(ini)
      varargout{1} = ini;
      n = min( size(ini,1) , Nout-1 );
      varargout(2:n+1) = ini(1:n,2);
   end
   return
end

%*********************************************************************
% Check for WildCards

ni = size(ini,1);
nk = size(key,1);

iwc = find( ~( sum( ( char(key(:,1)) == '*' ) , 2 ) == 0 ) );

if ~isempty(iwc)
   kwc =  key(iwc,1);
   nwc = size(kwc,1);
   for ii = nwc : -1 : 1
       jj = find( strwcmp( ini(:,1) , kwc{ii} ) );
       if ~isempty(jj)
          nj = size(jj,1);
          kk = iwc(ii) + [ -1  1 ];
          key = cat( 1 , key(1:kk(1),1) , ini(jj,1) , key(kk(2):end,1) );
       end
   end
end

%*********************************************************************
% Extract KeyWords

ni = size(ini,1);
nk = size(key,1);

key      = cat( 2 , key , cell(nk,4) );
key(:,2) = { [] };   % Value
key(:,3) = { '' };   % Class

ok = zeros(ni,1);
ik = zeros(nk,1);

for ii = nk : -1 : 1

    jj = find( strcmp( ini(:,1) , key{ii,1} ) );

    if ~isempty(jj)

       kk = cat( 1 , 1 , find(~ok(jj)) );

       if iwc  &  ( size(kk,1) > 2 )
          ll = find( jj(kk) == nk );
          if isempty(ll)
             jj = jj(kk(end));
          else
             jj = jj(kk(ll(end)));
          end
       else
          jj = jj(kk(end));
       end

       ok(jj) =  1;
       ik(ii) = jj;

    end

end

ii = find( ~( ik == 0 ) );

key(ii,:) = ini(ik(ii),:);


varargout{1} = key;
n = min( nk , Nout-1 );
varargout(2:n+1) = key(1:n,2);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = get_ini(str,k0,k1,cm,nl,msk);

% GET_INI   Get MarkerInitialisation

ini = cell(0,5);   %  { KeyWord  Value  Class  NDim  Size }

%---------------------------------------------------------------------
% Get MarkerIndize, Flag and Length

mark = { k0 k1 cm nl };

[im,fm,lm] = get_mark(str,mark,1);

if isempty(im)
   return
end

%---------------------------------------------------------------------
% NewLine before Start and End

nm = size(mark,2);

im = cat( 2 , 0  , im , size(str,2)+1 );
fm = cat( 2 , nm , fm , nm            );
  
%---------------------------------------------------------------------
% Check correct Statement of Marker

%       Flag FlagBefore True/False

chk = [  1    nm         1     % NewLine before KeyWordStart
         2    1          1     % KeyWordStart before KeyWordEnd
         3    3          0  ]; % No Comment before Comment


[im,fm] = chk_mark(im,fm,chk,1);

if isempty(im)
   return
end

if ~any( fm == 1 )
   return
end

%---------------------------------------------------------------------
% Get { KeyWord String  ''  EVAL_VarName }

ini = get_str(str,im,fm,lm,msk);

if isempty(ini)
   ini = cell(0,5);
   return
end

%---------------------------------------------------------------------
% Get { KeyWord Value Class }

ini = eval_key(ini);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [im,fm,lm] = get_mark(str,mark,ret);

% GET_MARK  Returns MarkerIndize, Flag and Length
%

nm = size(mark,2);

im = zeros( 1 , 0  );  % Index for Marker
fm = zeros( 1 , 0  );  % Flag  for Marker
lm = zeros( 1 , nm );  % Length of Marker

%---------------------------------------------------------------------
% Find Marker: Index & Flag

for ii = 1 : nm

    lm(ii) = size( mark{ii} , 2 );

    jj = findstr( str , mark{ii} );

    if ~isempty(jj)

       im = cat( 2 , im ,      jj );
       fm = cat( 2 , fm , ii+0*jj );
    
    elseif any( ii == ret )

       return
 
    end

end

%---------------------------------------------------------------------
% Sort

[im,si] = sort(im);

fm = fm(si);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [im,fm,cc] = chk_mark(im,fm,chk,mode);

% CHK_MARK  Check Correct Statement of Marker

% chk = [ Flag FlagBefore True/False ];

cc = ones(size(im));

for ii = 1 : size(chk,1)

    jj = find( fm == chk(ii,1) );

    if ~isempty(jj)

        ok = ( fm(jj-1) == chk(ii,2) );

        ok = ( 1 - chk(ii,3) ) + ( 2*chk(ii,3) - 1 ) * ok;

        jj = jj(find(~ok));

        cc(jj) = 0;

        if mode

          fm(jj) = [];
          im(jj) = [];

        end

    end

end

if mode
   return
end

jj = find(~cc);

im(jj) = [];
fm(jj) = [];


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = get_str(str,im,fm,lm,msk);

% GET_STR  Returns CellArray: { KeyWord String '' }

%---------------------------------------------------------------------
% Append KeyWordStart

im = cat( 2 , im , size(str,2)+1 );
fm = cat( 2 , fm , 1 );

%---------------------------------------------------------------------

ik0 = find( fm == 1 );  % KeyWordStart
icm = find( fm == 3 );  % Comment

%---------------------------------------------------------------------

nk = size(ik0,2) - 1;   % KeyWordNumber

ini = cell(nk,5);

ini(:) = { char(zeros(1,0)) };

ok = zeros(nk,1);

for ii = 1 : nk
  
    %----------------------------------------------------------
    % KeyWord
 
    i0 = ik0(ii);
    i1 = ik0(ii)+1;

    i00 = im(i0) + lm(fm(i0));
    i11 = im(i1) - 1;

    ini{ii,1} = mask( rmblank( str(i00:i11) , 2 ) , msk );

    ok(ii) = ~isempty(ini{ii,1});

    if ok(ii)

        i00 = im(i1) + lm(fm(i1)) * ( fm(i1) == 2 );

        i1  = ik0(ii+1);
        i11 = im(i1) - 1;

        s = str(i00:i11);

        %----------------------------------------------------------
        % Remove Comment for KeyWord

        ic = icm( find( ( i0 < icm ) & ( icm < i1 ) ) );

        if ~isempty(ic)
             i0  = im(ic) - i00 + 1;                              % Start
             i1  = im(ic+1) - im(ic) + ( im(ic-1) == im(ic)-1 );  % Lenght
             i1  = min( i1 , size(s,2)-i0+1 );
             %-----------------------------------------------------
             % Check first CommentLine for EVAL
               ind = ( 1 : i1(1) ) + i0(1) - 1;
               evl = rmblank( s(ind) , 2 );
               evl = sepname( evl , NaN , ' ' ); % { CommentMarker  'EVAL'  VarName }
               if prod(size(evl)) >= 3
                  if strcmp(evl{2},'EVAL')
                     ini{ii,4} = rmblank( evl{3} , 2 , char([39 34]) );
                  end
               end
             %-----------------------------------------------------
             ic  = grp2ind( i0 , i1 );
           s(ic) = [];
        end

        %----------------------------------------------------------
        % String for KeyWord

        ini{ii,2} = mask( rmblank(s,2) , msk );

    end
 
end

ini( find(~ok) , : ) = [];

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function evl_kwd_iii = eval_key(evl_kwd_iii)

% EVAL_KEY  Evaluate String from 2. Column of INI: { KeyWord String '' }
%
%   returns: { KeyWord  Value  Class }
%

evl_kwd_www = warnstat;
              warning('off');

evl_kwd_iii(:,3) = {'none'};   % class


for evl_kwd_jjj = 1 : size(evl_kwd_iii,1)

    evl_kwd_okk = 0;
    evl_kwd_vvv = [];

  %-----------------------------------------------------------
  % Check for EVAL
  %-----------------------------------------------------------
  if ~isempty( evl_kwd_iii{evl_kwd_jjj,4} )
  %-----------------------------------------------------------

     if ~isempty( evl_kwd_iii{evl_kwd_jjj,2} )

        if ~isempty( findstr( evl_kwd_iii{evl_kwd_jjj,2} , ...
                              evl_kwd_iii{evl_kwd_jjj,4} ) )

            evl_kwd_sss = cat( 2 , 'evl_kwd_vvv = ' , ...
                                    evl_kwd_iii{evl_kwd_jjj,4} , ';' );

            evl_kwd_sss = cat( 2 , evl_kwd_iii{evl_kwd_jjj,2} , char(10) , ...
                                   evl_kwd_sss );

            try
               eval(evl_kwd_sss);
               evl_kwd_okk = 1;
            catch
               evl_kwd_sss = strrep( evl_kwd_sss , char([39 10 39]) , char(10) );
               evl_kwd_sss = strrep( evl_kwd_sss , char([39 10   ]) , char(10) );
               evl_kwd_sss = strrep( evl_kwd_sss , char([   10 39]) , char(10) );
               evl_kwd_sss = evl_kwd_sss((1+strcmp(evl_kwd_sss(1),char(39))):end);
               try
                  eval(evl_kwd_sss);
                  evl_kwd_okk = 1;
               end
            end

        end

     end

  %-----------------------------------------------------------
  end
  %-----------------------------------------------------------

  %-----------------------------------------------------------
  if ~evl_kwd_okk
  %-----------------------------------------------------------

    for evl_kwd_rrr = {  ''  ';,' }

        evl_kwd_str = evl_kwd_iii{evl_kwd_jjj,2};

        if ~isempty(evl_kwd_rrr{1})
            evl_kwd_str = rmblank( evl_kwd_str , -2i , evl_kwd_rrr{1} );
        end

        for evl_kwd_kkk = [ 91  39  123 ]    %  []  ''  {}

            evl_kwd_sss = cat( 2 , char(evl_kwd_kkk) , ...
                                        evl_kwd_str , ...
                                   char(evl_kwd_kkk+2*(~(evl_kwd_kkk==39))) );

            % Check before warning if "''"
            % Warning: Out of range or non-integer values truncated during 
            %          conversion from double to character.

            lastwarn('');

            try

               evl_kwd_vvv = eval(evl_kwd_sss);

            catch

               warning('failed');

            end

            evl_kwd_okk = isempty(lastwarn);

            if  evl_kwd_okk  &  ischar(evl_kwd_vvv)  % Check valid Characters

                evl_kwd_okk = ( ( evl_kwd_vvv ==  9 ) |  ...
                                ( evl_kwd_vvv == 10 ) |  ...
                                ( evl_kwd_vvv == 13 ) |  ...
                        ( ( evl_kwd_vvv >=  28 )  &  ( evl_kwd_vvv <= 126 ) ) | ...
                        ( ( evl_kwd_vvv >= 160 )  &  ( evl_kwd_vvv <= 255 ) )       );

                evl_kwd_okk = (  evl_kwd_okk & ( evl_kwd_vvv == round(evl_kwd_vvv) ) );

                evl_kwd_okk = all( evl_kwd_okk(:) );

                if evl_kwd_okk & ( evl_kwd_kkk == 91 )
                % Character with "[]" ==> Check with "{}" for CELLSTR  !!!!!
                   try
                      evl_kwd_okk = iscellstr(eval(cat(2,'{',evl_kwd_str,'}')));
                   catch
                      evl_kwd_okk = 0;
                   end 
                end   

            end

            if evl_kwd_okk
               break
            end

        end

        if evl_kwd_okk
           break
        end

    end

  %-----------------------------------------------------------
  end
  %-----------------------------------------------------------

    if evl_kwd_okk

       evl_kwd_iii{evl_kwd_jjj,2} =       evl_kwd_vvv;
       evl_kwd_iii{evl_kwd_jjj,3} = class(evl_kwd_vvv);

       % Evaluate KeyWord
       try
          eval([ evl_kwd_iii{evl_kwd_jjj,1} ' = evl_kwd_vvv;' ]);
       end

    end

  evl_kwd_iii{evl_kwd_jjj,4} = ndims(evl_kwd_iii{evl_kwd_jjj,2});
  evl_kwd_iii{evl_kwd_jjj,5} = size(evl_kwd_iii{evl_kwd_jjj,2});

end

warning(evl_kwd_www);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = mask(str,msk,quota);

% MASK   Mask Characters in String

if nargin < 3
   quota = '';
end

for ii = 1 : size(msk,1)

    if ~isempty(msk{ii,1}) & ~isempty(msk{ii,2})

       str = strrep( str , cat(2,quota,msk{ii,1},quota) , msk{ii,2} );

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

m = 2*2^20; % 2 MBytes;  % MaximumFileSize [Bytes]
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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,cmp,wc,nn,cc] = strwcmp(str,cmp,wc,nn,cc)

% STRWCMP   Compare Strings, including WildCards
%
% OK = STRWCMP( STR , CMP , WildCard )
%
% STR: CharacterArray or CellStringArray to compare with Comp
%
% CMP: CharacterArray or CellStringArray with strings to compare with STR
%         strings can contains WildCards
%
% WildCard: specify WildCard to use, default: '*'
%
% OK : logical Array with same size of STR, contains 1 if 
%        any strings of CMP match string of STR
%
% Special Output:  [ok,cmp,wc,nn,cc] = STRWCMP( STR , CMP , [WildCard] )
%
%  to use in follwing statements: ok = STRWCMP( STR , cmp , wc , nn , cc );
%
%  which makes it a bit faster.
%
% see also: STRCMP, FINDSTR
%

Nin  = nargin;
Nout = nargout;

Msg = '';
 nl = char(10);

%***************************************************************
% Check Inputs

%---------------------------------------------------------------
% String

if Nin < 1
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
           'First Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% CompareString

if Nin < 2
   cmp = cell(0,0);
end

if ~( iscell(cmp) & isempty(cmp) )
   [ok,cmp] = chkcstr(cmp,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Second Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% WildCard

if Nin < 3
   wc = '*';
elseif ~( ischar(wc) & ( prod(size(wc)) == 1 ) )
   Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
         'WildCard must be a Single Character.' ];
end
  
%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%***************************************************************

si = size(str);

if ( isempty(str) | isempty(cmp) ) & ( Nout <= 1 ) 
   ok = zeros(si);
   return
end

cmp = cmp(:);

if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
   ok = ones(si);
   return
end

%***************************************************************
% Analyze CompareStrings

nc = size(cmp,1);

ok = ( Nin == 5 );

if ok
   ok = ( isequal(size(cc),[nc 1]) & iscell(cc) & ...
          isequal(size(nn),[nc 2]) & isnumeric(nn)    );
   if ok
      try
         ok = cat( 1 , cc{:} ); 
         ok = ( size(ok,2) == 3 );
      catch
         ok = 0;
      end
   end
end


%--------------------------------------------------
if ~ok
%--------------------------------------------------

  cc    = cell(nc,1);
  cc(:) = { zeros(0,3) };  % { [ Start End N  ] }

  nn = zeros(nc,2);        %   [ Ncmp  sum(N) ] 

  for ii = 1 : nc

    if ~isempty(cmp{ii})

       iwc = ( double(cmp{ii}) == double(wc) );

       if ~any( iwc )

          nn(ii,:) = size(cmp{ii},2);
          cc{ii}   = [ 1  nn(ii,:) ];

       else

         %--------------------------------------------
         % Remove Duplicate WildCards

         iwc = find( iwc );
         if ~( prod(size(iwc)) == 1 )
            jj = find( diff(iwc) == 1 );
            if ~isempty(jj)
               cmp{ii}(iwc(jj+1)) = [];
            end
         end

         %--------------------------------------------
         % Get Start End
     
         iwc = ( double(cmp{ii}) == double(wc) );
  
          n  = size(iwc,2);

          if ( n == 1 ) & ( iwc == 1 ) & ( Nout <= 1 ) % Wildcard only
             ok = ones(si);
             return
          end

          i0 = ~iwc(1);
          i1 = ~iwc(n);

         iwc = cat( 2 , ones(1,i0) , iwc , ones(1,i1) );

         iwc = find( iwc );

         iwc = iwc(:);

           n = size(iwc,1) - 1;

         cc{ii} = zeros(n,3);

         if n > 0      

            cc{ii}(:,[1 2]) = cat( 2 , iwc(1:n)+1 , iwc((1:n)+1)-1 ) - i0;

            cc{ii}(:,3) = cc{ii}(:,2) - cc{ii}(:,1) + 1;
 
         end

         nn(ii,:) = cat( 2 , size(cmp{ii},2) , sum(cc{ii}(:,3),1) );

       end

    end

  end

%--------------------------------------------------
end
%--------------------------------------------------

if ( Nout > 1 )

  if ( isempty(str) | isempty(cmp) )
     ok = zeros(si);
     return
  end

  if any(strcmp(cmp,wc))  % Wildcard only
     ok = ones(si);
     return
  end

end

%***************************************************************
% Compare

ok = zeros(si);

for ii = 1 : prod(si)
 
    s2 = size(str{ii},2);

    for jj = 1 : nc

        ok(ii) = ( ( s2 == 0 ) & ( nn(jj,1) == 0 ) );
       
        if ok(ii)
           break
        end
       
        if ( s2 >= nn(jj,2) ) & ~( nn(jj,1) == 0 )

           ok(ii) = compare(str{ii},cmp{jj},cc{jj});

           if ok(ii)
              break
           end
            
        end

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = compare( str , cmp , cc )

sc = size(cmp,2);

ok = 1;

for ii = 1 : size(cc,1) 

    s2 = size(str,2);

    ok = ( ok &  ( s2 >= cc(ii,3) ) );

    if ~ok
       break
    end

    ic  = ( cc(ii,1) : cc(ii,2) );

    i01 = ( cc(ii,[1 2]) == [ 1  sc ] );

    i0  = ( 1 : cc(ii,3) );
    i1  = s2 - cc(ii,3) + i0;

    i2  = cat( 2 , findstr( str , cmp(ic) ) , 0 );
    i2 = i2(1);
    i2 = ( i2 + cc(ii,3) - 1 ) * ~( i2 == 0 );

    i2  = s2       * (  i01(1) &  i01(2) & strcmp(str    ,cmp    ) ) + ...
          cc(ii,3) * (  i01(1) & ~i01(2) & strcmp(str(i0),cmp(ic)) ) + ...
          s2       * ( ~i01(1) &  i01(2) & strcmp(str(i1),cmp(ic)) ) + ...
          i2       * ( ~i01(1) & ~i01(2) );

    ok = ( ok & ~( i2 == 0 ) );

    if ~ok
       break
    end

    str = str( i2+1 : s2 );

end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

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
%    A positive complex Value for DIM, to removes Blanks only from Start,
%    A negative complex Value for DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 32  13  10  9 ];  % [ Space CR LF TAB ]
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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,nl] = catmsg(msg,str,n);

% CATMSG   Expands Message
%
% [ NewMsg , Separator ] = CATMSG( OldMsg , String , BlankNumber )
%
%  OldMsg       Original MessageString
%  String       CharacterArray of CellStringArray of Strings to append to OldMsg, 
%                 separated with NewLine and Blanks
%  BlankNumber  Number of Blanks, added to NewLine to separate Strings.
%                 Separator = cat( 2 , char(10) , char(32*ones(1,BlankNumber)) )
%
%  If the Value for BlankNumber has an Imaginary Part, the UpperCase Name of the 
%   calling M-File will added to OldMsg: 
%    OldMsg = cat( 2 , MFILENAME , ': ' , OldMsg )
%   The BlankNumber will expand by the Length of [ MFILENAME ': ' ]
%
% example:  
%
%  txt = { 'Size of X and Y must be agree.'
%          'Invalid Method.'                }
%
%  msg = catmsg( 'Invalid Inputs', txt )
%  msg = catmsg( 'Invalid Inputs', txt , 8 )
%  msg = catmsg( 'Invalid Inputs', txt , 1+i )
%

Nin = nargin;

err = '';
 nl = char(10);

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------
% Message

if Nin < 1
   msg = char(zeros(1,0));
   return
end

if isempty(msg)
   msg = char(zeros(1,0));
else
   if ~chkstr(msg,0)
      err = [ err  nl(1:(end*(~isempty(err)))) ...
              'Message must be a String.' ];
   end
end

%---------------------------------------------------------------
% String

if Nin < 2
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      err = [ err  nl(1:(end*(~isempty(err)))) ...
             'String must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% BlankNumber

if Nin < 3
   n = 0;
elseif ~( isnumeric(n) & ( prod(size(n)) == 1 ) )
      err = [ err  nl(1:(end*(~isempty(err)))) ...
              'BlankNumber must be a single Number.' ];
end
 
%---------------------------------------------------------------

if ~isempty(err)
   error(err)
end

%*********************************************************************

%---------------------------------------------------------------
% Prepend CALLER

if ( imag(n) == 0 )

   c = char(zeros(1,0));

else

   c = dbstack;

   if prod(size(c)) < 2
      c = 'Matlab';
   else
      c = sepname(c(2).name,0,filesep,0);  % Calling M-FileName
      c = upper(c);
     jj = find( double(c) == double('.') );
     if ~isempty(jj)
        c = c(1:max(jj)-1); 
     end
   end

   c = cat( 2 , c , ': ' );

   n = real(n) + size(c,2);

end

%---------------------------------------------------------------

nl = cat( 2 , nl , char( 32 * ones(1,ceil(n))) );

if ~isempty(str)
   str = strrep(str,char(10),nl);
   str = strhcat(str,nl);
end

if isempty(str) & isempty(msg)
   return
end

inl = ~( isempty(msg) | isempty(str) );

msg = cat( 2 , c , msg ,  nl(1:(end*inl)) , str );

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = sepname(name,dep,sep,mode);

% SEPNAME  separates Name
%
% String = SEPNAME( Name , Depth , Seperator ,Mode )
%
% Depth gives the number of Recursions
% 
% Use Depth = NaN  to get all parts of Name in String.
%  In this case String is a CellStringArray.
%
% Mode =  1  Recursion from Start
%        -1  Recursion from End
%         0  Recursion from End, 
%             special Handling for Class- and Private-Directories
%
% Defaults:  Depth     = 0
%            Seperator = FILESEP
%            Mode      = 0    ( 1  if Depth = NaN )
%


Nin = nargin;

if Nin < 1
   name = '';
end

if Nin < 2
   dep = 0;
end

dep_nan = isnan(dep);

if Nin < 3
   sep = filesep;
end

if Nin < 4
   mode = dep_nan;
end

%********************************************

if dep_nan
   str = cell(1,0);
else
   str = char(zeros(1,0));
end

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
% Take care for duplicate Seperators

if ~dep_nan | ( double(sep) == char(32) )

   is( find( is(:,1) > is(:,2) ) , : ) = [];

end

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------
if dep_nan
   
   ind = [ 1  ni ];

  flip = ~( mode == 1 );

   ind = ind( [ 1  2 ] + [ 1 -1 ]*flip );

   ind = ( ind(1) : 1-2*flip : ind(2) );

   is = is(ind,:);

   str = cell(1,ni);
  
   for ii = 1 : ni

       str{ii} = name( is(ii,1) : is(ii,2) );

   end

   return

end

%---------------------------------------------

ii = ni - 1 * ( ni > 1 );

nn = name( is(ii,1) : is(ii,2) );

ic = strcmp( nn(1) , '@'       );
ip = strcmp( nn    , 'private' );

id = 1 * ic + 2 * ip;

dep = dep + 1 + id * ( mode == 0 ) * ( ni > 1 );

dep = min( dep , ni );

%---------------------------------------------

is(1,1) = is(1,1) - 1 * ( is(1,1) > 1 );  % Start incl. Seperator

ind = ( 1 : ni-1 );

is(ind,2) = is(ind,2) + 1;                % End incl. Seperator

is(:,2) = is(:,2) - is(:,1) + 1;          % Length

%---------------------------------------------

ind = ( 1 : dep ) + ( ni - dep ) * (~( mode == 1 ));

is  = is(ind,:);

ind = grp2ind(is(:,1),is(:,2));

str = name(ind);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
