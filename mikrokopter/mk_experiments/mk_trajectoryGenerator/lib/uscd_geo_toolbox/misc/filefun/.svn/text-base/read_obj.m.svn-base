function [msg,ini] = read_obj(file,km,cm,sp,cl);

% READ_OBJ  Read ObjectFile
%
% [Msg,V] = READ_OBJ( ObjectFile , ObjectMarker , CommentMarker , ...
%                                  Seperator    , ContinueMarker       )
%
% Load ObjectNames and their Parameter-Value-Pairs from an ObjectFile.
% 
% The ObjectNames are marked in the File with the ObjectMarker.
%   A ObjectName ends normally with NL, CommentMarker or a additional Marker,
%    specified as second Value for ObjectMarker in a 2-Element CellString:
%
%    ObjectMarker = { ObjectNameStart ObjectNameEnd }
%
% The Comment's are marked in the File with the CommentMarker.
%   A Comment is valid until End of Line (NL).
%
% The Parameter-Value-Pairs, seperated by Seperator, follow the ObjectName.
%   A Value is valid until End of Line (NL), 
%     i.e. a Line can contain maximum ONE Parameter-Value-Pair.
%
%   Following Lines, beginning with the ContinueMarker, will added to the Value.
%
% Use "" around a Marker to mask them.
%
% Defaults:
%
%    ObjectMarker   = { '@'  ',' }
%    CommentMarker  =   '$$'
%    Seperator      =   '='
%    ContinueMarker =   '++'
%
%
%----------------------------------------------------------------------------
% OUTPUT
%
%  Msg contains the ErrorMessage. Msg is empty, if READ_OBJ works successful.
%
%  V   returns an [ N by 2 ]-CellArray for N Objects:
%
%     { OjectName  { Parameter Value} }
%
%
% see also:  LOAD_OBJ, STR2PAR, WRT_OBJ, LOADFILE
%
%



Nin  = nargin;
Nout = nargout-1;

msg = cell(0,1);

nl  = char(10);

ini = cell(0,2);  % { KeyWord  { Parameter Value } }

quota = '"';

ext = '.obj';  % DefaultExtension to use to expand File

%*********************************************************************
% Get Inputs

if Nin < 1
   msg = catmsg('Input File is missing.','',i);
   return
end

if Nin < 2
 km = { '@'  ',' };
end
if Nin < 3
 cm = '$$';
end
if Nin < 4
 sp = '=';
end
if Nin < 5
 cl = '++';
end

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------------
% File

if ~( chkstr(file) & ~isempty(file) )
   msg = cat( 1 , msg , {'Input File must be a nonempty String.'} );
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
              {' 2 KeyMarker defines { KeyWoordStart KeyWordEnd }.' } , ...
              {' KeyWordStart must be nonempty.' }  );
end

%---------------------------------------------------------------------
% CommentMarker

if ~chkstr(cm)
   msg = cat( 1 , msg , {'Input CommentMarker must be a String.'} );
end
%---------------------------------------------------------------------
% Seperator

if ~chkstr(sp)
   msg = cat( 1 , msg , {'Input Seperator must be a String.'} );
end

%---------------------------------------------------------------------
% ContinueMarker

if ~chkstr(cl)
   msg = cat( 1 , msg , {'Input ContinueMarker must be a String.'} );
end

%---------------------------------------------------------------------

if ~isempty(msg)
   msg = catmsg( 'Invalid Inputs,' , msg , i );
   return
end

%*********************************************************************
% Read File

[msg,str] = loadfile(file,'char',ext);

if ~isempty(msg)
   msg = catmsg( 'Invalid File.' , msg , i );
   if isempty(str)
      return
   end
   warning(msg);
end

msg = '';

if ( Nout == 0 )
   return
end

str = strrep(str,char([13 10]),nl);  % DOS: CRLF --> LF  
str = strrep(str,char(13),nl);       % MAC: CR   --> LF

%*********************************************************************
% Mask

msk = { k0  char(1)
        k1  char(2)
        cm  char(3)
        sp  char(4)
        cl  char(5)
      quota char(6)  };

str = mask(str,msk,quota);

msk = msk(:,[2 1]);  % Flip to remask

%*********************************************************************
% Initialisation for KeyWords
%
% { KeyWord  String }
%

ini = get_ini(str,k0,k1,cm,nl,msk);  


%*********************************************************************
% Get { KeyWord  { Parameter Value }  }

for ii = 1 : size(ini,1)

    ini{ii,2} = get_par(ini{ii,2},sp,cl,nl,msk);

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = get_ini(str,k0,k1,cm,nl,msk);

% GET_INI   Get MarkerInitialisation

ini = cell(0,2);   %  { KeyWord  String  }

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
% Get { KeyWord String }

ini = get_str(str,im,fm,lm,[1 3],0,msk,nl);


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = get_par(str,sp,cl,nl,msk)


ini = cell(0,2);  % { Parameter Value }

%---------------------------------------------------------------------
% Get MarkerIndize, Flag and Length

mark = { sp cl nl };

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

chk = [  1    nm         1     % NewLine before Seperator
         2    nm         1 ];  % NewLine before ContinueMarker


[im,fm] = chk_mark(im,fm,chk,0);

if isempty(im)
   return
end

if ~any( fm == 1 )
   return
end

%---------------------------------------------------------------------
% Get Parameter-Value

ini = get_str(str,im,fm,lm,[1 2],1,msk,nl);

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

function ini = get_str(str,im,fm,lm,fl,mode,msk,nl);

% GET_STR  Returns CellArray: { KeyWord String }
%
% im     MarkerIndex in str
% fm     MarkerFlag
% lm     MarkerLength
% fl     [ KeyWordFlag  CommentFlag ]
% mode   0  KeyWord   and Comment
%        1  Seperator and ContinueMarker
%

%---------------------------------------------------------------------
% Append KeyWordStart

im = cat( 2 , im , size(str,2)+1 );
fm = cat( 2 , fm , fl(1) );

%---------------------------------------------------------------------

ik0 = find( fm == fl(1) );  % KeyWordStart
icm = find( fm == fl(2) );  % Comment

%---------------------------------------------------------------------

nk = size(ik0,2) - 1;   % KeyWordNumber

ini = cell(nk,2);

ini(:) = { char(zeros(1,0)) };

ok = zeros(nk,1);

for ii = 1 : nk
  
    %----------------------------------------------------------
    % KeyWord
 
    i0 = ik0(ii) - mode;
    i1 = ik0(ii) - mode + 1;

    i00 = im(i0) + lm(fm(i0));
    i11 = im(i1) - 1;

    ini{ii,1} = mask( rmblank( str(i00:i11) , 2 ) , msk );

    ok(ii) = ~isempty(ini{ii,1});

    if ok(ii)

        %----------------------------------------------------------
        % Extract Full String for KeyWord

        i00 = im(i1) + lm(fm(i1)) * ( fm(i1) == 2 );

        i1  = ik0(ii+1) - mode;
        i11 = im(i1) - 1;

        s = str(i00:i11);

        %----------------------------------------------------------
        % Find CommentMarker for KeyWord

        ic = icm( find( ( i0+mode < icm ) & ( icm < i1 ) ) );

        if mode
           % Add Seperator
           ic = cat( 2 , ik0(ii) , ic );
        end

        %----------------------------------------------------------
        % Remove Comment for KeyWord
        % Extract from Seperator/ContinueMarker to next Marker==NewLine

        if ~isempty(ic)

             i0  = im(ic) - i00 + 1;   % Start
             i1  = im(ic+1) - im(ic);  % Lenght

             if mode
                jj = ( fm(ic) == fl(2) );  % ContinueMarker
                dd = lm(fm(ic)+1) - jj;
                i0 = i0 + dd;
                i1 = i1 - dd;
                s(i0(find(jj))) = nl;      % Continue --> NewLine
             else
                i1 = i1 + ( im(ic-1) == im(ic)-1 );
             end

             %%% i0 = i0 + lm(fm(ic)) * mode;
             %%% i1 = i1 + ( lm(fm(ic+1)) - lm(fm(ic)) ) *   mode + ...
             %%%           (    im(ic-1) ==   im(ic)-1 ) * (~mode);

             %%% i1(end) = i1(end) - lm(fm(ic(end)+1))*mode;

             ic  = grp2ind( i0 , i1 );

             if mode
                ic(find(ic>size(s,2))) = [];
                s = s(ic);
             else
                s(ic) = [];
             end

        end

        %----------------------------------------------------------
        % String for KeyWord

        ini{ii,2} = mask( rmblank(s,2) , msk );

    end
 
end

ini( find(~ok) , : ) = [];


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

