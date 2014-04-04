function [p,msg] = str2par(str,sp,cl);

% STR2PAR  Returns Property-Value-Pairs from TextInput
%
%  [ P , Msg ] = STR2PAR( String , Seperator , ContinueMarker )
%
% The Parameter-Value-Pairs are seperated by Seperator.
%   A Value is valid until End of Line (NL), 
%     i.e. a Line can contain maximum ONE Parameter-Value-Pair.
%
%   Following Lines, beginning with the ContinueMarker, will added to the Value.
%
% Use "" around a Marker to mask them.
%
% Defaults:
%
%    Seperator      =   '='
%    ContinueMarker =   '++'
%
% The returned 2-Row-CellStringArray "P" contains: { Parameter ; Value }
%
% ErrorMessages are returned in Msg.
%
%
% example: 
%
%  nl  = char(10);
%
%  str = [ 'Longitude = 57 34.45 W ' nl ...
%          'Latitude  = 67 54.2  N ' nl ...
%          'Depth     =  230 m '     nl ...
%          'Comment   = Test of '    nl ...
%          '          ++STR2PAR'     nl ...
%          'Profile   = 000'                 ];
%
%  p = str2par(str)
% 
%% p = 
%
%%    'Longitude'     'Latitude'      'Depth'    'Comment'      'Profile'
%%    '57 34.45 W'    '67 54.2  N'    '230 m'    [1x16 char]    '000'    
%
%


p   = cell(2,0);

msg = cell(0,1);

Nout = nargout;
Nin  = nargin;

%*********************************************************************

if Nin == 0
   msg = '';
   return
end

if Nin < 2
 sp = '=';
end
if Nin < 3
 cl = '++';
end

 nl = char(10);

quota = '"';

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------------
% String

if ~chkstr(str,0)
   msg = cat( 1 , msg , {'Input File must be a nonempty String.'} );
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

%*********************************************************************

if ~isempty(msg)
   msg = catmsg( 'Invalid Inputs,' , msg , i );
   if Nout < 2
      error(msg)
   end
   return
end

msg = '';

if isempty(str)
   return
end

str = strrep(str,char([13 10]),nl);  % DOS: CRLF --> LF  
str = strrep(str,char(13),nl);       % MAC: CR   --> LF

%*********************************************************************
% Mask

msk = { sp  char(1)
        cl  char(2)
      quota char(3)  };

str = mask(str,msk,quota);

msk = msk(:,[2 1]);  % Flip to remask

%*********************************************************************
% Get { Parameter Value }

p = get_par(str,sp,cl,nl,msk);

p = permute( p , [ 2  1 ] );


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

ini = get_str(str,im,fm,lm,[1 2],1,msk);

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

function ini = get_str(str,im,fm,lm,fl,mode,msk);

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

             ic  = grp2ind( i0 , i1 );

             if mode
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

