function s = bkmrk2struct(file)

% BKMSTRUCT  Converts NetscapeBookmarkFile to Structure
%
% S = BKMSTRUCT( BookmarkFile )
%
% default: BookmarkFile = '$HOME/.netscape/bookmarks.html'
%
% S = isdir:   0  |  1          % true for Folder
%      link:  URL | SubStruct
%      type:  LinkType { 'folder' | 'file' | 'http' | 'ftp' }
%      name:  LinkName
%
%------------------------------------------------------------------
% Example for expected BookmarkFileStructure
%
% <DL><p>
%     <DT><A HREF="URL" ... >LinkName</A>
%     <DT><H3 FOLDED ... >FolderName</H3>
%     <DL><p>
%         <DT><A HREF="URL" ... >LinkName</A>
%         <DT><H3 FOLDED ... >FolderName</H3>
%         <DL><p>
%             ...
%         </DL>
%         <DT><A HREF="URL" ... >LinkName</A>
%     </DL>
% </DL>
%
% A Link to an URL    starts with: <DT><A HREF="
% A Link to an Folder starts with: <DT>
%                     followed by: <DL> ... </DL>
%
%------------------------------------------------------------------
%
% see also: LOADFILE, GET_REF, GRP2IND, INSERT
%

if nargin < 1
   file = fullfile(getenv('HOME'),'.netscape/bookmarks.html');
end

if ~( exist(file,'file') == 2 )
    error(sprintf('File doesn''t exist: %s',file));
end

%************************************************************
% ReadFile

fprintf(1,'\nCall LOADFILE(''%s'') ... ',file);

[m,s] = loadfile(file);

if ~isempty(m)
    fprintf(1,'%s\n%s\n','error',m);
    return
end

fprintf(1,'%s\n\n','done');

%************************************************************
% Marker

f0 = '<DL>';       % FolderStart
f1 = '</DL>';      % FolderEnd

l0 = '<DT>';       % LinkStart

u0 = '<A HREF="';  % URL_Start
u1 = '"';          % URL_End

%------------------------------------------------------------
% Extract all Lines starting with FolderMark or LinkStart

s = get_line( s , { f0  f1  l0 } );

fprintf(1,'\n');

%------------------------------------------------------------

fprintf(1,'Convert to  CharArray ... ');

[m,s] = char2cell(s);

s = char(s);

fprintf(1,'%s\n','done');

%************************************************************
% Check for Marker

fprintf(1,'Check   for Marker');

n = size(s,1);

ok = zeros(n,1);

ok(findmark(s,f0)) =  1;   % FolderStart
ok(findmark(s,f1)) = -1;   % FolderEnd

ok(findmark(s,l0))      = 2; % Link (Name)
ok(findmark(s,[l0 u0])) = 3; % Link with URL

fprintf(1,'\n');

%************************************************************
% Check Folder-Start-End

fprintf(1,'Check FolderMarker');

[s,ok] = fold_chk(s,ok,f0,f1);

fprintf(1,'\n\n');

%************************************************************
% Recurse Structure

fprintf(1,'Get BookmarkStructure');

s = get_struct(s,ok,l0,u0,u1);

fprintf(1,'\n\n');

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ii,ok] = findmark(s,m)

if ~ischar(s)
    s = char(s);
end

nm = size(m,2);
if nm > size(s,2)
   return
end

n = size(s,1);

ok = ( ( s(:,1:nm) ==       m(ones(1,n),:)  ) | ...
       ( s(:,1:nm) == lower(m(ones(1,n),:)) )       );

ok = ( sum(ok,2) == nm );

ii = find(ok);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function s = get_line(s,l);

nl = char(10);

ok = zeros(size(s));

for r = l
    str = sprintf('%s',r{1});
    fprintf(1,'Extract Lines:   %s%s... ',str,char(32*ones(1,24-size(str,2))));
    for rr = { r{1} lower(r{1}) }
        [m,rf,jj] = get_ref(s,rr{1},nl);
        if ~isempty(m)
            fprintf(1,'%s\n%s\n','error',m);
            return
        end
        ok(grp2ind(jj(:,1),jj(:,2))) = 1;
    end
    fprintf(1,'%s\n','done');
end

s = s(find(ok));


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [s,ok] = fold_chk(s,ok,f0,f1);


if any( ok == 0 )
   jj = find( ~( ok == 0 ) );
    s =  s(jj,:);
   ok = ok(jj);
end

%************************************************************
% Check Folder-Start-End

%------------------------------------------------------------
% FolderEnd before first FolderStart

ii = find( ok ==  1 );
jj = find( ok == -1 );

if ~isempty(ii) & ~isempty(jj)
    ii = ( jj > ii(1) );
    if ~all(ii)
        jj = jj(find(ii));
         s =  s(jj,:);
        ok = ok(jj);
    end
end

%------------------------------------------------------------
% Insert FolderStart after Name (ok==2)

m = size(s,2);
n = size(s,1);

ii = zeros(n,1);
ii(1:(n-1)) = ( ( ok(1:(n-1)) == 2 ) & ~( ok(2:n) == 1 ) );
ii(n)       = ( ok(n) == 2 );

if any(ii)
   ii = find(ii);
   bl = char(32*ones(1,m));
   bl(1:size(f0,2)) = f0;
    s = insert(  s , ii , bl(ones(1,size(ii,1)),:) );
   ok = insert( ok , ii , 1 );
end

%------------------------------------------------------------
% Missing FolderEnd / FolderStart

ii = ( sum( ok == 1 ) - sum( ok == -1 ) );

if     ii > 0

       % Missing FolderEnd, Add FolderEnd
       bl = char(32*ones(1,m));
       bl(1:size(f1,2)) = f1;
        s = cat( 1 ,  s , bl(ones(ii,1),:) );
       ok = cat( 1 , ok , -1*ones(ii,1) );

elseif ii < 0

       % Missing FolderStart, Remove FolderEnd
       fn = cumsum(ok.*(abs(ok)==1));  % FolderNumber
       fn = fn .* ( ok == -1 );        % FolderEndNumber
       for jj = ( -1 : -1 : ii )
              kk  = find( ( fn == jj ) );
              kk  = kk(1);
           fn(kk) = NaN;          
       end
       jj = find(~isnan(fn));
        s =  s(jj,:);
       ok = ok(jj);

end

%------------------------------------------------------------
% Sort Start and End of Folders

fi = fold_ind(ok);  % [ StartIndex  EndIndex  Number ]

if isempty(fi)
   return
end

%------------------------------------------------------------
% Remove Folders without Name (ok == 2) before

ii = ones(size(ok));

jj = find( fi(:,1) == 1 );

ii(fi(jj,[1 2])) = 0;

jj = find( fi(:,1) > 1 );

jj = jj(find(~(ok(fi(jj,1)-1) == 2)));
if ~isempty(jj)
    ii(fi(jj,[1 2])) = 0;
end

if ~all(ii)
    ii = find(ii);
     s =  s(ii,:);
    ok = ok(ii);
end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fi,fn] = fold_ind(ok);

% FolderNumber 
fn = cumsum( ok .* ( abs(ok) == 1 ) ) + ( ok == -1 );

nf = sum( ok == 1 );  % Number of Folders

fi = zeros(nf,3);     % [ StartIndex  EndIndex  Number ]

if nf == 0
   fn = zeros(size(ok));
   return
end

fi(:,1) = find( ok == 1 );
fi(:,3) = fn(fi(:,1));

fk = fn .* ( ok == -1 );

for ii = 1 : max(fn)
       jj    = find( fi(:,3) == ii );
    fi(jj,2) = find( fk == ii );
end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = get_struct(s,mrk,l0,u0,u1)

c = struct( 'isdir' , {  0 } , ...
            'link'  , { '' } , ...
            'type'  , { '' } , ...
            'name'  , { '' }       );

if isempty(mrk)
   c = c(ones(0,1));
   return
end

[fi,fn] = fold_ind(mrk);

ii = find( fn == min(fn) );

n = size(ii,1);

c = c(ones(n,1));

ok = ones(n,1);

nl = size(l0,2) + 1;
 
for jj = 1 : n

  kk = ii(jj);

  b = rmblank(s(kk,nl:end));

  if ~isempty(b)
      nm = b;
      [m,r,k] = get_ref(nm,'<','>');
      for r = size(k,1) : -1 : 1
          nm(grp2ind(k(r,1),k(r,2))) = [];
      end
      c(jj).name = nm;
  end

  switch mrk(kk)

    %------------------------------------
    case 3
    % URL

      [m,r] = get_ref(b,u0,u1);
      if isempty(r)
         [m,r] = get_ref(b,lower(u0),u1);
      end
      if ~isempty(r) 
          c(jj).link = r{1,2};
          k = find( r{1,2} == ':' );
          if ~isempty(k)
              k = k(1) - 1;
              c(jj).type = r{1,2}(1:k);
          end
      end


    %------------------------------------
    case 2
    % Folder

       c(jj).isdir = 1;
       c(jj).link  = c(ones(0,1));
       c(jj).type  = 'folder';

       if ~isempty(fi)
           if any( fi(:,1) == kk+1 )
              ll = find( fi(:,1) == kk+1 );
              if fi(ll,2) > fi(ll,1)+1
                 ind = ( fi(ll,1)+1 : fi(ll,2)-1 );
                 c(jj).link = get_struct(s(ind,:),mrk(ind),l0,u0,u1);
              end
           end
       end
         
    %------------------------------------
    otherwise

        ok(ii) = 0;

  end

end

if ~all(ok)
    c = c(find(ok));
end
