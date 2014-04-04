function [file,msg,ok] = dbfile(ident,chk);

% DBFILE  Search for DataFiles
%
% [Files,Msg] = DBFILE( IDENT )
%
% IDENT is the identifer for DataFiles, a String or CellArray of Strings.
% Each String can be a individual FileName of a DataFile or 
%  a coded Identifer with the Separator by DBGET/DBSET,
%  default is the Colon-Character ":".
% 
%   Coded Ident: [Base]:[Name]:[Type]:[Numbers]
%                [Base]:[Name]:[Type]:[Pattern]
%
% Base, Name, Type and Numbers or Pattern can be Elements of a CellArray,
%  which can have up to 4 Columns. Use this form for numeric Values of Numbers.
%
%  Example:  { '[Base]:[Name]:[Type]'  [Numbers] }       2 Columns
%
%            { '[Base]' '[Name]' '[Type]'  [Numbers] }   4 Columns
%
%-----------------------------------------------------------------------------
%
% [Files,Msg,Ok] = DBFILE( IDENT , CHECK )
% 
% In case CHECK == 1 , DBFILE checks for individual Files only, 
%    Ok == 1 if they exist, 0 if not and NaN if invalid Path.
%
% Otherwise only existing Files are returned.
%
%-----------------------------------------------------------------------------
%
% DBFILE search or check for Files which full name match the syntax:
%
%  Root/[Base]/Name/Type/Name_???.###
%  Root/[Base]/Name/Name_???.###       (if empty Type)
%
%    where "???" matches the Numbers, 
%     the extension "###" should match Type
%     or any like "dat" "edt" "asc" etc. (see Extensions by DBSET)
% 
%  Root/[Base]/Name/Type/NamePattern
%  Root/[Base]/Name/NamePattern        (if empty Type)
%
%  Root/[Base]/Name/Pattern            (if FileSep in Pattern)
%
% Note: The WildCard "*" is allowed in Base, Name, Type and Pattern
%          for multiple matches. 
%       In case the string of Base contains the FileSeparator 
%          the Root-Directory is ignored!
%       In case the string for Pattern contains the FileSeperator
%          the Name is ignored in FileName:  Root/[Base]/Name/Typ/Pattern 
%
% Directories for Base, Name and Type are checked 
%  with all lowercase and all uppercase characters.
%
% The Root-Directory should be set by DBROOT.
% 
% The Separator for coding (":"), the Spacer between Name and Number ("_"), 
%  Digits for the Number and the default extensions to check 
%  can defined by DBSET.
%
%-----------------------------------------------------------------------------
% 
% see also: DBSCAN, DBROOT, DBSET, DBGET, FILESEP, FULLFILE
%


%-------------------------------------------------------
% Settings from Configuration

[dlm,ext,sep,dgt] = dbget('Separator','Extension','Spacer','Digits');

root = dbroot;

frm = cell(1,prod(size(dgt)));

for ii = 1 : size(frm,2)
    frm{ii} = sprintf('%%%.d.%.dd',dgt([ii ii]));  % Formats for FileNumbering
end

% NAME_###.EXT

c     = computer;
iswin = strcmp(upper(c(1:2)),'PC');

fs = filesep;  % Multiple allowed  !!!

%-------------------------------------------------------
% Check for RootPath

if chkstr(root,1)
   if any( root == dlm )
      root = sepname(root,NaN,dlm);
   end
end

[ok,root] = chkcstr(root);
if ok
   root(find(strcmp(root,''))) = [];
else
   root = {};
end

if isempty(root)
   root = {pwd};           % Use Current Directory
end
       
%-------------------------------------------------------
% Check ident

Nout = nargout;

file = {};
msg  = '';
ok   = [];

if isempty(ident)
   return
end

[msg,ident] = chkident(ident,dlm,{fs '.'});

if ~isempty(msg)
    if Nout < 2
       dbmsg(sprintf('%s\n',msg),1);
    end
    return
end

if isempty(ident)
   return
end

if nargin < 2
   chk = 0;
end

chk = isequal(chk,1);

if chk
   frm = frm(1);
   ext = ext(1);
end

%-------------------------------------------------------
% Check for Files

pp = pwd;

nid = size(ident,1);

ok = zeros(2,nid); % [ Ident ; FileOk ]

for id = 1 : nid

    f = ident{id,1};

    ok(1,id) = any( ( f == dlm ) | ( f == '*' ) | ~isempty(ident{id,2}) );

    if ~ok(1,id)

       [p,n,e] = fileparts(f);

       ok(2,id) = ~isempty(e);

       if ok(2,id)

          if isempty(p) & chk

             f = fullfile(pp,f);

             ok(2,id) = -1 + 2 * ( exist(f,'file') == 2 );

          else

             ok(2,id) = ( exist(f,'file') == 2 );

             if ok(2,id)
                ff = which(f);
                if ~isempty(ff)
                    f = ff;
                end
             end

             if ~( ok(2,id) | isempty(p) )
                 ok(2,id) = -1 * ( exist(p,'dir') == 7 );
             end

          end

       end

       ident{id,1} = {f};

    end

end

%-------------------------------------------------------
% Good Files

coded = ok(1,:);

ok = double(ok(2,:));

if chk
   gd = ~( ok == 0 );
   ok(find(ok== 0)) = NaN;
   ok(find(ok==-1)) = 0;   % Path exist
else
   gd = ( ok == 1 );
end

%-------------------------------------------------------
% No coded Ident

if ~any(coded)

    if ~any(gd)

        file = {};
          ok = [];
         msg = 'Can''t find Files.';

    else

       if ~( all(gd) | chk )
               gd =  find(gd);
            ident = ident(gd,:);
               ok =    ok(gd);
       end

       file = cat(1,ident{:,1});

       [h,si] = sort(file);
        h     = diff(double(char(h)),1,1);
        h     = ( sum(h==0,2) == size(h,2) ); 
       if any(h)
              h  = si( find(h) + 1 );
         file(h) = [];
           ok(h) = [];
             msg = 'Duplicate Files.';
       end

    end

    if ~isempty(msg) & ( Nout < 2 )
        dbmsg(sprintf('%s\n',msg),1);
    end

    ok = ok(:);

    return

end

%-------------------------------------------------------

file      = cell(nid,2);
file(:,1) = { cell(0,1) };
file(:,2) = { NaN };

if chk
   ind = ( 1 : nid );
else
   ind = find(gd);
end

for id = ind
    [h,ff] = chkcstr(ident{id,1});
    file{id,1} = ff;
    file{id,2} = ok(id);
end

msg = cell(0,1);

for id = find(coded)

    %------------------------------------------
    % Get Directories

    [pfd,str,bnt] = getpfad(root,ident{id,1},dlm,iswin,fs);

    m = '';

    if isempty(pfd)
       m = 'Can''t find Directory. ';
    end

    %------------------------------------------
    % Get File Numbers from String

    nn = cell(1,0);

    nr = [];   % Number from STR

    %%% if ( str(1) == dlm ), str = str(2:end); end

    if ~isempty(str)
        str = rmblank(str,2);
    end

    if isequal(str,':')

       nn = { [] }; 

    elseif ~isempty(str)

        mm = '';

        try
           nr = eval( str , cat(2,'[',str,']') );
        catch
           mm = 'Invalid Expression. ';
        end

        if ~isempty(nr)
            ok = isnumeric(nr);
            if ok
               nr = nr(:);
               ok = all( ( mod(nr,1) == 0 ) & ( nr >= 0 ) );
            end
            if ~ok
                mm = 'Invalid FileNumbers. ';
            end
        end

        if ~isempty(mm)
            % Check for Pattern: Letter, Number or WildCard , '_' ,  '.' 
            ok = ( isletter(str) | ( ( '0' <= str ) & ( str <= '9' ) ) | ...
                  ( str == '_' ) | ( str == '*' ) | ( str == '.' ) );
            if any(ok)
               for c = fs
                   ok = ( ok | ( str == c ) );
               end
            end
            if all(ok)
               nr = str;
               mm = '';
            else
                m = cat( 2 , m , mm );
            end
        end

        if isempty(mm)
           nn = {nr};
        end

    end

    %------------------------------------------
    % Append FileNumbers from numeric IdentCell

    nr = ident{id,2};

    if ~isempty(nr) | isempty(nn)
        nn = cat( 2 , nn , {nr(:)} );
    end

    %------------------------------------------------------------------
    % Get File Names for each Path 

    if isempty(m)

       file{id,1} = char(ones(0,1));
       file{id,2} =      ones(0,1);

       for nr = nn

           use_typ  = ~isempty(bnt{3});
           use_name =  ischar(nr{1});
           if use_name
              for c = fs
                  if any( nr{1} == c )
                     use_typ = -1;
                     break
                  end
              end
           end

           nd = size(pfd,1);
        
           ff = cell(nd,1);
           ff(:) = {cell(0,1)};

           for ii = 1 : nd

               for form = frm            
                   ff{ii} = getfile(pfd{ii},nr{1},use_typ,chk,ext,sep,form{1},iswin);
                   if ~isempty(ff{ii})
                       break
                   end
               end

           end

           ff = cat( 1 , ff{:} );
        
           if isempty(ff)
              if isempty(m), m = 'Can''t find Files.'; end
           else
              file{id,1} = cat( 1 , file{id,1} , ff);
              if chk & ( Nout == 3 )
                 nf = prod(size(ff));
                 ok = ones(nf,1);
                 for ii = 1 : nf
                     ok(ii) = ( exist(ff{ii},'file') == 2 );
                 end
                 file{id,2} = cat( 1 , file{id,2} , ok );
              end
           end

       end

       if isempty(file{id,1}) & chk ,  file{id,1} = ident{id,1}; end
       if isempty(file{id,2}) ,  file{id,2} = NaN; end

    end

    if ~isempty(m)
        m = sprintf('%s | %s',ident{id,1},m);
        msg = cat(1,msg,{m});
    end
    
end

  ok = cat( 1 , file{:,2} );
file = cat( 1 , file{:,1} );

if ~isempty(file)
    [h,si] = sort(file);
     h = diff(double(char(h)),1,1);
     h = ( sum(h==0,2) == size(h,2) );
     if any(h)
        h = si( find(h) + 1 );
        file(h) = [];
        if chk
           ok(h) = [];
        end
        msg = cat( 1 , msg , {'Duplicate Files.'} );
     end
end

if ~chk | isempty(file) 
    ok = [];
end
     
if isempty(msg)
   msg = '';
else
   msg = sprintf('%s\n',msg{:});
   if Nout < 2
      dbmsg(msg,1);
   end
end


%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function f = getfile(pfd,nr,use_typ,chk,ext,sep,frm,iswin)

% Get RODB-Files in Pfad

p = sepname(pfd,NaN,filesep,-1);

f = cell(0,1);

if prod(size(p)) < 2
   return
end

typ  = p{1};
name = p{2};

%---------------------------------------------------------------
% Check with Numbers and Extensions

n1 = ischar(nr);   % String

n0 = ( ~n1 & isempty(nr) );

if n1
   sep = '';
end

n  = max(1,size(nr,1));

nsp = size(sep,2);
nnr = size(sprintf(frm,0),2);

f = cell(n,1);

if ~use_typ
    name = typ;
end

nnm = size(name,2) * double( n0 );

nome = {name};

if iswin
   nome = lower(nome);
elseif use_typ == -1
   nome = { '' };
else
   for nm = { lower(name)  upper(name) }
       if ~strcmp( name , nm{1} )
           nome = cat( 2 , nome , nm );
       end
   end
end

if use_typ == 1
   ext = [ {typ}  ext ];
end

if n1
   if any( ( nr == '.' ) | ( nr == '*' ) )
      ext = {''};
   else
      ext = cat( 2 , ext ,{''} );
   end
end

if iswin
   ext = lower(ext);
else
   ext = cat( 1 , ext , lower(ext) , upper(ext) );
   ok      = ~strcmp( ext , ext([1 1 1],:) );
   ok(1,:) = 1;
   ext = ext(find(ok));
   ext = ext(:)';
end

for nm = nome

    for e = ext

        for ii = 1 : n
        
            if     n0
               if chk
                  m = frm;
               else
                  m = '*';
               end
            elseif n1
               m = nr;
            else
               m = sprintf(frm,nr(ii));
            end
       
            ff = sprintf('%s%s%s',nm{1},sep,m);

            if ~isempty(e{1})
                ff = sprintf('%s.%s',ff,e{1});
            end

            if chk
               f{ii} = {fullfile(pfd,ff)};
            else
               f{ii} = chkfile(pfd,ff,nnm,nsp,nnr,sep);
            end

        end
    
        ff = cat(1,f{:});
    
        if ~isempty(ff)
            break
        end

    end

    if ~isempty(ff)
        break
    end

end

f = ff;

%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function f = chkfile(pfd,ff,nnm,nsp,nnr,sep)

f = cell(0,1);

%----------------------------------------------
% Check for File

one = ~any( ff == '*' );

ff = fullfile( pfd , ff );

if one
   if exist(ff,'file') == 2
      f = {ff};
   end
   return
end

%----------------------------------------------
% Inquire for Files

pfd = fileparts(ff);

d  =  dir(ff);
ok = ~isempty(d);
 
if ok
   jj = ~( cat(1,d.isdir) | ( cat(1,d.bytes) == 0 ) );
   ok = any(jj);
   if ok
      if ~all(jj)
              jj = find(jj);
              d  = d(jj);
      end
   end
end

if ~ok
    return
end

f  = {d.name};
f  =  f(:);
n  =  size(f,1);

ok = ( nnm == 0 );

if ~ok
    f = sort(f);
end

ok = ok(ones(n,1));

nf = nnm + nsp + nnr;   % Length of Name without Extension

for ii = 1 : n
    [pn,ff,e] = fileparts(f{ii});
    if ~ok(ii)
        ok(ii) = ( size(ff,2) >= nf ) & ( size(e,2) > 1 );
        if ok(ii)
           ok(ii) = strcmp( ff(nnm+(1:nsp)) , sep  );
           if ok(ii)
              fn = ff( nnm+nsp+1 : end );
              ok(ii) = all( ( '0' <= fn ) & ( fn <= '9' ) );
          end
       end
    end
    if ok(ii)
       f{ii} = fullfile(pfd,[ff e]);
    end
end

if ~all(ok)
    if any(ok)

       ok = find(ok);
       f  = f(ok);
    else
       f = cell(0,1);
    end
end

%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [pfd,str,bnt] = getpfad(pfd,str,dlm,iswin,fs);

% Find SubDirectories for PFD,
%  matching Identifer of STR
%

bnt    = cell(1,3);
bnt(:) = {''};

if ~isempty(str)
    str = rmblank(str,2);
end

if isempty(str)
   return
end

%----------------------------------------------
% Split String into [ Base : Name : Type ]

if ~any( str == dlm )

    [p,n,e] = fileparts(str);

    str = [ n  e ];

    if ~isempty(p)
        for ii = [ 3  2 ]
            [p,n,e] = fileparts(p);
            bnt{ii} = [ n  e ];
            if isempty(p)
               break
            end
        end
        bnt{1} = p;
    end
    
else

    for ii = 1 : 3

        if ( str(1) == dlm )
           str = str(2:end);
        else
           [sub,str] = strtok(str,dlm);
           if ~isempty(sub)
               sub = rmblank(sub,2);
           end 
           bnt{ii} = sub; 
           if ~isempty(str)
               if ( str(1) == dlm )
                  str = str(2:end);
               end
           end
        end

        if isempty(str)
           break
        end

    end

end

if ~isempty(bnt{1})
    ok = ~isempty(fs);
    if ok
       for c = fs
           ok = any( bnt{1} == c );
           if ok
              break
           end
       end
    end
    if ok
       pfd = bnt(1);
       bnt{1} = '';
       [p,n,e] = fileparts(pfd{1});
       n = [ n  e ];
       if any( n == '*' )
          bnt{1} = n;
          pfd{1} = p;
       end
    end
end

%----------------------------------------------
% Check for single Directory in ROOT

if 0 %%% ~any( strcmp( bnt , '*' ) )

    bnt = { fullfile(bnt{:}) };
    
end

%----------------------------------------------
% Search for multiple Directories in ROOT

for sub = bnt
    
    name = sub;

    if ~isempty(sub{1})
        if iswin
           name = lower(name);
        else
           for nm = { lower(sub{1})  upper(sub{1}) }
               if ~strcmp( sub , nm{1} )
                   name = cat( 2 , name , nm );
               end
           end
        end
    end

    n = prod(size(pfd));

    ok  = zeros(n,1);

    for ii = 1 : n
        pfd{ii} =   chkdir(pfd{ii},name);
         ok(ii) = ~isempty(pfd{ii});
    end

    if ~all(ok)
        if all(~ok)
           pfd = {};
           break
        else
           ok = find(ok);
           pfd = pfd(ok);
        end
    end

    pfd = cat( 1 , pfd{:} );
    
end

%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function pfd = chkdir(base,name)

% Check for SubDirectories of Pfade

pfd = cell(0,1);

if isempty(name{1})
   pfd = {base};
   return
end

for sub = name

    d = fullfile(base,sub{1});

    if any( sub{1} == '*' )

       d = getdir(d);

       if ~isempty(d)
           p = fileparts(fullfile(base,sub{1}));
           for ii = 1 : size(d,1)
               d{ii} = fullfile(p,d{ii});
           end
           pfd = cat( 1 , pfd , d );
       end

    else

       if exist(d,'dir') == 7
          pfd = cat(1,pfd,{d});
       end

    end

end

%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = getdir(pfd)

% Get SubDirectoryNames of Pfad

d = dir(pfd);

if isempty(d)
   d = {};
   return
end

ok = cat(1,d.isdir);

d  = {d.name};
d  = d(:);

ok = ( ok & ~( strcmp(d,'.') | strcmp(d,'..') ) );

if     ~any(ok)
        d = {};
elseif ~all(ok)
        ok = find(ok);
        d = d(ok);
end

