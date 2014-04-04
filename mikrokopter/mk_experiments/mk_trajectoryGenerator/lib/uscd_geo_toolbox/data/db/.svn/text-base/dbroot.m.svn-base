function root = dbroot(pfad,opt)

% DBROOT set RootPath for DataBaseAccess
%
% DBROOT( Path , Option )
%
% Option = 'new'      Overwrites   old RootPath (default)
%          'begin'    Prepend  to  old RootPath
%          'end'      Append   to  old RootPath
%          'remove'   Removes from old RootPath
%
% DBROOT without any Input returns or displays the actual RootPath.
%
% see also: DBSET, DBGET
%

Nin  = nargin;
Nout = nargout;

dlm = dbget('Separator');

chk = ( Nin == 0 );

if  chk

    new = 0;
    app = 0;
    del = 0;

else

    if Nin < 2
       opt = 'new';
    end

    if ~chkstr(opt,1)
        error('Option must be a String.')
    elseif ~any(lower(opt(1)) == 'nber');
        error('Invalid Option.');
    else
        opt = lower(opt(1));
        new = ( opt == 'n' );
        del = ( opt == 'r' );
        app = any( opt == 'be' );
    end

end

%--------------------------------------------------------
% Check PathInput

if chk

   pfad = {};

else

   if isempty(pfad)  | isequal(pfad,{''})

      if app
         pfad = {};
      else
         root = {''};
         msg = dbset('Root',root);
         return
      end

   else

      [ok,pfad] = chkcstr(pfad);

      if ~ok
          error('Path must be a CharacterArray or CellArray of Strings.');
      end

      [pfad,pok] = chkpfd(pfad,dlm);

      msg = '';

      if isempty(pfad)
         msg = 'No valid Name for Path.';
      else
         if ~( all(pok) | del )
             jj = find(~pok);
             msg = sprintf('Invalid Path: %s\n',pfad{jj});
             if ~any(pok)
                 pfad = {};
             else
                 pfad = pfad(find(pok));
             end
         end
      end

      if ~isempty(msg)
          dbmsg(msg,1);
      end

   end

end

%--------------------------------------------------------
% Check Existing RootPath

if new

   root = {};

else

   root = dbget('Root');

   if isempty(root) | isequal(root,{''})
      
      root = {};

   elseif ~chkcstr(root,1)

      msg = 'Invalid RootPath.';

   else

      [root,rok] = chkpfd(root,dlm);

      msg = '';

      if isempty(root)
         msg = 'No valid Name for Root.';
      else
         if ~all(rok)
             jj = find(~rok);
             msg = sprintf('Invalid Root: %s\n',root{jj});
             if ~any(rok)
                 root = {};
             else
                 root = root(find(rok));
             end
         end
      end

      if ~isempty(msg)
          dbmsg(msg,1);
      end

   end

end

%--------------------------------------------------------
% Set RootPath

if del & ~( isempty(root) | isempty(pfad) )
   ok = zeros(size(root));
   for ii = 1 : size(root,1)
       ok(ii) = ~any( strcmp( pfad , root{ii} ) );
   end
   if any(~ok)
      if all(~ok)
         root = {};
      else
         root = root(find(ok));
      end
   end
elseif new
   root = pfad;
elseif app
   if strcmp(opt,'b')
      root = cat( 1 , pfad , root );
   else
      root = cat( 1 , root , pfad );
   end
end
       
if isempty(root)
   root = {''};
else

   % Remove equal Names:

   [h,si] = sort(root);

   h = diff(double(char(h)),1,1);     % Diff of sorted Strings

   h = ( sum(h==0,2) == size(h,2) );  % equal Strings

   h = cat( 1 , 0 , h );              % 1. String unique

   if any(h)

      h  = si(find(h));

      root(h) = [];

   end

end

msg = dbset('Root',{root});

if Nout == 0

   fprintf(1,'\n');
   fprintf(1,'   ''%s''\n',root{:});
   fprintf(1,'\n');

   clear root

end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [p,ok] = chkpfd(p,dlm);

% Check PathStrings with Delimiters and for exist

ok = [];

if isempty(p)
   return
end

p = p(:)';

n = size(p,2);

for ii = 1 : n
    if any( p{ii} == dlm )
       p{ii} = sepname(p{ii},NaN,dlm);
    else
       p{ii} = p(ii);
    end
end

p = cat( 2 , p{:} );

jj = strcmp(p,'');
if all(jj)
   p = {};
   return
elseif any(jj)
   n = n - sum(jj);
   p =  p(find(~jj));
end

p = p(:);

n = size(p,1);

ok = zeros(n,1);

p0 = cd;

for ii = 1 : n 
    ok(ii) = ( exist(p{ii},'dir') == 7 );
    if ok(ii)
       try
          cd(p{ii}), p{ii} = cd;
       catch
          p{ii} = '';
       end
       ok(ii) = ~isempty(p{ii});
    end
end

try, cd(p0), end