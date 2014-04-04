function translate(Pfad,str1,str2,varargin)

% TRANSLATE  translates String in Files to another 
%
% TRANSLATE(Path,str1,str2,'-r',FileSpec,LogID)
%
% Translate in all Files in Path matching FileSpec the
%   "str1" by "str2", works recursivly if Option '-r' is set.
%
% "str1" and "str2" can be CharacterArrays or CellStringArrays.
%
% Use CellStringArray for single ControlCharacters!
%
% TRANSLATE(File,str1,str2,'-r',LogID)
%
% The File to translate will defined by the Input(s) FileSpec.
%  It can be a String, starting with "." to define a Extension,
%         or a String, contains the WildCard "*".
%
% default: FileSpec = '.m'
%
% Multiple Inputs for FileSpecifications are allowed.
%
% LogID can be a FileIdentifer to log the Messages.
%
% default: LogID = 1   LogOutput to Matlab-CommandWindow
%

VarArg = varargin;
VarArg = VarArg(:);

Nin = size(VarArg,1);

nl = char(10);

%*******************************************************
% Check Inputs

if isempty(Pfad)
   Pfad = cd;
end

msg = '';
if ~chkstr(Pfad,1)
    msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
               'Pfad or File must be nonempty Strings.' );
end

[ok1,str1] = chkcstr(str1,0);
[ok2,str2] = chkcstr(str2,0);
if ~( ok1 & ok2 )
    msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
               'str1 and str2 must be Strings or CellArray of Strings.' );
end

str1 = str1(:);
str2 = str2(:);

n1 = size(str1,1);
n2 = size(str2,1);

if ~( ( size(str1,1) == size(str2,1) )  |  ( n1 == 1 )  |  ( n2 == 1 ) )
    msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
               'str1 and str2 must have the same Number of Strings.' );
end

if isempty(msg)

   if ( n1 == 1 ) & ~( n2 == 1 )
      str1 = str1(ones(n2,1));
   end

   if ( n2 == 1 ) & ~( n1 == 1 )
      str2 = str2(ones(n1,1));
   end

   nn = max(n1,n2);

end

%*******************************************************
% Check Options

recurse = 0;
ext     = cell(0,1);
lid     = 1;

for ii = 1 : Nin
  if chkstr(VarArg{ii},1)
     if strcmp(VarArg{ii}(1),'.');
        ext = cat( 1 , ext , {cat(2,'*',VarArg{ii})} );
     elseif any( VarArg{ii} == '*' )
        ext = cat( 1 , ext , VarArg(ii) );
     else 
        recurse = ( recurse | strcmp(VarArg{ii},'-r') );
     end
  elseif isnumeric(VarArg{ii})
     lid = VarArg{ii};
  end
end

if ~isempty(msg)
   logmsg(lid,msg,1)
   return
end

if isempty(ext)
   ext = {'*.m'};
end

%*******************************************************
% Check Pfad

if ~( exist(Pfad,'file') == 2 )
    if any(strcmp(Pfad,{ '.'  '*' }))
       Pfad = cd;
    elseif ~isempty(findstr(Pfad,'..'))
       p0 = cd;
       cd(Pfad);
       Pfad = cd;
       cd(p0);
    end
end

logmsg(lid,Pfad,0);

d = dir(Pfad);

if isempty(d)
  logmsg(lid,' ... can''t read File or Directory ',1);
  return
end

IsFile = ( (size(d,1)==1)  &  (d(1).isdir==0) );
if  IsFile
    if isempty(fileparts(d(1).name))
       d(1).name = fullfile(fileparts(Pfad),d(1).name);
    end
else
    if ~strcmp(Pfad(length(Pfad)),filesep)
        Pfad = [ Pfad filesep ];
    end
end

%*******************************************************

logmsg(lid,'',1)

%*******************************************************
% Find Files to Translate

ok = IsFile;
dd = d;

if ~ok
    is_file = ~cat(1,d.isdir);
       ok   = any(is_file);
    if ok
       is_file = find(is_file);
       d       = d(is_file);
       is_cmp  = strwcmp( {d.name} , ext );
          ok   = any(is_cmp);
       if ok
          is_cmp = find(is_cmp);
          d      = d(is_cmp);
       end
   end
end

%*******************************************************
% Translate Files

if ok

   for ii = 1 : size(d,1)

       name = d(ii).name;

       file = [ Pfad((1:(end*(~IsFile))))  name ];
       fid = fopen(file,'r');
       if fid == -1
          logmsg(lid,['   Error read '  file ],1);
       else
          bb = fread(fid,'char');
          fclose(fid);
          if any(bb==0)
             logmsg(lid,['   Invalid Character in '  file ],1);
          end
             bb = char(bb(:)');
             ok = zeros(nn,1);
             for jj = 1 : nn
                 ok(jj) = ~isempty( findstr( bb , str1{jj} ) );
                 if ok(jj)
                    bb = strrep(bb,str1{jj},str2{jj});
                 end
             end
             if any(ok)
                fid = fopen(file,'wt');
                if fid == -1
                   logmsg(lid,['   Error open to write '  file ],1);
                else
                   fwrite(fid,bb,'char');
                   fclose(fid);
                end
             end
       end 
       % fid ok
   end
   % ii
end
% ok

%*******************************************************
% Recurse

if recurse

  is_dir = find(cat(1,dd.isdir));

  for ii = is_dir(:)'

    if ~any(strcmp(dd(ii).name,{ '.'  '..' }))

      translate([Pfad dd(ii).name],str1,str2,'-r',ext{:});

    end

  end

end      

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function logmsg(lid,msg,nl);

if ~chkstr(msg,0)
   return
end

if isempty(lid)
   return
end

msg = cat( 2 , msg , char(10*ones(1,nl)) );
   
for ll = lid(:)'
    fprintf(ll,'%s',msg);
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
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
