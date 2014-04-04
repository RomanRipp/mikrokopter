function out = dbident(ids)

% DBIDENT Decode database variable keywords.
%
%  OUT = DBIDENT(IDENT)
%
%  IDENT: Variable array of DataBase-Variables
%
%          Character- or CellArray of Strings.
%          A single String can have multiple keywords
%            which are separated by the SeperatorCharacter,
%
% The SeparatorCharacter is defined by DBSET
%
%
% OutputStructure:
%
%   ID      : ID of keyword, NaN if not found
%   Name    : Name of keyword (CellString)
%   Format  : Format for keyword (CellString)
%   Length  : DataLength of keyword information
%   Type    : Header = 1 / Data = 2 
%   Value   : Value of keyword, separated from keyword
%                by the AssignmentCharacter
%   Units   : Units of Keyword and Values, separated from keyword
%              by parentheses, braces or brackets or any of
%               '|/\' if no assignment in statement. 
%
% The AssignmentCharacter is defined by DBSET
%
%  Names like "[name]_[number]" will originate to DB-Name
%
%  See DBINIT for more informations about identifers
%
% 

%********************************************
% Settings

global DBKEY DB_KEY 

if isempty(DBKEY)

   [nkey,nfrm,nlen,DBKEY] = dbinit;

   DB_KEY = double(upper(DBKEY(:,1:nkey)));

else

   [nkey,nfrm,nlen] = dbinit;

   % nkey: Length of KeyWord
   % nfrm: Length of FormatString
   % nlen: Length of LengthString

end

% Delimiters for multiple Identifer in a String

[dlm,sep] = dbget('Separator','Assignment'); 

out = struct( 'ID'     , {[]} , ...
              'Name'   , {{}} , ...
              'Format' , {{}} , ...
              'Length' , {[]} , ...
              'Type'   , {[]} , ...
              'Value'  , {{}} , ...
              'Units'  , {{}}       );

%********************************************
% Check Inputs

if nargin < 1
   ids = {};
end

if isempty(ids)
   return
end

%----------------------------------------------------------
% Check for single String, seperated by Delimiter

if chkstr(ids,1)
   if any( ids == dlm )
      ids = sepname(ids,NaN,dlm);
   end
end

%----------------------------------------------------------

[ok,ids] = chkcstr(ids);
if ~( ok | isempty(ids) )
    error('Identifer must be a CharacterArray or CellArray of Strings.');
end

if isempty(ids);
   return
end

%----------------------------------------------------------
% Check Identifers

val = ids(:);  % One Column

uni    = cell(size(val,1),2);
uni(:) = {''};

ids = val';    % One Row !!!

ref = cat( 1 , '|/\([{' , ...
               '|/\)]}'       );

ind = cat( 2 , {( 1 : size(ref,2) )} , ...
               {find( ~( ref(1,:) == ref(2,:) ) )} );

ok = zeros(size(ids));

for ii = find(~strcmp(ids,''))

    jj = ind{ 1 + any(ids{ii}==sep) };

    [txt,val{ii},uni(ii,:)] = splithead(ids{ii},sep,ref(:,jj));

    % Start with Letter !!!
    while ~isempty(txt)
           if isletter(txt(1))  
              break
           end
           txt = txt(2:end);
    end

    ok(ii) = ~isempty(txt);

    ids{ii} = txt;

end

if ~all(ok)
    if all(~ok)
       val = {};
       uni = {};
       return
    else
       ok  = find(ok);
       ids =  ids(ok);
       val =  val(ok);
       uni =  uni(ok,:);
    end
end

%----------------------------------------------------------

nid = size(ids,2);

id  = NaN * ones(nid,1);
frm = char(32*ones(nid,nfrm));
len = zeros(nid,1);

name = char(32*ones(nid,nkey));

ndb = size(DB_KEY,1);

for ii = 1 : nid

    txt = ids{ii};
    str = upper(txt);

    ok = 0; 
    cc = 0;

    while ( ok == 0 ) & (  cc < 3 )

        cc = cc + 1;

        if cc == 2
           str(find(str==' ')) = '_';
        elseif cc == 3
           str(find(str=='_')) = [];
        end

        s2 = size(str,2);

        ok = ( ( 0 < s2 ) & ( s2 <= nkey ) );
        
        if ok
           ok = double(str);
           ok = cat( 2 , ok , 32*ones(1,nkey-s2) );
           ok  = DB_KEY - ok(ones(ndb,1),:);
           ok  = cumprod( double( ok == 0 ) , 2 );
           ok  = sum(ok,2);
           jj  = ( ok == nkey );
           if any(jj)
              % All Characters match
              ok  = cumprod( double(~jj) , 1 );
              ok  = sum( ok , 1 ) + 1; 
           else
              % Min. Characters in String match
              ok  = cumprod( double( ok < s2 ) , 1 );
              ok  = sum( ok , 1 ) + 1; 
              ok  = ok * ( ok <= ndb ); % Take care for No match 
           end 
        end

    end
        
    if ( ok == 0 )

       jj = ~( ( ( '0' <= txt ) & ( txt <= '9' ) ) | ...
               ( ( 'A' <= txt ) & ( txt <= 'Z' ) ) | ...
               ( ( 'a' <= txt ) & ( txt <= 'z' ) )       );

       % Replace bad Characters by '_'
       if any(jj)
              jj  = find(jj);
          txt(jj) = ' ';
          txt = rmblank(txt,2);
          txt =  strrep(txt,' ','_');
       end

       % Remove duplicate '_'
       jj = ( txt == '_' );
       if any(jj)
          jj = find(jj);
          kk = ( diff(jj) == 1 );
          if any(kk)
             kk = find(kk) + 1;
             txt(jj(kk)) = [];
          end
       end

       % Check for Number at End

       nr = char(ones(1,0));
       while ( ( '0' <= txt(end) ) & ( txt(end) <= '9' ) )
              nr = cat( 2 , txt(end) , nr );
             txt = txt(1:(end-1));
             if isempty(txt)
                txt = nr;
                 nr = '';
                break
             end
       end

       if ~isempty(nr)

           sep = ( txt(end) == '_' );

           str = upper(txt(1:(end-double(sep))));
            s2 = size(str,2);

            ok = ( ( 0 < s2 ) & ( s2 <= nkey ) );
        
            if ok
               ok = double(str);
               ok = cat( 2 , ok , 32*ones(1,nkey-s2) );
               ok  = DB_KEY - ok(ones(ndb,1),:);
               ok  = cumprod( double( ok == 0 ) , 2 );
               ok  = sum(ok,2);
               jj  = ( ok == nkey );
               ok  = any(jj);
               if ok
                  % All Characters match
                  jj  = cumprod( double(~jj) , 1 );
                  jj  = sum( jj , 1 ) + 1; 

                  txt = rmblank(DBKEY(jj,1:nkey),2);
                  if sep
                     txt = cat(2,txt,'_');
                  end
                  txt = cat( 2 , txt , nr );
                
                  frm(ii,:)  = DBKEY(jj,nkey+(1:nfrm));
                  len(ii)    = sscanf(DBKEY(jj,nkey+nfrm+(1:nlen)),'%d',1);
                  id(ii)     = sscanf(DBKEY(jj,nkey+nfrm+nlen+1:end),'%d',1);

               end
            end

            if ~ok
                txt = cat( 2 , txt , nr );
            end
           
       end
          
       s2 = min(nkey,size(txt,2));

       name(ii,1:s2) = txt(1:s2);

    else

       name(ii,:) = DBKEY(ok,1:nkey);
       frm(ii,:)  = DBKEY(ok,nkey+(1:nfrm));

       len(ii)    = sscanf(DBKEY(ok,nkey+nfrm+(1:nlen)),'%d',1);

       id(ii)     = sscanf(DBKEY(ok,nkey+nfrm+nlen+1:end),'%d',1);

    end

end

typ = 1 + double( id > 0 );

typ(find(isnan(typ))) = 0;

out.ID     = id;
out.Name   = cellstr(name);
out.Format = cellstr(frm);
out.Length = len;
out.Type   = typ;
out.Value  = val;
out.Units  = uni;

