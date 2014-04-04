function v = str2cell(v)

% STR2CELL  Tries to Convert a String to CellArray
%
%  C = STR2CELL( String )
%
%  Try to build CellArray using EVAL(['{' String '}'])
%
%  In case of multiple Lines, separated by LF / CHAR(10)
%  a multiple Row CellArray for the nonempty Lines is returned.
%
% requires: RMBLANK
%

if ~chkstr(v,1) | isempty(v)
    return
end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove Ending Blanks, NewLines etc.

    v = strrep(v,char( 0),' ');
    v = strrep(v,char( 9),' ');
    v = strrep(v,char(13),' ');

    v = rmblank(v,2);

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove ending "," ";"

   while any( v(end) == ',;' )
      if size(v,2) == 1
         v = '';
      else
         v = v( 1 : (end-1) );
         v = rmblank(v,2-i);
      end
      if isempty(v)
         break
      end
   end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove surrounding "{}"

    if ~isempty(v)
        if v(1) == '{'
           if size(v,2) == 1
              v = '';
           else
              v = v( 2 : end );
              v = rmblank(v,2+i);
           end
        end
    end
    if ~isempty(v)
        if v(end) == '}'
           if size(v,2) == 1
              v = '';
           else
              v = v( 1 : (end-1) );
              v = rmblank(v,2-i);
           end
        end
    end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Split by Lines

ok = ( v == 10 );

if ~any(ok)
    return
end

[i0,lg] = ind2grp(find(~ok));

ng = size(i0,1);

vv =  cell(ng,1);
nn = zeros(ng,1);

for ii = 1 : ng

    jj = ( 1 : lg(ii) ) + i0(ii) - 1;

    vv{ii} = rmblank(v(jj),2);

    if ~isempty(vv{ii})

        vv{ii} = str2cell(v(jj));

        nn(ii) = iscell(vv{ii});

        if nn(ii)
           if isempty(vv{ii})
              nn(ii) = -1;
           else
              vv{ii} =      vv{ii}(:)';
              nn(ii) = size(vv{ii},2);
           end
        end

    end

end

%-----------------------------------------------------
% Get Good Lines

ok = ~( nn == 0 );

if ~any(ok)
    return
end

n  = sum(ok);
ok = find(ok);

nn = nn(ok);
vv = vv(ok);

m  = max(nn);

if m == -1
   v = cell(n,0);
   return
end

v = cell(n,m);

for ii = 1 : n
    if ~( nn(ii) == -1 )
        v(ii,1:nn(ii)) = vv{ii};
    end
end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [i0,l] = ind2grp(ii);

% IND2GRP  Built StartIndex and Length from IndexVector
%
% [ StartIndex , GroupLength ] = IND2GRP( Index )
%

i0 = zeros(0,1);
l  = zeros(0,1);

if isempty(ii);
   return
end

ii = ii(:);
n  = size(ii,1);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , n+1 );

l  = diff(i0,1,1);

i0 = ii(i0(1:end-1));

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

