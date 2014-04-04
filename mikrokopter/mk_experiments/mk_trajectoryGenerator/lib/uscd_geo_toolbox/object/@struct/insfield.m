function t = insfield(s,f1,f2)

% INSFIELD  Insert structure field.
%
%   S = INSFIELD(S,Field,{InsertFields}) inserts the fields
%        in the structure array S behind Field.  
%
%   S = INSFIELD(S,{InsertFields},Field) inserts the fields
%        in the structure array S before Field.  
%
%   Field is a CharacterString of a valid FieldName of S.
%   {InsertFields} is a CellStringArray of FieldNames to insert.
%
%   The changed structure is returned.
%
%   See also RMFIELD, RNFIELD, SETFIELD, GETFIELD, FIELDNAMES

%   C.Begler IFM-GEOMAR Kiel

if nargin < 3
   error('Not enough InputArguments.');
end

if ~isa(s,'struct'), error('S must be a structure array.'); end

ins = ( ischar(f1) & iscellstr(f2) ) - ( ischar(f2) & iscellstr(f1) );

if ins == 0
   error('Field must be a String, InsertFields a CellString.');
end

if ins == -1
   f3 = f2;
   f2 = f1;
   f1 = f3;
end

if isempty(f2)
   t = s;
   return
end


f0 = fieldnames(s); % Org
f3 = f0;            % New

ok = strcmp( f0 , f1 );

if ~any(ok)
    error(sprintf('Fieldname "%s" not exist.',f1));
end

ii = find(ok) - ( ins == -1 );
jj = ii + 1;

n = size(f3,1);
m = prod(size(f2));

ii = (  1 : ii );
jj = ( jj :  n );

f3 = cat( 1 , f3(ii) , f2(:) , f3(jj) );

ok = zeros(1,n+m);

ok(ii)   = ii;
ok(m+jj) = jj;

f = double(char(sort(f3)));
f = ( sum( diff(f,1,1) == 0 , 2 ) == size(f,2) );
if any(f)
   jj = find(f);
    f = sort(f3);
    m = sprintf(' "%s" ',f{jj});
   warning(sprintf('Duplicate FieldNames.\n%s',m));
end


% Now rename the fieldnames by copying the ones we want to keep.
t = [];
for ii = prod(size(s)) : -1 : 1,
    for jj = 1 : prod(size(f3)),
        if ok(jj)
           v = getfield(s,{ii},f0{ok(jj)});
        elseif any(strcmp(f3{jj},f0))
           v = getfield(s,{ii},f3{jj});
        else
           v = [];
        end
         t = setfield(t,{ii},f3{jj},v);
    end
end

if ~isempty(t)
    t = reshape(t,size(s));
end
