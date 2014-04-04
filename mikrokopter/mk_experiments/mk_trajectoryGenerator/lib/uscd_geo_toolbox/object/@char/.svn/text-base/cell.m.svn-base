function c = cell(s)

% CELL  Create cell array of strings from rows of character array
%
% CellArray = CELL(CharacterArray)
%
% CELL does NOT remove trailing whitespace characters from strings
%  of CharacterArray, like CELLSTR it does.
%
% see also: CHAR, CELLSTR
%

si = size(s);

m = si(2);

si(2) = min(m,1);

c = cell(si);

n = prod(si);

if n == 0 
   return
end

if ~( n == si(1) )

     s = permute(s,[2 1 (3:ndims(s))]);

     s = reshape(s,m,n);

     s = permute(s,[2 1]);

end

for ii = 1 : n
    c{ii} = s(ii,:);
end

