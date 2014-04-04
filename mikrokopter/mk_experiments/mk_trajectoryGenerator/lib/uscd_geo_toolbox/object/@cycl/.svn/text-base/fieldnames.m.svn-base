function [f,p,s] = fieldnames(obj);

% FIELDNAMES Returns FieldNames of an Object
%
% [ FieldNames , ParentClassName ] = FIELDNAMES( Object )
%
%  returns a CellStringArry of the FieldNames of the ObjectStructure
%   and the Name of ParentClass an Object.
%
% [ ... , ObjectStructure ] = FIELDNAMES( Object )
%
%  returns additionaly the ObjectStructure
%  

p = '';

if isempty(struct(obj))
   obj = feval(class(obj));
end

s = struct(obj);
s = s(1);

f = fieldnames(s);

f = f(:);

n = size(f,1);

e = getfield(s,f{n});

if isobject(e)

   if isa(obj,class(e))

      p = f{n};

      f = f( 1 : (n-1) );

   end

end