function display(obj,varargin)

% DISPLAY  Display for Object
%
% DISPLAY( Object )
%
%   Display the ObjectValue
%
% DISPLAY( Object , FieldNames )
%
%   Displays only the Fields from the CellStringArray FieldNames
%
% DISPLAY( Object , FieldName1 , FieldName2 , ... )
%
%
% calls: @class/disp.m
%

v = varargin;

name = inputname(1);

if isempty(v)

   eval([name ' = getfield(obj,''VALUE'')']);
 
   return

end

if prod(size(v)) == 1

    if isequal(v{1},class(obj))

       v = {};

    end

end

[txt,c,msg] = disp(obj,v{:});

if ~isempty(msg)
   error(msg)
end

nl = char(10);

   name = inputname(1);

   if ~isempty(name) & ( nargin == 1 )
     fprintf(nl)
     fprintf('%s = ',name);
     fprintf(nl);
   end

   fprintf(nl);
   fprintf(strrep(strrep(txt,'%','%%'),'\','\\'));
   fprintf(nl);

