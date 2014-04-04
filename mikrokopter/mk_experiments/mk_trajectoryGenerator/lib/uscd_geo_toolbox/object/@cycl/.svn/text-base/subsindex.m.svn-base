function v = subsindex(obj)

% SUBSINDEX  Subsindex for CYCL-Objects
%

obj = obj.VALUE;

ok = isnumeric(obj);
if ok & ~isempty(obj)
   ok = all( isfinite(obj) & ( obj >= 1 ) );
end

if ~ok
   error('Subscript indices must be finite numeric values, greater or equal 1.')
end

v = fix(double(obj)) - 1;    % ZERO-based !!!
