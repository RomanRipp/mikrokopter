function  v = obj2val(v)

% OBJ2VAL  Converts Object or Object-CellArrays to their Value
%

c0 = 'cycl'

switch class(v)

  case {c0}

     v = obj.VALUE;

  case 'cell'

    for ii = 1 : prod(size(v))
      if strcmp( class(v{ii}) , 'cycl' )
         v{ii} = v{ii}.VALUE;
      end
    end

end
