function obj = subsasgn(obj,s,val)

% SUBSASGN  Subsasgn for CYCL-Object
%
% 
% Note:  use OBJ.VALUE = ...  or  OBJ{:} = ... 
%        to redefine the ObjectValue
%

%---------------------------------------------------------
% Check for: obj{:} | obj{:}

if isequal( s(1) , struct('type',{'{}'},'subs',{{':'}}) ) | ...
   isequal( s(1) , struct('type',{'.' },'subs',{'VALUE'}) )
   s   = s(2:end);
end

if isempty(s)
   obj.VALUE = val;
   return
end

%---------------------------------------------------------
% Check for NOT  () | {}

if ~( isequal( s(1).type , '()' ) | ...
      isequal( s(1).type , '{}' )       )

  try
     obj.VALUE = subsasgn(obj.VALUE,s,val);
  catch
     error(lasterr)
  end

  return

end

%---------------------------------------------------------

  [msg,s(1).subs,n] = cyclind( size(obj.VALUE) , s(1).subs , 1 );

  if ~isempty(msg)
      error(msg)
  end

  nv = prod(size(val));

  if ~( ( prod(n) == nv )  |  ( any( nv == [ 1  0 ] ) ) )
      error('Subscripted assignment dimension mismatch.')
  end

  ns = prod(size(s));

  switch s(1).type

    %-------------------------------------------------------------------------------
    case '()'

      try
        if ns == 1
           obj.VALUE(s(1).subs{:}) = val;
        else
           obj.VALUE(s(1).subs{:}) = subsasgn(obj.VALUE(s(1).subs{:}),s(2:end),val);
        end
      catch
        error(lasterr)
      end

    %-------------------------------------------------------------------------------
    case '{}'

      try
        if ns == 1
           obj.VALUE{s(1).subs{:}} = val;
        else
           obj.VALUE{s(1).subs{:}} = subsasgn(obj.VALUE{s(1).subs{:}},s(2:end),val);
        end
      catch
        error(lasterr)
      end

  end

 