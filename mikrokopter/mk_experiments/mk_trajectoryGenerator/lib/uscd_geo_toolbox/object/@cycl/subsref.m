function obj = subsref(obj,s)

% SUBSREF  Subsref for Value of CYCL-Object
%
%
% Note:  OBJ.VALUE  or  OBJ{:}  returns the ObjectValue
%
%        multiple Arguments using with CellArrays doesn't works
%


obj = obj.VALUE;

%---------------------------------------------------------
% Check for: obj{:} | obj.VALUE

if isequal( s(1) , struct('type',{'{}'},'subs',{{':'}}) ) | ...
   isequal( s(1) , struct('type',{'.'},'subs',{'VALUE'}) )

   s   = s(2:end);

end

if isempty(s)
   return
end

%---------------------------------------------------------
% Check for: () | {}

if isequal( s(1).type , '()' )  | ...
   isequal( s(1).type , '{}' ) 

 
       [msg,s(1).subs] = cyclind( size(obj) , s(1).subs , 0 );

       if ~isempty(msg)
          error(msg)
       end

       try
          obj = subsref(obj,s(1));
       catch 
          error(lasterr)
       end

       s = s(2:end);

end

if isempty(s)
   return
end

%---------------------------------------------------------
% other

try
  obj = subsref(obj,s);
catch
  error(lasterr);
end


