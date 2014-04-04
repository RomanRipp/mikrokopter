function [str,msg] = getstr(obj,field,form);

% GETSTR Returns String for FieldValues of an Object
%
% [ String , Msg ] = GETSTR( Object )
%
%  Returns the Strings of all Fields of an Object
%
% [ String , Msg ] = GETSTR( Object , FieldNames )
% 
%  Returns the Strings for the specified FieldNames 
%
% [ String , Msg ] = GETSTR( Object , FieldNames , Format )
% 
%  Returns the Strings for the specified FieldNames, using the given Format.
%
% String, FieldNames, Format are CellStringArrays.
%
%
% Contains Object more then 1 Element, String is an MultiColumn-CellStringArray,
%  one Column for each Element.
%
%  calls: @class/history.m
%         @class/getfield.m
%         
%         val2str.m
%

  
str = cell(0,1);
msg = '';


 nl = char(10);

Nin = nargin;

if Nin == 0
   return
end


Nout = nargout;

%*************************************************************
% Check Inputs

%-------------------------------------------------------------
% FieldNames

if Nin < 2

  [par,field] = history(obj);

   field = cat(1,field{:});

else

  [ok,field] = chkcstr(field,0);

  if ok
     field = field(:);
  else
     msg = 'Input "FieldNames" must be a String or CellStringArray.';
  end

end
   
nf = size(field,1);

%-------------------------------------------------------------

if ~isempty(msg)
  
  if Nout < 2
     error(msg);
  end

  return

end


%*************************************************************

n = max(prod(size(struct(obj))),1);

str    = cell(nf,n);
str(:) = { '' };

if isempty(struct(obj))
   return
end

ok = zeros(nf,1);

form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

for jj = 1 : n

  for ii = 1 : nf

    command = 'GETFIELD';

    [val,msg1] = getfield(obj(jj),field{ii});

    if isempty(msg1)

%-------------------------------------------------------------
% Use following example for multiple ObjectFields
%
%       switch field{ii}
%
%        %----------------------------------------------------
%        case  'class'
%
%          txt = cat(2,'''',val,'''');
%
%        %----------------------------------------------------
%        case 'size'
%
%          txt = cat(2,'[',sprintf(' %.0f ',val),']');
%
%        %----------------------------------------------------
%        case  { 'created'  'modified' }
%
%          val = val(:);
%            s = size(val,1);
%          val = cat( 1 , val , zeros(6-s,1) );
%
%          val = datenum(val(1),val(2),val(3),val(4),val(5),val(6));
%
%          txt = datestr(val,0);
%
%        %----------------------------------------------------
%        otherwise
%-------------------------------------------------------------

          txt = displ(val);
 
%-------------------------------------------------------------
%       end
%       % switch field
%-------------------------------------------------------------
          
    end

    ok(ii) = isempty(msg1);

    if ok(ii)

      str{ii,jj} = txt;

    else

      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
         'Error using ' ,  command  , ' for Field "' , field{ii} , '": ' , msg1 );
 

    end

  end
  % ii = 1 : nf

end
% jj = 1 : n

  
if ~isempty(msg)  &  ( Nout < 2 )
    error(msg);
end
