
function [obj,msg] = setfield(obj,field,val)

% SETFIELD Set a Value to a Field of an Object
%
% [Object,Msg] = SETFIELD( Object , FieldName , Value )
%


msg = '';
 nl = char(10);

Nout = nargout;

%---------------------------------------------------
% Check Input field

if ~chkstr(field,1)

    msg = 'FieldName must be a String.';

    if Nout < 2
       error(msg);
    end

    return

end


%---------------------------------------------------
% Get ObjectField-  and ParentFieldName
 
[f,par] = fieldnames(obj);

%---------------------------------------------------

is_field  = any( strcmp( field , f ) );
is_parent = ~isempty(par);
 
if ~( is_field | is_parent )

    msg = cat(2,'Invalid FieldName "',field,'".');

    if Nout < 2
       error(msg);
    end

    return

end

%---------------------------------------------------
% SUBSREF

if is_field

  sub = field;

else

  sub = par;

end


command = cat( 2 , 'obj(ii).' , sub , ' = vv;' );


for ii = 1 : prod(size(struct(obj)))  

  msg1 = '';

  try
    
     if is_field

         vv = val;

     else

         [vv,msg1] = setfield( getfield(obj(ii),par) , field );

     end

     if isempty(msg1)
         eval(command)
     end

  catch
  
     msg1 = lasterr;
 
  end

  if ~isempty(msg1)

     msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , msg1 );

  end

end

if ~isempty(msg)

   if  ( Nout < 2 )
       error(msg);
   end

   return

end


