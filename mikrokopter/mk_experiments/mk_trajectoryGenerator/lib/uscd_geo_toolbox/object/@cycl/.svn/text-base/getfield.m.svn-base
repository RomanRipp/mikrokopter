function [val,msg] = getfield(obj,field)

% GETFIELD  Returns specific FieldValue of Object
%
% [Value,Msg] = GETFIELD( Object , FieldName )
%
% calls: @class/subsref.m
%        @class/fieldnames.m
%        @class/getfield.m
%

val  = [];
msg  = '';

 nl = char(10);

Nin  = nargin;
Nout = nargout;

msg0 = [ upper(mfilename) ': ' ];

%---------------------------------------------------
% Check Input obj

if ~( prod(size(struct(obj))) == 1 )
  
   msg = [ msg0 'Object must contain 1 Element.' ];

end


%---------------------------------------------------
% Check Input field

if Nin < 2

    msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
               msg0 , 'FieldName is missing.' );

else

  if ~chkstr(field,1)

      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 msg0 , 'FieldName must be a String.' );

  end

end

%---------------------------------------------------

if ~isempty(msg)

    if Nout < 2
       error(msg);
    end

    return

end

%---------------------------------------------------

if strcmp( field , class(obj) );

   val = obj;

   return

end

%---------------------------------------------------
% Get FieldName, ParentName, ObjectStructure

[f,p,obj] = fieldnames( obj );


%---------------------------------------------------
if any( strcmp( field , f ) )

% Get Field Directly

  s = struct( 'type' , { '.' } , ...
              'subs' , { field }      );

  val = subsref( obj , s );

%---------------------------------------------------
elseif isempty(p)

% No Parent

  msg = cat( 2 , msg0 , 'Invalid FieldName "' , field , '".' );

else

% Try Parent

  obj = getfield(obj,p);  % ParentObject

  [val,msg] = getfield( obj , field );

end


if ~isempty(msg) & ( Nout < 2 )

    error(msg);

end

