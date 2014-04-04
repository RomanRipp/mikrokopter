function [msg,v] = subref(v,s)

% SUBREF  single SUBSREF for ObjectValue
%
% calls:   @class/private/cyclind.m
%


  msg = '';

  if any( strcmp( s.type , { '()'  '{}' } ) )
     s.subs = cyclind( size(v) , s.subs );
  end

  try
     v = subsref(v,s);
  catch
     msg = lasterr;
  end

  return

  c = class(v);
  msg0 = cat(2,'Invalid Syntax "',s.type,'" for class "',c,'".');

  switch s.type

     %--------------------------------------------
     case  '()'

       ind = cyclind( size(v) , s.subs );

       v = v(ind{:});

     %--------------------------------------------
     case  '.'

       is_obj = isobject(v);
       is_str = isstruct(v);

       if ~( is_str | is_obj )
          msg = msg0;
          return
       end

       %-----------------------------------------------
       % Try Object

       if isobject(v)
          try
             v = subsref(v,s);
             is_str = 0;
          catch
             v = struct(v);
             is_str = 1;
          end
       end

       %-----------------------------------------------
       % Try Structure

       if is_str

          if any( strcmp( s.subs , fieldnames(v) ) );

             try 
                  v = getfield( v , s.subs );
                  msg = '';
             catch
                  msg = lasterr;
             end

          else

             msg = cat(2,'Invalid Fieldname "',s.subs,'".');

          end   

       end

     %--------------------------------------------
     case '{}'

       if strcmp( class(v) , 'cell' )

            ind = cyclind( size(v) , s.subs );

            v = v(ind{:});

            if isempty(v)
 
               v = [];

            elseif prod(size(v)) == 1

               v = v{1};

            end

       else

           msg = msg0;

       end

     %--------------------------------------------
     otherwise
 
       msg = msg0;

  end
