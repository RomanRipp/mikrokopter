function [val,msg] = defval(typ,si,val);

% DEFVAL   Returns DefaultValue for Type and Size
%
% [val,msg] = DEFVAL(typ,si,val);
%
% The 3. Input "val" is used, if the "typ" 'struct' is requested,
%   to create a Structure with new Size from existing Structure "val".
%
% typ: 'numeric'
%      'cellstr'
%      'cellvec'
%      'string'
%      'object'
%      <class>
%      ''
%

Nin = nargin;

if Nin < 1
   typ = '';
end

if Nin < 2
   si = [ 0  0 ];
end

if Nin < 3
   val = [];
end


msg = '';


si(find(isnan(si))) = 0;
 

switch typ

  case { 'char'  'string' }
 
          val = char( 32 * ones(si) );

  case { 'double'  'single'  'numeric'  '' }

          val = NaN * ones(si);

  case {  'int8'     'int16'    'int32' ...
         'uint8'    'uint16'   'uint32' ...
         'logical'  'sparse'        }

          val = feval( typ ,  zeros(si) );

  case 'cell'

          val = cell(si);

  case 'cellstr'

          val = cell(si);

          if ~any( si == 0 )
              val(:) = { '' };
          end

  case 'cellvec'

          val = cell(si);

          if ~any( si == 0 )
              val(:) = { [] };
          end

  case 'struct'

       if ~strcmp(class(val),'struct');

           val = struct( 'NaN' , { [] } );

       end

       ns = prod(size(si));

      ind = cell(ns,1);

      for ii = 1 : ns
          ind{ii} = ones(1,si(ii));
      end
 
      val = val(ind{:});

  otherwise

     try

        val = feval(typ,si);

     catch

        msg = lasterr;

     end

end
