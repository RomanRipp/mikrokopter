function [txt,c,msg] = disp(obj,varargin)

% DISP  Returns Text of Object to Display
%
% [ Text , CellString , Msg ] = DISP( Object )
%
%   Displays the full Object
%
%
% [ Text , CellString , Msg ] = DISP( Object , FieldNames )
%
%   Displays only the Fields from the CellStringArray FieldNames
%
%
% [ Text , CellString , Msg ] = DISP( Object , FieldName1 , FieldName2 , ... )
%
% calls: @class/history.m
%        @class/getstr.m
%

txt = char(0,0);
  c = cell(0,5);
msg = '';
  
Nout = nargout;

  Nin = nargin;

show_all = ( Nin == 1);

si = size(struct(obj));
 n = prod(si);

%---------------------------------------------------------
% Get Strings

if show_all

   [par,field] = history(obj);

    field = cat(1,field{:});

else

   if Nin == 2
      field = varargin{1};
   else
      field = varargin;
   end

   [ok,field] = chkcstr(field,0);

   if ok
      field = field(:);
   else

      msg =  'Fieldnames must be Strings.';

      if Nout < 2
         error(msg);
      end

      return

   end 

end


[str,msg] = getstr(obj,field);

if ~isempty(msg)

    if Nout < 3
       error(msg);
    end

    return

end

%---------------------------------------------------------
% Replace NewLine-Character

 str = strrep(str,char(10),'\\');

%---------------------------------------------------------
% Make ColumnVectors, Separate with Empty Lines

is_single = ( n <= 1 );

np     = 3 - 2 * is_single;

ns     = size(str,2);

sep    = cell(np,ns);
sep(:) = { '' };

field = field(:,1);

   nf = size(field,1);


fld = cat( 1 , sep , field(:,ones(1,ns)) );
fld = fld(:);

str = cat( 1 , sep , str );
str = str(:);


%---------------------------------------------------------
% Add Class and Parent

if show_all

  si = sprintf('%.0f by ',si);
  si = si(1:(end-4));

  str = cat( 1 , {par{end}} , {strhcat(par(1:(end-1)),':')} , {si} , str );

  fld = cat( 1 , {'Class'} , {'Parent'} , {'Size'} , fld );

end


fld = char(fld);

fld( find( double(fld) == 32 ) ) = char(10);

fld = strrep( cellstr(fld) , char(10) , char(32) );

%---------------------------------------------------------
% Build Text

nl = char(10);

nt = size(fld,1);

c = cell(nt,5);

c(:,1) = { '   ' };
c(:,2) = fld;
c(:,3) = { '  =  ' };
c(:,4) = str;
c(:,5) = { nl };

di = nf+np;
i0 = 0 + 3 * show_all;

for ii = 1 : np
  i1 = i0 + ii;
  c( i1:di:nt , [1 2 3] ) = { '' };
end

for ii = 1 : n*(~is_single)
  i1 = i0 + 2 + (ii-1) * di;
  c( i1 , 1 ) = { sprintf( '%% --- %3.0f of %.0f  --- %%' , [ii n] ) };
end

if ~show_all
   c(1,:) = [];
end

txt = permute(c,[2 1]);

txt = cat(2,txt{:});
