function str = var2mstr(val,bl)

% VAR2MSTR   Converts Variables to String with Matlab expression
%
% String = VAR2MSTR( V )
%
%  Converts 2-dimensional MatlabVariables 
%   to their String-Expression, using by EVAL 
%             
% Variables could be Numeric, Char, Cell- or StructArrays.
%
% String = VAR2MSTR( V , N )
%
%  real(N) = Number of Blanks around the Seperator
%  imag(N) = Number of Blanks around a numeric
%
%  real(N) < 0  ==> No "," - Separator along 2. Dimension
%
% default: N = 1
%
 
if nargin < 2
   bl = [];
end

if isempty(bl)
   bl = 1 + 0*i;
end

sp =     real(bl(1));   % Blank arround Seperator
nb = abs(imag(bl(1)));  % Blank arround Numbers

nb = char(32*ones(1,nb));

   % Usage to convert Number x to char:
   % vk = floor(log(abs(x))/log(10))+1;
   % vk = vk*(vk>0);
   % sprintf(sprintf(kform,nk+vk),x);

 % Seperator for Dimensions
   sep0 = char(32*ones(1,abs(sp)));

   sep1 = [ sep0  ';'  sep0 ];

   sep0 = char(32*ones(1,abs(sp)));
   if sp < 0
      sep2 = sep0;
   else
      sep2 = [ sep0  ','  sep0 ];
   end

 % Make 2-Dimensional
   if isobject(val)
      si = size( struct(val) );
   else
      si = size(val);
   end

  val = reshape(val,si(1),prod(si(2:end)));

   si = size(val);
   si = si * (~isempty(val));

try

   str = get_str(val,si,sep0,sep1,sep2,nb,bl);

catch

   str = -1;

end

if isequal(str,-1)

   % Build String from Size and class

   cl = class(val);

   ff = [ '[]' ; '{}' ];

   ff = ff(1+strcmp(cl,'cell'),:);

   str = sprintf('%.0fx',si);

   str = cat(2,ff(1),str(1:end-1),' ',cl,ff(2));

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = get_str(val,si,sep0,sep1,sep2,nb,bl)

 str = -1;

 %**************************************************
 % Char
 if ischar(val) 

   if si(1) <= 1
      str = [ ''''  val '''' ];
   else
      str = [ '['  sep0 ]; 
      for ii = 1 : si(1)
        str = [ str '''' val(ii,:) '''' ...
                sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  ']' ];
   end 

 %**************************************************
 % Numeric
 elseif isnumeric(val)

   if ~strcmp(class(val),'double')
       val = double(val);
   end

   clear eps
   nk    = floor(abs(log(eps)/log(10)))-1 ;
   kform = [ nb  '%%.%.0fg'  nb ];

   if prod(si) == 0
      str = '[]';
   elseif prod(si) == 1
      str = numstr(val,nb,nk,kform);
   else
      str = [ '['  sep0 ]; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
            str = [ str  numstr(val(ii,jj),nb,nk,kform)  sep2(1:(end*(jj~=si(2)))) ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  ']' ];
   end 

 %**************************************************
 % CellArray
 elseif iscell(val)

   if prod(si) == 0

         str = '{}';

   else

      str = [ '{'  sep0 ] ; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(val{ii,jj},bl)           ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  '}' ];
 
   end

 %**************************************************
 % StructArray
 elseif isstruct(val)

    fnames = fieldnames(val);
    fnames = fnames(:);
    nf     = size(fnames,1);
   
    fsep = ',';

    str = [ 'struct(' ];
  
    for ff = 1 : size(fnames,1);

      str = [ str  var2mstr(fnames{ff},bl) fsep '{' ]; 

      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(getfield(val,{ii,jj},fnames{ff})) ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end

      str = [ str '}' fsep(1:(end*(ff~=nf)))  ];

    end

    str = [ str ')' ];

 %**************************************************
 % StructArray
 elseif isobject(val)

    str = disp(val);

 end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = numstr(val,nb,nk,kform)

vi = imag(val);

if ~( vi == 0 )

    vr = real(val);

    sr = '';
    if ~( vr == 0 )
        sr = numstr( vr  , '' , nk , kform );
    end

    add_plus = ( ~isempty(sr) & ( ( vi > 0 ) | isnan(vi) ) );

    if abs(vi) == 1
       si = 'i';
       if vi < 0
           si = [ '-' si ];
       elseif add_plus
           si = [ '+' si ];
       end
    else
       si = numstr( vi , '' , nk , kform );
       if add_plus
           si = [ '+' si ];
       end
       si = [ si '*i' ];
    end

    str = [ nb sr si nb ];

    return

end

if val == 0
   str = '0';
elseif isnan(val)
   str = 'NaN';
elseif isinf(val)
   str = 'Inf';
   if sign(val) == -1
      str = [ '-' str ];
   end
elseif ( abs(mod(val,pi)) < 1e3*eps ) & ( val < 1/(1e3*eps) )
   val = round(val/pi);
   str = 'pi';
   if ~( val == 1 )
       if val == -1
          str = '-pi';
       else
          str = [ numstr(val,'',nk,kform) '*pi' ];
       end
   end
else
   vk = floor(log(abs(val))/log(10))+1;
   vk = vk*(vk>0);
  str = sprintf(sprintf(kform,nk+vk),val);
end
        
str = [ nb  str  nb ];
