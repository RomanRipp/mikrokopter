function c = num2char(c,r)

% NUM2CHAR   Converts numeric Array's to CharacterArrays
%
%  C = NUM2CHAR( Array , Replace );
%
%  default: Replace == ' ';  % SpaceCharacter
%
%  The values of Array will transformed into Interval [ 0 .. 255 ]
%  Bad Characters will replaced.
%
%  Good Characters are:  9, 10, 13, [ 28 .. 126 ], [ 160 .. 255 ]
%

 
Nin = nargin;
Msg = '';
 nl = char(10);

if Nin == 0
   c = char(zeros(0,0));
   return
end

%---------------------------------------------------------------

if ~isnumeric(c) | ischar(c)
   Msg = 'Input must be a Numeric- or Character-Array.';
end

%---------------------------------------------------------------

if Nin < 2
   r = 32;
elseif isempty(r)
   r = [];
else
   if ~( ( isnumeric(r) | ischar(r) ) & ( prod(size(r)) == 1 ) )
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
              'Value for REPLACE must be empty or a single Numeric or Character.' ];
   end
end

%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%---------------------------------------------------------------

if isempty(c)
   c = char( zeros(size(c)) );
   return
end


  if ~isa(c,'double')
     c = double(c);
  end

c = floor( c - 256 * floor( c / 256 ) );


if ~isempty(r)

   if ~isa(r,'double')
      r = double(r);
   end

   r = floor( r - 256 * floor( r / 256 ) );

end


c( find( ~( ( c ==  9 ) |  ...
            ( c == 10 ) |  ...
            ( c == 13 ) |  ...
            (  28 <= c  &   c <= 126 ) | ...
            ( 160 <= c  &   c <= 255 )        ) ) ) = r;


c = char(c);