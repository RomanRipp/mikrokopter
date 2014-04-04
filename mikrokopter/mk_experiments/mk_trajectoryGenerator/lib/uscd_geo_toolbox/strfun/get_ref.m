function [Msg,out,ind] = get_ref(str,s1,s2,c1,c2)

% GET_REF  Returns a Reference from a String
%
% [Msg,Reference] = GET_REF( String , S1 , S2 , C1 , C2 );
%
%  S1:  StartMarker  Characters which marks the Begin of Reference
%  S2:    EndMarker  Characters which marks the End   of Reference
%
%  C1:  Characters in Reference to remove,
%        default: [ CR NL TAB ]
%  C2:  Characters closed to C1 to remove,
%        default: [ Space '-' ]
%
%  Reference = { FullString  ReferenceString }  2-Column Cell-Array
%
%
% [Msg,Reference,Index] = GET_REF( ... )
%
%  Returns Index = [ StartIndex Length ] of FullString.
%
%  To Remove the References from String, use: 
%     
%    String( grp2ind( Index(:,1) , Index(:,2) ) ) = [];
%  

Msg = '';

out =  cell(0,2);
ind = zeros(0,2);

nl = char(10);

Msg0 = 'GET_REF: ';

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);

%---------------------------------------------------------

if ~( ischar(str) &  ( prod(size(str)) == size(str,2) ) )
  Msg = '1. Input must be a String.';
end
 
if ~( ischar(s1) &  ( prod(size(s1)) == size(s1,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end
 
if ~( ischar(s2) &  ( prod(size(s2)) == size(s2,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end

%---------------------------------------------------------

if nargin < 4

  c1 = [ 13  10  9 ];

else

  if ischar(c1)
    c1 = double(c1);
  end
  ok = isnumeric(c1);
  if ok & ~isempty(c1)
    c1 = c1(:)';
    ok = all( ( mod(c1,1) == 0 )  & ...
              ( c1 >= 0 ) & isfinite(c1)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if nargin < 5

  c2 = [ double('-')  32  ];

else

  if ischar(c2)
    c2 = double(c2);
  end
  ok = isnumeric(c2);
  if ok & ~isempty(c2)
    c2 = c2(:)';
    ok = all( ( mod(c2,1) == 0 )  & ...
              ( c2 >= 0 ) & isfinite(c2)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
end

if isempty(str)  |  isempty(s1)  |  isempty(s2)  |  ...
   ~isempty(Msg)
 
  return

end


%**********************************************************

n = size(str,2);

%--------------------------------
% Start of Reference at End of s1

i1 = findstr(str,s1);

if isempty(i1)
  return
end

i1 = i1 + size(s1,2) - 1;
i1( find( i1 > n ) ) = [];

if isempty(i1)
  return
end
 
%--------------------------------
% End of Reference at Begin of s2

i2 = findstr(str,s2);

if isempty(i2)
  return
end

i2( find( i2 < i1(1) ) ) = [];

if isempty(i2)
  return
end

%----------------------------

ok = zeros(1,n);  % Start/End
ii = zeros(1,n);  % Index

ok(i2) = 2;  % End   of Reference
ok(i1) = 1;  % Start of Reference

ii(i2) = i2;
ii(i1) = i1;

jj = find(ok);

ok = ok(jj);
ii = ii(jj);

%----------------------------
% Following Start-End

ok = find( abs(diff(ok)) == 1 );

if isempty(ok)
   return
end

%----------------------------
% Start's

ok = ok(1:2:end);

n = prod(size(ok));

out      = cell(n,2);
out(:,1) = { [ s1  s2 ] };
out(:,2) = { '' };

ind = zeros(n,2);
ind(:,1) = ii(ok+0)'-size(s1,2)+1;
ind(:,2) = ii(ok+1)'+size(s2,2) - ind(:,1);

for jj = 1 : n

  if ( ii(ok(jj)+1) > 1+ii(ok(jj)) )
  % Reference not empty

    % String in Reference
 
    sr = str( ii(ok(jj))+1 : ii(ok(jj)+1)-1 );


    % Full String

    out{jj,1} = [ s1  sr  s2 ];


    % Remove bad Characters from Reference

    sr = rmblank( sr , 2 );

    b1 = zeros( 1 , size(sr,2) );
    for r1 = c1
      b1 = ( b1  |  ( double(sr) == r1 ) );
    end

    b2 = zeros( 1 , size(sr,2) );
    for r2 = c2
      b2 = ( b2  |  ( double(sr) == r2 ) );
    end
  

    i0 = find( diff(cat(2,0,( b1 | b2 )  )) ==  1 );  % Start of Group
    i1 = find( diff(cat(2,  ( b1 | b2 ),0)) == -1 );  % End   of Group
 
     b = 2 * b1 + 1 * b2;
    cb = cumsum(b,2);

    lg = i1 - i0 + 1;              % Length of Group
    sg = cb(i1) - cb(i0) + b(i0);  % Sum    of Group         

    ib = find( sg > lg );          % Groups incl. b1

    lg = lg(ib);
    i0 = i0(ib);

    % Build IndexVector to remove

     b = grp2ind(i0,lg);

     sr(b) = [];
 
     out{jj,2} = sr;

  end
  % Reference not empty

end
% jj

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);