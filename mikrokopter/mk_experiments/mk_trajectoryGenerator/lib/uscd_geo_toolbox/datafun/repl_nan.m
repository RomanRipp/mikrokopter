function   Cneu = repl_nan( Corg , nxy , set_nan )

% REPL_NAN  replace NaN-Values by mean of surrounding non-NaN
%
% C = REPL_NAN( Corg , N , [Set_NaN] )
%
% N = Nxy | [ Nx  Ny ] Depth to replace NaN-Values
%
% Not replaced NaN-Values will set to Set_NaN
%
% see also: ZEROSIZE
%

Nin = nargin;

if Nin < 2 ,
   nxy = [ 1  1 ]; 
else
  jj = find(~finite(nxy)); nxy(jj)=0*jj;
  if     isempty(nxy) 
   nxy = [ 1  1 ];
  elseif length(nxy) == 1
   nxy = [ nxy  nxy ];
  else
   nxy = nxy(1:2);
  end
end

if Nin < 3
   set_nan = NaN;
end

is_nan = find( isnan(Corg) );

Cneu         = Corg      ;
Cneu(is_nan) = 0 * is_nan;

si = size(Corg);

quot = zeros( si );

% Zahler
zxy = [0 0];

for ii = 1 : max(nxy) 

 for dimen = 1 : 2
 % 1 .. Y-Richtung
 % 2 .. X-Richtung

    if zxy(3-dimen) < nxy(3-dimen)

        N = si(dimen);
        M = si(3-dimen);

      ind = [ 1 : N-1 ]';

      for versch = [ 0   1 ]
      % 0 suche drueber   1 suche drunter

      not_nan = find( [ zeros(versch,M) ;   ... 
                       xor( isnan( Corg(ind+(~versch),:) +         ...
                                   Corg(ind+  versch,:)   ) ,     ...
                            isnan( Corg(ind+(~versch),:)   )   );   ...
                         zeros(~versch,M) ]     ... 
                     == 1 ) ; 

       Cneu(not_nan) = Cneu(not_nan) + Corg( not_nan - (2*versch-1) );
       quot(not_nan) = quot(not_nan) + 1;
      end
      % versch
 
    end
    % ZaehlerCheck 

    Corg = Corg';
    Cneu = Cneu';
    quot = quot';
 
    % Zaehler um 1 erhoehen
     zxy(3-dimen) = zxy(3-dimen)+1;

   end
   % dimen

   not_nan = find( quot ~= 0 );

   if isempty(not_nan)
      break
   end

   Corg(not_nan) = Cneu(not_nan) ./ quot(not_nan);
   Cneu(not_nan) = Corg(not_nan)                 ;

   quot = 0*quot;
  
end
% n

if 1 %%% nargin == 3

       is_nan  = find( isnan(Corg) );
  Cneu(is_nan) = set_nan;

end

