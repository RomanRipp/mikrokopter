function p = combine(n,k)

% COMBINE  returns Combinations of Elements
%
% P = COMBINE( N , K )  returns all combinations 
%                        of K elements on N choices
%
% P = [ M by K ]  IndexVector,  M = N! / ( K! * (N-K)! )
%
% for large numbers of N, COMBINE works faster then 
%
%       NCHOOSEK( 1:N , K ) (SpecializedFunctions)
%   or  COMBNTNS( 1:N , K ) (MappingToolbox)
%
%------------------------------------------------------
%                  
% see also:  PAIRS, NCHOOSEK, COMBNTNS
%

%***************************************************************
% Check Inputs

Nin = nargin;

if Nin == 0
   error('Not enough InputArguments.');
end

if Nin < 2
   k = [];
end

msg = cell(0,1);
v   = {  n    k  };
l   = { 'N'  'K' };

for ii = 1 : Nin
    ok = ( isnumeric(v{ii}) & ( prod(size(v{ii})) == 1 ) );
    if ok
       ok = ( ( v{ii} == round(v{ii}) ) & ( v{ii} >= 0 ) );
    end
    if ~ok
        msg = cat(1,msg,{sprintf('%s must be a positive Integer.',l{ii})});
    end
end

if isempty(msg)
   if Nin == 2
      if k > n 
         msg = {'K must be in Intervall [ 0 .. N ].'};
      end
   else
      k = min(n,1);
   end
end

if ~isempty(msg)
    error(sprintf('%s\n',msg{:}));
end

p = comb(n,k);

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function p = comb(n,k)

if any( k == [ 0  1  n ] )

   if     k == 0
      p = zeros(1,0);
   elseif k == 1
      p = ( 1 : n )';
   else
      p = ( 1 : k );
   end   

   return

end

%------------------------------------------------------

m = nok(n,k);

p = zeros(m,k);

z = 0;

for ii = 1 : n-1
 
     m = nok(n-ii,k-1);
    jj = z + (1:m);

    p(jj,1) = ii;
    p(jj,2:k) = comb(n-ii,k-1) + ii;

     z = z + m;

end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function n = nok(n,k)

n = prod(n-k+1:n) / prod(1:k);
