function y = noise(n,scl,x,mode);

% NOISE  generate noises with periodicity
%
% Y = NOISE( N , Scale , [Start] , [Mode] )
%
%  Scale = Intervall/TimeScale
%
%  Scale = [ 1 by 1 ] or M-Element Vector
%  Start = [ 1 by 1 ] or P-Element Vector
%          not given ==> RandomValue
%
%  Y = N-Element Vector or M-P-N-Matrix
%
%  Mode = 0 | 1
%
%  Mode = 1 ==> equal Noise per TimeStep
%

Nin = nargin;

%************************************************************
% Check Inputs

if Nin < 2
   error('Not enough InputArguments.');
end

msg = cell(0,1);

%------------------------------------------------------------

ok = ( prod(size(n)) == 1 );
if ok
   ok = ( ( mod(n,1) == 0 ) & ( n > 0 ) );
end
if ~ok
    msg = cat( 1 , msg , {'N must be a positive Integer.'} );
end

%------------------------------------------------------------
% Scale

sc = size(scl); pc = prod(sc);

ok = ( ~( pc == 0 )  & ( pc == max(sc) ) );
if ~ok
    msg = cat( 1 , msg , {'Scale must be Scalar or nonempty Vector.'} );
end

%------------------------------------------------------------
% StartValue

if Nin < 3
   x = [];
end

if isempty(x)
   sx = size(x);
   sx = sx + ( sx == 0 );
   x = randn(sx);
end

sx = size(x); px = prod(sx);

ok = ( ~( px == 0 )  & ( px == max(sx) ) );
if ~ok
    msg = cat( 1 , msg , {'StartValue must be Scalar or nonempty Vector.'} );
end

%------------------------------------------------------------
% Mode for Random on X

if Nin < 4
   mode = 0;
else
   mode = isequal(mode,1);
end

%------------------------------------------------------------

if ~isempty(msg);
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%************************************************************

if any( scl < 0 )
   msg = cat( 1 , msg , {'Negative Scale.'} );
end

if any( abs(scl) > 1 )
   msg = cat( 1 , msg , {'ABS(Scale) exeeds ONE.'} );
end

if any( abs(x) > 2*pi )
   msg = cat( 1 , msg , {'ABS(StartValue) exeeds 2*PI.'} );
end

if ~isempty(msg);
    msg = sprintf('%s\n',msg{:});
    warning(sprintf('Strange Inputs.\n%s',msg));
end

%************************************************************
% Check Dimensions, Flip to [ X by SCL by N ]

vx = ( px > 1 );
vc = ( pc > 1 );

if vx
   x = x(:);
end

if vc
   scl = scl(:)';
end

if vx | vc

   di     = [ min(find(sx==px))  min(find(sc==pc)) ]; % Dimension for [ X SCL ]
   di     = di .* [ vx  vc ];
   di(2)  = di(2) + ( di(2) == di(1) );           % Raise SCL if equal to X

   dn     = ones( 1 , max(di)+1 );

   ii     = find(di);                             % Valid Dimensions
   di     = di(ii);

   dn(di) = 0;
   dn     = min(find(dn));                        % Dimension for N

   di     = [ di dn ];
   ii     = [ ii 03 ];
   
   id     = ones( 1 , max(max(di),3) );
   pm     = 0 * id;
   pm(di) = ii;
   id(ii) = 00;

   while any(id)
            ii = min(find( id));
            jj = min(find(~pm));
         pm(jj) = ii;
         id(ii) = 00;
   end

else

   pm = [ 3  1  2 ];  % N in 1. Dimension

end

%************************************************************

ox = ones(1,px);
oc = ones(1,pc);

if mode
   si = [  1  1 ];
   or = 1;
else
   si = [ px  1 ];
   or = oc;
end

scl = 1 - scl;
sql = sqrt( 1 - scl.^2 );

scl = scl(ox,:);
sql = sql(ox,:);

y   = zeros(px,pc,n);

y(:,:,1) = x(:,oc) ./ sql;

for ii = 2 : n
    y(:,:,ii) = y(:,:,ii-1) .* scl + randn(si) * or;
end

y = y .* sql( : , : , ones(1,n) );

%************************************************************
%  Flip back

if ~all( pm == ( 1 : size(pm,2) ) )
    y = permute(y,pm);
end
