function [spind, pout] = pspike(pin,tol,wlen)
% PSPIKE  markiere einzelne Spikes und interpoliere ggf.
%
%        I = PSPIKE(P,TOL,WLEN) sucht einzelne Spikes in P durch Vergleich
%        zwischen dem Wert des mittleren Elementes und dem Medainwert in
%        einem gleitenden Suchfenster der L"ange WLEN. Ist die Differenz
%        gr"osser als TOL, so wird der Index des mittleren Elementes in
%        I aufgenommen.
%
%        I = PSPIKE(P) setzt TOL = 2.0 und WLEN = 3 ein.
%
%        [I,POUT] = PSPIKE(...) interpoliert P(I) und gibt den interpolierten
%        Vektor als POUT aus.
%


pin = pin(:);
if nargin < 3,
	wlen = 3;	% Fensterl"ange
	tol = 2.0;	% Toleranz
elseif nargin < 2,
	tol = 2.0;	% Toleranz
end;
if rem(wlen,2)==0,
	error('WLEN must be odd.');
end;

pwork = pin(:)';
pend = pin(length(pin));
count = 0;
spind = [];
N = (wlen-1)/2;

% erzeuge eine Matrix, in der pwork jeweils um 1 verschoben ist, so dass
% man 'median' vektorisiert darauf loslassen kann:

Pakt = zeros(wlen,wlen-1+length(pwork));
sizP = size(Pakt);

% Jetzt die "Schw"anze" auf pend setzen, damit nicht versehentlich der letzte
% Wert in der Nachbarschaft von Nullen als Spike interpretiert wird:

Pakt(:,(sizP(2)-2*N+1):sizP(2)) = zeros(sizP(1),2*N)+pend;

for n=1:wlen,
    Pakt(n,(wlen-n)+1:(wlen-n)+length(pwork)) = pwork;
end;

isn = any(isnan(Pakt(:)));

while 1,
	Mid = Pakt(N+1,:);		% mittl. Werte
        if isn
           Med = mednan(Pakt);
        else
	   Med = median(Pakt);
        end
	ii = find(abs(Mid-Med)>tol)-1;
	if isempty(ii),
		break;
	elseif ii <= 0,
		break;
	else
	ii = ii(1)-N+1;			% immer nur den ersten gefundenen Spike mitnehmen
	end;
	count = count+1;	% Verkuerzung beruecksichtigen
	% kleiner Kunstgriff: Der als Spike identifizierte Wert pwork(ii)
	% soll rausgeschmissen werden; um das durch Manipulation von Pakt 
	% zu erreichen, muessen die entsprechenden Indizes berechnet werden:
	zind = 0:wlen-1;
	ind1 = ii + (wlen-1);							% (wlen-1) ist die Zahl der Nullen in der 1. Zeile
	ind2 = (length(pwork)) + (wlen-2);		% in der 2. Zeile steht eine Null weniger
	idx  = ind1 + zind*ind2;

	Pakt=Pakt';
	Pakt(idx)=[]; 
	Pakt=reshape(Pakt,sizP(2)-count,sizP(1));
	Pakt=Pakt';
	
	pwork(ii) = [];	% diesen ignorieren und wieder von vorne suchen
	spind = [spind; ii+count-1]; 	% Index merken (+count, weil pwork 
											% ja jedesmal um 1 kuerzer wird
end;

spind = union1(sort(spind));
% disp(['pspike: marked ',num2str(length(spind)),' values as spikes.']);

if nargout > 1,
	if ~isempty(spind),	
%		if any(diff(spind)<2),
%			beep(3);
%			disp('pspike: warning: two or more adjacent spikes indicated!')
%			disp('        Interpolation may not work!')
%		end;
		pout = pin;
		spind = union1(spind);
%		disp('pspike: interpolating...')
		iv = spind-1;	% alle Punkte vor einem Spike
		in = spind+1;	% alle Punkte nach einem Spike
		iv(find(diff(iv)==1)+1)=[];	% alle direkt aufeinanderfolgenden Indizes
		in(find(diff(in)==1)  )=[];	% wegwerfen (sind selber Spikes)
		is = union1(iv,in);	% ein Satz von Indizes, der die Spikes umgibt
		pout(spind) = interp1(is,pin(is),spind);
	else
		pout = pin;
	end;
end;


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function u = union1(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9)

% UNION1 Form the union of 'sets'.
%   U = UNION1(X,Y) forms the union of the elements in X and Y,
%   removing repeated elements. X and Y may be any size, U will be
%   a sorted row vector. 
%   The function accepts 0 to 10 input arguments,
%   for example U=UNION1(X) or U=UNION1(X,Y,Z). 

% Author: Robert Piche, Tampere University of Technology, Finland
%         (piche@butler.cc.tut.fi)    
% Last revised: Jan 18, 1994

x = [];
for k=1:nargin
  eval(['x=[x x' num2str(k-1) '(:).''];' ]);
end

u = sort(x); 
n = length(u);
if n > 1
   u(find([0 ~(u(2:n)-u(1:n-1))])) = [];
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function y = mednan(x,dim)

% MEDNAN   Median value, ignores NaN
%
%   For vectors, MEDNAN(X) is the median value of the elements in X.
%   For matrices, MEDNAN(X) is a row vector containing the median
%   value of each column.  For N-D arrays, MEDNAN(X) is the median
%   value of the elements along the first non-singleton dimension
%   of X.
%
%   MEDNAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then mednan(X,1) is [1.5 2.5 3.5] and median(X,2) is [1
%                                                         4]
%
%   See also MEDIAN, MEAN, STD, MIN, MAX, COV.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.11 $  $Date: 1997/11/21 23:23:56 $

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end
if isempty(x), y = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);


% Run over Columns

y  = NaN * ones(1,size(x,2));

for ii = 1 : size(x,2)

  jj = ~isnan(x(:,ii));
  n  = sum(jj,1);
  r  = rem(n,2);  
  r  = r - 1 * ( n == 0 );
   
       jj  = find(jj);

  if     r == 0  % Odd number of elements along DIM

     y(ii) = (x(jj(n/2),ii) + x(jj(n/2+1),ii))/2;

  elseif r == 1  % Even number of elements along DIM

     y(ii) = x(jj((n+1)/2),ii);     

  end

end

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);
