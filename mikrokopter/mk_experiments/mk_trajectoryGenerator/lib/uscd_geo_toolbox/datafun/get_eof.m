function [good,ev,ew,amp,fac,cv]=geteof(data,w,filtpar,gap_par);

% GET_EOF calculates eigen-vectors, -values and the covariance matrix
%
% function [good,ev,ew,amp,fac,covmat]=geteof(data,z,[filtpar],[gap_par]);
%
% calculates the eigen-vectors, -values and the covariance matrix
% of a given data set 'data'
%
% input  : data   		- the data set (no data : NaN)
%          w                    - Weight for  1. Dimension of Data
%	   [filtpar]	[0]	- =0 : normal EOFs  
%				  >0 : EOFs of filtered cov-matrix,
%		                       provides you with orthogonal
%				       smooth modes (not the 'best' ones!!)
%				       filtpar is used as input for
%				       'filt2d'
%	   [gap_par]	[100]	- if a data-set is gappy, only
%				  'gap_par'-percent of the maximum number of
%				  eigenvectors is used
%
% output : good                 - vector of "good" data-rows, which used for EOF's  
%          ev     		- 'vertical' eigen-vectors of 'data'
%          ew     		- corresponding eigen-values
%	   amp   		- amplitudes of ev
%          fac    		- number of realisations of pairs
%	   covmat 		- covariance matrix
%
% version 1.1.3		last change 14.03.1996

% Gerd Krahmann, IfM Kiel, Feb 1993
%
% Ref.: U.Send's  EOF notes

% added filtering mode 		G.Krahmann, Sep 1994
% changed filtering mode	G.Krahmann, Jan 1995
% changed output order		G.Krahmann, Sep 1995
% transponing for speed		G.Krahmann, Feb 1996, removed
% removed problem with amps	
% and introduced 'gap_par'	G.Krahmann, Mar 1996
% new MatlabCode                Ch.Begler,  Apr 2001
% remove transponing            Ch.Begler,  Apr 2001
% added weight mode             Ch.Begler,  Apr 2001


Nin  = nargin;
Nout = nargout;

%---------------------------------------------------
% Check Inputs

if ~( isnumeric(data)  &  ~isempty(data) &  ...
   ( size(data,1)*size(data,2) == prod(size(data)) ) )
  error('"data" must be a 2-dimensional , nonempty Numeric.');
end


if Nin < 2
  w = [];
elseif ~( isnumeric(w)  & ...
          ( isempty(w) |  ( size(data,1) == prod(size(w)) ) ) ) 
  error('"w" must matching Length of 1. Dimension of "data" or empty.');
end

if Nin < 3
  filtpar = [];
end

if Nin < 4
  gap_par = [];
end



if isempty(filtpar)
  filtpar = NaN;
end

if isempty(gap_par)
  gap_par = 100;
end

 
s1 = size(data,1);
s2 = size(data,2);

good = zeros(s1,1);

%-----------------------------------------
% Weight

if isempty(w)

 w = ones(size(data,1),1);

end

 
%-----------------------------------------
% Check to Calculate Amplitude

is_ampl = ( Nout > 3 );


%-----------------------------------------
% Look for NaN-Columns

 col_ok = ( sum(isnan(data),1) < s1(1) );

 data = data( : , find(col_ok) );

 if isempty(data)

     ev =  NaN*ones(s1,s1);
     ew =  NaN*ones(s1,1);
     amp = zeros(s2,s1);
 
     return

 end

 s2 = size(data,2);

  
%-----------------------------------------
% Look for NaN-Rows

% Remove Rows with full NaN 
%  from Begin and End in 1. Dimension

nn = isnan(data);  % IsNaN

n0 =  1 + min( sum(cumprod(nn,1),1) );

n1 = s1 - min( sum(cumprod(nn(s1:-1:1,:),1),1) );


data = data(n0:n1,:);
   
  nn =   nn(n0:n1,:);

   w =    w(n0:n1,:);

good(n0:n1) = 1;

%-----------------------------------------

is_nan  = any( nn(:) );

if is_nan

    if 0
      h1  = ~isnan(data);
      fac = zeros(s1,s1);
      for ii = 1 : s1
        for jj = 1 : s1
          fac(ii,jj) = sum( h1(ii,:) .* h1(jj,:) ) ;
        end
      end
    end

    % NaN --> 0

    data(find(nn)) = 0;

 
     nn = ~nn;          % NOT NaN

    fac = nn * nn';

else

    fac = s2;

end

 
%-----------------------------------------
% calculate covariance matrix

cv = data * data';
cv = cv ./ ( fac + ( fac == 0 ) );

%-----------------------------------------
% filter covariance matrix 

if ( filtpar < 1 )
   cv = filt2d( cv , 2 , filtpar );
end


%-----------------------------------------
% calculate Singular Value Decomposition

 w = sqrt(w);

dw = diag(w);

[u,s,v] = svd(dw*cv*dw);

ev = diag(1./w) * u;
ew = diag(s);

ew = ew(:)';

s3 = size(ev,1);
o1 = ones(s3,1);

f = sqrt( sum((ev.^2),1) / s3 );

ev = ev ./ ( o1 *  f );
ew = ew .* (  f .^ 2 );


%-----------------------------------------
% project to get amplitudes

if is_ampl

    amp = zeros(s1,size(col_ok,2));
    cc  = find(col_ok);

  if is_nan 

    for ii = 1 : s2

       ok = find( nn(:,ii) );
      lok = ceil( size(ok,1) * gap_par/100 );

      ca  =  ev(ok,:)' * ev(ok,:);
      ca1 = inv( ca(1:lok,1:lok));

      amp(1:lok,cc(ii)) = ca1 * ( ev(ok,1:lok)' * data(ok,ii) );

    end

  else

    amp(:,cc) = ( ev' * data );

  end

     amp = amp';

end    

%-----------------------------------------
% Fill up removed NaN-Rows

if  ( n0 > 1 )  |  ( n1 < s1 )

   s0 = n1 - n0 + 1;  % Size of ev, ew

  %--------------------------------------- 
  ev1 = zeros(s1,s1);

  % Diagonal from (n1+1,n1+1)  -->  (s1,s1)
    i1 = ( n1+1 : s1 );
    i2 = i1;

    ev1( i1+(i2-1)*s1 ) = 1;

  % Diagonal from (1,s0+1)  -->  (n0-1,s0+n0-1)
    i1 = ( 1 : n0-1 );
    i2 =  i1 + s0;

    ev1( i1+(i2-1)*s1 ) = 1;

  %--------------------------------------- 

    ev1( n0 : n1 , 1 : s0 ) = ev;

    ev = ev1;


  %--------------------------------------- 

    ew1       = zeros(s1,1);
    ew1(1:s0) = ew;
    ew        = ew1;

end


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = filt2d(x,n,f)

%FILT2D Applies 2-D filter to matrix by filtering rows and columns.
%  Y = FILT2D(X,N,F) applies 1-D Butterworth filter to all columns and then
%  to all rows. Since these filters are linear, order of filtering does not
%  matter and symmetry of matrices is preserved.
%
%  n+1 is the length of the filter f the cutoff frequency, in units of
%  Nyquist (1/2 sample frequency).
%  To filter out point-to-point ripples, n=2, f=.5-.75 suffices
%  The columnwise filter routine is MATLAB's filtfilt, which preserves phase
%  and works well at ends of interval.

%  addapted from 'filt2d.m' by U. Send, IfM Kiel
%  C. Mertens, IfM Kiel
%  $Revision$ $Date$

[b,a] = butter(n,f);
[m,n] = size(x);
y = NaN*ones(m,n);
for j= 1:n
  y(:,j) = filtfilt(b,a,x(:,j));
end
for i= 1:m
  y(i,:) = (filtfilt(b,a,(y(i,:))'))';
end

