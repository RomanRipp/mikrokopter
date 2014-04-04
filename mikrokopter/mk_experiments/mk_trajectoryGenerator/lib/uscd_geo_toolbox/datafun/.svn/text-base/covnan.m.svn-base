function [x,c] = covnan(x,flag)

% COVNAN   Covariance matrix and Correlation Coefficients, ignores NaN
%
% [ C , CC ] = COVNAN( X , Flag )
%
%   If X is a vector, C is the variance.  
%   If X is a matrice, where each row is an observation, 
%                        and each column a variable,
%                     C is the covariance matrix.
%  
%           DIAG(COVNAN(X))  is a vector of variances for each column, 
%      SQRT(DIAG(COVNAN(X))) is a vector of standard deviations.
%
%   Flag = 0  normalizes by (N-1) where N is the number of observations.  
%            This makes COVNAN(X) the best unbiased estimate of the covariance matrix 
%             if the observations are from a normal distribution.
%            (default)
%
%   Flag = 1  normalizes by N and produces the second moment matrix 
%              of the observations about their mean.
%
%  CC is a matrix of correlation coefficients for the variables
%
%  CC(i,j) = C(i,j) / SQRT( C(i,i) * C(j,j) )
%
%
% See also: COV, CORRCOEF, STD, VAR
%

Nout = nargout;

%*****************************************
% Check Inputs

if nargin < 1
   error('Not enough InputArguments.')
end

if nargin < 2
   flag = 0;
end

%-----------------------------------------
% Check X

if ~( isnumeric(x)  &  ~isempty(x) &  ...
      ( size(x,1)*size(x,2) == prod(size(x)) ) )
  error('X must be a 2-dimensional , nonempty Numeric.');
end

% RowVector --> ColumnVector

if ( size(x,1) == 1 )  & ( size(x,2) > 1 );
   x = permute(x,[2 1]);
end

m = size(x,1);

%*****************************************
% Check for NaN's

is_nan = any(isnan(x(:)));

%-----------------------------------------
if is_nan

   %--------------------------------------
   % Remove Rows with full NaN 
   %  from Begin and End in 1. Dimension

   n  = double(isnan(x));  % True for NaN

   i0 = 1 + min( sum(cumprod(n,1),1) );

   i1 = m - min( sum(cumprod(n(m:-1:1,:),1),1) );

   x = x(i0:i1,:);
   
   n = n(i0:i1,:);

   m = size(x,1);

   is_nan = any(n(:));

end

%-----------------------------------------
if is_nan

   %--------------------------------------
   % Mask NaN's with ZERO

    % NaN --> ZERO

    x(find(n)) = 0;

      n = ~n;          % True for NOT NaN

     sc = n' * n;      % Scale for Normalisation

     sc = sc + ( sc == 0 ); 

     ok = n;           % True for NOT NaN

    %-------------------------------------
    % Number of NOT-NaN-Elements per Column

    n = sum(n,1);

    n = n + ( n == 0 );

%-----------------------------------------
else

     ok = 1;           % Always TRUE

      n = m;           % Number of Observations

     sc = m;           % Scale for Normalisation

end

%*****************************************

%-----------------------------------------
% Special Case 

if m == 1
   x = 0;
   c = NaN;
   return
end

%-----------------------------------------
% Remove Mean

x = x - ones(m,1) * ( sum(x,1) ./ ( n + ( n == 0 ) ) );

x = x .* ok;   % NaN --> ZERO

%-----------------------------------------
% Build Matrice

x = x' * x;

%-----------------------------------------
% Normalize

sc = sc - isequal( flag , 0 );

sc = sc + ( sc == 0 );

x   = x ./ sc;

%-----------------------------------------

if Nout < 2
   return
end

%-----------------------------------------
% CorrCoeff

c = diag(x);  % Variance

c = sqrt( c * c' );

c(find(c==0)) = NaN;

c = x ./ c;

