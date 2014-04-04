function [C,varargout] = cor4flat(Mode,C,varargin)

% COR4FLAT  Correct Matrices to use SURFACE with Flat shading
%
%  [ C , X , Y , Z ] = COR4FLAT( Mode , C , X , Y , Z );
%
%  C = CData, X = XData, Y ...
%
%   Mode = 'in'   C will set to mean of Vertex
%                   No change in Size of Matrice
%
%   Mode = 'out'  X, Y, Z will be interpolate to Mean of Vertex
%                  and 1 More at the Borders, 
%
%
% Inputs:  C = M x N x K
%          X = M x N  or  1 x N   or   1 x 1 
%          Y = M x N  or  M x 1   or   1 x 1  
%          Z = M x N  or  M x 1   or   1 x 1 
% 
% Ouputs:  C,X,Y,Z = M   x  N     for  mode == 'in'
%          C,X,Y,Z = M+1 x  N+1   for  mode == 'out'
%
%
%  after call of COR4FLAT use
%
%   h = surface( 'xdata'    , X      , ...
%                'ydata'    , Y      , ...
%                'zdata'    , Z      , ...
%                'cdata'    , C      , ...
%                'facecolor', 'flat' , ...
%                'edgecolor', 'none'        );
%
%   Warning:  Don't use  FaceColor Interp   !!!!!!
%


Nin  = nargin  - 2;
Nout = nargout - 1;

varargout = cell(1,Nout);

if Nin < 0
  error('Inputs Mode and C are missing.');
elseif Nin > 3
  error('Too many input arguments.')
end

nn = min(3,Nout);

s1 = size(C,1);
s2 = size(C,2);

%******************************************************
% Check X Y Z

% Defaults for X Y Z
 def = { ( 1 : s2 ) 
         ( 1 : s1 )
         0           };

for ii = 1 : nn

 if Nin < ii
    v = def{ii};
 else
    v = varargin{ii};
 end

 if any( size(v) == 1 )
    v = v(:);
    if  all( size(v) == 1 )
        v =  v * ones(s1,s2); 
    elseif ( size(v,1) == s1 )  &  any( ii == [2 3] )
        v =  v * ones(1,s2);
    elseif ( size(v,1) == s2 )  &  ( ii == 1 )
        v =  ones(s1,1) * v' ; 
    end
 end

  if ~isequal( [s1 s2] , size(v) )
    error(sprintf('Size of %.0f. Input must match Size of C.',ii));
  end

  varargout{ii} = v;

end

%******************************************************

cl = class(C);

if ~strcmp(cl,'double');
    C = double(C);
end


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Length > 300000  ==>  Use parts off Length == 250000 == ( 500 x 500 )  

st1 = (s1+1) + ( 500 - (s1+1) ) * ( s1*s2 > 3*1e5 );
st2 = (s2+1) + ( 500 - (s2+1) ) * ( s1*s2 > 3*1e5 );


i1 = [ ( 1 : st1 : s1 )  s1+1 ];
i2 = [ ( 1 : st2 : s2 )  s2+1 ];



%******************************************************
if strcmp(Mode,'in')
%******************************************************

  for ii1 = 1 : size(i1,2)-1
    for ii2 = 1 : size(i2,2)-1

       ind1 = ( i1(ii1) : i1(ii1+1)-2 );
       ind2 = ( i2(ii2) : i2(ii2+1)-2 );

      C(ind1,ind2,:) = ( C( ind1+0 , ind2+0 , : ) + ...
                         C( ind1+0 , ind2+1 , : ) + ...
                         C( ind1+1 , ind2+0 , : ) + ...
                         C( ind1+1 , ind2+1 , : )  ) / 4;

    end
  end

%******************************************************
else
%******************************************************

  sc = size(C);

  sc1    = sc;
  sc1(1) = 1;

  sc2    = sc;
  sc2(1) = sc2(1)+1;
  sc2(2) = 1;

  C = cat( 1 , C , ones(sc1) );
  C = cat( 2 , C , ones(sc2) );
  
  C(:,s2+1,:) = C(:,s2,:);
  C(s1+1,:,:) = C(s1,:,:);

  for ii = 1 : nn

    X0 =  varargout{ii};
    X1 =  NaN*ones(s1+1,s2+1);

     % Inner

     for ii1 = 1 : size(i1,2)-1
       for ii2 = 1 : size(i2,2)-1

          ind1 = ( i1(ii1) : (i1(ii1+1)-1) );
          ind2 = ( i2(ii2) : (i2(ii2+1)-1) );

          ind1 = ind1( 1 : ( end - ( ind1(end) == s1 ) ) );
          ind2 = ind2( 1 : ( end - ( ind2(end) == s2 ) ) );
[ind1 NaN ind2 ]

         X1(ind1+1,ind2+1) = ( X0( ind1+0 , ind2+0 ) + ...
                               X0( ind1+0 , ind2+1 ) + ...
                               X0( ind1+1 , ind2+0 ) + ...
                               X0( ind1+1 , ind2+1 )  ) / 4;

       end
     end

     % Border's

     X1(2:s1,   1) = ( X0( 1:s1-1 ,  1 ) + ...
                       X0( 2:s1-0 ,  1 )  ) / 2;
     X1(2:s1,s2+1) = ( X0( 1:s1-1 , s2 ) + ...
                       X0( 2:s1-0 , s2 )  ) / 2;

     X1(   1,2:s2) = ( X0(  1 , 1:s2-1) + ...
                       X0(  1 , 2:s2-0)  ) / 2;
     X1(s1+1,2:s2) = ( X0( s1 , 1:s2-1) + ...
                       X0( s1 , 2:s2-0)  ) / 2;

     % EdgePoints

     X1([ 1 s1+1 ((s2+1)-1)*(s1+1)+1 (s1+1)*(s2+1)] ) = ...
     X0([ 1 s1+0 ((s2+0)-1)*(s1+0)+1 (s1+0)*(s2+0)] );
  
     varargout{ii} = X1;                         

  end

%******************************************************
end
%******************************************************

if ~strcmp(cl,'double');
    C = feval(cl,C);
end


                