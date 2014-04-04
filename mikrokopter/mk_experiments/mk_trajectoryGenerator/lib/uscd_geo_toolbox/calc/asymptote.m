function y = asymptote(x,varargin);

% ASYMPTOTE Asymptotic Approximation
%
% Y = ASYMPTOTE( X , X0 , Mode , Parameter )
%
% Y( X < X0 ) = 0;
% Y( X > X0 ) = X - X0
%
% Mode = 'cos' | 'log' | 'hyp'
%
% Parameter:  finite larger ZERO
%
% Either X or Parameter must be a Skalar or both a Vector.
%

Nin = nargin;

x0 = [];
md = '';
a  = [];

%******************************************************
% Check Inputs

msg = {};

if Nin == 0
   error('Not enough Input Arguments.');
end

sx = size(x); px = prod(sx); cx = class(x);

if ~strcmp(cx,'double')
    msg = cat( 1 , msg , {'X must be double.'} );
end

vx = ( px == max(sx) );  % True for Vector X

for ii = 1 : Nin-1

    m = '';

    v = varargin{ii};

    s = size(v); p = prod(s); c = class(v);

    ok = ( ~isempty(v) & strcmp(c,'char') & ...
           ( p == s(2) ) &  isempty(md)       );
    if ok
       md = v;
    elseif ~isempty(v) & strcmp(c,'double') 
       ok = ( isempty(x0) & ( p == 1 ) );
       if ok
          x0 = v;
       elseif isempty(a) 
          ok = ( ( px == 1 ) | ( p == 1 ) | ...
                 ( vx & ( p == max(s) ) ) );
          if ~ok
              m = 'Either X or Parameter must be a Skalar or both a Vector.';
          else
             ok = all( isfinite(v(:)) & ( v(:) > 0 ) );
             if ~ok
                 m = 'Parameter must be finite larger ZERO';
             end
          end
          if ok
             a = v;
             if ( ( px > 1 ) & ( p > 1 ) & ~isequal(sx,s) )
                [x,a] = meshgrid(x,a);
             elseif ( ( px == 1 ) & ( p > 1 ) )
                x = x + zeros(size(a));
             end
          end
       end
    end

    if ~ok
        str = sprintf('Invalid or unrecognized %.0f. Input',ii+1);
        if ~isempty(m)
            str = sprintf('%s: %s',str,m);
        end
        str = sprintf('%s.',str);
        msg = cat( 1 , msg , {str} );
    end

end

%---------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs\n%s',msg);
    error(msg)
    return
end

%******************************************************

if isempty(x0); x0 =  0;  end
if isempty(md); md = 'c'; end
if isempty(a);  a  =  1;  end

if isempty(x)
   y = x;
   return
end

if isinf(x0)
    if x0 > 0
       y = 0;
    else
       y = inf;
    end
    y = y * x;
    return
end

m = lower(md(1));

% COS: 2*a / LOG: log(10)*a
% a = a * ( 1 + ( m == 'c' ) + ( log(10) - 1 ) * ( m == 'l' ) );


x = a .* ( x - x0 );

switch m

   case 'c'

     d  = ( 1 - pi/2 ); 
     y  = ( 1 - cos( x - d ) ); 

     ii = find( x < d ); y(ii) = 0;
     ii = find( x > 1 ); y(ii) = x(ii);

   case 'l'

     y = log( exp(x) + 1 );

   case 'h'

     y = ( sqrt( 1 + (x+1).^2 ) - 1 );  % Hyperbel

     ii = find( x < -1 ); y(ii) = 0;

   otherwise

     error(sprintf('Invalid Mode "%s".',md));

end

y = y ./ a;

  
