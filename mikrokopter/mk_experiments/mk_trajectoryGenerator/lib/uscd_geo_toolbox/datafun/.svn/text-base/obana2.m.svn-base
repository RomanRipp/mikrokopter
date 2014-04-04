function y=obana(x,posx,posy,grdx,grdy,r1x,r2x,r1y,r2y);

% OBANA3 fast data gridding using Gaussian weights. 
%
%   ZI = obana2(Z,X,Y,XI,YI,xrad,xcut,yrad,ycut)
%
% influence and cut-off radii may be given as vectors or matrices
%
%  uses :	gauss.m sumnan.m
%
% version 1.1.0		last change 19.06.1997

% M. Visbeck
%	changed to variable influence and cut-off radii
%	G.Krahmann, IfM Kiel, Jun 1997
1
% make isotropic if no y values are specified
Nin = nargin;

if Nin < 8, r1y=r1x; end
if Nin < 9, r2y=r2x; end


% blow up influence and cut-off radii
if size(r1x,2)==1
  r1x=r1x*ones(1,size(grdx,2));
end
if size(r2x,2)==1
  r2x=r2x*ones(1,size(grdx,2));
end
if size(r1x,1)==1
  r1x=ones(size(grdx,1),1)*r1x;
end
if size(r2x,1)==1
  r2x=ones(size(grdx,1),1)*r2x;
end
if size(r1y,2)==1
  r1y=r1y*ones(1,size(grdx,2));
end
if size(r2y,2)==1
  r2y=r2y*ones(1,size(grdx,2));
end
if size(r1y,1)==1
  r1y=ones(size(grdx,1),1)*r1y;
end
if size(r2y,1)==1
  r2y=ones(size(grdx,1),1)*r2y;
end

r1x=r1x(:)';
r2x=r2x(:)';
r1y=r1y(:)';
r2y=r2y(:)';

% setup input and target vector positions

si = size(grdx);

posx = posx(:)';
posy = posy(:)';

grdx = grdx(:)';
grdy = grdy(:)';

x = x(:);

m = size(grdx,2);

% reset output

y  = NaN * zeros(1,m);

% loop over each output value

pc = 0;

fprintf(1,'\nOBANA: ');

for ii = 1 : m

  % display process
  if ( ii/m > pc ) 
     fprintf(1,'\rOBANA: %3.0f%%',pc*100);
     pc = pc+0.1; 
  end

  % positions difference
  dx = grdx(ii) - posx;
  dy = grdy(ii) - posy;

  % norm with cutoff radius
  jj = sqrt( (dx/r2x(ii)).^2 + (dy/r2y(ii)).^2 );

  % select only values within cutoff radius
  jj = find( abs(jj) < 1 );

  if ~isempty(jj)

    % norm with inflence radius
    d = sqrt( (dx(jj)/r1x(ii)).^2 + (dy(jj)/r1y(ii)).^2 );

    % get factors using a gauss distribution
    d = gauss(d,1);

    % sum up values
    s = sumnan(d);

    if s > 0
       y(ii) = ( d * x(jj) ) / s;
    end

  end

end

fprintf(1,'\n\n');

% reshape to GridSize

y = reshape(y,si);
