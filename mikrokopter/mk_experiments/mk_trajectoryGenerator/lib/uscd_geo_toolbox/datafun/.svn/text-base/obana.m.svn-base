function y=obana(x,posx,posy,gridx,gridy,r1x,r2x,r1y,r2y);

% OBANA simple objective mapping
%
% ZI = OBANA( Z , Y , X , YI , XI , R1x,R2x,R1y,R2y);
%
% Z:  values
% X, Y  positions of each value (nonuniform)
%
% XI, YI  target positions  of each value  (grid)
%
% r1: influence radius (in x and y direction)
% r2: cutoff radius    (in x and y direction)
%
% see also: OBANA3, GRIDDATA, OBJMAP
%
% M. Visbeck

% make isotropic if no y values are specified
if nargin<9, r1y=r1x;
 if nargin<8, r2y=r2x; end
end

% setup complex target vector positions
i=sqrt(-1);
[ly,lx]=size(gridx);
gridv=gridx(:)'+i*gridy(:)';


% setup complex input vector positions
xv=x(:)';
ld=length(xv);
posv=posx(:)'+i*posy(:)';

% reset sums
nn=zeros(1,ly*lx);
yv=zeros(1,ly*lx);

% loop over each input value
t0=clock;
iper=0;
for j=1:ld

 per=j./ld.*100;
 if floor(per)>=25 & iper < 1
  disp(['25%  ', sprintf('%6.1f min',etime(clock,t0)./60)])
  iper=1;
 elseif floor(per)>=50 & iper <2
  disp(['50%  ', sprintf('%6.1f min',etime(clock,t0)./60)])
  iper=2;
 elseif floor(per)>=75 & iper <3
  disp(['75%  ', sprintf('%6.1f min',etime(clock,t0)./60)])
  iper=3;
 end

 % positons difference
 dp=gridv-posv(j);

 % norm with cutoff radius
 cr=real(dp)/r2x+i*imag(dp)/r2y;

 % select only values in cutoff radius
 ix=find( abs(cr)<1);
 if length(ix)>0

  % norm with inflence radius
  ir=real(dp(ix))/r1x+i*imag(dp(ix))/r1y;

  % get factors using a gauss distribution
  n=gauss(abs(ir),1);

  % sum up weights and values
  nn(ix)=nn(ix)+n;
  yv(ix)=yv(ix)+n.*xv(j);

 end
end

% devide by sum of weights
ii=find(nn>0);
yv(ii)=yv(ii)./nn(ii);
ij=find(nn==0);
yv(ij)=yv(ij).*nan;

% reshape to taget position size
y=reshape(yv,ly,lx);













