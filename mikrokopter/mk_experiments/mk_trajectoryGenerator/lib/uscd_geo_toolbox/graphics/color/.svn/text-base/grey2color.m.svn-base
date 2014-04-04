function m = grey2color(g,c);

% GREY2COLOR  Transforms a GreyMap to a ColorBase
%
% ColorMap = GREY2COLOR( GreyMap , BaseColor )
%
% GreyMap   = UINT8 or in Range [ 0 .. 1 ]
%             [ N by 1 ] GreyScale, 
%             [ N by 3 ] RGB-Map, 
%                        Grey = 0.299*R + 0.587*G + 0.114*B
%
% BaseColor = UINT8 or in Range [ 0 .. 1 ]
%             [ 1 by 3 ] RGB-Tripel
%          or Matlab ColorSpecifier or XRGB-ColorName
%
% ColorMap  = [ N by 3 ] ColorMap based on BaseColor
%
% see also: RGB2GREY, COLSPEC, XRGB, XRGBROWSE
%
% Example:
%
%% MATLABROOT/toolbox/matlab/demos/durer.mat
%%
%% X              648 x 509    2638656   double
%% caption          2 x 28         112   char
%% map            128 x 3         3072   double
% 
% load durer
%
% figure('colormap',map)
% image(X,'cdatamapping','direct')
% axis image 
%
% figure('colormap',grey2color(map(:,1),[0 0 0.5]))
% image(X,'cdatamapping','direct')
% axis image 
%
%%---------------------------------------------------
%
% load clown
%
% figure('colormap',map)
% image(X,'cdatamapping','direct')
% axis image 
%
% figure('colormap',grey2color(map,'purple'))
% image(X,'cdatamapping','direct')
% axis image 
%


m = zeros(0,3);

%********************************************************************
% Check Inputs

msg = {};

%---------------------------------------------
% Check GrayMap

if isempty(g)
   return
end


cl = class(g);
ok = ( ( ndims(g) == 2 ) & any( size(g,2) == [ 1  3 ] ) & ...
        any(strcmp(cl,{'double'  'unit8'})) );
if ok
   if strcmp(cl,'uint8')
      g = double(g)/255;
   else
     ok = all( abs(g(:)-0.5) <= 0.5 );
   end
   if  ok & ( size(g,2) == 3 )
       g = g * [ 0.299 ; 0.587 ; 0.114 ];
   end
end

if ~ok
    msg = cat(1,msg,{'GrayMap must be UINT8 or in Range of 0 .. 1.'});
end

%---------------------------------------------
% Check Color

if ischar(c)

   try
      c = colspec(c,1);
      if isempty(c)
         msg = cat(1,msg,{'Unknown Color by COLSPEC.'});
      end
   catch
      msg = cat(1,msg,{sprintf('Error get Color by COLSPEC.\n%s',lasterr)});
   end

else

   cl = class(c);
   ok = ( isequal(size(c),[1 3]) & any(strcmp(cl,{'double'  'unit8'})) );
   if ok
      if strcmp(cl,'uint8')
         c = double(c)/255;
      else
        ok = all( abs(c-0.5) <= 0.5 );
      end
   end

   if ~ok
       msg = cat(1,msg,{'Numeric Color must be a RGB-Tripel in UINT8 or in Range of 0 .. 1.'});
   end

end

%---------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%********************************************************************

c = rgb2hsv(c);

n = size(g,1);

m = zeros(n,3);

m(:,1) = c(:,1);  % Hue

m(:,2) = c(2) * ( 1 - g );   % Sat

m(:,3) = c(3) + ( 1 - c(3) ) * g;

m = hsv2rgb(m);

