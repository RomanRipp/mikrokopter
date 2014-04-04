function [fig,hl,ll] = plot_gebco(xy,i0,lg,lv,mm);

% PLOT_GEBCO  Plots GEBCO-IsoLines, read by READ_GEBCO
%
% [Fig,H,L] = PLOT_GEBCO( XY , I0 , LG , LV , [MXY] )
%
%   XY = [ Lon  Lat ], NaN-separated Segments
%   I0 = StartIndex of Segment in XY
%   LG = Length of Segment
%   LV = Level of Segment
%  MXY = [ MX  MY ] of Segment (optional)
% 
%  Fig = FigureHandle
%    H = [ LineHandle  [MeanPointHandle] ],  [ N by 1|2 ]
%    L =   Levels (unique)                   [ N by 1   ]
%
 
Nin = nargin;

ll = sort(lv);
ll(find(diff(ll,1,1)==0)+1) = [];

nc = size(ll,1);

cl = terrain(3*nc);
cl = cl(1:nc,:);
cl = cl(nc:-1:1,:);

cl(find(ll==0),:) = 0;

cm = rgb2hsv(cl);
cm(:,1) = cm(:,1) + 1/3;
cm(:,1) = cm(:,1) - floor(cm(:,1));
cm = hsv2rgb(cm);


hl = zeros(nc,1+(Nin==5));

figure('menubar','none','toolbar','none')
hold on

for il = 1 : nc

    ii = find( lv == ll(il) );
    jj = grp2ind(i0(ii),lg(ii));

    hl(il,1) = plot(xy(jj,1),xy(jj,2),'color',cl(il,:));

    if Nin == 5

       hl(il,2) = plot(mm(ii,1),mm(ii,2), ...
                       'linestyle','none', ...
                       'marker' , '.' , ...
                       'markersize' , 10 , ...
                       'color' , cm(il,:)      );
    end
end

try

   mercator(gca);

   axe_zoom('new',gca,'mercator'); 
   axe_zoom('on',gca);

end

