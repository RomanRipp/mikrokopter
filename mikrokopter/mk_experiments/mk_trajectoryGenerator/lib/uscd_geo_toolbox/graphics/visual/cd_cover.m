function cd_cover(file,image1,txt,image2)

% CD_COVER  CoverImage for a CompactDisk
%
% CD_COVER(output_file,image,txt,image_back)
%

if nargin == 1
  image1 = [file '.jpg'];
  image2 = [file '.jpg'];
  txt = [file '.txt'];
end


[IF,CMAP] = imread(image1,image1(end-2:end));

if isempty(CMAP);
  CMAP = [ 0 0 0 
           1 1 1 ];
end

figure(1)

clf

[IF,CMAP] = imread(image1,image1(end-2:end));

set(gcf, 'PaperUnits', 'centimeters', ...
         'PaperOrientation', 'Landscape', ...
         'PaperPosition', [0 0 29.6774 20.984], ...
         'colormap' , CMAP )


image([0 .5], [0 1], IF)

fid1=fopen(txt);
tt=fscanf(fid1,'%c');
fclose(fid1);
cc = find(tt==char(10));

tit = tt(1:cc(1));
son = tt;
% son(1:cc(1)) = [];  % Remove first Title-Line from Text *** !!! ***

text(.55,.15,son, 'VerticalAlignment', 'top')

set(gca, 'Units', 'centimeters', ...
         'Position', [2.85 4.5 24 12], ...
         'XTick', [], ...
         'YTick', [], ...
         'Box', 'on', ...
         'XLim',[0 1])

wysiwyg
print(1,'-dpsc2', [file '.ps1'])

figure(2)

clf
hold on

set(gcf, 'PaperUnits', 'centimeters', ...
         'PaperType', 'A4', ...
         'PaperOrientation', 'Landscape', ...
         'PaperPosition', [0 0 29.6774 20.984])

if exist('image2')

  IF = imread(image2,image2(end-2:end));
  image([0 13.8], [0 1], IF)
  set(gca,'YDir', 'Reverse')

else

  set(gca,'YDir', 'Reverse')

end

plot([0 0],[0 1],'--k')
plot([13.8 13.8],[0 1],'--k')

text(0,.95, tit, 'Rotation', 90, 'FontWeight', 'bold')
text(13.9,.05, tit, 'Rotation', -90, 'FontWeight', 'bold')

text(6.9,.15,tit, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'EraseMode', 'xor')
text(6.9,.25,son, 'HorizontalAlignment', 'center',  'VerticalAlignment', 'top', 'EraseMode', 'xor')

set(gca, 'Units', 'centimeters', ...
         'Position', [2.85 4.5 15.1 11.8], ...
         'XTick', [], ...
         'YTick', [], ...
         'Box', 'on', ...
         'XLim',[-0.6 14.4], ...
         'YLim', [0 1])

wysiwyg
print(2,'-dpsc2', [file '.ps2'])