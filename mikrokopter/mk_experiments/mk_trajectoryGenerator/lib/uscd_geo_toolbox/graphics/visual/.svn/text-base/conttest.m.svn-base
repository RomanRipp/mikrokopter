% CONTTEST  Example for CONTFILL with 2 ColorMaps in a figure

 fig = figure;

cmaps = { gray(128)  jet(128) };

axe = zeros(2,1);

for ii = 1 : 2

    cm = cmaps{ii};        % Die ColorMap

    axe = subplot(2,1,ii);

    [cs,h] = contfill(rand(10,10));  % h sind die Handles zu den Patches

    % caxis([0.2 0.8])

    hold on % spaetestens hier !!!

    set(axe,'layer','top');  % Sieht besser aus

    %-----------------------------------------------------------------
    % Scaled Colors ---> True Colors

    clm = get(axe,'clim');   % das ColorLimit, damit in die ColorMap skalieren

    xc = linspace(clm(1),clm(2),size(cm,1))';

    for hh = h(:)'
        cc = get(hh,'cdata');
        if ~isnan(cc)
            cc = min(max(cc,clm(1)),clm(2)); % Check for Out of Range
            cc = interp1(xc,cm,cc);          % Color [ R G B ]
            set(hh,'facecolor',cc)
        end
    end

    %-----------------------------------------------------------------
    % ColorBar

    axc = colorbar('horiz');           % Die ColorBarAxe

    hc = findobj(axc,'type','image');  % Der bunte Streifen im ColorBar

    cc = get(hc,'cdata');   % ColorIndize oder FarbWerte zum skalieren

    if strcmp(get(hc,'cdatamapping'),'scaled')

       clm = get(axc,'clim');
       xc = linspace(clm(1),clm(2),size(cm,1))';
       cc = min(max(cc,clm(1)),clm(2)); % Check for Out of Range
       prm = ( size(cc,1) == 1 );       % ColorBar horiz
       cc = interp1(xc,cm,cc);    
       cc = permute(cc,[1 3 2]);
       if prm
          cc = permute(cc,[2 1 3]);
       end

    else   %%% CDataMapping indexed

       nc = size(cm,1);
       cc = cat(3,cc,cc+nc,cc+2*nc);
       cc = cm(cc);

    end

    % Nu ist cc eine 3D-Matrix mit RGB-FarbWerten   

    set(hc,'cdata',cc); 


    drawnow


end     
