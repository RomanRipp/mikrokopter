function [strlat,strlon,strlat2,strlon2] = pos2str(pos)
%POS2STR Converts a decimal position into a position string.
%  [STRLAT,STRLON,STRLAT2,STRLON2] = POS2STR(POS)
%
%  Input: - pos                  : decimal position [lat,lon] in degrees
% 
%  Output: - strlat              : position string of latitude ## ##.### N
%          - strlon              : position string of longitude ### ##.### E
%          - strlat2             : position string of latitude ## ##.##N
%          - strlon2             : position string of longitude ### ##.##E

%  Gerd Krahmann, IfM Kiel, Jun 1995, last change 22.04.1996
%  C. Mertens, IfM Kiel, minor changes
%  F. Morsdorf, IfM Kiel, changed input 
%  $Revision: 1.2 $ $Date: 1999/8/02 09:05:40 $

[m,n] = size(pos);
estrlon = '';estrlat = '';estrlon2 = '';estrlat2 ='';
for i = 1:m
latd=fix(pos(i,1));
latm=(pos(i,1)-latd)*60;
lond=fix(pos(i,2));
lonm=(pos(i,2)-lond)*60;
if latd<0 | latm<0
  ns='S';
  nsv=-1;
else
  ns='N';
  nsv=1;
end
if lond<0 | lonm<0
  ew='W';
  ewv=-1;
else
  ew='E';
  ewv=1;
end

% form latitude string
if abs(latd)<10
  strlat=['0',int2str(latd*nsv),' '];
else
  strlat=[int2str(latd*nsv),' '];
end
if abs(latm)<10
  if latm==0
    strlat=[strlat,'00.000'];
  else
    strlat=[strlat,'0',sprintf('%5.3f',latm*nsv)];
  end
else
  strlat=[strlat,sprintf('%6.3f',latm*nsv)];
end
if isempty(find(strlat=='.'))
  strlat=[strlat,'.'];
end
l=length(strlat);
while l<9
  strlat=[strlat,'0'];
  l=length(strlat);
end
strlat = [strlat,' ',ns];
estrlat = strvcat(estrlat,strlat);
strlat = estrlat;
% form longitude string
if abs(lond)<10
  strlon=['00',int2str(lond*ewv),' '];
elseif abs(lond)<100
  strlon=['0',int2str(lond*ewv),' '];
elseif abs(lond)>=100
  strlon=[int2str(lond*ewv),ew];
end
if abs(lonm)<10
  if lonm==0
    strlon=[strlon,'00.000'];
  else
    strlon=[strlon,'0',sprintf('%5.3f',lonm*ewv)];
  end
else
  strlon=[strlon,sprintf('%6.3f',lonm*ewv)];
end
if isempty(find(strlon=='.'))
  strlon=[strlon,'.'];
end
l=length(strlon);
while l<10
  strlon=[strlon,'0'];
  l=length(strlon);
end
strlon = [strlon,' ',ew];
estrlon = strvcat(estrlon,strlon);
strlon = estrlon;
% make second format
strlat2=[strlat(1:2),setstr(176),strlat([4:8,3])];
strlon2=[strlon(1:3),setstr(176),strlon([5:9,4])];
if strlat2(1)=='0'
  strlat2(1)=' ';
end
if strlon2(1)=='0'
  strlon2(1)=' ';
end
if strlon2(2)=='0'
  strlon2(2)=' ';
end
estrlon2 = strvcat(estrlon2,strlon2);
strlon2 = estrlon2;
estrlat2 = strvcat(estrlat2,strlat2);
strlat2 = estrlat2;

end
