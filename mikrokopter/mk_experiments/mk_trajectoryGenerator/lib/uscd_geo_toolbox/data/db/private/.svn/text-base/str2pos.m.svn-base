function [pos,degs,mins] = str2pos(posstr)
%STR2POS Converts a position string into a decimal position.
%  [POS,DEGS,MINS] = STR2POS('posstr')
%
%  Input: - posstr              : string with position information
%                                 may look like
%                                       43.1234
%                                   or  43N07
% 
%  Output: - pos                : decimal position [degrees]
%          - degs               : full degrees
%          - mins               : decimal minutes
%                                 degs and mins are negative for W or S

%  Gerd Krahmann, IfM Kiel, Jan 1995, last change 31.01.1995
%  C. Mertens, IfM Kiel, minor changes
%  $Revision: 1.1 $ $Date: 1996/12/03 09:07:56 $

% check which format is used
posstr = upper(posstr);
good=find( (posstr=='N') | (posstr=='S') | (posstr=='E') | (posstr=='W') );
if isempty(good)
  pos=str2num(posstr);
  degs=fix(pos);
  mins=(pos-degs)*60;
else
  degs=str2num(posstr(1:good-1));
  if length(degs) > 1
    mins = degs(2);
    degs = degs(1);
  else
    mins=str2num(posstr(good+1:length(posstr)));
  end
  if isempty(mins)
    mins=0;
  end
  pos=degs+mins/60;
  if ( (posstr(good)=='S') | (posstr(good)=='W') )
    degs=-degs;
    mins=-mins;
    pos=-pos;
  end
end

