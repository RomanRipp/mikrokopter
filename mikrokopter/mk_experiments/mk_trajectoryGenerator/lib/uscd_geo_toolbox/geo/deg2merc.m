function  y = deg2merc(y,ext)

 %  function  y = deg2merc(y,ext)

 if ext 
  y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;
 end
