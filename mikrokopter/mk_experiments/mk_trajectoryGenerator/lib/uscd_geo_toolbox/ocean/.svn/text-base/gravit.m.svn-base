function g = gravit(lat);

% GRAVIT  Gravitation Acceleration from Latitude
%
% GRAVIT( Latitude )
%
%-------------------------------------------------------------------------------
%
% gravit        - Calcul de l'acceleration de gravite en fonction
%                 de la latitude en degres decimaux. (formule GRS-80)
%
%-------------------------------------------------------------------------------
% Version:
% -------
%  1.01 Criation (d'aprhs gravi8, chaine hydro)          14/06/94 F. Gaillard
% 
%-------------------------------------------------------------------------------

lat   = lat * pi/180;

lat2  = sin( 2 * lat );
lat   = sin( 1 * lat );

ge     =  9.780318;
b2     =  0.530244e-02;
b4     = -0.585000e-05;

g = ge * ( 1 + b2 * lat.*lat + b4 * lat2.*lat2 );
