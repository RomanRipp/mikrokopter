function [fs,ssi,ppi,is_win,fn] = bestfont

% BESTFONT  Returns an optimal FontSize for this computer
%
% [ FontSize , ScreenSize , ScreenPixelsPerInch ] = BESTFONT
%
% [ ... , IsWin , FixedWidthFontName ] = BESTFONT
%
% returns true for PCWIN-System and the FixedWidthFontname
%

is_win = strcmp( upper(computer) , 'PCWIN' );

uni = get(0,'units');       
      set(0,'units','pixels')
ssi = get(0,'ScreenSize');  
      set(0,'units',uni);
          
ppi = get(0,'ScreenPixelsPerInch');

is_tall = -1 + ( ssi(4) >=  480 ) + ...
               ( ssi(4) >=  600 ) + ...
             1*( ssi(4) >= 1050 );

fs =  8 + 2 * is_tall - 1 * is_win;

if isunix
   fs = fs + ( 14 - fs ) * strcmp(get(0,'terminalprotocol'),'none');
end

if is_win
   fn = 'courier';
else
   fn = { 'arrial'  get(0,'fixedwidthfontname') };
   fn = fn{ 1 + ( ssi(4) >= 1050 ) } ;
end
