function h = menupos(action,axe,h);

% MENUPOS creates a Menu and Show X-Y_Position
%
% MenuHandle = MENUPOS('new',axe,Form)
%
%   creates Menu and call MENUPOS('on',axe,MenuHandle)
%
%
%  MENUPOS('on'    ,axe,MenuHandle)    
%  MENUPOS('off'   ,axe,MenuHandle) 
%  MENUPOS('delete',axe,MenuHandle) 
%
% internal for FigureWindowButtonMotionFcn:
%
%    MENUPOS('move'  ,axe,MenuHandle)
%
%

 Nin = nargin;
 if Nin < 1
   action = 'new';
 end
 if Nin < 2
   axe = gca;
 end

 
 fig = get(axe,'parent');


switch action

case 'new'


  if Nin < 3
    form = '%10.3f';
  else
    form = h;
  end
  
  h = uimenu('parent'  , fig , ...
             'userdata', form      );

  menupos('on',axe,h);

case { 'on'  'off' }

   clear eps
   Hform = sprintf( '%%.%.0ff' , ceil(abs(log(eps)/log(10)))+1 );

     Axe = sprintf(Hform,axe);
     H   = sprintf(Hform,h);

   % Command 
   comm = ['menupos(''move'','  Axe ','  H ');'];   

   % Old MotionFcn
   CB = get(fig,'windowbuttonmotionfcn');

   % Remove existing Command from MotionFcn
   ii = findstr(CB,comm);
   if ~isempty(ii)
     ii = (-sort(-ii(:)))';
     for iii = ii 
         CB((1:size(comm,2))+iii-1) = [];
     end
   end


   if strcmp(action,'on')
   % Add Command to MotionFcn
      if ~isempty(CB)
        if ~any(strcmp(CB(size(CB,2)),{';' ','}))
          CB = [ CB  ',' ];
        end
      end
      CB = [ CB  comm  ];
   end

   set(fig,'windowbuttonmotionfcn', CB );

   
case 'move'

  cp = get(axe,'currentpoint');
  
  set(h,'label',sprintf( get(h,'userdata') , ...
                           cp(1,[1 2]) ) );

case 'delete'

  menupos('off',axe,h);

  delete(h)

end
   
   