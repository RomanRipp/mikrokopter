function [Msg , H , CB_down ] = rectmove( action   , Hin , sets , ...
                                         CB_motion , CB_up , ...
                                         color,linewidth,linestyle)

% RECTMOVE  draws and move a Rectangle
%
% [Msg, RectHandle] =  RECTMOVE( Action , Handle , Sets )
%
% Action      Handle      Sets
%
% 'new'       AxeHandle  [ CentreX  CentreY  ExtX  ExtY  Rotation , ...
%                           Xmin  Xmax  Ymin  Ymax  , ...
%                          1 , ...     % enable Moving of Centre
%                          2 , ...     % enable Rotation
%                          3       ]   % enable Changing of Extension
%                        []  takes defaultvalues for Sets !!!
%
% For action 'new'  more Options are:
% 
%   RECTMOVE( ... , ButtonUpFcn , ButtonMotionFcn , ...
%                   Color , LineWidth , LineStyle )
%
% [Msg, RectHandle, CB_down ] = RECTMOVE( ... )
%
%  returns ths the ButtonDownFcn.
%
%
% The action 'newp'  specifies a PatchObject (default),
% the action 'newl'            a  LineObject for the RectAngle.
%
% To get the Position and Rotation:
%
%  [Msg, V, ... ] = RECTMOVE( 'get' , RectHandle )
%
%  with V = [ Xmin  Xmax  Ymin  Ymax  ...
%             CentreX  CentreY  ExtX  ExtY  Rotation ] .
%
%  The Rotaion is measured to the X-Axis, anticlockwise positiv.
%
%
% The UserData of the RectHandle contains:
%   { set  CB_down CB_motion CB_up  ...  }.
%
%   where set = [ ...
% 
%  1       2       3    4    5   6      7      
%  CentreX CentreY ExtX ExtY Phi Xmouse Ymouse  ...
%
%  8    9    10   11   12      13     14
%  Xmin Xmax Ymin Ymax is_move is_rot is_ext    ];
%
%
%
% Internal Actions:
%
%    'down',RectHandle                         ButtonDownEvent
%    'move',RectHandle                         Moves Rect
%    'rot' ,RectHandle                         Rotates Rect 
%    'ext' ,RectHandle, is_bord                Changes Extension
%    'up'  ,RectHandle                         ButtonUpEvent
%    'set' ,RectHandle, [ CentreX  CentreY  ExtX  ExtY  Rotation , ...
%                          Xmin  Xmax  Ymin  Ymax , ...
%                         is_move is_rot is_ext         ]
%

Nin = nargin;

if Nin < 1
   action = 'new';
end
if Nin < 2
    Hin = gca;
end

Msg        = '';
H          = [];
CB_down   = '';


% Format for Handle
clear eps
form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );


if ~any(strcmp(action,{ 'new'  'newl'   'newp' }))

  ok = ( isnumeric(Hin) & ( prod(size(Hin)) == 1 ) );
  if ok
     ok = ishandle(Hin);
     if ok
        ud = get(Hin,'userdata');
        ok = ( iscell(ud) & ( prod(size(ud)) >= 3 ) );
        if ok
           ok = ( isnumeric(ud{1}) & ...
                  ( size(ud{1},1) == 1 ) &  ( size(ud{1},2) == 14 )      & ...
                  ischar(ud{2}) & ( prod(size(ud{2})) == size(ud{2},2) ) & ...
                  ischar(ud{3}) & ( prod(size(ud{2})) == size(ud{2},2) ) & ...
                  ischar(ud{4}) & ( prod(size(ud{4})) == size(ud{4},2) )  );
        end
     end
  end

  if ok
    set0       = ud{1};
    CB_down    = ud{2};
    CB_motion  = ud{3};
    CB_up      = ud{4};
    ud0        = ud(5:length(ud));
  else
    Msg = ' RECTMOVE: Invalid Handle.';
    return
  end

end


switch action

 %************************************************************************
 case { 'new'  'newl'   'newp' }

   xl = get(Hin,'xlim'); yl = get(Hin,'ylim');
 
   set(Hin,'xlim',xl,'ylim',yl)

   if Nin < 3,  sets = [] ;     end 
   if isempty(sets)  
                pos = [ mean(xl) mean(yl) 0.5*diff(xl) 0.5*diff(yl)  0 ];
               sets = [  pos  xl yl 1  2  3 ] ;     end 
   if Nin < 4, CB_motion = '';                      end
   if Nin < 5,     CB_up = '';                      end
   if Nin < 6,     color = [ 250   19   64 ] / 255; end
   if Nin < 7, linewidth =  2 ;                     end
   if Nin < 8, linestyle = '-';                     end

   sets = sets(:)';

   set1 = [ NaN NaN 0 0 0 xl yl  ];

   sets = [ sets   set1( length(sets)+1 : 9 ) ];

   if length(sets) == 9
    sets = [ sets [ 1  2  3 ] ];
   end

    mre = [ any(sets(10:length(sets))==1) , ...
            any(sets(10:length(sets))==2) , ...
            any(sets(10:length(sets))==3)       ];     % [ Move Rot Ext ]


   % Place 6 .. 7  are for old CurrentPoint of Mouse !!!!!!!!!!!!!!
   sets = [ sets(1:5) sets(1:2) sets(6:9) mre ];

% sets = [ ...
%  1       2       3    4    5   6      7      8    9    10   11   12      13     14
%  CentreX CentreY ExtX ExtY Phi Xmouse Ymouse Xmin Xmax Ymin Ymax is_move is_rot is_ext

if ~strcmp( action , 'newl' )
       
   H =  patch( ...
     'xdata'     , nan*ones(5,1) , ...
     'ydata'     , nan*ones(5,1) , ...
     'facecolor' ,'none', ...
     'edgecolor' , color , ...
     'erasemode' , 'xor' , ...
     'linestyle' , linestyle , ...
     'linewidth' , linewidth , ...
     'clipping'  , 'off', ...
     'parent'       , Hin ); 

else

    H =  line( ...
     'xdata'     , nan*ones(5,1) , ...
     'ydata'     , nan*ones(5,1) , ...
     'color'     , color         , ...
     'marker'    , 'none'        , ...
     'markerfacecolor' , 'none' , ...
     'markeredgecolor' , 'none' , ...
     'erasemode' , 'xor' , ...
     'linestyle' , linestyle , ...
     'linewidth' , linewidth , ...
     'clipping'  , 'off', ...
     'parent'       , Hin ); 

end


   CB_down = sprintf( [ 'rectmove(''down'','  form ');'] , H );

   set(H,'userdata'      , { sets CB_down CB_motion CB_up } , ...
         'buttondownfcn' , CB_down)
 
   rectmove('set',H,sets);

 %************************************************************************
 case 'set'

  sets = sets(:)';

  sets(size(sets,2)+1:size(set0,2)) = set0(size(sets,2)+1:size(set0,2));
 
  sets(5) =  sets(5) - 180 * ceil( (sets(5)-90) / 180 ) ;

  lims = [ sets(8:9)'  sets(10:11)' ];  % [ Xlims  Ylims ]

  phi = -sets(5)*pi/180;   

 
  T = [  cos(phi)    sin(phi)  ; ...
        -sin(phi)    cos(phi)        ];


     % Rect 
       
      xr = sets(3) * [ -1  1  1 -1 -1 ] * 0.5;
      yr = sets(4) * [ -1 -1  1  1 -1 ] * 0.5;
 
      xy = [ T * [ xr ; yr ] ] + sets(1:2)'*ones(1,5);
 

   if all(( ( ( sets( 8) <= xy(1,:) ) + ...
              ( sets( 9) >= xy(1,:) ) + ...
              ( sets(10) <= xy(2,:) ) + ...
              ( sets(11) >= xy(2,:) )       ) == 4 ))

     
     set(Hin,'xdata',xy(1,:)', ...
             'ydata',xy(2,:)', ...
             'userdata',[ {sets CB_down CB_motion CB_up}  ud0 ] );

   end

 %************************************************************************
 case 'get'


% set0 = [ ...
%  1       2       3    4    5   6      7      8    9    10   11   12      13     14
%  CentreX CentreY ExtX ExtY Phi Xmouse Ymouse Xmin Xmax Ymin Ymax is_move is_rot is_ext 


  set0(5) =  set0(5) - 180 * ceil( (set0(5)-90) / 180 ) ;

  phi = -set0(5)*pi/180;  

 
  T = [  cos(phi)   sin(phi)  ; ...
        -sin(phi)   cos(phi)        ];


     % Rect 
       
      xr = set0(3) * [ -1  1  1 -1 ] * 0.5;
      yr = set0(4) * [ -1 -1  1  1 ] * 0.5;
 
      xy = [ T * [ xr ; yr ] ] + set0(1:2)'*ones(1,4);
  

      H = [ min(xy(1,:))  max(xy(1,:))  ...
            min(xy(2,:))  max(xy(2,:))  set0(1:5) ];


%************************************************************************
otherwise


 axe = get(Hin,'parent');
 fig = get(axe,'parent');

        cp = get(axe,'currentpoint');  
 selection = get(fig,'selectiontype');

  cp = cp(1,1:2);

% set0 = [ ...

%  1       2       3    4    5   6      7      8    9    10   11   12      13     14
%  CentreX CentreY ExtX ExtY Phi Xmouse Ymouse Xmin Xmax Ymin Ymax is_move is_rot is_ext 

 phi = -set0(5) * pi/180;   
 
      % Normale in RectExtension
       T = [  cos(phi)    sin(phi)  ; ...
             -sin(phi)    cos(phi)        ];  % [ nn1 nn2 ]

  
 switch action

  %******************************************************
  case 'down'

   is_bord = [ 0  0 ];
 
   move_act = '';

   % alt & edge == rot
   % alt & line == move
   % normal     == ext


   % Search if clicked on Edge

     % Length of one Pixel in Units 'Data'
     upp = ppunit(axe);

     % Rotation of Vector ( Centre --> CurrentPoint ) to phi=0
     % InversRotationMatrix = Transp.RotationMatrix
       d0 = T' * [cp(1:2)-set0(1:2)]';

      is_bord = ( abs( set0(3:4)/2 - abs(d0') ) <= 5*upp ) .* ...
                  sign(d0') ;
     
 
  switch upper(get(Hin,'type'))

    case 'LINE'

     if strcmp(selection,'alt')  &  [ set0(12) | set0(13) ]

       move_act = 'move';

       if all(is_bord) & set0(13)
        move_act = 'rot';  
       end

     elseif strcmp(selection,'normal') & set0(14) 

       move_act = 'ext'; 

     end

    case 'PATCH'
 
     if strcmp(selection,'normal')  &  ( set0(12) | set0(14) )

       move_act = 'move';

       if any(is_bord) & set0(14)
        move_act = 'ext';  
       end

     elseif strcmp(selection,'alt') & set0(13) 
 
       move_act = 'rot'; 

     end

   end


  % Store old WindowButtonMotion/UpFcn
  %  as  String for a Vector of AsciiCodes

  mot_fcn0 = get(fig,'windowbuttonmotionfcn');
   up_fcn0 = get(fig,'windowbuttonupfcn'    );

  mot_fcn = [ '['  sprintf('%4.0f',abs(mot_fcn0)) ']' ];
   up_fcn = [ '['  sprintf('%4.0f',abs( up_fcn0)) ']' ];


  set( fig , 'windowbuttonmotionfcn'    , ...
    [ sprintf(['rectmove(''' move_act ''','  form  ',[%1.0f %1.0f]);' ] , ...
            [Hin,is_bord])  CB_motion ] , ...
              'windowbuttonupfcn'       , ...
       sprintf(['rectmove(''up'','  form  ',' ...
                mot_fcn ',' up_fcn ');'       ] , Hin ) )


 %******************************************************
 case  'move'


  set0(1:2) = set0(1:2) + ( cp(1:2) - set0(6:7) );


 %******************************************************
 case  'ext'

     % sets == is_bord % !!!  [ <1,0,-1>  <1,0,-1> ] 

     % Rotation of Vector ( Centre --> CurrentPoint ) to phi=0
     % InversRotationMatrix = Transp.RotationMatrix
       d0 = T' * [cp(1:2)-set0(1:2)]';
   
     % opposite Sides
     opp = [  set0(1:2)-sets.*set0(3:4)/2 ; ...
              set0(1:2)+d0'     ];

     set1 = [ mean(opp)  abs(diff(opp)) ] .* [ (sets~=0) (sets~=0) ] + ...
              set0(1:4) .* [ (sets==0) (sets==0) ];
       
     set0(1:4) = [ set0(1:2) + ( T * (set1(1:2)-set0(1:2))' )' set1([3 4]) ];


 %******************************************************
 case  'rot'
    
      phi =  180/pi * ( phi - ( ...
               atan2(cp(2)-set0(2),(cp(1)-set0(1))) - ...
               atan2(set0(7)-set0(2),(set0(6)-set0(1)))       ) );

  set0(5) = -(phi-180*ceil(phi/180));

 
 %******************************************************
 case 'up'

   % Set old WindowButtonMotion/UpFcn back

   set(fig,'windowbuttonmotionfcn' , char(sets) ),  
   set(fig,'windowbuttonupfcn'     , char(CB_motion)   ),

   eval( CB_up , 'disp('' RECTMOVE: Error in ButtonUpFunction.'')' )
   


 end


  set0(6:7) = cp(1:2);

  rectmove('set',Hin,set0);
     
end 


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel  and AxesPixelSize
%
%  [ UnitPerPixel, PixelSize ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
 

fig = get(axe,'parent');

 
fig_uni = get(fig,'units');
set(fig,'units','pixels')
figpos = get(fig,'position');
set(fig,'units',fig_uni);

axe_uni = get(axe,'units');
set(axe,'units','normalized')
axepos = get(axe,'position');
set(axe,'units',axe_uni);

dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));


axesize = axepos.*figpos([3 4 3 4]);

mode   = get(axe,'dataaspectratiomode');


if strcmp(mode,'manual')
   aspect = get(axe,'dataaspectratio');

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * dx/dy * aspect(2)/aspect(1);
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * dy/dx * aspect(1)/aspect(2);

   ii = find( sum( ( pos <= [1;1]*axesize([3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy ] ./ axesize([3 4]) ; % UnitPerPixel
