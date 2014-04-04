function [Msg , H ] =  linemove(action,Hin,in1,in2,...
                                CB_motion,CB_up, ...
                                varargin)

% LINEMOVE  draws and move a Rectangle
%
% [Msg, RectHandle] =  LINEMOVE( Action , Handle , Data , Set )
%
% Action      Handle      Data          Set
%
% 'new'       AxeHandle  [XData YData]  [XParameter YParameter]
%
%
% For action 'new'  more Options are:
% 
%   linemove( ... , ButtonMotionFcn , ButtonUpFcn , ...
%                   LineProperty , PropertyValue , ... )
%
% The UserData of the RectHandle contains:
%   { set  data  CB_down CB_motion CB_up  ...  }.
%
%   where set = [ XParameter YParameter  ;
%                 XDown      YDown       ;
%                 CurrentPointX CurrentPointY   ]  % [ 3 by 2 ]
% 
%         data = [ XOriginal ... ; 
%                  YOriginal ... ; 
%                  XDataLast ... ; 
%                  YDataLast ...   ];              % 4 Row's 
%
%
% Internal Actions:
%
%    'down',RectHandle                         ButtonDownEvent
%    'move',RectHandle                         MotionFcn 
%    'up'  ,RectHandle                         UpFcn 
%

Nin = nargin;

if Nin < 1
   action = 'new';
end
if Nin < 2
    Hin = gca;
end

Msg       = '';
H         = [];


Msg0 = 'LINEMOVE: ';

nl = char(10);

% Format for Handle
clear eps
form = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );


if ~( ischar(action)  &  ~isempty(action)  & ...
      ( prod(size(action)) == size(action,2) ) )

  Msg = [ Msg0  'Action must be a String.' ];

  return

end


action = upper(action);

if ~strcmp(action,'NEW')
   
   ok = ishandle(Hin);
   if ok
      ok = ( strcmp( get( Hin , 'type' ) , 'line' )  &  ...
             strcmp( get( Hin , 'tag'  ) , 'LineMoveObject' ) );
   end
   
   if ~ok
     Msg = [ Msg0  'Invalid Handle.' ];
     return
   end

   ud = get(Hin,'userdata');

end


switch action

 %****************************************************
 case 'NEW'


    if Nin < 3, data = NaN*ones(2,1); 
    else,       data = in1;
    end
    if Nin < 4, par = NaN;
    else,       par = in2;
    end 
    
   if Nin < 5, CB_motion = '';                      end
   if Nin < 6,     CB_up = '';                      end
               

    H =  line( ...
     'xdata'     , nan , ...
     'ydata'     , nan , ...
     'color'     , [ 250   19   64 ] / 255 , ...
     'marker'    , 'none'        , ...
     'markerfacecolor' , 'none' , ...
     'markeredgecolor' , 'none' , ...
     'erasemode' , 'xor' , ...
     'linestyle' , '-'   , ...
     'linewidth' , 1     , ...
     'clipping'  , 'off' , ...
     'parent'    , Hin   , ...
     varargin{:}         , ...
      'tag' , 'LineMoveObject'  );


   CB_down = sprintf( [ 'linemove(''down'','  form ');'] , H );

   set(H,'userdata'      , { NaN*ones(3,2) NaN*ones(4,1) ...
                             CB_down CB_motion CB_up } , ...
         'buttondownfcn' , CB_down)
      
     
      Msg1 = linemove('par',H,par);
      Msg2 = linemove('data',H,data);
      
      if ~isempty(Msg1) | ~isempty(Msg2)
         
         delete(H)
         H = [];
         
         Msg = [ Msg1 nl(1:(end*(~isempty(Msg1)))) Msg2 ];
         
      end
     
 %****************************************************
 case 'PAR'
    
    par = in1;
    
    ok = isnumeric(par);
    if ok
      ok =  ~( ( size(par,1) ~= 1 )   |  ( size(par,2) ~= 2 )  | ...
               any( real(par) < 0 )   | ...
               any( isinf(par)    )       );
    end
          
    if ~ok
      Msg = [ Msg0 'Parameter must be [ 1 by 2 ] Vector,' ...
                   ' containing positive, finite Elements or NaN ' ];
      return
    end
    
    par1 = real(par);
    par2 = imag(par);
    
    for ii = 1 : 2
       if isnan(par1(ii))
          par2(ii) = NaN;
       end
          
       if par2(ii) < 0
          if ~isnan(par1(3-ii))
             par2(ii) = par1(3-ii);
          else
             par2(ii) = 0;
          end
       end
    end
    
          
       
    ud{1}(1,:) = par1 + i*par2; 
    
    set(Hin,'userdata',ud);
    
    
 %****************************************************
 case 'DATA'
    
    data = in1;
    
    ok = isnumeric(data);
    if ok
      ok =  ( ( ndims(data)==2 ) & any( size(data) == 2 ) & ...
                 ~any(isinf(data)) );
    end
          
    if ~ok
      Msg = [ Msg0 'Data must be Matrice with 2 Columns or Row''s,' ...
                   ' containing finite Elements or NaN ' ];
      return
    end
    
    if ( size(data,1) ~= 2 )  &  ( size(data,2) == 2 )
       data = data';
    end
    
    ud{2} = data([1 2 1 2],:); 
    
    set(Hin,'xdata',data(1,:), ...
            'ydata',data(2,:), ...
            'userdata',ud);
         
         
     
 %****************************************************
 case 'SETDEV'
    
    m   = in1;   % Target
    dev = in2;   % Deviation
    
    par1 = real(ud{1}(1,:));
    par2 = imag(ud{1}(1,:));
    
       ok = isnumeric(dev);
       if ok
          ok =  ~( ( size(dev,1) ~= 1 )  |  ( size(dev,2) ~= 2 )  | ...
                   any(isinf(dev)) );
       end
          
       if ~ok
         Msg = [ Msg0  'Deviation must be [ 1 by 2 ] Vector,' ...
                       ' containing finite Elements or NaN ' ];
         return
       end
       
    if ( sum(isnan(par1)) == 2 )  |  ( sum( dev == 0 ) == 2 )
       return
    end
    
    if ischar(m)
       if strcmp(in1,'absolute')
          xy0 = ud{2}([3 4],:); % [ XDataLast YDataLast ];
          m   = ud{1}( 2   ,:); % [ XDown     YDown     ];
       else
          xy0 = [ get(Hin,'xdata') ; get(Hin,'ydata') ];
          m   = ud{1}(3,:);    % [ cpx_old  cpy_old ]
       end
    else
       
       ok = isnumeric(m);
       if ok
         ok =  ~( ( size(m,1) ~= 1 ) | ( size(m,2) ~= 2 ) | ...
                    any(isinf(m)) );
       end
          
       if ~ok
         Msg = [ 'Target must be [ 1 by 2 ] Vector,' ...
                 ' containing finite Elements or NaN ' ];
         return
       end
      
          xy0 = [ get(Hin,'xdata')  get(Hin,'ydata') ];
          jj = find(isnan(m));
          if length(jj) == 2
             return
          elseif length(jj) == 1
             m(jj) = interp1(xy0(3-jj),xy0(jj),m(3-jj));
          end
          ud{1}([2 3],:) = ones(2,1)*m([1 2]);
          ud{2}([3 4],:)   = xy0;
          set(Hin,'userdata',ud);
       
    end
     
    xy1 = xy0;
    
    N1 = ones(1,size(xy0,2));
    
    pot1                 = -0.5*( (xy0-(m')*N1) ./ ((par1')*N1) ).^2;
              is_zero1   = ( pot1 < (-1e2) ) + isnan(pot1);
    pot1(find(is_zero1)) = 0;
    
    fak1 = exp(pot1) - is_zero1 ; 

    
    par2                 =  ( (par2') * N1 );
              is_inf     =  ( par2 == 0 );
    par2(find(is_inf))   =  NaN;
    
    pot2                 = -0.5*( (xy0([2 1],:)-[m']*N1) ./ par2 ).^2;
              is_zero2   = ( pot2 < (-1e2) ) + isnan(pot2);
    pot2(find(is_zero2)) = 0;     
    
    fak2 = exp(pot2) - is_zero2 + is_inf;
    
    xy1([2 1],:) = xy0([2 1],:) + ((dev([2 1])')*N1).*fak1.*fak2;
          
     
    set(Hin,'xdata',xy1(1,:), ...
            'ydata',xy1(2,:)       );
    
         
 %****************************************************
 case 'DOWN'
    
     axe = get(Hin,'parent');
     fig = get(axe,'parent');
     
     selection = get(fig,'selectiontype');
     cp        = get(axe,'currentpoint');
     
     if strcmp(selection,'alt')
        linemove('undo',Hin);
        return
     end  
     if ~strcmp(get(fig,'selectiontype'),'normal')
       return
     end
   
    
   % Store old WindowButtonMotion/UpFcn
   %  as  String for a Vector of AsciiCodes

    mot_fcn0 = get(fig,'windowbuttonmotionfcn');
     up_fcn0 = get(fig,'windowbuttonupfcn'    );

    mot_fcn = [ '['  sprintf('%4.0f',abs(mot_fcn0)) ']' ];
     up_fcn = [ '['  sprintf('%4.0f',abs( up_fcn0)) ']' ];


   set( fig , ...
        'windowbuttonmotionfcn' , ...
          [ sprintf(['linemove(''move'','  form  ');' ] , Hin ) ...
                      ud{4} ]  , ...
        'windowbuttonupfcn' , ...
          sprintf(['linemove(''up'','  form  ',' ...
                     mot_fcn ',' up_fcn ');'       ] , Hin ) )

   ud{1}([2 3],:) = ones(2,1)*cp(1,[1 2]);  % [ xdown ydown ; cpx cpy ]
   ud{2}([3 4],:) = [ get(Hin,'xdata') ; get(Hin,'ydata') ]; % [XDataLast YDataLast ]
               
   set( Hin , 'erasemode' , 'xor' , ...
              'userdata'  , ud         );

   
 %****************************************************
 case  'MOVE'
   
    cp = get( get(Hin,'parent') , 'currentpoint' );
    cp = cp(1,[1 2]);
    
    modes = { 'absolute'  'iterative' };
    
    Mode = 1;
    
    linemove('setdev',Hin,modes{Mode},cp-ud{1}(Mode+1,:));
 
    ud{1}(3,:) = cp;
    
    set(Hin,'userdata',ud);
    
 
 %****************************************************
 case 'UP'
    
     axe = get(Hin,'parent');
     fig = get(axe,'parent');

   % Set old WindowButtonMotion/UpFcn back
     
   set(fig,'windowbuttonmotionfcn' , char(in1) ),  
   set(fig,'windowbuttonupfcn'     , char(in2)   ),
  
   eval( ud{5} , 'disp('' LINEMOVE: Error in ButtonUpFunction.'')' )

   
 %****************************************************
 case 'UNDO'

    set(Hin,'xdata',ud{2}(3,:), ...
            'ydata',ud{2}(4,:)      );

   eval( ud{5} , 'disp('' LINEMOVE: Error in ButtonUpFunction.'')' )
        
 %****************************************************
 case 'RESET'
    
    ud{2}([3 4],:) = ud{2}([1 2],:);
    
    set(Hin,'xdata',ud{2}(1,:), ...
            'ydata',ud{2}(2,:), ...
            'userdata',ud );
    
 end
