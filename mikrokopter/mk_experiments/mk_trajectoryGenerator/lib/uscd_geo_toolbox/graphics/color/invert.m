function msg = invert(h,varargin);

% INVERT  Inverts and Rotates the Color of Object
%
% Msg = INVERT( Handle , f , 'r' , 'u' , NoTag )
%
%  f = rotate + i*invert
%
%  rotate = [ -1 .. 1 ]
%  invert =  +1          invert Color after Rotate
%            -1          invert Color before Rotate
% 
%  'r'       Recurse trough Children
%  'u'       Undo    the procedure (reverse)
%
%  NoTag     CellStringArray which Tags, which Objects are not to invert.
%             The Wildcard '*' can be used at Begin and/or End.
%
% The Color of Light-Objects will not inverted !!!
%

msg = '';

%----------------------------------------------------------
% Check Input

if nargin < 1
   msg = 'Input Handle is missing.';
   return 
end

ok = ( strcmp(class(h),'double')  &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   msg = 'First Input must be a single GraficHandle.';
   return
end

%----------------------------------------------------------
% Get Options

recurse = 0;
invert  = 0;
rotate  = 0;
undo    = 0; 

notag   = cell(0,1);

Vin = varargin;
Nin = length(Vin);

for ii = 1 : Nin

    recurse = ( recurse | isequal(Vin{ii},'r') );
    undo    = ( undo    | isequal(Vin{ii},'u') );

    if strcmp( class(Vin{ii}) , 'double' )  &  ...
       ( prod(size(Vin{ii})) == 1 )

       rotate = real(Vin{ii});
       invert = imag(Vin{ii});

    end

    if iscellstr(Vin{ii})
       notag = cat(1,notag,Vin{ii}(:));
    end

end

%----------------------------------------------------------

if ( mod(rotate,1) == 0 )  &  ( invert == 0 )
  % Nothing to do
   return
end


%----------------------------------------------------------
% Run

f = 1 - 2 * undo;

invert = f * invert;
rotate = f * rotate;

shh = get( 0 , 'showhiddenhandles' );
      set( 0 , 'showhiddenhandles' , 'on' );


setobj(h,invert,rotate,notag,recurse);


set( 0 , 'showhiddenhandles' , shh );
    
    

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setobj(h,invert,rotate,notag,recurse)


[p,t] = getprop(h,notag);

if isempty(p)
   return
end

invert = invert * ~strcmp(t,'light');

if strcmp(t,'axes')
% Get LabelColors

   tl = cat( 1 , get(h,'xlabel') , ...
                 get(h,'ylabel') , ...
                 get(h,'zlabel') , ...
                 get(h,'title')        );

   cl = get( tl , 'color' );

end


for pp = p(:)'

    v = get(h,pp{1});

    d = 2 + strcmp( pp{1} , 'cdata' );  % Dimension to Control

    if strcmp( pp{1} , 'pointershapecdata' )  &  invert

       v = 3 - v ;

    elseif  ( ( size(v,d) == 3 ) & ~isempty(v)  &  ...
          any(strcmp(class(v),{ 'double' 'uint8' }))       );

        is8 = strcmp(class(v),'uint8');

        if is8
           v = double(v) / 255;
        end

        v = modify(v,invert,rotate);

        v = min(max(v,0),1);

        if is8
           v = uint8( v * 255 );
        end

    end

    set( h , pp{1} , v );


end


if strcmp(t,'axes') 
% Set LabelColors back !!!

  for ii = 1 : size(tl,1)

     set(tl(ii),'color',cl{ii});

     if ~recurse
        setobj(h,invert,rotate,notag,0);
     end
  end

end


if ~recurse
    return
end

%----------------------------------------------------
% Recurse

ch = get(h,'children');

if isempty(ch)
   return
end

for h = ch(:)'
    setobj(h,invert,rotate,notag,recurse);
end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function v = modify(v,invert,rotate)

% Modify Colors


mode = { 'rotate'  'invert' };

mode = mode( [ 1  2 ] + [ 1  -1 ] * ( invert < 0 ) );  

for mm = mode

  switch mm{1}

   case 'invert'

     if ~( invert == 0 )
         v = 1 - v;
     end

   case 'rotate'

     if ~( rotate == 0 )

       v = rgb2hsv(v);

       if size(v,3) == 3

         v(:,:,1) = v(:,:,1) + rotate;
         v(:,:,1) = v(:,:,1) - floor(v(:,:,1));
 
       else

         v(:,1) = v(:,1) + rotate;
         v(:,1) = v(:,1) - floor(v(:,1));

       end

       v = hsv2rgb(v);
    
     end
 
  end

end

%  c = [  ( ( 1 - br ) * c + br ) ; ...
%         c*br                          ];


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [p,t] = getprop(h,notag)


p = cell(0,1);

t = get(h,'type');

tag = get(h,'tag');

nt = size(tag,2);


for tt = notag(:)'

    if isempty(tt{1}) & isempty(tag)
       return
    end

    if ~isempty(tt{1}) & ~isempty(tag)

       n2 = size(tt{1},2);

       wc = ( double(tt{1}) == double('*') );

       if wc(1) & wc(n2)

          ret = ~isempty( findstr(tag,tt{1}(2:n2-1)) );

       elseif wc(1)

          n   = min(nt,n2-1);
          ret = strcmp(tt{1}(2:n2),tag(nt-n+1:nt));

       elseif wc(n2)

          n   = min(nt,n2-1);
          ret = strcmp(tt{1}(1:n2-1),tag(1:n));
          
       else

          ret = strcmp(tt{1},tag);

       end

       if ret   
          return
       end

    end

end

switch t

 case 'figure'      

         p =  { 'color'       
                'colormap'
                'pointershapecdata'  };

 case 'axes'

         p =  { 'color'
                'colororder'
                'xcolor'
                'ycolor'
                'zcolor'          };

 case 'uicontrol'
 

         p =  { 'cdata'
                'foregroundcolor'
                'backgroundcolor' };

 case 'image'

         p =  { 'cdata' };


 case 'light'

         p =  { 'color' };

 case 'line'

         p =  { 'color'
                'markerfacecolor'
                'markeredgecolor'  };

 case 'patch'

         p =  { 'cdata'                % 3. Dimension == 3
                'facevertexcdata'      % 2. Dimension == 3
                'facecolor'
                'edgecolor'  
                'markerfacecolor'
                'markeredgecolor'  };

 case 'rectangle'

         p =  { 'facecolor'
                'edgecolor'        };

 case 'surface'

         p =  { 'cdata'                % 3. Dimension == 3
                'facecolor'
                'edgecolor'  
                'markerfacecolor'
                'markeredgecolor'  };

 case 'text'

         p =  { 'color' };  


end
