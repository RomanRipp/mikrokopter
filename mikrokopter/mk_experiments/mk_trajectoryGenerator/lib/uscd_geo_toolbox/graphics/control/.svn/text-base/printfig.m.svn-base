function  [Msg,file] = printfig(fig,varargin);

% PRINTFIG  Print Figure
%
% [ Msg , file ] = PRINTFIG( FigureHandle , Device , Option );
%
%  The Input Option must be follow the Input Device!
%
%  empty Device ==> opens PRINTDLG(fig)
%  
%-----------------------------------------------------------------------
%  Devices:   
%
%  Send to Windows-ClipBoard  (Windows only)
%
%    meta     % Send figure to clipboard in Metafile format
%    bitmap   % Send figure to clipboard in bitmap format
%
%
%  Print to PostScript-File
%
%    ps       % Level 2 PostScript for color printers
%    eps      % Encapsulated Level 2 Color PostScript
%
%   Option:  cmyk  % Use CMYK colors instead of RGB
% 
%  Save Image-Filr, use GETFRAME, IMWRITE
%
%    bmp      % Windows Bitmap (BMP)
%    jpg      % Joint Photographic Experts Group (JPEG)
%    tif      % Tagged Image File Format (TIFF)
%    png      % Portable Network Graphics (PNG)
%    pcx      % Windows Paintbrush (PCX)
%
%  Options for jpg :   Q     % Quality, a number between 0 and 100  
%
%-----------------------------------------------------------------------
% Additional Inputs
%  
%  for using PRINTDLG (empty Device), MetaFile or PostScript are:
%
%  ... , AxesHandle , NewFontSize ...
%
% Scales all FontSize of Axes- or TextObjects in Figure to new NewFontSize,
%   refering to FontSize of AxesHandle. 
% This works NOT for Objects with FontUnits 'normalized'!
% The Units for NewFontSize are 'points'.
% 
% The Input NewFontSize must follow the Input AxesHandle!
%
%  using PRINTDLG (empty Device)
%  
%    PRINTFIG( FigureHandle , AxesHandle , NewFontSize )
%    PRINTFIG( FigureHandle , AxesHandle )
%
%  PostScript ( Device =  'ps'  |  'eps' )
%
%    PRINTFIG( FigureHandle , Device , Option , AxeHandle , NewFontSize );
%    PRINTFIG( FigureHandle , Device , Option , AxeHandle );
%
%    PRINTFIG( FigureHandle , AxesHandle , NewFontSize , Device , Option ); 
%
%  MetaFile
%
%    PRINTFIG( FigureHandle , AxesHandle , NewFontSize , Device ); 
%
%
% default: NewFontSize = 11  (points)
%


Msg  = '';
file = '';

Msg0 = 'PRINTFIG: ';

 nl = char(10);

VarIn = varargin;
VarIn = VarIn(:);


%---------------------------------------------------
% allowed Devices

dv0 = { 'meta'  'bitmap'  'ps'   'eps'   ...
        'bmp'   'jpg'     'tif'  'png'  'pcx'  };

is_win = strcmp( upper(computer) , 'PCWIN' );

dv0 = dv0( 1+2*(~is_win) : end );

%---------------------------------------------------
% Get Parameter from VarIn

if nargin < 1
  fig = get(0,'currentfigure');
  if isempty(fig)
     return
  end
end

 dev = '';
 opt = [];

 axe = [];
 fs  = 11;

MsgV = '';

if ~isempty(VarIn)

  Nin = size(VarIn,1);
  ok  = zeros(Nin,1);

  for ii = 1 : Nin

    if ~ok(ii)
 
      val = VarIn{ii};

       jj = ii + ( 0 : ( ii < Nin ) );

      %----------------------
      % Check for AxesHandle
  
      if ( isnumeric(val)  & (  prod(size(val)) == 1 ) );
         if ishandle(val);
            ok(jj) = strcmp(get(val,'type'),'axes');
            if ok(ii)
               axe = val;
               if ii < Nin
                  fs = VarIn{ii+1};
               end
            end
         end
      end

      %----------------------
      % Check for Device

      if ~ok(ii)
          ok(jj) = ( ischar(val)  &  ( prod(size(val)) == size(val,2) ) );
          if ok(ii)
             dev = val;
             if ii < Nin
                opt = VarIn{ii+1};
             end
          end
      end

    end
    % ~ok(ii)

  end
  % ii  

  if any(~ok)
     MsgV = 'Can''t identify all Inputs.';
  end

end
% ~isempty(VarIn)

%***************************************************
% Check Parameter

%---------------------------------------------------
% FigureHandle

ok = ( isnumeric(fig)  & (  prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp(get(fig,'type'),'figure');
   end
end

if ~ok
  Msg = 'FigureHandle required.';
end


%---------------------------------------------------
% Device

ok = ( any( strcmp( dev , dv0 ) ) | isempty(dev) );
if ~ok
    Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Device must be any of: '  strhcat(dv0,', ')  '.'  ];
end


%---------------------------------------------------
% Option

if ok & ~isempty(opt)

   switch dev

     %----------------------------------------------
     case { 'ps'  'eps' }

       if ~isequal(opt,'cmyk')
         Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
                 'Option for Postcript must be ''cmyk''.' ];    
       end

     %----------------------------------------------
     case 'jpg'

        ok = ( isnumeric(opt)  & (  prod(size(opt)) == 1 ) );

        if ok
           ok = all( ( 0 <= opt )  &  ( opt <= 100 ) );
        end

        if ~ok
          Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
                  'Option for JPEG must be a Number between 0 and 100.' ];
        end

   end
end

%---------------------------------------------------
% FontSize

if ~isempty(axe)
  ok = ( isnumeric(fs)  &  ~isempty(fs)  &  ...
          ( prod(size(fs)) == 1 ) );
  if ok
     ok = ( isfinite(fs)  & ( fs > 0 ) );
     if ~ok
        Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Value for NewFontSize must be a single Number larger ZERO.' ]; 
     end  
  end
end


%---------------------------------------------------

if ~isempty(Msg) | ~isempty(MsgV)
    Msg = [ Msg0  Msg  ...
            nl( 1 : (end*( ~isempty(Msg) & ~isempty(MsgV) )) ) ...
                  MsgV  ];
    return
end
 

%***************************************************
% Store ResizeFcn

  ResizeFcn = get( fig , 'ResizeFcn' );

              set(  fig , 'ResizeFcn' , '' );


%***************************************************
% PRINDLG

if isempty(dev)

   if ~isempty(axe)
   % Scale FontSize

      hf = repfont(fig,axe,fs);  % [ Handle OriginalFontSize ]

   end

   try
     f = printdlg(fig);
     if ~strcmp( upper(computer) , 'PCWIN' )
        waitfor(f);
     end
   catch
     Msg = lasterr;
   end

   if ~isempty(Msg)
       Msg = [ Msg0  'Error call PRINTDLG.'  nl Msg ];
   end

   set( fig , 'ResizeFcn' , ResizeFcn );
 
   if ~isempty(axe)
   % Reset FontSize

      repfont(fig,hf);

   end

   return

end

%***************************************************
% Set to Windows-ClipBoard

if any( strcmp( dev , { 'meta'  'bitmap' } ) )

   if ~isempty(axe)  &  strcmp(dev,'meta')
   % Scale FontSize

      hf = repfont(fig,axe,fs);  % [ Handle OriginalFontSize ]

   end

   dev = cat( 2 , '-d' , dev );
   f   = cat( 2 , '-f' , epsstr(fig) );

   try
      print( f , dev );
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)
       Msg = [ Msg0  'Error call PRINT( '  dev ' ).'  nl Msg ];
   end

   set( fig , 'ResizeFcn' , ResizeFcn );
 
   if ~isempty(axe)  &  strcmp(dev,'meta')
   % Reset FontSize

      repfont(fig,hf);

   end

   return

end


%***************************************************
% Get FileName

 ext = cat( 2 , '*.' , dev );

 name = get(fig,'name');
 if isempty(name)
    name = epsstr(fig);
 end

 name = cat( 2 , 'Save Figure ' , name , ' to ... ' );

 [file,pfad]= uiputfile( ext , name );

 if isequal(pfad,0) | isempty(file)

    file = '';

    set( fig , 'ResizeFcn' , ResizeFcn );

    return

 end

 file = cat( 2 , pfad , file );

 if ( exist(file,'dir') == 7 )

   Msg = 'Selected Filename is a Directory.';

   file = '';

   set( fig , 'ResizeFcn' , ResizeFcn );

   return
 
 end

%***************************************************
% PostScript

if any( strcmp( dev ,  { 'ps'  'eps' } ) )


   if ~isempty(axe)
   % Scale FontSize

      hf = repfont(fig,axe,fs);  % [ Handle OriginalFontSize ]

   end


   Option = cell(1,0);

   if ~isempty(opt)
      if isequal(opt,'cmyk')
         Option = { '-cmyk' };
      end
   end

   dev = cat( 2 , '-d' , dev , 'c2'  );  % Level2 Color
   f   = cat( 2 , '-f' , epsstr(fig) );

   PrintIn = cat( 2 , { f } , { dev } , Option , { file } );

   try
      print( PrintIn{:} );
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)
       Msg = [ Msg0  'Error call PRINT( '  ...
                      strhcat(PrintIn,' , ')  ' ).'  nl Msg ];
       file = '';
   end

   set( fig , 'ResizeFcn' , ResizeFcn );


   if ~isempty(axe)
   % Reset FontSize

      repfont(fig,hf);

   end


   return


end

%***************************************************
% Image

%---------------------------------------------------
% JPEG Option

  Option = cell(1,0);
  if strcmp( dev , 'jpg' )  &  ~isempty(opt) 
     Option = { 'Quality'  opt };
  end

  Option = cat( 2 , { file } , { dev } , Option );

%---------------------------------------------------
% Capture Figure

  figure(fig);

  ws = get(fig,'windowstyle');
       set(fig,'windowstyle','modal');

  drawnow

  c = getframe(fig);

       set(fig,'windowstyle',ws);

  if isempty( c.colormap )
     c = cat( 2 , { c.cdata } , Option );
  else
     c = cat( 2 , { c.cdata } , { c.colormap } , Option );
  end

   try
      imwrite( c{:} );
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)
       Msg = [ Msg0  'Error call IMWRITE( C , ''' file  ''' , ' ...
                      ''''  dev  ''' ).'  nl Msg ];
       file = '';
   end

   set( fig , 'ResizeFcn' , ResizeFcn );


%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hf = repfont(fig,axe,f1)

% REPFONT  Scale all Axes- and TextFontSizes of Figure
%


%------------------------------------------
% Reset

if nargin < 3
  
   hf = axe;

   if ~isempty(hf)
 
     for ii = 1 : size(hf,1)

       set( hf(ii,1) , 'fontsize' , hf(ii,2) );

     end

   end

   return

end

%------------------------------------------
% Get Axes- and TextHandles

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

hf = cat( 1 , findobj(fig,'type','axes') , ...
              findobj(fig,'type','text')       );

      set(0,'showhiddenhandles',shh);

if isempty(hf)
   
   hf = zeros(0,2);
 
   return

end

%------------------------------------------
% ScaleValue

fu = get(axe,'fontunits');
     set(axe,'fontunits','points');

sc = f1 / get(axe,'fontsize');

     set(axe,'fontunits',fu);

%------------------------------------------
% Scale

n = size(hf,1);

hf = hf(:,[1 1]);  % [ Handle  FontSize ]
              
ok = zeros(n,1);   % Ok for NOT 'normalized'

leg = 'legend';

for ii = 1 : n

   fu = get(hf(ii,1),'fontunits');
   
   is_leg = ~isempty( findstr( lower(get(hf(ii,1),'tag')) , leg ) );
               
   ok(ii) = ~( strcmp( fu , 'normalized' )  |  is_leg );

   if ok(ii)

     % Original Value
      hf(ii,2) = get(hf(ii,1),'fontsize'); % Original Value

      set(hf(ii,1),'fontunits','points');

      set(hf(ii,1),'fontsize',get(hf(ii,1),'fontsize')*sc);

      set(hf(ii,1),'fontunits',fu);

   end  

end

hf = hf(find(ok),:);

%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ h , form ] = epsstr( h );

% EPSSTR  Transform Number exact into String,
%          
% using Matlab's floating point relative accuracy.
%
%  Form = EPSSTR;   
%    returns Format for using with SPRINTF
%
%  [ String , Form ] = EPSSTR( Number ); 
%    returns exact String for Number
%
%  Form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 )
%


form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

if nargin < 1

  h = form;

  return

end


if ~isnumeric(h)  |  ( prod(size(h)) > 1 )
 error('Handle must be a single Numeric.');
end


  h = sprintf(form,h);
