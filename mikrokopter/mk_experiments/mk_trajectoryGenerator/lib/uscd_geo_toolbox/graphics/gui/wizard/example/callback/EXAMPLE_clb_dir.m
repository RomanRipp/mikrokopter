function dir_sel(par,group,action,clb,varargin)

% DIR_SEL Callback for Directory-Section (Ordner)
%

ud = get( par , 'userdata' );

ch = getfield( ud.Children , group );


switch upper(action)

%*****************************************************
case 'LIST'

 switch upper(group)

  %--------------------------------------
  case 'SELECT'
  % PopUp

    val = get(ch.List,'value');
    str = get(ch.List,'string');
    lud = get(ch.List,'userdata');

    switch lud(val) 

     %-----------------------------
     case  1 
     % Go UP

       OldPfad = str{val};
       NewPfad = str{val};

       
     %-----------------------------
     case  0
     % Same Directory, go one DOWN

       val_new = val + 1 * ( val < size(str,1) );
     
       val_new = val_new - 1 * ( lud(val_new) == -2 ) ...
                             * (     val_new  >   1 );
       
       OldPfad = str{ val     };
       NewPfad = str{ val_new };

     
     %-----------------------------
     case -1
     % Go DOWN
       
       NewPfad = str{ val };
       OldPfad = str{ val - 1 * ( val > 1 ) };

     %-----------------------------
     case  -2
     % Change Disk, Windows only

       OldPfad = '';
       NewPfad = str{val};

       
     %-----------------------------
     otherwise

        return
      
    end

  %--------------------------------------
  case 'DIR'
  % DirectoryListBox

     nr  = get(ch.List,'value');

     lud = get( ch.List , 'userdata' ); % { IsDir Name }

     name = lud{nr,2};

 
     %------------------------------------------------
     % File | Dir ==> Get Full Name

     if ~( lud{nr,1} == -1 )

       val = get(ud.Children.Select.List,'value');
       str = get(ud.Children.Select.List,'string');

       pfad = str{val};

       fs = filesep;
       fs = fs(1:(end*(~strcmp(pfad(end),fs))));

       file = cat( 2 , pfad , fs , name );

     end

     %------------------------------------------------
     % Dir ==> Set ToolTip of UP-Button

     if ( lud{nr,1} == 1 )

       set( ud.Children.Select.Up , 'enable' , 'on' , ...
                             'tooltipstring' , file         );

     else

       set( ud.Children.Select.Up , 'enable' , 'off' , ...
                             'tooltipstring' , ''          );

     end

     %-----------------------------------------------
     % Check SelectionType
     
     fig = get(ch.List,'parent');

     if ~strcmp( get(fig,'selectiontype') , 'open' )
       
         return

     end

     %------------------------------------------------
     % File ==> display in FileInfo
     if ( lud{nr,1} == 0 )

       if isempty(name)
          return
       end

       EXAMPLE_clb_dir(par,'File','NewFile',clb,file);
 
       return

     end


     %------------------------------------------------
     % Get Parent Directory

     str = get(ud.Children.Select.List,'string');

     if ( lud{nr,1} == -1 )
     % Go DOWN
     
       lud = get(ud.Children.Select.List,'userdata');

       i0 = find( lud == 0 );

       if isempty(i0)
          return
       end

       OldPfad = str{ i0 };
       NewPfad = str{ i0 + 1 * ( i0 < size(str,1) ) };

       set(ud.Children.Select.List,'value',i0);

     else
     % Go UP
     
       val = get(ud.Children.Select.List,'value');

       NewPfad = str{val};

       fs = filesep;
       fs = fs(1:(end*(~strcmp(NewPfad(end),fs))));

       OldPfad = cat( 2 , NewPfad , fs , name );

     end

  %--------------------------------------
  otherwise

     return

  end


   EXAMPLE_clb_dir(par,'Select','NEWDIR',clb,NewPfad,OldPfad);


%*****************************************************
case 'UP'

  % Go Up to selected Directory of DirectoryContentsList
  
  nr  = get(ud.Children.Dir.List,'value');
  lud = get(ud.Children.Dir.List,'userdata'); % { IsDir Name }

  if ~( lud{nr,1} == 1 )
     return
  end

  name = lud{nr,2};
  
  val = get(ch.List,'value');
  str = get(ch.List,'string');

  NewPfad = str{val};

       fs = filesep;
       fs = fs(1:(end*(~strcmp(NewPfad(end),fs))));

  OldPfad = cat( 2 , NewPfad , fs , name );
  
  EXAMPLE_clb_dir(par,group,'NEWDIR',clb,NewPfad,OldPfad);


%*****************************************************
case 'DOWN'

  val = get(ch.List,'value');
  str = get(ch.List,'string');
  lud = get(ch.List,'userdata');

   i0 = find( lud == 0 );

   if isempty(i0)
      return
   end

   i1 = i0 + 1 * ( i0 < size(str,1) );
     
   i1 = i1 - 1 * ( lud(i1) == -2 ) * ( i1 > 1 );
        
   OldPfad = str{ i0 };
   NewPfad = str{ i1 };

   set(ch.List,'value',i0);

  EXAMPLE_clb_dir(par,group,'NEWDIR',clb,NewPfad,OldPfad);


%*****************************************************
case 'NEW'

  set(ch.List,'visible','off')

  set(ch.Edit,'visible','on' , ...
              'string' , ''        );

    hh = [ ch.Up ch.Down ch.New ch.Delete ];
  sets = cell(1,size(hh,2));
  for ii = 1 : size(hh,2)
      sets{ii} = get(hh(ii),'enable');
  end

  set( ch.New , 'userdata' , { hh sets } );


 %---------------------------------------
 % ParentDirectory

  str = get(ch.List,'string');
  val = get(ch.List,'value');

 set( ch.Edit , 'tooltipstring' , ...
       cat( 2 , get(ch.Edit,'userdata') , str{val} , '  ' ) );
 
%*****************************************************
case 'EDIT'

  nl = char(10);

  name = get( ch.Edit , 'string' );

  sets = get( ch.New , 'userdata' );

  for ii = 1 : size(sets{1},2)
     set( sets{1}(ii) , 'enable' , sets{2}{ii} );
  end


  set(ch.Edit,'visible','off' , ...
              'string' , ''        );

  set(ch.List,'visible','on')

 %---------------------------------------
 % ParentDirectory

  str = get(ch.List,'string');
  val = get(ch.List,'value');

  pfad = str{val};

 %---------------------------------------
 % Check Input

  [Msg,name,dirname] = chk_name(name,pfad,'','','dir');

  if ~isempty(dirname) & ~isempty(name)
     if ( exist( dirname , 'dir' ) == 7 )
        [pfad,n] = fileparts(dirname);
        EXAMPLE_clb_dir(par,group,'NEWDIR',clb,pfad,dirname);
        return
     end
  end

  if isempty(dirname)
     if ~isempty(Msg)
        main_dlg(par,Msg,'Invalid Input','warn');
     end  
    return
  end

 %---------------------------------------
 % Create Directory

  msg = [ 'Create Directory   ' dirname ];

  EXAMPLE_clb_msg( ud.Root , msg , 'new' );

 
  command = cat(2,'mkdir ',dirname);
 
  if isunix
    [status,Msg] = unix(command);
  elseif strcmp(computer,'PCWIN') 
    [status,Msg] =  dos(command);
  else
    command = [ '! ' command ];
    eval(command,'status = -1; Msg = lasterr;')
  end

      
  if isempty(dir(dirname))

    EXAMPLE_clb_msg( ud.Root , 'Error' , 'append' );

    Msg = ['Can''t create Directory.' ...
             nl nl dirname nl  ...
             nl(1:(end*(~isempty(Msg)))) Msg ];

     warndlg(Msg,'Warning','warn'); 

     return

  end



  EXAMPLE_clb_dir(par,group,'NEWDIR',clb,pfad,dirname);


%*****************************************************
case 'DELETE'

    val = get(ch.List,'value');
    str = get(ch.List,'string');
    lud = get(ch.List,'userdata');

    if ~( lud(val) == 1 )
       return
    end

    msg = [ 'Remove Directory   ' str{val} ];

    EXAMPLE_clb_msg( ud.Root , msg , 'new' );

    [Msg,ok] = remdir( str{val} , 0 );

    if ok == -1 

        EXAMPLE_clb_msg( ud.Root , 'cancel' , 'append' );

        return

    end

    if ~isempty(Msg)

        EXAMPLE_clb_msg( ud.Root , 'Error' , 'append' );

        warndlg(Msg,'Warning','warn');

    end

    if exist( str{val} , 'dir' ) == 7
       return
    end

    if any( lud == 0 )
       EXAMPLE_clb_dir( par , group , 'DOWN' , clb );
       return
    end


    ok = ( size(lud,1) > val );
    if ok
       ok = ( lud(val+1) == 1 );
    end

    val = val - 1 + 2*ok;
    
    val = max(val,1);
    val = min(val,size(str,1));    

    OldPfad = str{val};

    i0 = find( lud == 0 );
    if isempty(i0)
       NewPfad = str{val};
    else
       NewPfad = str{i0}; 
    end


   EXAMPLE_clb_dir(par,'Select','NEWDIR',clb,NewPfad,OldPfad);

    
%*****************************************************
case 'NEWDIR'

   is_win = strcmp( upper(computer) , 'PCWIN' );
   
   if isempty(varargin)
   % call by MENU_SORT, update

      str = get(ch.List,'string');
      val = get(ch.List,'value');

      NewPfad = str{val};
      OldPfad = str{val};

   else

      NewPfad = varargin{1};
      OldPfad = varargin{2};

   end

   if is_win
   % Problems with UpperCase Names

     NewPfad = lower( NewPfad );
     OldPfad = lower( OldPfad );
     
   end
   
       fig = get( ch.List , 'parent' );

   pointer = get( fig , 'pointer' );
             set( fig , 'pointer' , 'watch' );


   ok = ( strcmp(NewPfad,OldPfad)  &  ~isempty(varargin) );

   if ok

      str = get(ch.List,'string');
      val = get(ch.List,'value');
      lud = get(ch.List,'userdata');

      ok  = ~isempty(str);
      if ok
         ok = ~isempty(str{val});
      end

   end

   if ~ok
   % New History

     msg = [ 'Get History of   ' NewPfad ];

     EXAMPLE_clb_msg( ud.Root , msg , 'new' );

     c = dirhist( NewPfad , {ud.SortMode} );

     str = c(:,1);
     lud = cat(1,c{:,2});

     if is_win
        str = lower(str);
     end
    
     if isempty(OldPfad)

       val = find( strcmp( str , NewPfad ) );

       if isempty(val)
          val = 1;
       end
      
     else
     % Search for best match  of NewPfad in OldPfad

        n  = size(str,1);

       cmp = double( char( cat( 1 , str , {OldPfad} ) ) );

       cmp = double( cmp(1:n,:) == ones(n,1)*cmp(n+1,:) );

       cmp = sum( cumprod(cmp,2) , 2 );

       [cmp,val] = max( cmp );
                 
     end

     set( ch.List , 'value'    , val , ...
                    'string'   , str , ...
                    'userdata' , lud       );


   end

   tip = '';

   on = ( lud(val) == 1 );
 
   sets = { 'off'  'on' };

   sets = sets{ 1 + on };

   i0 = find( lud == 0 );

   if ~isempty(i0) & strcmp( sets , 'on' )
      tip = str{i0(1)};
   end

   set( ch.Down , 'enable' , sets , ...
           'tooltipstring' , tip        );
 
   EXAMPLE_clb_dir( par , 'Dir' , 'Contents' , clb );

   %------------------------------------------

   NewDir = str{val};
   OrgDir = ud.Directory;

   sets = { 'off'  'on' };
   sets = sets{ 1 + ( lud(val) == 1 ) };

   set( ch.Delete , 'enable' , sets );
   
   %------------------------------------------
   % Search for best match  of NewDir in OrgDir
   %  only if DOWN

   if size(NewDir,2) < size(OrgDir,2)

     fs = filesep;
     fs = fs(1:(end*(~strcmp(NewDir(end),fs))));

     lud = get( ud.Children.Dir.List , 'userdata' );

     is_dir = find( ( cat(1,lud{:,1}) == 1 ) );

     lud = lud( is_dir , : );
     
     lud(:,1) = { cat( 2 , NewDir , fs ) };

     lud = permute( lud , [ 2  1 ] );

     [msg,lud] = char2cell( strhcat( lud , '' , 2 ) );

     if isempty(msg)

     % Search for best match  of OrgDir in [ NewDir  DirNames ]

        n  = size(lud,1);

       cmp = double( char( cat( 1 , lud , {OrgDir} ) ) );

       cmp = double( cmp(1:n,:) == ones(n,1)*cmp(n+1,:) );

       cmp = sum( cumprod(cmp,2) , 2 );

       [cmp,val] = max( cmp );

       if cmp > ( size(NewDir,2) + size(fs,2) )

         set( ud.Children.Dir.List , 'value' , is_dir(val) );

         get( ud.Children.Dir.List , 'listboxtop' );  % a Trick

         % EXAMPLE_clb_dir( par , 'DIR' , 'LIST' , clb );

          set( ud.Children.Select.Up , 'enable' , 'on' , ...
                                'tooltipstring' , lud{val}        );

       end

     end

   end

   ud.Directory = NewDir;

   set( par , 'userdata' , ud );

   %------------------------------------------

   if clb

     % Call Functions  % new_dir( ud.Parent , NewDir )

   end

   set( fig , 'pointer' , pointer );

%*****************************************************
case 'CONTENTS'

% Get new DirectoryContents

  val = get(ud.Children.Select.List,'value');
  str = get(ud.Children.Select.List,'string');
 

  if isempty(str)

     l   = { '' };
     lud = { NaN  '' };

     pfad = '';

  else

     pfad = str{val};

     msg = [ 'Get Contents of   ' pfad ];

     EXAMPLE_clb_msg( ud.Root , msg , 'new' );

     [c,l] = dircont( pfad , { ud.SortMode } ); 

     if isempty(l)

        EXAMPLE_clb_msg( ud.Root , 'can''t read Directory' , 'append' );

        l   = { ' ..' };
        lud = { -1  '..' };

     else

        lud = c(:,[1 4]); % { IsDir  Name }

     end

  end

  set( ch.List , 'value'    , 1   , ...
                 'string'   , l   , ...
             'listboxtop'   , 1   , ...
                 'userdata' , lud       ); 

  set( ch.Text , 'string' , ...
        cat( 2 , get(ch.Text,'userdata') , pfad ) );


  val  = get(ch.List,'value');

  sets = { 'off'  'on' };

  sets = sets{ 1 + ( lud{val,1} == 1 ) };

    fs = filesep;
    fs = fs(1:(end*(~strcmp(pfad(end),fs))));

   tip = cat( 2 , pfad , fs , lud{val,2} );

   tip = tip( 1 : (end*strcmp(sets,'on')) );


  set( ud.Children.Select.Up , 'enable' , sets , ...
                        'tooltipstring' , tip        );


%*****************************************************
case 'NEWFILE'


  file = varargin{1};

   nf  = size(file,2);
   n   = min( nf , 4 );

   ext = lower( file( nf-n+1 : nf ) );

  ImgExt = { '.jpg'  '.jpeg'  '.tif'  '.tiff' ...
             '.bmp'  '.png'   '.pcx'  '.xpm'  ...
             '.ps'   '.psc'   '.eps'  '.epsc'      };

  par = get( par , 'parent' );

  pointer = get( par , 'pointer' );
            set( par , 'pointer' , 'watch' );

  switch lower(ext) 

    %--------------------------------------------
    case ImgExt

      msg = [ 'Read Image from File   ' file ];

      EXAMPLE_clb_msg( ud.Root , msg , 'new' );

      fig = '';
      txt = '';

      try
         [Msg,fig] = showimg(file);
      catch
          Msg = lasterr;
      end

      if ~isempty(Msg)

        Msg = [ 'Error call SHOWIMG.'  char([10 10])  Msg ];

        EXAMPLE_clb_msg( ud.Root , 'can''t read Image' , 'append' );

      end


      str = { ' *** Image *** '
              ''
              cat(2,' ',txt)       };


    %--------------------------------------------
    case '.mat'
 
      [pfad,name,ext] = fileparts(file);


      [ Msg , txt ] = whosfile(file);

      if ~isempty(Msg)

         EXAMPLE_clb_msg( ud.Root , 'Error' , 'append' );
   
      else

         pre = cat( 2 , 'Contents of MAT-File: ' , name , ext );

         str = cat( 1 , {''} , {''} , {pre} , {''} , {''} , txt , {''}  );

         str{2} = char( double('=') * ones( 1 , size(char(str),2) ) );
         str{4} = char( double('=') * ones( 1 , size(char(str),2) ) );

      end

    %--------------------------------------------
    otherwise

      msg = [ 'Read File   ' file ];

      EXAMPLE_clb_msg( ud.Root , msg , 'new' );

      [Msg,str] = readfile(file);

      if ~isempty(Msg)

        EXAMPLE_clb_msg( ud.Root , 'can''t read File' , 'append' );

      end

  end
  % switch

  set( par , 'pointer' , pointer );

  if ~isempty(Msg)

       str = cat( 1 , {' *** Can''t read File! ***'} , {''} );

       [Msg,txt] = char2cell(Msg);

       if isempty(Msg)

          str = cat( 1 , str , txt );

       end

  end

  d = dir(file);

  txt = cat( 2 , get(ch.Text,'userdata') , file , ...
                 '   ' , d.date );

  set( ch.Text , 'string' , txt );


  set( ch.List , 'value'  , 1   , ...
                 'string' , str , ...
             'listboxtop' , 1       );

end
