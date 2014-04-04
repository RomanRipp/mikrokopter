function [msg,varargout] = cnfcdf(varargin)

% CNF_CDF  GUI to create a Configuration for gridded NetCDF-Database
%           to using LOAD_CDF  
%
% UserInterface to select Dimension- and VariableNames for
%  4 Dimensions ( X  Y  Z  Time ).
% Additional Variables could be masked by a Keyword.
% The Configuration could be saved in an AsciiFile..
% 
% InputOptions are:
% 
% [MSG,Fig]   = CNF_CDF( CDF_FILENAME )
%
% [MSG,Fig,V] = CNF_CDF( CONFIG_FILENAME )
%
%
% [MSG,V,DIM,VAR,ATT,TXT] = CNF_CDF( CONFIG_FILENAME , 'check' )
%
% checks, if the Dimension- and VariableNames in the ConfigFile
%  are valid for the NetCDF-File, which FileName stored in
%   the ConfigFile under: # FileName
%
% If the NetCDF-File is ok, but the Data in the ConfigFile are
%   not valid, V contains the NetCDF-File-Name.
%
% If  all Data are valid, V is an CellArray with 
%   4 Columns :{ FileName DimName DimLength VarName }
%   4 Rows:  Dimension X; Y; Z; T .
%  The FileName is the First CELL-Element ( V{1,1} ).
%  The Second Cell-Element ( V{2,1} ) contains an
%    [ n by 2 ] - CellArray, contains such other Keywords in the
%   first Column which Values are VariableNames (second Column)
%     in the NetCDF-File.
%   
% For the following OutPut's see LOOK_CDF.
%
% To REDEFINE the NetCDF-FileName of the ConfigFile use:
%
%    CNF_CDF( CONFIG_FILENAME , 'check' , NetCDF_FileName )
%
% 
%-------------------------------------------------------------
%
% See also: LOOK_CDF, LOAD_CDF, LOAD_ASC, WRT_ASC
%
%-------------------------------------------------------------
%
% Requested M-Files (should be located in PRIVATE-Directory):
%
% look_cdf
% load_asc
% wrt_asc
% mkgui
% stdgui
% msg_logo
% msg_list
% txt_logo
% tag_frame
% tab_list
% readfile
% char2cell
% get_ref
% findtag
% rmblank
%

Nout = nargout;

varargout = cell(1,Nout-1);

nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

tag = upper(fcn);

msg0 = sprintf('%s: ',tag);

[msg,fig,action,file,hndl,varin] = checkin(tag,varargin{:});

if ~isempty(msg)
    msg = [ msg0  msg ];
    if Nout == 0
       error(msg);
    end
    return
end

dlg = {};
out = {};

%*************************************************************

msg0 = sprintf('%s(%s): ',tag,action);

switch upper(action)

%*************************************************************
case 'CHECK'

   [msg,inf,dim,var,att,txt] = chkcnf(file,Nout,varin);

   out = { inf dim var att txt };

%*************************************************************
case 'NEW'


     %----------------------------------------------
     % Check for File

     if ~isempty(file)

         [msg,dim,var,att,txt] = look_cdf(file);
         is_cdf = isempty(msg);

         is_cnf = ~is_cdf;
         if is_cnf
            [m,inf,dim,var,att,txt] = chkcnf(file);
            is_cnf = isempty(m);
            if ~is_cnf
                is_cdf = chkstr(inf{1,1},1);
                if is_cdf
                   file = inf{1,1};
                elseif chkcstr(inf{1,1},1);
                   msg = m;
                end
            end
         end

        if ~( is_cdf | is_cnf )
            msg = sprintf('%sFile must be a valid Config- or NetCDF-File.\n%s', ...
                           msg0,msg);
            if Nout == 0
               error(msg);
            end
            return
        end

     end

     %----------------------------------------------

     [msg,fig] = mkgui(gui_config(fcn),'NetCDF-Config');

     if ~isempty(msg)
         try, delete(fig), end
         msg = sprintf('%sError call MKGUI.\n%s',msg0,msg);
         if Nout == 0
            error(msg);
         end
         return
     end

     ud = get(fig,'userdata');

     %------------------------------------------------------
     % Modify DataVariableEdit

     hc = get(ud.Children.Frame,'userdata');
     hc = hc.Children.Control;
     
     tud = get(hc.Data.List,'userdata');

     set(tud.edit{1}(3),'style','popupmenu','callback','');

     set([hc.CNF.Edit hc.CDF.Edit],'horizontalalignment','left');

     set([hc.L.Dim hc.L.Var],'horizontalalignment','center');

     %------------------------------------------------------
     % InfoFigures

     f1 = newfig('NetCDF-Info',sprintf('%s_CDF',tag));
     f2 = newfig('Config-Info',sprintf('%s_CNF',tag));

     hm = ud.Children.Menu;

     set([hm.Load hm.Save],'enable','off');

     hm = get(hm.Info,'userdata');
     hm = hm.Children;

     set(hm.CDF,'userdata',f1);
     set(hm.CNF,'userdata',f2);

     %------------------------------------------------------

     form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

     f = sprintf(form,fig);

     %------------------------------------------------------
     % New callback of TAB_LIST-Controls, 
     %  to use "private/tab_list.m" in Callback !!!

     ht = findtag(fig,'TAB*','type','uicontrol');

     for hh = ht(:)'
         cb = get(hh,'callback'); 
         if ~isempty(cb)
             cb = sprintf('%s(%s,''EVAL_CB'',''%s',fcn,f,strrep(cb,'(',''','));
             set(hh,'callback',cb);
         end
     end

     %------------------------------------------------------

     CloseFcn  = sprintf('%s(%s,''Close'');',fcn,f);
     DeleteFcn = sprintf('%s(%s,''Delete'');',fcn,f);

     ud.Figures = [ fig  f1  f2 ];

     set( fig , 'CloseRequestFcn' , CloseFcn  , ...
                      'DeleteFcn' , DeleteFcn , ...
                            'Tag' , tag       , ...
                        'visible' , 'on'      , ...
                       'userdata' , ud             );


     d = struct( 'String' , { {} } , ...
                 'Rang'   , { [] } , ...
                 'Length' , { [] }       );
     v = struct( 'String' , { {} } , ...
                 'Dim'    , { [] } , ...
                 'Rang'   , { [] }       );

     setappdata(fig,'Config',struct('Dim',{d},'Var',{v}));
     setappdata(fig,'Modify',0);
     setappdata(fig,'File','');

     %-------------------------------------------------------

     if isempty(file)

        out = { fig };

     else

        if ~isempty(msg)
            msgbox(msg,'Invalid ConfigFile','warn');
        end

        if is_cnf
           msg = setcnf(fig,file,inf,dim,var,txt);
           fcn = 'setcnf';
           out = { fig inf dim var att txt };
        else
           msg = setcdf(fig,file,dim,var,txt);
           fcn = 'setcdf';
           out = { fig dim var att txt };
        end

        if ~isempty(msg)
            msg = sprintf('Error call %s(''%s'').\n%s',upper(fcn),file,msg);
            dlg = { msg 'Error' 'warn' };
            msg = [ msg0 msg ];
        else
           sets = { 'on'  'off' };
           set(ud.Children.Menu.Save,'enable',sets{1+isempty(getappdata(fig,'File'))});
        end

     end


%*************************************************************
case 'CNF'

   ud = get(fig,'userdata');

   hc = get(ud.Children.Frame,'userdata');
   hc = hc.Children.Control.CNF;

   old = get(hc.Edit,'userdata');

   sets = { 'on'  'off' };

   switch upper(varin{1})

     %----------------------------------------------
     case 'EDIT'

       file = rmblank(char(get(hc.Edit,'string')));

       if isempty(file)
          set( hc.Edit , 'string' , old );
          return
       end

       str = chkfile(file,'cnf');

       if ~isempty(str)
           file = str;
       else
           [p,f,e] = fileparts(file);

           ok = ~isempty(p);
           if ok
              ok = ~( exist(p,'dir') == 7 );
           end
           if ~ok
               p = pwd;
           end

           if isempty(e)
              e = get(hc.Browse,'userdata');
           end

           file = cat( 2 , fullfile(p,f) , e );
       end

     %----------------------------------------------
     case 'BROWSE'

        ext = cat(2,'*',get(hc.Browse,'userdata'));

       [str,pfad] = uigetfile(ext,'Select a ConfigFile');

       if isequal(str,0)
          return
       end

       file = [ pfad  str ];

   end

   set( hc.Edit , 'string'   , file , ...
                  'userdata' , file        );

   set(ud.Children.Menu.Load,'enable',sets{1+isempty(str)});

   set(ud.Children.Menu.Save,'enable',sets{1+isempty(getappdata(fig,'File'))});

   setappdata(fig,'NewConfigFile',1);

   out = {file};

%*************************************************************
case 'CDF'

  ud = get(fig,'userdata');

  hc = get(ud.Children.Frame,'userdata');
  hc = hc.Children.Control.CDF;

  old = get(hc.Edit,'userdata');

   sets = { 'on'  'off' };

  switch upper(varin{1})

     %----------------------------------------------
     case 'EDIT'

       file = rmblank(char(get(hc.Edit,'string')));

       if isempty(file)
          set(hc.Edit,'string',old)
          return
       end

       str = chkfile(file,'cdf');

       if isempty(str)
          msg = sprintf('File "%s" not found.',file);
       elseif isequal(str,old)
          set(hc.Edit,'string',old)
          return
       else
          file = str;
       end

       if isempty(msg)
          if ~change(fig)
              set(hc.Edit,'string',old)
              return
          else      
              msg = setcdf(fig,file);
          end
       end

     %----------------------------------------------
     case 'BROWSE'

        if ~change(fig)
            return
        end

       ext = cat(2,'*',get(hc.Browse,'userdata'));

       [file,pfad] = uigetfile(ext,'Select a NetCDF-File');

       if isequal(file,0)
          return
       end

       file = [ pfad  file ];

       if isequal(file,old)
          return
       end

       msg = setcdf(fig,file);

   end

   set(ud.Children.Menu.Save,'enable',sets{1+isempty(getappdata(fig,'File'))});

   if ~isempty(msg)
       set(hc.Edit,'string',old)
   else
       out = {file};
   end

   dlg = msg;

%*************************************************************
case { 'X' 'Y' 'Z' 'T' }

     switch upper(varin{1})

       case 'DIM'

            setvar(fig,action);

       case 'VAR'

            editvar(fig,action);
     
     end

%*************************************************************
case 'DATA'

  switch upper(varin{1})

     %-----------------------------------------------
     case 'TABDELETE'

        setappdata(fig,'Modify',1),

     %-----------------------------------------------
     case 'TABSETQUEST'

        str = varin{3};  % { nrow  Key   ''  }
       
          nr = str{1};
         str = rmblank(str{2},2);

        if isempty(str)
           msgbox([ 'Wrong Input.' ] , 'Error' ,'warn' );
           msg = 0;
           return
        end

        ud = get(fig,'userdata');
        hc = get(ud.Children.Frame,'userdata');
        hc = hc.Children.Control;
     
        tud = get(hc.Data.List,'userdata');

        % PopupMenu
        val = get(tud.edit{1}(3),'value');
        var = get(tud.edit{1}(3),'string');
        usd = get(tud.edit{1}(3),'userdata');

        if ischar(var)
           var = cellstr(var);
        end

        var = var{val};

        nr = min(nr,get(hc.Data.List,'max')+1);

        % New Output for TabList
          msg = { nr  [' ' str]  [' ' var] };
          setappdata(fig,'Modify',1);   % IS_CHANGE

          return

     %-----------------------------------------------
     case { 'TABSET'  'TABINS' }

        ud = get(fig,'userdata');
        hc = get(ud.Children.Frame,'userdata');
        hc = hc.Children.Control;
     
        tud = get(hc.Data.List,'userdata');

        % PopupMenu
        val = get(tud.edit{1}(3),'value');
        usd = get(tud.edit{1}(3),'userdata');
 
        dlg = tab_list(hc.Data.List,'userdata',{ varin{3}{1} usd(val) });

   end

%*************************************************************
case 'INFO'

    fig = get(hndl,'userdata');

    fig = getfield(fig.Children,varin{1});

    fig = get(fig,'userdata');

    set(fig,'visible','on');
   
    figure(fig);

%*************************************************************
case 'EVAL_CB'

    if isempty(varin)
       return
    end

    if ~chkstr(varin{1},1)
      
        dlg = 'First Option for EVAL_CallBack must be a FunctionName.';

        msg = [ msg0  dlg ];

    else

        try

            feval(varin{:});

        catch

            dlg = sprintf('Error call FEVAL(''%s'').\n%s',varin{1},lasterr);

            msg = [ msg0  dlg ];

        end

    end

%*************************************************************
case 'MENU'

 switch upper(varin{1})

  %*************************************************************
  case 'LOAD'

  ud = get(fig,'userdata');

  hc = get(ud.Children.Frame,'userdata');
  hc = hc.Children.Control.CNF;

  if ~change(fig)
      return
  end
 
  file = rmblank(get(hc.Edit,'string'),2); 

  if isempty(file)
     [msg,file] = cnf_cdf(fig,'CNF','Browse');
     if ~isempty(msg) | isempty(file)
         return
     end
  end

  msg = setcnf(fig,file);

  if isempty(msg)
     setappdata(fig,'Modify',0);
  else
     dlg = msg;
  end

  %*************************************************************
  case 'SAVE'

  ud = get(fig,'userdata');

  hc = get(ud.Children.Frame,'userdata');
  hc = hc.Children.Control;

  [fk,key,km,cm] = keydef;  % key = { #Dim  #Var }

  cnf = getappdata(fig,'Config');
  % { Dim.String Dim.Rang Dim.Length Var.String Var.Dim Var.Rang }

  %-----------------------------------------------------
  % Get Dimension and VariableSelection

  sk = size(key);

  var    = cell(sk+[0 1]);  % { #Dim  #Var  #Len }
  var(:) = {''};

  for ii = 1 : sk(1)

      h = getfield(hc,key{ii,1}(1));

      for jj = 1 : sk(2)

          hh = getfield(h,key{ii,jj}(2:end));

         val = get(hh,'Value');
         str = get(hh,'String');

         if ~strcmp(str{val},'none')
             var{ii,jj} = str{val};
             if ( jj == 1 )
                var{ii,3} = cnf.Dim.Length(val);
             end
         end

      end
         
  end

  c = var(find(~strcmp(var(:,1),'')),1);
  c = double(char(sort(c)));
  c = ( sum(diff(c,1,1)==0,2) == size(c,2) );

  if any(c)
     msg = 'Select different Dimensions.';
  end


  %-----------------------------------------------------
  % Get  Other Variables in Tabular

  tud = get(hc.Data.List,'userdata');

  okey = {};
  ovar = tud.userdata;

  no = min(size(ovar,1),get(hc.Data.List,'max'));

  if no == 0
     ovar = {};
  else
     okey = cellstr(rmblank(char(tud.string{1}),2));
     ovar = ovar(1:no);

        c = okey(find(~strcmp(okey,'')),1);
        c = double(char(sort(c)));
        c = ( sum(diff(c,1,1)==0,2) == size(c,2) );

        if any(c)
           msg = [ msg nl(1:(end*(~isempty(msg)))) ...
                   'Select different KeyWords for DataVariables.' ];
        end
  end

        v = cat(1,key(:),okey);
        v = sort(v);
        c = double(char(v));
        c = ( sum(diff(c,1,1)==0,2) == size(c,2) );
        if any(c)
           c = find(c);
          jj = cat( 1 , 1 , find( diff(c,1,1) > 1 )+1 );
          jj = c(jj);
           msg = [ msg nl(1:(end*(~isempty(msg)))) ...
                   'KeyWords are not allowed for DataVariables:' nl ...
                    strhcat(v(jj),', ') ];
        end

  if ~isempty(msg)
      msgbox(msg,'Invalid Configuration','warn');
      return
  end

  %-----------------------------------------------------

  file = rmblank(char(get(hc.CNF.Edit,'string')),2);

  if ~isempty(file)

      if ( exist(file,'file') == 2 )

          ok = questdlg('File exist!  Overwrite?' , ...
                         'Stop','Yes','No','Cancel','Yes'); 

          switch ok
                 case 'Yes'
                       is_change = 0;   
                 case 'Cancel'  
                       msg = 'cancel';
                       return;
                 case 'No'
                       file = '';
          end

      end

  end

  if isempty(file)

     ext = cat(2,'*',get(hc.CNF.Browse,'userdata'));
     [file,pfad] = uiputfile(ext,'Save Configuration as');
     if isequal(file,0)
        msg = 'cancel';
        return
     end

     file = [pfad,file];

     set(hc.CNF.Edit,'string',file,'userdata',file);
 
  end

  %-----------------------------------------------------

  key = permute(key,[2 1]);
  key = key(:);

  var(:,1) = num2cell(var(:,[1 3]),2);  % { #Dim Len }
  var(:,2) = num2cell(var(:,2),2);
  var = permute(var(:,[1 2]),[2 1]);
  var = var(:);

  %----------------------------------------------------
  % Check if NetCDF-File is in SearchPath

  ncfile = getappdata(fig,'File');
  [p,f,e] = fileparts(ncfile);
  if isequal(which([f e]),ncfile)
     ncfile = [f e];
  end

  %----------------------------------------------------

  setmsg(fig,sprintf('Write CONFIG-File "%s".',file),'new');

  msg = wrtcnf( file , fk , key , okey , ncfile , var , ovar , km , cm );


  if ~isempty(msg)
      setmsg(sprintf('error\n%s',msg),'append'); 
      return
  end 

  setappdata(fig,'Modify',0);
  setappdata(fig,'NewFileConfig',0);

  setinfo(fig,file);


  %*************************************************************
  case 'QUIT'

    cnf_cdf(fig,'Close');

end

%*************************************************************
case 'CLOSE'

   if ~change(fig)
       return
   end

   cnf_cdf(fig,'Delete');


%*************************************************************
case 'DELETE'

    ud = get(fig,'userdata');

    try
       for f = ud.Figures
           try, delete(f); end
       end
    catch
       delete(fig);
    end

%**************************************************************
otherwise

    dlg = sprintf('Invalid action "%s".',action);

    msg = [ msg0 dlg ];

end

%**************************************************************

n = min(size(out,2),Nout-1);

varargout(1:n) = out(1:n);


if ~isempty(dlg) & ~isempty(fig);

    if chkstr(dlg)
       dlg = { dlg 'Error' 'warn' };
    end
    msgbox(dlg{:});
end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,inf,dim,var,att,txt] = chkcnf(file,Nout,v);

if nargin < 2
   Nout = 10;
end

if nargin < 3
   v = {};
end


inf = {};
dim = {};
var = {};
att = {};
txt = {};

msg = '';

nl = char(10);

 %-----------------------------------------------------------

 [fk,key,km,cm] = keydef;  % key = { Dim  Var }

 nk = size(key,1);
 nb = size(key,2);

 inf      = cell(nk,4);

 inf(2,1) = { cell(0,2) };  % { Key VarName }

 inf(:,3) = { NaN };        % DimLength

 inf(:,[2 4]) = { '' };     % { Dim Var }

%-----------------------------------------------------------
% Check for FileName
 
if ~isempty(v)
    v = v{1};
    if ~chkstr(v,1)
        v = '';
    end
end

if ~isempty(v)

    ncfile = v;

    m = 'Invalid Input for NetCDF-FileName';

else

 %-----------------------------------------------------------
 % Load FileName

 [ msg , v ] = load_asc(file,fk,km,cm);

 if ~isempty(msg)
     return
 end

 if isempty(v{1,3});
    msg = 'Didn''t find';
 else
    v = v(find(~strcmp(v(:,3),'none')),:);
    if isempty(v)
       msg = 'Didn''t find valid';
    else
       ncfile = v{1,2};
       if isempty(ncfile)
          msg = 'Empty';
       end
    end
 end

 if ~isempty(msg)
     msg = sprintf('%s KeyWord "%s" in ConfigFile "%s".',msg,fk,file);
     return
 end

 m = sprintf('Invalid NetCDF-FileName in Keyword "%s".',fk);

end

inf{1,1} = {ncfile};   % CellString !!!

%------------------------------------------------------------
% Check NetCDF-File

 if     Nout >= 6
    [msg,dim,var,att,txt] = look_cdf(ncfile);
 elseif Nout == 5
    [msg,dim,var,att] = look_cdf(ncfile);
 else
    [msg,dim,var] = look_cdf(ncfile);
 end
     
 % dim = { DimName  DimLength }
 % var = { VarName VarType Ndim [dim] Nattr }
 
 if ~isempty(msg)
     msg = sprintf('%s\nError call LOOK_CDF(%s).\n%s',m,ncfile,msg);
     return 
 end 

 f = which(ncfile);
 if ~isempty(f)
     ncfile = f;
 end

 inf{1,1} = ncfile;  % String !!!

 %------------------------------------------------------------
 % Check for Dim- and Var-KeyWords

 [ msg , v ] = load_asc(file,key(:),km,cm);
 
 if ~isempty(msg)
     return
 end

 v = permute(v,[2 1]);
 v = reshape(v,size(v,1),nk,nb);
 v = permute(v,[2 1 3]);

 ok = ~( strcmp(v(:,3,:),'none') | strcmp(v(:,3,:),'') );
 ok = permute(ok,[1 3 2]);

 if ~any(ok(:,1))
     return
 end

 %----------------------------------------
 % Check for ("#DIM" NOT found) AND  ("#VAR" found)

 bad = ( ok(:,2) & ~ok(:,1) );

 if any(bad)
    bad = find(bad);
    msg = sprintf('Empty Dimensions for Variables: %s.',strhcat(key(bad,2),', '));
 end

 %----------------------------------------------------------
 % Check DimensionNames 

  ind = find(ok(:,1));

  for ii = ind(:)'

      val = v{ii,2,1};

      if ~isempty(val)
          if iscell(val)
             val = val{1};
          end
      end

      if chkstr(val,1)
         jj = strcmp(dim(:,1),val);
         if any(jj)
            inf{ii,3} = dim{find(jj),2};
            inf{ii,2} = val;
         else
            inf{ii,3} = -1;    % Not found
         end
      elseif ~isempty(val)
         inf{ii,3} = 0;
      end

  end
 
  len = cat(1,inf{:,3});

  for bad = [ 0  -1 ]
      jj = ( len == bad );
      if any(jj)
         jj = find(jj);
         str = strhcat(v(jj,1,1),', ');
         if bad == 0
            str = sprintf('DimensionNames must be Strings: %s.',str);
         else
            str = sprintf('Invalid DimensionNames for %s.',str);
         end
         msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , str );
         inf(jj,3) = { NaN };
          ok(jj,1) = 0; 
      end
  end

  %------------------------------------------------------------
  % Check VariableNames and VariableDimensions

  vok = NaN*zeros(nk,1);

  ok(:,2) = ( ok(:,1) & ok(:,2) );

  ind = find(ok(:,2));

  for ii = ind(:)'

      val = v{ii,2,2};

      if ~isempty(val)
          if iscell(val)
             val = val{1};
          end
      end

      if ~chkstr(val,1)
          vok(ii) = 0;
      else
          jj = strcmp(var(:,1),val);
          if ~any(jj)
              vok(ii) = -1;
          else
             jj = find(jj);
             if ~( var{jj,3} == 1 )
                 vok(ii) = -2;
             else
                 if ~strcmp(dim{var{jj,4}+1,1},v{ii,2,1})
                     vok(ii) = -3;
                 else
                     if strcmp(var{jj,2},'char')
                        vok(ii) = -4;
                     else
                        inf{ii,4} = val;
                     end
                 end
             end
          end
      end

  end

  for bad = [ -1 -2 -3 ]
      jj = ( vok == bad );
      if any(jj)
         jj = find(jj);
         str = strhcat(v(jj,1,2),', ');
         if     bad == -1
                str = sprintf('Invalid VariableNames for %s.',str);
         elseif bad == -2
                str = sprintf('Multiple Dimensions for Variables: %s.',str);
         elseif bad == -3
                str = sprintf('Invalid DimensionNames for Variables: %s.',str);
         elseif bad == -4
                str = sprintf('Invalid Type "char" for Variables: %s.',str);
         end
         msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , str );
      end
  end


  %-----------------------------------------
  % Search for Other VariableKeys

  [m,k] = load_asc(file,'*',km,cm);

  k = k(strcmp(k(:,3),'char'),:);
  
  if isempty(k)
     return
  end

  n = size(k,1);

  ok = zeros(n,1);

  key = key(:);

  for ii = 1 : n
      
      if ~any(strcmp(k{ii,1},key))
          if chkstr(k{ii,2},1)
             ok(ii) = any(strcmp(k{ii,2},var(:,1)));
          end
      end

  end
       
  if any(ok)
     ok = find(ok);
      k = k(ok,[1 2]);
      % Check for duplicate KeyWords, take last
      [c,si] = sort(k(:,1));
      c = double(char(c));
      c = ( sum( ( diff(c,1,1) == 0 ) , 2 ) == size(c,2) );
      if any(c)
         c = find(c);     % Leave last
         k(si(c),:) = [];
      end
      inf{2,1} = k;
  end



%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  msg = wrtcnf( file , fkey , gkey , dkey , ...
                               fvar , gvar , dvar , km , cm );

msg = '';

fid = fopen(file,'wt');

if fid == -1
   msg = sprintf('Can''t open File "%s" for writing.',file);
   return
end

%***************************************************
% Basic Comment

cl = '*********************************************';

dt = clock;
dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
dt = datestr(dt,0);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

c = { 'Configuration for gridded NetCDF-DataBase'
      sprintf('Matlab: %s by %s',dt,upper(fcn))
      'see also: LOAD_CDF'  };

fprintf(fid,'%s%s\n',cm,cl);

for cc = c(:)'
    fprintf(fid,'%s %s\n',cm,cc{1});
end

fprintf(fid,'\n');

%***************************************************
% FileName

fprintf(fid,'%s %s %s ''%s''\n\n',km{1},fkey,km{2},fvar);

%***************************************************
% DataVariables

if ~isempty(dvar)

    fprintf(fid,'%s%s\n',cm,cl);
    fprintf(fid,'%s %s\n\n',cm,'DataVariables');

    s2 = size(char(dkey),2);

    for ii = 1 : size(dkey,1)
        bl = char(32*ones(s2-size(dkey{ii},2)));
        fprintf(fid,'%s %s%s %s ''%s''\n',km{1},dkey{ii},bl,km{2},dvar{ii});
    end

    fprintf(fid,'\n');

end


%***************************************************
% Grid

fprintf(fid,'%s%s\n',cm,cl);
fprintf(fid,'%s %s\n\n',cm,'Dimensions of Grid');

frm = { '''%s'''   '''%s''  %.0f' };

nl ={ '\n'  '\n\n' };

for ii = 1 : size(gkey,1)

    ff = frm{ 1 + ( size(gvar{ii},2) == 2 ) };
    nn = nl{2-mod(ii,2)};
    
    fprintf(fid,'%s %s %s ',km{1},gkey{ii},km{2});
    fprintf(fid,ff,gvar{ii}{:});
    fprintf(fid,nn);

end

%***************************************************

fclose(fid);


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = change(fig);

 ok = ~getappdata(fig,'Modify');

 if ~ok
   
    ok = questdlg(['Unsaved Changes!  Save them now?'] , ...
                   'Stop','Yes','No','Cancel','Yes'); 

    switch ok
       case 'Yes'
            msg = cnf_cdf(fig,'Menu','Save'); 
             ok = isempty(msg);
       case 'No'  
             ok = 1;
       case 'Cancel'
             ok = 0;
    end

 end    

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = setcnf(fig,file,inf,dim,var,txt);

Nin = nargin;

msg = '';

ud = get(fig,'userdata');

hc = get(ud.Children.Frame,'userdata');
hc = hc.Children.Control;


if Nin < 3
   setmsg(fig,sprintf('Call CHKCNF(''%s'')',file));
   [msg,inf,dim,var,att,txt] = chkcnf(file);
   if ~isempty(msg)
       setmsg(fig,'error','append');
       if chkstr(inf{1,1},1)              % NetCDF-File ok
          set(hc.CDF.Edit,'string',inf);
          msg1 = setcdf(fig,inf{1,1},dim,var,txt);
          if ~isempty(msg1)
              msg = sprintf('%s\n%s',msg,msg1);
          end
       end
       return
   end 
end

msg = setcdf(fig,inf{1,1},dim,var,txt,0);
if ~isempty(msg)
    return
end

%***************************************************************

  dim_str = inf(:,2);

  jj = all( char(dim_str) == 32 , 2 );
  if any(jj)
     jj = find(jj);
     dim_str(jj,:) = { 'none' };
  end

  var_str = inf(:,4);

  jj = all( char(var_str) == 32 , 2 );
  if any(jj)
     jj = find(jj);
     var_str(jj,:) = { 'none' };
  end

  % Check, if DimensionNames and VariableName
  %  right set by SETCDF
  % And Set DimensionPopUpValue and VariablePopUp

  id = 'XYZT';

  for ii = 1 : size(id,2)

       h = getfield(hc,id(ii));

       % Value in DimensionPopUp
       jj = find( strcmp( dim_str{ii} , get(h.Dim,'string') ) );

       set(h.Dim,'value',jj)

       % Set VariablePopup

       old = get(h.Dim,'userdata'); 
             set(h.Dim,'userdata',0);  % OldValue

       setvar(fig,id(ii));

       set(h.Dim,'userdata',old)

  end
  % ii
 

  % Set if OtherVariables right in CDF-File
  %  Set Tabular

  if ~isempty(inf{2,1})

      tud = get(hc.Data.List,'userdata');

      okey = inf{2,1}(:,1);
      ovar = inf{2,1}(:,2);

      usd = get(tud.edit{1}(3),'userdata');
      str = get(tud.edit{1}(3),'string');

      for ii = 1 : length(ovar)

          % Value in VariablePopUp
          jj = strcmp(ovar{ii},usd);

          if any(jj)
             jj = find(jj);
             tab_list(hc.Data.List,'string',{ ii [' ' okey{ii}] [' ' str{jj}] });
             tab_list(hc.Data.List,'userdata',{ ii  usd(jj) });
          end

      end
      % ii

  end
  % ~isempty(inf{2,1})

  set(hc.CNF.Edit,'string',file,'userdata',file);

  msg = setinfo(fig,file);


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = setcdf(fig,file,dim,var,txt,mode);

Nin = nargin;

msg = '';

if Nin < 6
   mode = 1;
end

mode = isequal(mode,1);

ud = get(fig,'userdata');

hc = get(ud.Children.Frame,'userdata');
hc = hc.Children.Control;

old = get(hc.CDF.Edit,'userdata');

if Nin < 3
   setmsg(fig,sprintf('Call LOOK_CDF(''%s'')',file));
   [msg,dim,var,att,txt] = look_cdf(file);
   if ~isempty(msg)
       setmsg(fig,'error','append');
       set(hc.CDF.Edit,'string',old);
       return
   end
end

%-------------------------------------------------------------
%  Set OtherVariablesPopUp in Tabular

tud = get(hc.Data.List,'userdata');

  nv = size(var,1);

  [h,si] = sort( -cat(1,var{:,3}) );

  str = cell(nv,1);

  for ii = 1 : nv

      vdim = strhcat(dim(var{si(ii),4}+1),',');

      str{ii} = sprintf('%s(%s)',var{si(ii),1},vdim);

  end

  set(tud.edit{1}(3),'string', str , ...
                     'value' ,  1  , ...
                  'userdata' , var(si,1)  );

  % Check, if old Entrys in Tabular ok

  if ~isempty(tud.string{1})

      tab_list(hc.Data.List,'delete',( 1 : size(tud.string{1},1)) );

      if mode

         jj = find( ~all( tud.string{1} == 32 , 2 ) );

         okey = cellstr(tud.string{1}(jj,:));
         ovar = cellstr(tud.string{2}(jj,:));

         no = size(okey,1);
         ok = zeros(no,1);
         for ii = 1 : no
             cmp = strcmp(rmblank(ovar{ii},2),str);
             if any(cmp)
                ok(ii) = find(cmp);
             end
         end

         if ~all( ok == 0 )
             ok = find(ok);
             tab_list( hc.Data.List , 'string' , { okey(ok)  ovar(ok) } );
             tab_list( hc.Data.List , 'userdata' , { ok  var(si(ok),1) } );
         end

      end

  end


%-------------------------------------------------------------
% Get Variables with single Dimensions

  jj = ( cat(1,var{:,3}) == 1 );    % Ndim == 1
 
  if ~any(jj)
      msg = 'No Variable with single Dimension found.';
      var_str = { '' };
      var_dim = [];
  else
      jj = find(jj);
      var_str = var(jj,1);
      var_dim = cat(1,var{jj,4});
  end

  nv = size(var(jj,1),1);

%-------------------------------------------------------------
% Sort Dimensions
 
  dim_str = dim(:,1);   

  nd = size(dim,1);  

 id = 'XYZT';

 dim_mat = lower(char(dim_str));

 dim_rang = ( 1 : nd )';

 rang = { 'x' 1  ; 'lon' 1  ; ...
          'y' 2  ; 'lat' 2  ; ...
          'z' 3  ; 'dep' 3  ; ...
          't' 4  ; 'tim' 4  };    %  { Name  Rang  ID }

 nr = 4;

 for ii = 1 : size(rang,1);

     jj = strmatch(rang{ii,1},dim_mat);

     if ~isempty(jj)
         dim_rang(jj) = rang{ii,2};
     end

 end

 for ii = 1 : nd
 
     dim_rang(ii) = dim_rang(ii) + nr * (~any(strcmp(dim_str{ii},var_str)));

 end

%-------------------------------------------------------------
% Sort Variables

 var_mat = lower(char(dim_str));

 var_rang = dim_rang( var_dim+1 );

%-------------------------------------------------------------
% Set DimensionPopupMenus

 for rr = 1 : nr
  
     jj = cat( 1 , find( dim_rang == rr ) , ...
                   find( dim_rang == rr+nr )  );

     if isempty(jj)
        jj = nd+1;   % 'none'
     end

     h = getfield(hc,id(rr));

     set( h.Dim , 'string' , [ dim_str ; {'none'} ] , ...
                   'value' , jj(1) , ...
                'userdata' , 0      )

 end

 dim_len = cat(1,dim{:,2});


     d = struct( 'String' , { dim_str  } , ...
                 'Rang'   , { dim_rang } , ...
                 'Length' , { dim_len  }       );
     v = struct( 'String' , { var_str  } , ...
                 'Dim'    , { var_dim  } , ...
                 'Rang'   , { var_rang }       );

     c = struct('Dim',{d},'Var',{v});

     setappdata(fig,'Config',c);
 
%-------------------------------------------------------------
% Set VariablePopupMenus

 setvar(fig,'XYZT',c);

%-------------------------------------------------------------
% Set InfoWindow

 full = which(file);

 if ~isempty(full)
     file = full;
 end

 set(ud.Figures(2),'name',file)

 set(get(ud.Figures(2),'userdata'),'string',txt)
  
 set(hc.CDF.Edit,'string',file,'userdata',file);

 setappdata(fig,'File',file)

 setappdata(fig,'Modify',1)  % IS_CHANGE


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setvar(fig,id,cnf)

 if nargin < 2
    id = 'XYZT'
 end

 if nargin < 3
    cnf = getappdata(fig,'Config');
 end

 is_change = getappdata(fig,'Modify');

ud = get(fig,'userdata');

hc = get(ud.Children.Frame,'userdata');
hc = hc.Children.Control;

 for ii = id(:)'

     h = getfield(hc,ii);

   old = get(h.Dim,'userdata');
   val = get(h.Dim,'value');

   is_change = ( is_change  |  ~( val == old ) );
   
   set(h.Dim,'userdata',val);

   if  ~( val == old )

     dim = val-1;
    
     vv = find( cnf.Var.Dim == dim );       % var_dim == dim

     if isempty(vv)
      s_i = 1;
     else
      [hilf,s_i] = sort( cnf.Var.Rang(vv) );  % sort(var_rang)
     end

     set(h.Var,'string', [ cnf.Var.String(vv) ; {'none'} ] , ...
               'value' , s_i(1)          );
   end

 end

 setappdata(fig,'Modify',is_change);


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function editvar(fig,id)

 if nargin < 2
    return
 end

 is_change = getappdata(fig,'Modify');

ud = get(fig,'userdata');

hc = get(ud.Children.Frame,'userdata');
hc = hc.Children.Control;

h = getfield(hc,id);

   old = get(h.Var,'userdata');
   val = get(h.Var,'value');

   is_change = ( is_change  |  ( val ~= old ) );
   
   set(h.Var,'userdata',val);

 setappdata(fig,'Modify',is_change);


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = setinfo(fig,file);

setmsg(fig,sprintf('Call READFILE(''%s'')',file),'new');

[msg,txt] = readfile(file);

if ~isempty(msg)
    setmsg(fig,'error','append');
    return
end

ud = get(fig,'userdata');

 set(ud.Figures(3),'name',file)

 set(get(ud.Figures(3),'userdata'),'string',txt)

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setmsg(fig,msg,mode,varargin)

if nargin < 2
   return
end

if nargin < 3
   mode = 'new';
end

ud = get(fig,'userdata');
ud = get(ud.Children.Frame,'userdata');

h = ud.Children.Logo.Message;

msg_list(h,'Message',msg,mode,varargin{:});

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fk,key,km,cm] = keydef

fk = 'FileName';

        base = { 'Dim'
                 'Var'   };     nb = size(base,1);
         pre = 'XYZT';          np = size(pre,2);
     
      key        = cell( np , nb );

      for ii = 1:np
       for jj = 1:nb
        key{ii,jj} = cat( 2 , pre(ii) , base{jj} );
       end
      end 

      km = { '#'  ':' };  % KeyMarker
      cm =   '%';         % CommentMarker

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function file = chkfile(str,mode);

if nargin < 2
   mode = '';
end

switch mode
 case 'cdf'
      ext = { '.cdf'  '.nc'  '' };
 case 'cnf'
      ext = { '.cnf'  '' };
 otherwise
      ext = { '.cnf'  '.cdf'  '.nc'  '' };
end

   file = '';
   for e = ext
       file = cat(2,str,e{1});
       if exist(file,'file') == 2
          f = which(file);
          if ~isempty(f)
              file = f;
          end
          break
       end
   end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,fig,action,file,hndl,v] = checkin(tag,varargin)

Nin = nargin - 1;

msg    = '';
fig    = [];
action = '';
file   = '';
hndl   = [];
v      = cell(1,0);

if Nin == 0
   action = 'new';
   return
end

%**************************************************
% Check for FileName or Handle

val = varargin{1};

is_str = chkstr(val,1);

if is_str

   file = chkfile(val);

   if isempty(file)
      msg = sprintf('File "%s" not found.',val);
      return
   end

else

   ok = ( ( prod(size(val)) == 1 ) & strcmp(class(val),'double') );
   if ok
      ok = ishandle(val);
      if ok
         hndl = val;
         [m,fig] = recpar(hndl,'Root');
         ok = isempty(m);
         if ok
            fig = fig(1);
            ok = strcmp(get(fig,'tag'),tag);
         end
         if ~ok
             msg = sprintf('HandleInput must orign from %s-Figure.',tag);
             return
         end
      else
         msg = 'Numeric Input must be a Handle.';
         return
      end
   else
      msg = 'Input must be a FileName or Handle.';
      return
   end

end

%**************************************************
% Check for Action

if Nin < 2
   if ~isempty(hndl)
       msg = 'Input ACTION is missing.';
   else
       action = 'new';
   end
   return
end

val = varargin{2};

if ~chkstr(val,1)
    msg = 'Input ACTION must be a String.';
    return
end

action = val;

v = varargin(3:Nin);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = newfig(name,tag);

% NEWFIG   Create New ListBoxFigure


  %--------------------------------------
  scr_uni = get(0,'units');      set(0,'units','pixels')
  scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);
  
  ppi    = get(0,'ScreenPixelsPerInch');

  %--------------------------------------

   is_tall = -1 + ( scr_si(4) >=  480 ) + ...
                  ( scr_si(4) >=  600 ) + ...
                  ( scr_si(4) >= 1024 );

   is_win = strcmp( upper(computer) , 'PCWIN' );

   fontsize =  8 + 2 * is_tall - 2 * is_win;

 
   if is_win
      fontname = 'courier';
   else
      fontname = { 'arrial'  get(0,'fixedwidthfontname') };
      fontname = fontname{ 1 + ( scr_si(4) >= 1050 ) };
   end

  %--------------------------------------
  
  fs = floor( 25/18 * ppi/100 * fontsize );  % Points --> Pixels

  ww = 80;                   % Character
  hh = 50;
                          
  fig00 = ceil([ 0.55*ww  1.2*hh ] * fs);                         
  fig11 = floor([ 1/2  2/3 ].*scr_si(3:4));
  
  figpos = NaN*ones(1,4);
  
  figpos(3:4)= fig00 + ( fig11 - fig00 ) .* ( fig11 < fig00 );
  
  voffs = max( 60 , min( ceil(1/6*scr_si(4)) , 80 ) );
  hoffs = 20;

  figpos(1) = scr_si(3)-hoffs-figpos(3);
  figpos(2) = scr_si(4)-voffs-figpos(4);
  

 fig  = figure('position'   , figpos , ...
               'numbertitle', 'off'  , ...
               'menubar'    , 'none' , ...
               'toolbar'    , 'none' , ...
               'name'       , name   , ...
               'tag'        , tag    , ...
               'createfcn'  , ''     , ...
               'visible'    , 'off'  , ...
          'handlevisibility','callback' );


  hl = uicontrol( 'parent'    , fig        , ...
                  'style'     , 'listbox'  , ...
            'backgroundcolor' , [1 1 1]    , ...
            'foregroundcolor' , [0 0 0]    , ...
                  'units'     ,'normalized', ...
                  'position'  , [0 0 1 1] , ...
                  'min'       , 0         , ...
                  'max'       , 2         , ...
                  'fontunits' , 'points'  , ...
                  'fontsize'  , fontsize  , ...
                  'fontname'  , fontname  , ...
                  'string'    , ''        , ...
                  'tag'       , 'LIST'    , ...
        'horizontalalignment' , 'left'           );


    ParentCB = 'get(gcbo,''parent'')';
     CloseCB = sprintf('set(%s,''visible'',''off'');',ParentCB);
   
   hc = uimenu( 'parent'      , fig     , ...
                'accelerator' , 'C'     , ...
                'label'       , 'Close' , ...
                'tag'         , 'CLOSE' , ...
                'callback'    , CloseCB       );

   set( fig , 'userdata' , hl , 'CloseRequestFcn' , CloseCB );

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Config = gui_config(fcn)

Config = struct( 'Menu'    , { gui_menu(fcn)    } , ...
                 'Logo'    , { gui_logo         } , ...
                 'Control' , { gui_control(fcn) }        );


%**************************************************************
function Config = gui_menu(fcn)

info = struct( 'CDF' , { { 'NetCDF-File'  0  fcn  [] } } , ...
               'CNF' , { { 'Config-File'  0  fcn  [] } }  );

Config = struct( 'Info' , { { 'Info'  0  info  [] } } , ...
                 'Save' , { { 'Save'  0  fcn   [] } } , ...
                 'Load' , { { 'Load'  0  fcn   [] } } , ...
                'Blank' , { { '   '   0  ''    [] } } , ...
                 'Quit' , { { 'Quit'  0  fcn   [] } }       );

%**************************************************************
function Config = gui_logo(fcn)

cl0 = [ 0.0  0.0   0.0 ];  % LogoFore
cl1 = [ 0.95  1.0   1.0 ];  % LogoBack


logo_txt = { 'NetCDF' 
             'Configuration'
             'for gridded DataBase' };

Config = { logo_txt  cl0 cl1 };

%**************************************************************
function Config = gui_control(fcn)

cl = [ 0.95  1.0   1.0 ];

wwt = 10;  % FileText
wwl = 40;  % FileEdit
wwb = 8;   % Browse

wwd = 4;   % DimText
wwp = 25;  % PopUp

Config = cell(10,2);

%-------------------------------------------------------------------------------
% File

Text1 = stdgui('Type' , 'text' , 'Width' , wwt, 'String' , 'Config-File' );
Text2 = stdgui('Type' , 'text' , 'Width' , wwt, 'String' , 'NetCDF-File' );

Edit1 = stdgui( 'Type' , 'edit' , 'Width' , wwl , 'Color' , cl , 'UserData' , '' , ...
                'CBFcn' , fcn );
Edit2 = stdgui( 'Type' , 'edit' , 'Width' , wwl , 'Color' , cl , 'UserData' , '' , ...
                'CBFcn' , fcn );

Brws1 = stdgui( 'Type' , 'pushbutton' , 'Style' , 'ediT' , 'Width' , wwb , 'Color' , cl , ...
                'String' , 'Browse' , 'UserData' , '.cnf' , 'CBFcn' , fcn );
Brws2 = stdgui( 'Type' , 'pushbutton' , 'Style' , 'ediT' , 'Width' , wwb , 'Color' , cl , ...
                'String' , 'Browse' , 'UserData' , '.cdf' , 'CBFcn' , fcn );

File1 = struct( 'Text' , { { 1  Text1 } } , ...
                'Edit' , { { 2  Edit1 } } , ...
              'Browse' , { { 3  Brws1 } }       );

File2 = struct( 'Text' , { { 1  Text2 } } , ...
                'Edit' , { { 2  Edit2 } } , ...
              'Browse' , { { 3  Brws2 } }       );


Config(1:3,:) = { 'CNF'    { File1 } 
                  'CDF'    { File2 } 
                  'Sep1'   {  NaN  }        };
                 
%-------------------------------------------------------------------------------
% Select

pre = { '' 'X' 'Y' 'Z' 'Time' };
suf = { 'Txt'  'Dimension'  'Variable' };

t = stdgui('Type','text','Width',wwd);
p = stdgui('Type','popupmenu','Width',wwp,'Color',cl,'CBFcn',fcn,'UserData',0);

np = size(pre,2);
ns = size(suf,2);

ini = cell(2,ns);

for ii = 1 : np
    
    for jj = 1 : ns
        if any( [ ii  jj ] == 1 )
           c = t;
           if ii == 1
              if jj > 1
                 c.Width  = p.Width+2;
                 c.String = suf{jj};
                 c.Style  = 'Control';
              end
           else
              c.String = pre{ii};
           end
        else
           c = p;
        end

        ini(1,jj) = { suf{jj}(1:3) };
        ini(2,jj) = { { { jj  c } } };

    end

    if ii == 1
       lab = 'L';
    else
       lab = pre{ii}(1);
    end

    Config{ii+3,1} = lab;

    Config{ii+3,2} = { struct(ini{:}) };

end

Config(9,:) = { 'Sep2'   { NaN }  };

%-------------------------------------------------------------------------------
% Table
  
 opt = {  ...
          ...  
  'ntab'  ,    2 , ...
  'nrow'  ,    6  , ...
  'ncol'  ,  [ 30  44 ]  , ...
  'title' ,  { 'KeyWord'  'NetCDF-Variable'}  , ...
  'form'  ,  { ''   ''    } , ...
  'vedit' ,  { 'on' 'on'  } , ...
  'cedit' ,  'on'       , ...
  'cset'  ,  'on'       , ...
  'cins'  ,  'on'       , ...
  'cdel'  ,  'on'       , ...
  'creset',  'on'       , ...
  'cundo' ,  'off'      , ...
  'credo' ,  'off'      , ...
  'min'   ,   0         , ...
  'max'   ,   2         , ...
 'msgfcn' ,  fcn               };

Text = stdgui('Type' , 'text' , 'String' , 'Select DataVariables' );

List = stdgui('Type' , 'tab_list' , 'Color' , cl , 'CBFcn' , fcn, 'Option' , opt);

List = struct( 'Text' , { { 1+0*i  Text } } , ...
               'List' , { { 1+1*i  List } }       );

%**************************************************************

Config(10,:) = { 'Data' , { List } };

Config = permute(Config,[2 1]);

Config = struct(Config{:});

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,t] = recpar(h,typ);

% RECPAR  returns ParentHistory of Handle
%
% [Msg,HandleHist,TypeHist] = RECPAR( Handle , StopType )
%
%  recurse trough the Parents unto ( ParentType == StopType )
%
%  StopType starts with UpperCase, returns History excluding
%     Handle with ( ParentType == StopType )
%
%  default: StopType = 'Root'  (recurse unto 'figure')
%
%  HandleHist(end) == Handle
%    TypeHist(end) == HandleType
%

Msg = '';
 c  = zeros(0,1);
 t  =  cell(0,1);

if nargin < 1
   Msg = 'Input Handle is missing.';
   return
end

if nargin < 2
   typ = 'Root';
end

%-----------------------------------------------

if isempty(h)
   return
end

ok = ( isnumeric(h) &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   Msg = 'First Input must be a Single Handle.';
   return
end

if ~( ischar(typ) & ~isempty(typ) & ...
      ( prod(size(typ)) == size(typ,2) ) )
   Msg = 'Type must be a String';
   return
end

%-----------------------------------------------

c = h;
t = { get(h,'type') };

z = 1;

t0 = lower(typ);

while ~( c(1) == 0 )  &  ( ~strcmp(t{1},t0) | ( z == 1 ) )

   z = z + 1;

   c = cat( 1 ,         get(c(1),'parent')  , c );

   t = cat( 1 , { lower(get(c(1),'type')) } , t );

end


if strcmp( t{1} , t0 )

  n = 1 + strcmp( typ(1) , upper(typ(1)) );

  c = c(n:z);
  t = t(n:z);

else

   Msg = [ 'Handle has no Parents with Type '''  typ '''.' ]; 

end

%-----------------------------------------------
   
%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


