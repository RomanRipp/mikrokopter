function [msg,file,pre] = write_cnf(ini,file,name,cnf,prop,unit)

% WRITE_CNF  Writes ConfigFile or returns HelpText for Configuration
%
% [Msg,File] = WRITE_CNF( INI , FileName , Name )   Writes ConfigFile
%
% [Msg,Text] = WRITE_CNF( INI , '-' )        Returns HelpText
% [Msg,Help] = WRITE_CNF( INI , '%' )        Returns HelpText with "%"
%
% 
% INI = { FieldName  [{Default}]  [Type]  [Size]  [Optional]  [Comment] }
%
%     FieldName  String
%     Default    DefaultValue
%     Type       String
%     Size       RowVector with Integer or NaN's
%     Optional   CellArray of optional Values
%     Comment    String or CellArray of Strings for Description of Field
%
%
% WRITE_CNF( ... , CNF ) use the Values of the Structure CNF as DefaultValues.
%
% WRITE_CNF( ... , Properties , Units ) appends the Units
%      to the DefaultValues of Properties in the HelpText.
% 
% see also: LOAD_CNF, CHECK_CNF
%

Nin = nargin;
msg = '';
pre = '';

%-----------------------------------------------

if Nin < 2, file = ''; end
if Nin < 3, name = ''; end
if Nin < 4, cnf  = []; end
if Nin < 5, prop = {}; end
if Nin < 6, unit = {}; end

if Nin < 1
  msg = 'Not enough Input Arguments.';
  return
end

if isempty(ini)
   return
end

msg = cell(0,1);

%--------------------------------------------------------------------
% INI

n2 = 6;

s2 = size(ini,2);

ok = ( iscell(ini) & ~isempty(ini) & ( ndims(ini) == 2 )  & ( s2 > 1 ) );

is_val = 1;

if ~ok

    msg = cat( 1 ,  msg , {'INI must be a CellArray: { FieldName [Default] Type [Size] [Optional] [Help] }.'} );

else

    s2 = size(ini,2);

    ok = ( ( s2 < n2 ) & chkcstr(ini(:,2)) );  % 2. Column is Type
    if ok & ( s2 >= 3 )
       ok = ~chkcstr(ini(:,3));                % 3. Column is not Type
    end
    if ok
       ini = ini(:,[1 1 (2:s2)]);              % Insert DefaultValue
       ini(:,2) = {[]};
       s2 = s2 + 1;
       is_val = 0;
    end

    if     s2 > n2

           ini = ini(:,1:n2);

    elseif s2 < n2

           ini = ini(:,[ (1:s2) ones(1,n2-s2) ]);

           if s2 < 2, ini(:,2) = { [] };        end  % Default
           if s2 < 3, ini(:,3) = { 'none' };    end  % Type
           if s2 < 4, ini(:,4) = { [NaN NaN] }; end  % Size
           if s2 < 5, ini(:,5) = { cell(0,1) }; end  % Optional
           if s2 < 6, ini(:,6) = { '' };        end  % Help

    end

    ok = chkcstr(ini(:,1));
    if ok
       ok = ~any(strcmp(ini(:,1),''));
    end

    if ~ok
        msg = cat( 1 ,  msg , {'Elements of 1. Column of INI must be nonempty Strings (FieldNames).'} );
    end

    if ~chkcstr( ini(:,3) )
        msg = cat( 1 ,  msg , {'Elements of 3. Column of INI must be Strings (Type).'} );
    end

    try, v  = cat(2,ini{:,4}); catch, v = zeros(2,1); end

    ok = ( isnumeric(v) & ( ( prod(size(v)) == size(v,2) ) | isempty(v) ) );
    if ok
       ok = all( ( mod(v,1) == 0 ) | isnan(v) );
    end

    if ~ok
        msg = cat( 1 ,  msg , {'Elements of 4. Column of INI (Size) must be RowVectors with Integers or NaNs.'} );
    end

end

%--------------------------------------------------------------------
% File

if ~isempty(file)
    if ~chkstr(file)
        msg = cat( 1 ,  msg , {'FileName must be a String.'} );
    end
end

hlp = ( isequal(file,'-') | isequal(file,'%') );

%--------------------------------------------------------------------

if ( Nin == 4 ) & chkcstr(name) & chkcstr(cnf)
    prop = name;
    unit = cnf;
    name = '';
    cnf  = [];
end

if ( Nin == 5 ) & chkcstr(cnf)
    unit = prop;
    prop = cnf;
    if isstruct(name)
       cnf  = name;
       name = '';
    else
        cnf = [];
    end
end

%--------------------------------------------------------------------

if ~( chkstr(name) | isempty(name) )
    msg = cat( 1 , msg , {'Name must be a String.'} );
end

if ~isempty(cnf)
    if ~( isstruct(cnf)  &  ( prod(size(cnf)) == 1 ) ) 
        msg = cat( 1 , msg , {'CNF must be a Structure with Length 1.'} );
    end
end

if ~isempty(prop)
    [ok,prop] = chkcstr(prop);
    if ~ok
        msg = cat( 1 , msg , {'Properties must be a CellArray of Strings.'} );
    end
end

if ~isempty(unit)
    [ok,unit] = chkcstr(unit);
    if ~ok
        msg = cat( 1 , msg , {'Units must be a CellArray of Strings.'} );
    elseif ~isempty(prop)
        np = prod(size(prop));
        nu = prod(size(unit));
        if ~( ( nu == 1 ) | ( nu > np ) )
            msg = cat( 1 , msg , {'Size of Units must match Size of Properties.'} );
        elseif nu == 1 
            unit = unit(ones(1,np));
        end
    end
end

%--------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    ok  = NaN * zeros(size(ini,1),1);
    return
end

msg = '';

%-----------------------------------------------------------------
% Get File

if isempty(file)

   [file,pfad] = uiputfile( '*.ini' , 'Save ConfigFile ... ' );

   if isequal(file,0) | isequal(pfad,0) | isempty(file)
      file = '';
      return
   end

   file = fullfile(pfad,file);

end
 
%-----------------------------------------------------------------
% Get Text

cm = '';
nl = '';

if hlp
    cm = '';
    nl = '';
else
    cm = '%';
    nl = char(10);
end

pre = cm;

% Insert NewLine in OptComment before ColorStrings
red = '''red''';
rep = sprintf('\n%% %s',red);

ni = size(ini,1);
if hlp
   nn = zeros(ni,5); % Size of [ Parameter TypeSize Optional Comment Value ]
end

fld = {};
is_cnf = ~isempty(cnf);

if is_cnf
   fld = fieldnames(cnf);
end

is_val = ( is_val | is_cnf );

is_prp = ~( isempty(prop) | isempty(unit) );

for ii = 1 : ni

    val = ini{ii,2};   % Default Value
    typ = ini{ii,3};   % Type and Size as Comment
    siz = ini{ii,4};
    opt = ini{ii,5};   % Optional as Comment
    cmt = ini{ii,6};   % Comment

    ini(ii,[2 3 4 5]) = {''};

    %---------------------------------------------
    % Size+Type as Comment

    if chkstr(typ,1) & ~isequal(typ,'none')
       if ~isempty(siz) & isnumeric(siz)
           if ~all(isnan(siz))
               siz = var2mstr(siz(:)',1);
               siz(find(siz==',')) = [];
               ini{ii,2} = [ strrep(siz,'NaN','-')  ' ' ];
           end
       end
       ini{ii,2} = sprintf('%s %s %s%s',cm,ini{ii,2},typ,nl);
    end

    %---------------------------------------------
    % Optional as Comment

    if ~isempty(opt)
        [ok,opt] = chkcstr(opt);
        if ok 
           if ~hlp
               opt = var2mstr(opt(:)',-1);
               if strcmp(opt([1 2]),'{ ');
                  opt = opt(3:end);
               end
               if strcmp(opt([end-1 end]),' }');
                  opt = opt(1:end-2);
               end
               opt = strrep(opt,red,rep);
               opt = sprintf('%s %s%s',cm,opt,nl);
           else
               opt = opt(:);
           end
           ini{ii,3} = opt;
        elseif hlp
           ini{ii,3} = {};
        end
    end

    %---------------------------------------------
    % Comment

    if ~isempty(cmt)
        [ok,cmt] = chkcstr(cmt);
        if ok
           if ~hlp
               frm = sprintf('%s%s S%s',cm,cm,nl);
               frm = strrep(frm,'%','%%');
               frm = strrep(frm,'S','%s');
               cmt = sprintf(frm,cmt{:});
           else
               cmt = cmt(:);
           end
           ini{ii,4} = cmt;
        elseif hlp
           ini{ii,4} = {};
        end
    end

    %---------------------------------------------
    % Value

    if is_val

       ok = is_cnf;
       if ok
          ok = any(strcmp(fld,ini{ii,1}));
       end
       if ok
          val = getfield( cnf , ini{ii,1} );
       elseif iscell(val) & ( prod(size(val)) == 1 )
          val = val{1};
       end

       str = var2mstr(val,1);
       
       if isnumeric(val)
          if ( str(1) == '[' ) & ~isequal(str,'[]')
             if ~hlp
                 str = strrep(str,'[',' ');
                 str = strrep(str,']',' ');
             end
             str = strrep(str,',','');
             str = strrep(str,';',[nl ' ']);
          elseif ( prod(size(val)) == 1 )
             if is_prp
                ok = strcmp(ini{ii,1},prop);
                if any(ok) & hlp
                   ok = find(ok); uni = unit{ok(1)};
                   if imag(val) == 0
                      str = [ str  uni ];
                   else
                      str = strrep('+',[uni '+']);
                      str = strrep('-',[uni '-']);
                   end
                end
             end
             if ~( str(1) == ' ' )
                 str = [ '  '  str ];
             end
          end
          if hlp 
             if any( str == 10 ) 
                str = sepname(str,NaN,char(10));
             else
                str = {str};
             end
          end   
       elseif hlp
          str = {str};
       elseif chkcstr(val)
          str = strrep(str,'''#''','''"#"''');  % Mask KeyWord
          str = strrep(str,'''%''','''"%"''');  % Mask Comment
       end
       
       ini{ii,5} = str;

    end

    if hlp
       for jj = 1 : 5
           nn(ii,jj) = size(ini{ii,jj},1);
       end
    end

end

ini = ini(:,1:5);

%-----------------------------------------------------------------
% HelpText

if hlp

   bl1 = char(32*ones(1,2));
   bl2 = char(32*ones(1,3));
   nl = char(10);

   if file == '%'
      pre = '% ';
   else
      pre = '  ';
   end

   name = char(ini(:,1));

   n1 = size(char(name),2);
   n2 = max( 32 , size(char(ini(:,2)),2) );

   def = '#default: ';

   for ii = 1 : ni

       %-----------------------------------
       % TypeSize+Optional

       tsopt = cell(0,1);
       if ~( nn(ii,2) == 0 )
           tsopt = cat( 1 , tsopt , ini(ii,2) );
       end

       opt = ini{ii,3};
       jj  = 0;
       while jj < nn(ii,3)
             s2 = 0;
             str = '';
             while ( s2 == 0 ) | ( s2+size(opt{jj+1},2)+3 < n2 )
                    jj = jj + 1;
                   str = sprintf('%s ''%s''',str,opt{jj});
                   str = str( ( 1 + ( s2 == 0 ) ) : end );
                   s2  = size(str,2);
                   if jj == nn(ii,3), break, end
             end
             tsopt = cat( 1 , tsopt , {str} );
             if jj == nn(ii,3), break, end
       end

       ini{ii,2} = char(tsopt);

       %-----------------------------------
       % Comment + Value

       if ~isempty(ini{ii,5})
           ini{ii,5} = cat( 2 , ones(nn(ii,5),1)*def , char(ini{ii,5}) );
           ini{ii,5} = cellstr(ini{ii,5});
       else
           ini{ii,5} = cell(0,1);
       end

       if isempty(ini{ii,4})
          ini{ii,3} = ini{ii,5};
       else
          ini{ii,3} = cat( 1 ,ini{ii,4} , ini{ii,5} );
       end

       ini{ii,3} = char(ini{ii,3});
 
       %-----------------------------------
       % Append BlankRows

       nn(ii,2) = size(ini{ii,2},1);
       nn(ii,3) = size(ini{ii,3},1);

       n = max(nn(ii,[1 2 3]));

       ini{ii,1} = cat( 1 , name(ii,:) , char(32*ones(n-nn(ii,1),n1)) );
       ini{ii,2} = cat( 1 ,  ini{ii,2} , char(32*ones(n-nn(ii,2),size(ini{ii,2},2))) );
       ini{ii,3} = cat( 1 ,  ini{ii,3} , char(32*ones(n-nn(ii,3),size(ini{ii,3},2))) );

   end

   ini = insert(ini(:,1:3),(1:ni),{''},1);

   ni = size(char(ini(:,1)),1);
   oi = ones(ni,1);

   file = cat( 2 , char(oi*pre) , char(ini(:,1)) , ...
                   char(oi*bl1) , char(ini(:,2)) , ...
                   char(oi*bl2) , char(ini(:,3))  );

   file = cellstr(file);

   pr = { 'Parameter        [Size]  Type / Values            Comment'
          '                 Optional' };

   bl  = cat( 2 , pre , char(32*ones(1,size(char(pr),2)+3)) );
   bl(find(bl==' ')) = '-';
   bl(1) = pre(1);

   pr = pr(:,[1 1]);
   pr(:,1) = {pre};
   pr = permute(pr,[2 1]);
   pr = sprintf('%s%s\n',pr{:});

   pr = sprintf('%s\n%s%s\n%s',bl,pr,bl,pre);

   file = sprintf('%s\n',pr,file{:});

   return

end

%-----------------------------------------------------------------
% Write File

txt = ' ConfigurationFile';
if ~isempty(name)
    txt = sprintf('%s for %s',txt,name);
end
lin = char( '-' * ones(size(txt)+[0 3]) );
txt = sprintf('%s%s\n%s%s\n%s%s\n\n',cm,lin,cm,txt,cm,lin);

ini = permute(ini,[2 1]);

frm = '# %s\n%s%s%s\n%s\n\n';

fid = fopen(file,'wt');
if fid == -1
   msg = sprintf('Can''t open ConfigFile "%s" for writing.',file);
   return
end

fprintf(fid,'%s',txt);

fprintf(fid,frm,ini{:});

fclose(fid);


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

