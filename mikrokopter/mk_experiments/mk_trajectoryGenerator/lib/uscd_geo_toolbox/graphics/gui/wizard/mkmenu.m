function [Msg,HM] = mkmenu(par,group,cnf,tag);

% MKMENU  Creates UIMenus from Structure
%
% MKMENU( ParentHandle ,  Group , Config , Tag );
%
%
%  FieldValues of Config for UIMenu-Definition:
%
%    { Label  SeparatorOn  CallbackFcn  UserData }
%    {          0 | 1      SubStructure          }
%
%  FieldValues for UIContextMenu-Definition:
%
%    { NaN  MenuStructure  CallbackFcn  UserData }
%


%-----------------------------------------------------------

Msg = '';
HM  = []; 

nl = char(10);

if isempty(cnf)
  return
end


%************************************************************


field = fieldnames(cnf);
field = field(:)';

    n = size(field,2);

%------------------------------------------------------------
% HandleStructure

HM      = cell( 2 , n );
HM(1,:) = field;
HM(2,:) = { {[]} };

HM      = struct( HM{:} );

HC      = zeros( n , 1 );

%------------------------------------------------------------

sets = { 'off'  'on' };

if ~isempty(group) & ~isempty(tag)
   tag = cat( 2 , tag , '.' , group );
elseif isempty(tag)
   tag = group;
end

for ii = 1 : n

  cc = getfield( cnf , field{ii} );


  CB    = '';
  usd   = [];

  if ischar( cc{3} )  &  ~isempty(cc{3})
       CB = [ cc{3} '(gcbo,'''  group  ''',''' field{ii} ''',1);'  ];
  end

  if prod(size(cc)) >= 4
     usd = cc{4};
  end

  tag1 = cat( 2 , tag , '.' , field{ii} );

  if ischar(cc{1})

    %----------------------------------
    % UIMenu

    enb = sets{ 1 + (~isempty(rmblank(cc{1},2))) };
    sep = sets{ 1 + isequal(cc{2},1) };
 
    HC(ii) = uimenu( ...
       'parent'          , par      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'label'           , cc{1}    , ...
       'tag'             , tag1     , ...
       'enable'          , enb      , ...
       'separator'       , sep      , ...
       'checked'         , 'off'    , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );

     ch = cc{3};

   else

     %----------------------------------
     % UIContextMenu

     HC(ii) = uicontextmenu( ...
       'parent'          , par      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'tag'             , tag1     , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );
 
     ch = cc{2};

   end 

  if isstruct(ch)
  % Children

    try
       [Msg,ch] = mkmenu(HC(ii),field{ii},ch,tag);
    catch
        Msg = lasterr;
    end

    if ~isempty(Msg) 
       Msg = [ tag  ': Error call MKMENU( '  field{ii}  ' ).' ...
               nl Msg ];
       return
    end

  else

    ch = [];

  end

  setappdata( HC(ii) , 'Children' , ch );

  HM = setfield( HM , field{ii} , HC(ii) );

end

