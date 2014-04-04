function [ Msg, KeyWords, Comments, V ] = cnf_chk(CNFile,CNFName,Kmark,Cmark)

% CNF_CHK  First Check for ConfigurationFiles, using LOADASC
%
% Msg = CNF_CHK( ConfigFile )
%
% Returns ErrorMessages in an 3-Row CellArray for:
%
%  Msg = { [ 0 | 1 ]  msg(Using LOADASC(...,'@LOOK')   to get KeyWords)
%          [ 0 | 1 ]  msg(Using LOADASC(...,KeyWords)  to get Data)
%          [ 0 | 1 ]  msg(Check KeyWords and DataErrors (Nkey by 1 Cell))   }
%
%  >> find( cat(1,~Msg{:,1}) )  gives the Location of Error's
%
%
% use
%
%     [ Msg , KeyWords , Comments, V ] = CNF_CHK( ConfigFile )
%
%   to get the KeyWords, Commentaries and the Data in an CellArray
%    (see LOADASC)
%
% use [ ... ] =  CNF_CHK( ConfigFile , ConfigName , KeyMarker , CommentMarker ) 
%
%   to give optional an Name for the Configuration, to use in the Messages,
%       and the Marker, default:
%
%     KeyMarker = '#';  CommentMarker = '%';
%
%
%  see also:  LOADASC  WRITE_ASC
%


Nin = nargin;


nl = char(10);;


Msg        = cell(3,2);
Msg(:,1)   = {1};
Msg(1:2,2) = {''};
Msg{ 3 ,2} = cell(0,0);



KeyWords = cell(0,0);
       V = cell(0,0);


if Nin < 1
 Msg{1} = 'Error using CNF_CHK, not enough InputArguments.';
 return
end

if Nin < 3
 CNFName = '';
end

if Nin < 3
 Kmark = '#';
end

if Nin < 4
 Cmark = '%';
end


if ~isempty(CNFName)
 CNFName = [ ' for ' CNFName ];
end

  

%--------------------------------------------------------------------
% 1. Read Config


[ Msg{1,2}, KeyWords , Comments ] = loadasc(CNFile,'@LOOK',Kmark,Cmark);
  
   Msg{1,1} = isempty(Msg{1,2});

   if ~Msg{1,1}
    Msg{1,2} = ['Invalid ConfigFile'  CNFName  ': ' CNFile  nl ...
                 Msg{1,2} ];
    return
   end
  
[Msg{2,2} , V ] = loadasc(CNFile,KeyWords,Kmark,Cmark);

   Msg{2,1} = isempty(Msg{2,2});

   if ~Msg{2,1}
     Msg{2,2} = ['Invalid Data in ConfigFile'  CNFName  ': ' CNFile ...
                  nl  Msg{2,2}  ];
    return
   end

  

Msg{3,2}    = cell(size(V,1),1);
Msg{3,2}(:) = {''};

  
for ii = 1 : size(V,1)


  msg = first_chk(V{ii,1},KeyWords{ii},Kmark);

    Msg{3,1} = ( Msg{3,1}  &  isempty(msg) );

  if ~isempty(msg)
    Msg{3,2}{ii} = msg;
  end


end
% ii




%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = first_chk(V,KeyWord,Kmark)

% msg = first_chk(V)
%
% First Check for KeyWordData, loaded by LOAD_ASC
%  Check for Duplicate Keyword  and Read_KeyWordError's
%
%  msg is the Message for an DialogBox
%

nl = char(10);


  msg = '';

  % Check for Duplicate Keyword

  if size(V) > 1
   
    msg = [ 'Duplicate Keyword '  Kmark  KeyWord   ];

  end
 

  % Check for Read-Errors

  if ~isempty(V{:,3})
     
     msg = [ msg  nl(1:(end*(~isempty(msg)))) V{:,3} ];

  end
