function [ Msg , fig , H , T , V ] = b_create(pfad,file,Name,Logo);

% B_CREATE  Creates BrowserFigure
%
% [ Msg , BrowserFigure ] = B_CREATE( HelpPath , MatFile )
%
%  creates a BrowserFigure with the Items, specified
%   in the ItemStructure of the Variable "V", which is
%     in the MAT-File "<HelpPath>/<MatFile>", 
%  or in the ItemStructure of the "BROWSER.INI"-Files,
%    which should exist DirectoryStructure of "HelpPath". 
%
%
% [ Msg , BrowserFigure ] = B_CREATE( HelpPath , MatFile , Name , Logo )
%
%   specifies the FigureName and the Logo of the BrowserFigure.
%
%
% [ Msg , BrowserFigure , IH , IT , V ] = B_CREATE( HelpPath , ... )
%
%   returns the Handles (IH), and the Tags (IT) of the 
%    Items, and the ItemStructure V.
%
%  Msg contains ErrorMessages, empty if all ok.
% 
%  


Msg = '';
fig = [];
H   = [];
T   = {};
V   = [];

Msg0 = 'B_CREATE: ';

nl = char(10);

Nin = nargin;

%-----------------------------------------------------------

if Nin < 1
  Msg = [ Msg0 'Input HelpPath is missing.' ];
  return
end

if Nin < 2
  file = '';
end

if Nin < 3
  Name = 'Browser';
end

if Nin < 4
  Logo = { '@'  'Online'  'Documentation' };
end

%-----------------------------------------------------------

if ~( ischar(pfad)  &  ~isempty(pfad)  &  ...
      ( prod(size(pfad)) == size(pfad,2) ) ) 

  Msg = [ Msg0  'HelpPath must be a nonempty String.' ];

  return

end

if ~( ischar(file)  &  ...
      ( prod(size(file)) == size(file,2) ) ) 

  Msg = [ Msg0  'MatFile must be a String.' ];

  return

end


%------------------------------------------------------
% Append FileSperator at pfad 

 fs   = filesep;
 pfad = cat( 2 , pfad , ...
                 fs(1:(end*(~strcmp(pfad(end),fs)))) );

%------------------------------------------------------
% Try to load an Exisiting ItemStructure 
%  from MAT-File in pfad

if ~isempty(file)

  % Try to load an Exisiting ItemStructure 
  %  from MAT-File in pfad

   file = cat( 2 , pfad , file );

   try
     v = whos('-file',file);
   catch
     v = [];
   end

   if ~isempty(v)

     if any( strcmp( cellstr(str2mat(v.name)) , 'V' ) )
  
        fprintf(nl)
        fprintf([Msg0 'Try Structure from '  ...
                      strrep(file,'\','\\')  '  ...  ' ]);

        load(file,'V');  % ItemStrucure

       [Msg,fig,H,T,V] = browser('New',Name,Logo,V);

       if isempty(Msg)

           fprintf([ 'ok'  nl ]);

           return

       end

           fprintf([ 'Error'  nl ]);

     end

  end

end
%------------------------------------------------------

% MAT-File doesn't exist or contains invalid ItemStructure

  fprintf(nl)

  fprintf([Msg0 'Read Structure of '  ...
                strrep(pfad,'\','\\')  '  ...  ' ]);

  [Msg,fig,H,T,V] = browser('New',Name,Logo,pfad);
 
  if isempty(Msg)
  % successfull

    fprintf([ 'ok'  nl ]);

    %--------------------------------------------
    % Save Structure into MAT-File
    if ~isempty(file)

        fprintf([Msg0 'Save Structure "V" into '  ...
                    strrep(file,'\','\\')  '  ...  ' ]);
  
        try

           save(file,'V','-mat');

           fprintf([ 'ok'  nl ]);

        catch

           fprintf([ 'Error'  nl ]);

        end  
        
    end
    %--------------------------------------------

    if nargout < 2

       set(fig,'visible','on');

    end

  else

     fprintf([ 'Error'  nl ]);

     Msg = [ Msg0 'Error create BrowserFigure.' ...
                   nl  Msg ];

  end



