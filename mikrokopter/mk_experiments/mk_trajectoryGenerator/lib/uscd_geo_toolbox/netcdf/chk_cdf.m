function [Msg,FileOk,IsCDF,Dim1,Dim,Var,Att,Txt] = chk_cdf(DefaultFile,Label,...
                                             Files,Labels,MsgCall)

% CHK_CDF   Checks Variables and Dimensions of  NetCDF-Files
%
% Checks, if all Dimensions and Variables with their Dimensions 
%   from a DefaultFile are also contained in NetCDF-Files.
%
% [ Msg , FileOk, IsCDF, DimFiles ] = CHK_CDF( DefaultFile, Label, ...
%                                  Files, Labels, MsgCall )
%
% Inputs:
%  DefaultFile is the FileName of the DefaultFile
%  Label  is a Label for the DefaultFile for Messages
%
%  Files is a CellStringArray of Files to check [ NFiles by 1 ]
%  Labels are Labels for the Files for Messages, 
%         [ 1 by 1 ] or [ NFiles by 1 ],  CellStringArray
%
%  MsgCall is the Command to evaluate for a Message, 
%          contains the Variable:  "msg" = { Text  'new'|'append' }
%
%  If Files are empty or not given, 
%    only the DefaultFile will be checked.
%
%
% OutPuts:
%  Msg      Message, if any Error in Files Occured,
%              empty, if all ok !!!!!!
%  FileOk   [ NFiles by 1 ], contains 1 if File ok, 0 if not ok.
%  IsCDF    [ NFiles by 1 ], contains 1 if File is a NetCDF-File, 0 if not.
% 
%  DimFiles  CellArray for Dimensions of the Files, by LOOK_CDF,
%             with NFiles as Length for 3. Dimension.
%  
% Other Outputs:  [ ... , Var, Att, Text ] = CHK_CDF( ... )
%
%  gives the Output from 
%       [ ... , Dim, Var, Att, Text ] = LOOK_CDF( DefaultFile ) 
%
%
% See also:  LOOK_CDF  CAT_CDF
%

Nin = nargin;

Msg      = '';
FileOk   = [];
IsCDF    = [];
DimFiles = {};

Dim = {};
Var = {};
Att = {};
Txt = '';


Msg0 ='CHK_CDF: ';

nl = char(10);

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);
nl2 = char([ 10 32*ones(1,nm0+2) ]);
nl4 = char([ 10 32*ones(1,nm0+4) ]);


%-------------------------------------------------
% Check Inputs

if Nin < 5
 MsgCall = '';
end
if Nin < 4
 Labels = 'File';
end
if Nin < 3
 Files = {};
end
if Nin < 2
 Label = 'File';
end



if ischar(Files)
 if isempty(Files)
  Files = {};
 else
  Files = cellstr(Files);
 end
end
Files = Files(:);

if iscellstr(Label)
 Label = Label{1};
end

 if ischar(Labels)
   Labels = cellstr(Labels);
 end
 Labels = Labels(:);

 NF = size(Files,1);

 NL = size(Labels,1);

if NL < NF
  Labels = [ Labels ; Labels( NL*ones(NF-NL,1) ) ];
end

Dim1 = cell(0,2,NF+[~NF]);

FileOk = zeros(NF,1);
IsCDF  = zeros(NF,1);


%*****************************************************
% Read DefaultFile

  msg = {[ 'Read data from default ' Label ': ' ...
            DefaultFile ]  'new' };
  eval(MsgCall,'0;');

  [Msg,Dim,Var,Att,Txt] = look_cdf(DefaultFile);

  if ~isempty(Msg)
   Msg = [Msg0 'Error Read default '  Label  ': ' ...
                DefaultFile nl Msg  ];
   msg = { 'Error'  'append' };
   eval(MsgCall,'0;');
   return
  end

  NDim = size(Dim,1);
  NVar = size(Var,1);
 
  Dim1 = cell(NDim,size(Dim,2),NF+[~NF]);
  Dim1(:,:,:) = Dim(:,:,ones(size(Dim1,3),1));

  Dim1(:,2,:) = { NaN };



%**********************************************************
% Check Files with DefaultFile

  for ii = 1 : NF      

     msg = {[ 'Read data from '  Labels{ii} ': ' ...
               Files{ii} ]  'new' };
     eval(MsgCall,'0;');
 
    [msg,dim1,var1] = look_cdf(Files{ii});

    if ~isempty(msg)
       Msg = [Msg  nl0(1:(end*(~isempty(Msg)))) ...
              'Error Read '  Labels{ii} ': ' Files{ii} nl msg  ];
       msg = { 'Error'  'append' };
       eval(MsgCall,'0;');
    else

      IsCDF(ii)  = 1;

      msg = '';
      
      %----------------------------------------------------------------
      % Check Dimensions

      dim_ind = zeros(NDim,1);
      for jj = 1 : NDim
        % Compare DimensionNames
        if any(strcmp(Dim{jj,1},dim1(:,1)))
          dim_ind(jj) = find(strcmp(Dim{jj,1},dim1(:,1)));
        end
      end

      if any(~dim_ind)
         jj = find(~dim_ind);  
         msg = [msg  nl0(1:(end*(~isempty(msg)))) ...
                'Missing Dimensions in '  Labels{ii}  ': '  Files{ii} ...
                nl4 strhcat(Dim(jj,1),', ',4,nl4)    ];
      end
 

      %--------------------------------------------------------------  
      if isempty(msg)
      % Dimensions ok ==>  Check Variables

        var_ind = zeros(NVar,3);   %  [ VarNameOk VarDimOk VarTypeOk]
        for jj = 1 : NVar
          % Compare VariableNames
          if any(strcmp(Var{jj,1},var1(:,1)))

            var_ind(jj,1) = find(strcmp(Var{jj,1},var1(:,1)));

            % Compare VariableDimensions
            var_ind(jj,2) = isequal( Dim( Var{jj,4}+1 , 1 ) , ...
                                    dim1( var1{var_ind(jj,1),4}+1 , 1 ) );

            % Compare VariableType
            var_ind(jj,3) = isequal( Var{jj,2} , var1{var_ind(jj,1),2} );

          end
        end

        if any( ~var_ind(:,1) )
           jj = find( ~var_ind(:,1) );  
           msg = [msg  nl0(1:(end*(~isempty(msg)))) ...
                  'Missing Variables in '  Labels{ii}  ': '  Files{ii} ...
                  nl4 strhcat(Var(jj,1),', ',4,nl4)    ];
        end

        if any( var_ind(:,1) & ~var_ind(:,2) )
           jj = find( var_ind(:,1) & ~var_ind(:,2) );  
           msg = [msg  nl0(1:(end*(~isempty(msg)))) ...
                  'Invalid Dimensions for Variables in '  Labels{ii}  ': '  Files{ii} ...
                  nl4 strhcat(Var(jj,1),', ',4,nl4)    ];
        end

        if any( var_ind(:,1) & ~var_ind(:,3) )
           jj = find( var_ind(:,1) & ~var_ind(:,3) );  
           msg = [msg  nl0(1:(end*(~isempty(msg)))) ...
                  'Invalid Type for Variables in '  Labels{ii}  ': '  Files{ii} ...
                  nl4 strhcat(Var(jj,1),', ',4,nl4)    ];
        end

      end
      % isempty(msg),   Check Variables
      %--------------------------------------------------------------  

  
      if isempty(msg)
           FileOk(ii) = 1;
         Dim1(:,:,ii) = dim1(dim_ind,:);
      else
          Msg = [Msg  nl0(1:(end*(~isempty(Msg)))) msg ];
          msg = { 'Invalid File'  'append' };
          eval(MsgCall,'0;');        
      end

    end
    % ~isempty(msg)

  end
  % ii = 1 : NF


 if ~isempty(Msg)
  Msg = [ Msg0  'Error check Files' nl0 Msg ];
 end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ~( ischar(str)  |  iscellstr(str) )
   error('StringArray must be a CharArray or CellStringArray.');
end

if iscellstr(str)
  str = char(str);
end

str = double(str);

    jj    = find( sum( ( str == 32 ) , 2 ) == size(str,2) );
str(jj,:) = [];
 
str = cellstr(char(str));


str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


