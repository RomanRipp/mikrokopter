function [c,cmap] = catxpm(inform,outfile,BlankColor)

% CATXPM   Concatenates Multiple Frames from XPM-Files into MovieMatrix
%
% [ C , ColorMap ] = CATXPM( XPM_Frame_File_Format )
%
%  C  indexed ColorMatrice, refers to ColorMap, Frames in 3. Dimension
%
%  C = [ Ny  by  Nx  by  NFrames ]
%    
%  XPM_Frame_File_Format  Format for multiple FramesFiles, using sprintf
%
%  The FileNames of the XPM-Files to read will build like:
%
%    XPM_Frame_File = sprintf( XPM_Frame_File_Format , FrameNumber );
%
%    By default, CATXPM starts with FrameNumber = 0,
%    The FrameNumber will raised by 1.
%
%    If the specified XPM_Frame_File is not found, the Loop breaks.
%
%  To give the explicit Numbers of Frames to read, use:
%
%    CATXPM( { XPM_Frame_File_Format  FrameNumbers } )
%
%    where FrameNumbers is a Vector.
%
%  To Write the Frames into a single XPM-File, using WRTXPM, use:
%
%    CATXPM( XPM_Frame_File_Format , XPM_Out_File , BackGroundColor )
%
%    The multiple Frames are concatenated in the  XPM_Out_File, 
%     use READXPM to read this File with multiple Frames.
%
%    If the Input BackGroundColor is given, this Color will added to
%     the ColorMap of the XPM_Out_File as recommended BackGroundColor
%     to use for Color None (NaN-Values of C).
%     This last Color in the ColorMap is masked with Blank-Characters.
%
%  C = CATXPM( ... ) returns a True-ColorMatrice, Frames in 4. Dimension.
%
%  C = [ Ny  by  Nx  by  3  by  NFrames ]
%
%--------------------------------------------------------------------------
%
%  see also: READXPM  WRTXPM
%
%--------------------------------------------------------------------------
%
% example:
%
%% You work under Linux/Unix, ImageMagic and Gimp are installed.
%%
%% you have an animated GIF with N Frames: anim.gif
%%
%% make shure, that the animation is UnOptimized, use GIMP:
%%
%% > gimp anim.gif
%%   RightClick on Image, select <Fileter>-<Animation>-<Animation UnOptimize>
%%   RightClick on new Image, select <File>-<Save As ...>
%%    and save the UnOptimized Image as GIF, select Save as <Animation>
%%
%% extract the Frames, use the ImageMagick-Tool CONVERT:
%%
%% > convert anim.gif anim.xpm
%%
%% now you have N XPM-Files: anim.xpm.0, anim.xpm.1, ... anim.xpm.N
%%
%
%  [c,cmap] = catxpm( 'anim.xpm.%.0f' );
%
%  s = size(c);
%
%  figure('units'    , 'pixels' , ...
%         'position' , [ 100  100  s(2)  s(1) ] , ...
%         'colormap' , cmap   , ...
%         'menubar'  , 'none'        );
%
%  axes('units'    , 'normalized'             , ...
%       'position' , [ 0  0  1  1 ]           , ...
%       'xlim'     , [ 1  s(2) ] + 0.5*[-1 1] , ...
%       'ylim'     , [ 1  s(1) ] + 0.5*[-1 1] , ...
%       'ydir'     , 'reverse'  , ...
%       'visible'  , 'off'      , ...
%       'nextplot' , 'add'              );
%
%  h = image( 'cdata'        , c(:,:,1) , ...
%             'cdatamapping' , 'direct' , ...
%             'erasemode'    , 'none'         );
%  while 1
%   for ii = 1 : s(3)
%       set(h,'cdata',c(:,:,ii));
%       drawnow
%       pause(0.1);
%   end
%  end
%
%% Or use:
%
% catxpm( 'anim.xpm.%.0f' , 'anim.xpm' , [ 0 0 0 ]);
%
% [Msg,c,cmap,ccm] = readxpm( 'anim.xpm' );
% 
% if all( double(ccm(end,:)) == 32 )     % BackGroundColor
%    c(find(isnan(c))) = size(cmap,1);
% end
%



nl = char(10);

Nin  = nargin;
Nout = nargout;

%-------------------------------------------------------
% Check for: { Format Nr }

nr = [];

ok = ( iscell(inform) & ( prod(size(inform)) >= 2 ) );
if ok
       nr = inform{2};
   inform = inform{1};
end

% Check FormatString

ok = ( ischar(inform) & ~isempty(inform)  & ...
       ( prod(size(inform)) == size(inform,2) )    );

if ~ok
   error('Format must be a String.');
end

%-------------------------------------------------------
% Check for OutputFile

if Nin < 2

 outfile = '';

else

  ok = ( ischar(outfile) & ~isempty(outfile)  & ...
         ( prod(size(outfile)) == size(outfile,2) )    );

  if ~ok
     error('OutputFilename must be a String.');
  end

end

%-------------------------------------------------------
% Check for BlankColor

if Nin < 3

 BlankColor = zeros(0,3);

else

  ok = ( isnumeric(BlankColor)  & isequal(size(BlankColor),[1 3]) );

  if ok
     if isa(BlankColor,'uint8')
        BlankColor = double(BlankColor) / 255;
     end
     if any(isnan(BlankColor))
        BlankColor = NaN*ones(1,3);
     else
        ok = all( ( 0 <= BlankColor ) & ( BlankColor <= 1 ) );
     end
  end

  if ~ok
     error('BlankColor must be a RGB Color-Triple with Values between 0 and 1.');
  end

end


%-------------------------------------------------
% AutoNumbering

is_nr = ~isempty(nr);

if ~is_nr

   nr = ( 0 : 1001 );

end

%*******************************************************
% Read Data, Loop over Files

c    = [];
cmap = zeros(0,3);

for ii = nr(:)'

    file = sprintf(inform,ii);

    fprintf(nl);
    fprintf([ 'Try ' file ' ... ']);

    if ( exist(file,'file') == 2 ) 

       [Msg,c,cmap] = addxpm(c,cmap,file);

       if isempty(Msg)
          fprintf('ok');
       else
          fprintf('error');
          fprintf(nl);
          fprintf(nl);
          fprintf(Msg);
          fprintf(nl);
       end

    else

       fprintf('not found');

       if ( ii >= 1 ) & ~is_nr
          fprintf(', break');
          fprintf(nl)
          break
       end

    end

    fprintf(nl)

end

fprintf(nl);
   
%*******************************************************
% Write OutputFile

if ~isempty(outfile)

  %---------------------------------------------------
  % Cat BlankColor if all ok

  if ~isempty(BlankColor) 
      cmap = cat( 1 , cmap , BlankColor );
  end

  %---------------------------------------------------
  % Call WRTXPM: [Msg,CharacterColorMap,ColorMapText] = WRTXPM( ... )

  fprintf('Write XPM-File: %s  ... ',outfile);

  [Msg,ccm,ctxt] = wrtxpm(c,cmap,outfile);

  if isempty(Msg)

     fprintf('ok');

  else

     fprintf('error');
     fprintf(nl);
     fprintf(nl);
     fprintf(Msg);
 
  end

     fprintf(nl);
     fprintf(nl);

  %---------------------------------------------------
  % Set BlankColor if all ok

  if ~isempty(BlankColor) & isempty(Msg)

     cpp = size(ccm,2);
     ind = ( 1 : cpp ) + 1;

     ct = ctxt;

     ct(end,ind) = char(32*ones(1,cpp));

     ct = cellstr(ct);

     comm = cat( 2 , '/* Recommended BackgroundColor */' , nl );

     ct = cat( 1 , ct(1:end-1) , {comm} , ct(end) );

     ct = strhcat(ct,'',size(ct,1)+1);

     ctxt = permute(ctxt,[2 1]);
     ctxt = char(ctxt(:)');

     fprintf('Set BackgroundColor to XPM-File: %s ... ',outfile);

     msg = '';

     try

       fid = fopen(outfile,'r');
       bb = fread(fid,'char');
       fclose(fid);

       bb = strrep( char(bb(:)') , ctxt , ct );

       fid = fopen(outfile,'wt');
       fprintf(fid,bb);
       fclose(fid);

     catch

       msg = lasterr;

     end

     if ~isempty(msg)

        fprintf('ok');

     else


       fprintf('error, file may be corrupt')
       fprintf(nl);
       fprintf(nl);
       fprintf(msg);
 
     end

     fprintf(nl);
     fprintf(nl);


   end

end

%------------------------------------------------------------
% RGB-Image

if Nout < 2

 cmap = [ cmap ; NaN*ones(1,3) ];

 c(find(isnan(c))) = size(cmap,1);

 c = permute(c,[1 2 4 3]);

 n = size(cmap,1);

 c(:,:,2,:) = c(:,:,1,:)+1*n;
 c(:,:,3,:) = c(:,:,1,:)+2*n;

 c = cmap(c);


end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,cmap] = addxpm(c,cmap,file)

 [Msg,cc,cm] = readxpm(file);

 if ~isempty(Msg)
    return
 end


 %--------------------------------------------------------
 % New Image

 if isempty(c)

    c = cc;
    cmap = cm;

    return

 end

 %--------------------------------------------------------
 % Check ImageSize

 s0 = size(c);
 s1 = size(cc);

 if ~isequal(s0([1 2]),s1)
    Msg = sprintf('Invalid Size of Image: %.0f x %.0f; must be: %.0f x %.0f' , ...
                  s1(2),s1(1),s0(2),s0(1));
    return
 end

 %--------------------------------------------------------
 % Append New Image

 c = cat( 3 , c , cc );

 %--------------------------------------------------------
 % Compare ColorMap

 if isequal(cm,cmap)

    return

 end

  
 s0 = size(cmap,1);
 s1 = size(cm,1);

 ok = zeros(s1,2);

 for ii = 1 : s1

     eq = ( sum( cmap == cm(ii*ones(s0,1),:) , 2 ) == 3 );

     ok(ii,1) = any(eq);

     if ok(ii,1)

        ok(ii,2) = min(find(eq));

     else

        ok(ii,2) = max(max(ok(:,2)),s0)+1;

     end

 end

 %------------------------------------------------
 % Append New Colors

 if any(~ok(:,1))

    jj = find(~ok(:,1));
    cmap = cat( 1 , cmap , cm(jj,:) );

 end

 %------------------------------------------------
 % Check ColorOrder

 if isequal( ok(:,2) , cumsum(ones(s1,1),1) )

    % No Change

    return

 end

 %------------------------------------------------
 % New ColorOrder of cc

 n3 = size(c,3);

 cc = NaN * cc;

 for ii = 1 : s1

     cc( find( c(:,:,n3) == ii ) ) = ok(ii,2);
 
 end

 c(:,:,n3) = cc;

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function translate(Pfad,str1,str2,varargin)

% TRANSLATE  translates String in Files to another 
%
% translate(Pfad,str1,str2,'-r','.ext')
% translate(File,str1,str2,'-r','.ext')
%

VarArg = varargin;
VarArg = VarArg(:);

Nin = size(VarArg,1);

nl = char(10);


recurse = 0;
ext     = '.m';

for ii = 1 : Nin
  if strcmp(VarArg{ii}(1),'.');
    ext = VarArg{ii};
  end
  recurse = ( recurse | strcmp(VarArg{ii},'-r') );
end

if isempty(Pfad)
 Pfad = cd;
end

if any(strcmp(Pfad,{ '..' }))
  return
end

if any(strcmp(Pfad,{ '.'  '*' }))
  Pfad = cd;
end

fprintf(Pfad);

d = dir(Pfad);

if isempty(d)
  fprintf([' ... can''t read File or Directory ' nl ]);
  return
end

IsFile = ( (size(d,1)==1)  &  (d(1).isdir==0) );
if ~IsFile
  if ~strcmp(Pfad(length(Pfad)),filesep)
    Pfad = [ Pfad filesep ];
  end
end


fprintf(nl)

is_file = find(~cat(1,d.isdir));

for ii = is_file(:)'

  name = d(ii).name;

  ok = IsFile;

  if ~ok

    % Check for m-File
    kk = findstr(name,ext);
  
    ok = ( length(kk) == 1 );
 
    if ok

       ok = (  kk == length(name)-length(ext)+1 );
 
    end

  end

  if ok

      file = [ Pfad((1:(end*(~IsFile))))  name ];

      fid = fopen(file,'r');
      if fid == -1
        fprintf(['   Error read '  file nl ]);
      else
        bb = fread(fid,'char');
        fclose(fid);
        if any(bb==0)
           fprintf(['   Invalid Character in '  file nl ]);
        else
          bb = strrep(char(bb(:)'),str1,str2);
          bb = strrep(char(bb(:)'),'%','%%');
          bb = strrep(char(bb(:)'),'\','\\');
          fid = fopen(file,'wt');
          if fid == -1
            fprintf(['   Error open to write '  file nl ]);
          else
            fprintf(fid,bb);
            fclose(fid);
          end
        end
      end 

  end

end

if recurse

  is_dir = find(cat(1,d.isdir));

  for ii = is_dir(:)'

    if ~any(strcmp(d(ii).name,{ '.'  '..' }))

      translate([Pfad d(ii).name],str1,str2,'-r',ext);

    end

  end

end      

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
 n = 10;
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

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


