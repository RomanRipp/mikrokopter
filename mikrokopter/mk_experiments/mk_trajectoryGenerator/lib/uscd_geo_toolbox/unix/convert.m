function w = convert(varargin)

% CONVERT  Calls the ImageMagicTool CONVERT
%
%  convert [ options ... ] file [ file... ] file
%  convert [ options ... ] @filelist        file
%
% CONVERT without an Input lists the Options
%
% Msg = CONVERT( ... )
%
% Returns the standard Output in Msg
%

command = ['convert ' sprintf('%s ',varargin{:}) ];

fprintf(1,'\nUNIX: %s\n\n',command);

%%% msg = lib_kim(0);

[s,w] = unix(command);

%%% msg = lib_kim(1);

if nargout == 0
   if ~isempty(w)
       fprintf(1,'\n%s\n\n',w);
   end
   clear w
end
