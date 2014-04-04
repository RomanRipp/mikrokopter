function str = ger2eng(str)

% GER2ENG   Converts German MonthNotation to English
%
% EnglishString = GER2ENG( GermanString );
%
%  mär --> mar
%  mai --> may
%  okt --> oct
%  dez --> dec
%
% In other case you get problems to use DATEVEC, DATENUM
%

if nargin < 0
   str = '';
   return
end

if ~( iscellstr(str) | ischar(str) );
   error('Input must be a CharArray or CellStringArray.');
end


is_char = ischar(str);

if is_char
   str = cellstr(str);
end

rep = { ...
  'mär'  'mar'
  'mai'  'may'
  'okt'  'oct'
  'dez'  'dec'  };


for ii = 1 : prod(size(str))

    s = str{ii};

    if ~isempty(s)  

        s = lower(s); 

        for jj = 1 : size(rep,1)

            kk = findstr( s , rep{jj,1} );

            if ~isempty(kk)

               n = size(rep{jj,1},2);

               for ll = kk(:)'

                   ind = ll + (1:n) - 1;     

                    nn = find( ~( double( s(ind) ) == double( str{ii}(ind) ) ) );

                    str{ii}(ind) = rep{jj,2};

                    str{ii}(ind(nn)) = upper(str{ii}(ind(nn)));

               end
               % ll == kk(findstr)
            end
            % ~isempty(findstr) 
        end
        % jj = rep
    end
    % ~isempty(str{ii})
end
%ii

if is_char
   str = char(str);
end