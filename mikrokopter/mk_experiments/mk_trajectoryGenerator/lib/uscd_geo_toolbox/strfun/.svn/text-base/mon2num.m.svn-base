function str = mon2num(str,mon,cc);

% MON2NUM  Replace a MonthString in  String by a Number
%
% STR = MON2NUM(STR)
%
% Replace in the String STR all english and german 
%  abbrevations (3 characters) for months by their Number.
%
% STR = MON2NUM( STR , MonthString , ReplaceCharacter )
%
% MonthString define the strings to replace by the Number
% MonthString can be CharacterArray or CellArray wich 
%  elements are Strings or CellArray of Strings for multiple
%  matches per Number. MON2NUM without any Input returns the
%  default value for MonthString:
%
%  MonthString = MON2NUM
%
% ReplaceCharacters are characters surrounding the strings 
%  from MonthString which will replaced by a Blank.
%  default: '-./+*^\:'   mathematical Operators
%
%
% see also: FINDSTR, STRREP
%

Nin = nargin;

%-------------------------------------------------

if Nin > 1
   if isempty(str)
      return
   elseif ~chkstr(str)
      error('String required.');
   end
end

%-------------------------------------------------

if Nin < 2

   mon = {  { 'JAN'       }
            { 'FEB'       }
            { 'MAR' 'MÄR' }
            { 'APR'       }
            { 'MAY' 'MAI' } 
            { 'JUN'       }
            { 'JUL'       }
            { 'AUG'       }
            { 'SEP'       }
            { 'OCT' 'OKT' }
            { 'NOV'       } 
            { 'DEC' 'DEZ' }  };

else

   if ischar(mon)
      mon = cellstr(mon);
   end

   if ~iscell(mon)
       error('MonthString must be a CharacterArray or CellArray.');
   end

   mon = mon(:);

   n  = size(mon,1);

   for ii = 1 : size(mon,1)
       [ok,mon{ii}] = chkcstr(mon{ii});
       if ~ok
           error('Elements of MonthString must be Strings or CellArrays of Strings.');
       end
       mon{ii} = upper(mon{ii}(:)');  % Row
   end

end

%-------------------------------------------------

if Nin < 3
   cc = '-./+*^\:';
elseif ~isempty(cc)
    if ~ischar(cc)
        error('Invalid Input for ReplaceCharacters.')
    end
    cc = cc(:)';
end

%-------------------------------------------------

if Nin < 1
   if nargout == 0
      fprintf(1,'\n');
      for ii = 1 : size(mon,1)
          len = size(mon{ii}{1},2);
          frm = sprintf('%%%.0f.%.0fd',len,len);
          rep = sprintf(frm,ii);
          fprintf(1,'   %s : {',rep);
          fprintf(1,' ''%s'' ',mon{ii}{:});
          fprintf(1,'}\n');
      end
      fprintf(1,'\n');
      clear str
   else
      str = mon;
   end
   return
end

%*************************************************

% UpperCase

upp = upper(str);

siz = size(str,2);

for ii = 1 : size(mon,1)

    for mm = mon{ii}

        len = max(1,size(mm{1},2));
        frm = sprintf('%%%.0f.%.0fd',len,len);
        rep = sprintf(frm,ii);
        ins = ~( size(rep,2) == len );

        if len == siz

           if strcmp( upp , mm{1} )
              str = rep;
              if ins
                 upp = upper(str);
                 siz = size(str,2);
              end
           end
  
        elseif ( 0 < len ) & ( len < siz )

           kk = findstr(upp,mm{1});
           if ~isempty(kk)
               for jj = kk(end:-1:1)
                   ll =  jj + [ -1  len ];
                   ll =  ll( (1+(ll(1)<1)) : (2-(ll(2)>siz)) );
                   if ~any(isletter(str(ll)));
                       if ~isempty(cc)
                           for c = cc
                               ok = ( str(ll) == c );
                               if any(ok)
                                   ok = find(ok);
                                   str(ll(ok)) = ' ';
                               end
                           end
                       end
                       if ins
                          str = cat( 2 , str(1:ll(1)) , rep , str(ll(2):end) );
                       else
                          str(jj-1+(1:len)) = rep;
                       end
                   end
               end
               if ins
                  upp = upper(str);
                  siz = size(str,2);
               end
           end

        end

    end

end

