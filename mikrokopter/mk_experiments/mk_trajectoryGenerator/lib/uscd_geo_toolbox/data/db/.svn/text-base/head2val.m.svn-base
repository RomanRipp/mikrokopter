function [val,str,msg] = head2val(val,frm,len,vnan,dlm);

% HEAD2VAL  Converts between HeaderValue and HeaderString
%
% [Value,String,Msg] = HEAD2VAL( Value  , Format )
% [Value,String,Msg] = HEAD2VAL( String , Format )
%
% Format = '%s'                  % String
%          '%f' '%g' '%e' '%d'   % [ 1 by Length ]
%          'lon'  'lat'  'pos'   % Lon / Lat / [ Lat Lon ]
%          'date' 'time'         % Date & Time
%          'list'                % CellString / Delimiter
%
%          'data'                % Data, any numeric
%
% HEAD2VAL( ... , Length , Dummies , [Delimiter] )
%
% Length = Length of Vector or ASCII-Code for Delimiter if CellString
%
% Dummies = Vector of possible Values for NaN
%
% In case of CellString the Length can be the ASCII-Code for the Delimiter
%  or the Delimiter is the 5. Input: Value = HEAD2VAL( ... , Delimiter);
%
% see also: STR2VAL
%

msg = '';
str = '';

Nin = nargin;

if Nin < 1
   val = [];
   str = '';
   return
end

if Nin < 2
   frm = '';
end

if Nin < 3
   len = 0;
elseif ~( mod(len,1) == 0 )
   error('Length must be a integer.')
end

if Nin < 4
   vnan = [];
end

if Nin < 5
   dlm = [];
end

if isempty(frm)
   frm = '%s';
end

acc  = 1e-10;  % Accuracy

%**********************************************
switch frm

%**********************************************
case '%s'

   if ischar(val)
      val = rmblank(val,2);
      str = val;
   else
      [val,str,frm,msg] = str2val(val,'char');
   end

%**********************************************
case { '%d' '%g' '%f' '%e' 'data' }

   if len == 0
      len = inf;
   end

   uselen = isfinite(len);

   siz = [ 1  len ];

   if isempty(val)
      if uselen
         val = NaN * ones(siz);
      else
         return
      end
   end

   switch frm
     case '%d'
          frm = '%.0f';
     case '%e'
          frm = '%.6e';
     otherwise
          frm = '%.7g';
   end

   [val,str,frm,msg] = str2val(val,frm,vnan,acc);

   if uselen & ~isequal(size(val),siz)
      if isempty(val)
         val = NaN*ones(siz);
         str = '';
      else
         val = val(:)';
         szv = size(val,2);
         if szv > len
            val = val(1:len);
         elseif szv < len
            val = cat( 2 , val , NaN*ones(1,len-szv) )
         end
      end
   end

%**********************************************
case 'date'

  %%% NaN if empty ???

   form = 'date0D';   % YYYY/MM/DD

  [val,str,frm,msg] = str2val(val,sprintf('date0D',frm),vnan,acc);

  % Decimal Day if NonZeros Hours !!!

   val(len) = str2val(val(len:end),'time1');  

   val = val(1:len);

%**********************************************
case 'time'

  %%% NaN if empty ???

   form = 'time2';
   if ( len == 3 )
      form = 'time3';
   end

  [val,str,frm,msg] = str2val(val,form,vnan,acc);

  flv = floor(val);

  switch len

    case 2 % [ Hour  Minute ]

        val = [ flv  floor(60*(val-flv)) ];

    case 3  % [ Hour Minute Second ]

        mnt = 60*(val-flv);
        flm = floor(mnt);

        val = [ flv  flm floor(60*(mnt-flm)) ];

  end  

%**********************************************
case { 'lat' 'lon' 'pos' }        % Lat / Lon

   if isempty(val)
      val = NaN;
      return
   end

   [val,str,frm,msg] = str2val(val,sprintf('G%s2',frm),vnan,acc);

%**********************************************
case { 'list' 'cell' }
     % len == Seperator

   if isempty(dlm)
      if len == 0
         dlm = ':;,';
      else
         len = abs(len);
         len = len + 99 * ( len < 32 );  % !!!
         dlm = char(len);
      end
   end

   for c = dlm
       if any( str == c )
          dlm = c;
          break
       end
   end

   [val,str,frm,msg] = str2val(val,{dlm});

   if strcmp(frm,'unit')
      for ii = 1 : prod(size(val))
          val{ii} = rmblank(val{ii},2);
      end
      return
   end

   for ii = 1 : prod(size(val))
       val{ii} = splithead(val{ii});
   end

%**********************************************
otherwise

  if ( len < 0 )
     [val,str,msg] = head2val(val,'cell',len,vnan,dlm);
  else
      error(sprintf('Invalid Format "%s".',frm))
  end

end

