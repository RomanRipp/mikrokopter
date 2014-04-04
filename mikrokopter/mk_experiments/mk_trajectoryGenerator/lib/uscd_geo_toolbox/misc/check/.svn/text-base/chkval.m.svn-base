function [val,str,form,msg] = chkval(val,typ,si,form,opt,rec);

% CHKVAL   Check Value with Type, Size and Format
%
% [val,str,form,msg] = CHKVAL(val,typ,si,form,opt,rec);
%
% calls:  str2val.m
%          defval.m
%        classchk.m
%
% typ: 'numeric'
%      'cellstr'
%      'string'
%      'object'
%      <class>
%      ''
%

Nin = nargin;

if Nin < 4
   msg = 'Not enough InputArguments.';
   return
end

if Nin < 5
   opt = cell(0,2);
end

if Nin < 6
   rec = 1;   % Recurse
end

msg = '';
str = '';

rec = isequal( rec , 1 );

%****************************************************
% Check with Optional Values

if ~isempty(opt)

   [val,str,ff,msg] = str2val(val,opt);

   if isempty(msg);
      form = ff;
      return
   end

end

%****************************************************

is_obj = ( ischar(form)  &  ~isempty(form) );
if is_obj
   is_obj = strcmp(form(1),'@');
end

if ~is_obj
    typ = lower(typ);
end

if ischar(val) & strcmp(typ,'cellstr')
   val = cellstr(val);
end


%-------------------------------------
% Convert with Format

if is_obj & strcmp( class(val) , form(2:(end-strcmp(form(end),'@'))) )

   [val,str,ff,msg] = str2val(val,'@');

else

   [val,str,form,msg] = str2val(val,form);

end

if ~isempty(msg)     % Invalid Value
    return
end


%------------------------------------------------
% Check for Empty Value

if isempty(val); 

   if ~isempty(opt)

      [val,str,ff,msg] = str2val(opt{1,1},opt(1,:));

      if isempty(msg)
         form = ff;
         return
      end

   end

   [val,msg] = defval(typ,si,val);

   [val,str,form,msg] = str2val(val,form);

   return
      
end

%****************************************************
% Check Type

cl = class(val);

[ok,val] = classchk(val,typ);

if ~ok
    msg = sprintf( 'Value must be of Type "%s".' , typ );
    return
end

%****************************************************
% Check for changed Type, 
%  call CHKVAL again, without Recurse !!!

if  ~strcmp( class(val) , val ) & rec
     [val,str,form,msg] = chkval(val,typ,si,form,opt,0);  
end

%****************************************************
% Check Size

 nd = size(si,2);

 ok = ( ndims(val) == nd );

 if ok

    ind = cell(nd,1);
       
    for jj = 1 : nd
                  
       if isnan(si(jj))
          ind{jj} = ( 1 : size(val,jj) );
       else
          ind{jj} = ( 1 : si(jj) );
       end

    end  

    for jj = 1 : nd

       if isnan(si(jj))
          ok = ( ok & ( size(val,jj) >= 1 ) );
       else
          ok = ( ok & ( size(val,jj) == si(jj) ) );
       end

    end

    if ok
       val = val(ind{:});
    end

 end

 if ~ok

    si = sprintf('%.0f x ',si);
    si = si( 1 : end-3 );

    msg = sprintf( 'Value must have Size: [ %s ].' , si );

 end

