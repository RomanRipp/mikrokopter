function str = split_str(ht,str,ww,offs)

%  SPLIT_STR  Split a Strings into multiple Lines
%
%  String1 = SPLIT_STR( TextHandle , String0 , MaximumWidth , OffsetString )
%
%  Split String0 into Strings with maximum Width of MaximumWidth,
%   measured in Units of TextHandle, at Characters like  ' ' or  ','  or  ';' .
%   using the TextHandle to determine the Width.
%
%   String0   CharArray (String)
%   String1   CellStringArray  (Splitted String)
%

str = { str };

ok = 0;

while ~ok

   str1 = str{end};

   % Remove Right Blanks
        ii   = find( abs(str1) == 32 );
%        jj1  = find( ii == ( 1 : length(ii) ) );
        jj2  = find( ii == ( -length(ii)+1 : 0 )+size(str1,2) );
   str1(jj2) = [];
%   str1(jj1) = [];
 
   if isempty(str1)
     ok = 1;
     str(end) = [];
   else

     set(ht,'string',str1);
     ext  = get(ht,'extent');

     ok = ( ext(3) < ww );

     if ~ok

      % Locations to Split
       ss = find( ( str1 == ' ' )  |  ...
                  ( str1 == ',' )  |  ...
                  ( str1 == ';' )         );

       jj = find( ss <= size(str1,2)*ww/ext(3) );

       if ~isempty(jj)
          jj = jj( find( ss(jj) > 1 ) );
       end

       if isempty(jj)
          ss = ( 1 : size(str1,2) );
          jj = find( ss < size(str1,2)*ww/ext(3) );
       else
          ss = [ ( 1 : ss(jj(1))-1 )  ss ];
          jj = ( 1 : size(ss,2) );
       end

       kk  = size(jj,2)+1;
       ok1 = 0;
       while ~ok1  &  ( kk > 1 )
         kk = ( kk - 1 ); 
         set(ht,'string',str1(1:ss(jj(kk))));
         ext = get(ht,'extent');
         ok1 = ( ext(3) < ww );
       end

       str2 = [ offs  str1( jj(kk)+1 : end ) ];             
       str1 = str1(1:jj(kk));

       str = [ str(1:end-1) ;  { str1 }  ;  { str2 }  ];
       str = str( 1 : end-isempty(str2) );
  
     end 
     % ~ok

   end

end
% while   