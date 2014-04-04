function [dat,s_i,out] = sort_grp(dat,sort_col,ang_col,grp_col,out_col,chk);

% SORT_GRP  Sort a Matrix by Rows and in Order of specified Columns  
%
% [dat,s_i,out] = sort_grp(dat,sort_col,ang_col,grp_col,out_col,chk);
%
%  Sort "dat" via Columns of "sort_col"
%  
%  Transforms Angles of Columns of "ang_col" to:
%
%      [    0 .. 360 )  for Sort
%      ( -180 .. 180 ]  for OutPut
%
%  Extract Groups of "dat" at same Values of Columns "grp_col"
%  Returns Groups in "out" at Columns of "out_col"  
%
%  Check for same Values in Group at Column "chk"
%


Nin  = nargin;
Nout = nargout;

s_i = [];
out = [];

if isempty(dat)
   return
end

if Nin < 2
   sort_col = ( 1 : size(dat,2) );
end

if Nin < 3
   ang_col = [];
end

if Nin < 4
   grp_col = [];
end


if Nin < 5
    out_col = grp_col;
end

if Nin < 6
    chk = [];
end


%***********************************
% Sort via Winkel & Radius

 %---------------------------------
 % Transform Winkel --> [ 0 .. 360 )

 if ~isempty(ang_col)
   dat(:,ang_col) = dat(:,ang_col) - 360 * floor( dat(:,ang_col) / 360 );
 end

 %---------------------------------
 % Sort

 n = size(dat,1);

 s_i = ( 1 : n )';

 for ii = sort_col(end:-1:1)
    [h,ind] = sort(dat(s_i,ii));
       s_i  = s_i(ind);
 end

 %---------------------------------
 % Transform Winkel --> ( -180 .. 180 )

 if ~isempty(ang_col)
   dat(:,ang_col) = dat(:,ang_col) - 360 * ceil( ( dat(:,ang_col) - 180 ) / 360 );
 end


 dat = dat(s_i,:);

 if isempty( grp_col )  | isempty(out_col)

   return

 end

%***********************************



 out = dat(:,grp_col);  % Planet:  Radius  Winkel  


%***********************************
% Groups of Planet

 ok = ones( n , 1 );

 bad = (  sum( ( diff(out,1,1) == 0 ) , 2 ) == size(out,2) );

 ok( find(bad) + 1 ) = 0;


 is_grp = find(ok);  % Start of  Group

 grp = cumsum(ok);   % GroupIndex for s_i --> same Planet

 %-----------------------------------------------------
 % Check for NaN-Values in chk-Column
 % Same Value for each Group

 if ~isempty(chk)

   bad = isnan(dat(:,chk));

   if all( bad )

      dat(:,chk) = 0;

   else
      
     md = median( dat( find(~bad) , chk ) );

     for ii = 1 : max( grp );

       jj = find( grp == ii );

       kk = find( ~isnan( dat(jj,chk) ) );

       if isempty(kk)
          dat(jj,chk) = md;
       else
          dat(jj,chk) = median( dat(jj(kk),chk) );
       end

     end

   end

 end

 out = dat( is_grp , out_col );
