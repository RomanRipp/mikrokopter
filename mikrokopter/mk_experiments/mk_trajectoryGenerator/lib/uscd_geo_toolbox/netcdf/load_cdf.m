function [msg,dat,varargout] = load_cdf(file,varargin)

% LOAD_CDF  Read Data from NetCDF-File, using a ConfigFile
%
% [MSG, V ] = LOAD_CDF( CONFIG_FILENAME ,  Options )
%
% The ConfigFile, created by CNF_CDF, contains Informations 
%  about the NetCDF-FileName and  the Names for 4 Dimensions
%  and correspondending Variables you can specify 
%  for the NetCDF-File.
%  The 4 Dimensions are interpreted as X, Y, Z and Time.
%
% Options are:
%
%  To give the Range of the Dimensions (see NCMEX):
%
%   '#dim' , [ Start Count Stride ]  |  [ Stride ] ,
%
%   where # can be "x","y","z" or "t"
%
%
%  To give the Range of the Variables:
%
%    '#var' , V ,     
%
%   where # can be "x","y","z" or "t"  and  V is a Vector,
%   specified the Range.
%   A Range specified by Variables has priority over a Range
%   by Dimensions, if a DimensionRange is given additional to
%   a VariableRange, only the StrideValue is used to read
%   the Data from the NetCDF-File.
%
%  To give the Data to extract: 
%
%    'data' , < NAMES > ,
%
%   where NAMES is an Cell- or CharArray of Strings for the
%    Name of the KeyWord in the ConfigFile, describing the
%    Variable in the NetCDF-File.
%         
%   The Data, matching the Range, will be interpolated.
%   The MatlabDimension of the Data is [ Y X Z T ... ],
%    thats independed of the Dimension in the CDF-File.
%
%  To give a FillValue for Data matching
%   <missing_value> or <_FillValue> :
%
%   'fill', FillValue ,
%
%   where FillValue is the Value, to set Data with 
%    <missing_value> or <_FillValue>, default is NaN.
%   FillValue = [] meens no change of such values.
%
%  To REDEFINE the NetCDF-FileName, defined with the KeyWord 
%   "#FileName" in the ConfigFile :
%
%   'file', NetCDF_FileName , 
%
%  With the Option  'orig'  the Data will NOT interpolated,
%   to the VariableRange, you get them as Original.
%
% [MSG, V , DIM , VAR , ATT ] = LOAD_CDF( ... )
%
%   returns the Information about the Dimensions, Variables, 
%     Attributes from LOOK_CDF.
%
% See also:  LOOK_CDF  CNF_CDF  NCMEX  LOAD_ASC
%


Nin  = nargin;
Nout = nargout - 2;


nl = char(10);

msg = '';
dat = [];

varargout = cell(1,Nout);

msg0 = 'LOAD_CDF: ';
if Nin < 1
 msg = [ msg0 'Not enough InputArguments.']; 
 return
end 

%**************************************************************************
% Read InputArguments (varargin)

pre =  cellstr(['xyzt']');
suf = { 'dim' ; 'var' };

N = size(pre,1);

[msg,cells,var_req,zz,is_origin,miss_fill,ncfile] = checkin(varargin,pre,suf);

if ~( ischar(file) & ~isempty(file) & ( prod(size(file)) == size(file,2) ) )
    msg = cat(1,{'CONFIG_FILENAME must be a nonempty String.'},msg); 
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('%s Invalid Inputs.\n%s',msg0,msg);
    return
end
 
%**************************************************************************
% Check CONFIGFILE  &   CDF-Inquire

[msg,cnf,dim,var,att,txt,look_msg] = cnf_cdf(file,'check',ncfile);

%  cnf = { CDF_File DimName DimLength VarName }

if ~isempty(msg)

    [m,cnf,dim,var,att,txt,look_msg] = cnf_cdf([file '.cnf'],'check');

    if ~isempty(m)
        msg = sprintf('%s Invalid ConfigFile "%s".\n%s',msg0,file,msg);
        return
    end

end

%----------------------------------------------------------------

out = { cnf dim var att txt };

n = min(Nout,size(out,2));

varargout(1:n) = out(1:n);

msg = cell(0,1);

%**************************************************************************
%----------------------------------------------------------------
% CDF-Inquire allready done by "look_cdf(ncfile)" in  cnf_cdf(file,'check')

ncfile = cnf{1,1};

ndim = size(dim,1);
nvar = size(var,1);

%----------------------------------------------------------------

dims = cells(:,1);
vars = cells(:,2);

if ~( zz == 0 )
    var_req = cellstr(char(var_req(1:zz)));
else
    var_req = cell(0,0);  
end

% DimensionDefault for all Dimensions: [ start  count  stride ] 
dimd      =  num2cell(ones(ndim,1)*[ 0  1  1  ],2);  

%---------------------------------------------------------------------
% Get Numbers for requested Dimensions and DimensionVariables

dim_nr = NaN * ones(N,1);
var_nr = NaN * ones(N,1);

is_interp = 0;

for ii = 1 : N

 % Interpolate Data only if Variable requested
 is_interp = ( is_interp | (~isempty(vars{ii})) );


 if ~isempty(cnf{ii,2})
     dim_nr(ii) = find(strcmp(cnf{ii,2},dim(:,1))) - 1;
     % Default for the BasicDimensions: Full Size
     dimd{dim_nr(ii)+1,1} = [0 dim{dim_nr(ii)+1,2} 1];    
 elseif dims{ii}{2}
     % Dimension in VarArgin requested
     msg = cat(1,msg,{sprintf('Dimension %s%s not defined in "%s".', ...
                              upper(pre{ii}),upper(suf{1}),file)});
     dims{ii}{2} = 0;
 end

 if ~isempty(cnf{ii,4})
     var_nr(ii) = find(strcmp(cnf{ii,4},var(:,1))) - 1;
 elseif ~isempty(vars{ii})
     % Variable in VarArgin requested
     msg = cat(1,msg,{sprintf('DimensionVariable %s%s not defined in "%s".', ...
                              upper(pre{ii}),upper(suf{1}),file)});
     vars{ii} = [];
 end

end


%---------------------------------------------------------------------
% Get Numbers for requested DataVariables

dat_nr = [];

if ~isempty(var_req)  &  isempty(cnf{2,1})

    msg = cat(1,msg,{sprintf('No other Variables defined in %s.',file)});

elseif isempty(var_req) & ~isempty(cnf{2,1})
   
   var_req = cell(size(cnf{2,1},1),1);
   var_req(:) = {''};
   
end

if ~isempty(var_req)

     zz = size(var_req,1);

 dat_nr = NaN*ones(zz,1);

 for ii = 1:zz

  val = var_req{ii};
  if ~isempty(val)


    jj = find( strcmp( val , cnf{2,1}(:,1) ) );
    if isempty(jj)
        msg = cat(1,msg,{sprintf('Requested DataVariable %s is not defined in %s.', ...
                                 val,file)});
    else
      nn = find( strcmp( cnf{2,1}{jj(1),2} , var(:,1) ) );
      if isempty(nn)
         msg = cat(1,msg,{sprintf('Requested CDF-Variable %s for DataVariable %s doesn''t exist in %s', ...
                           cnf{2,1}{jj(1),2},val,ncfile)});
      else
         dat_nr(ii) = nn(1) - 1 ;
      end
    end

  end
  % ~isempty(val)

 end
 % ii

end


if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('%sInvalid request in Inputs.\n%s',msg0,msg);
    return
end

msg = '';

%-----------------------------------------------------------------------

[fid,stat] = ncmex('open',ncfile,'nowrite');

if ~( ( fid > 0 ) & ( stat == 0 ) )
   ncfile1 = which(ncfile);
   [fid,stat] = ncmex('open',ncfile1,'nowrite');
   if ~( ( fid > 0 ) & ( stat == 0 ) )
      msg = sprintf('%sError open File %s as NetCDF-File.',msg0,ncfile);
      return
   end
end


%-----------------------------------------------------------------------
% Load DimVariables and calculate [ start count stride ]
%

val01 = zeros*ones(N,2);

is_flip   = zeros(N,1);   % Set, if FlipDimension for
                          %       monotonicly Increasing
is_transf = zeros(N,1);   % Set, if LongitudeTransformation
is_intv   = zeros(N,1);


for ii = 1:N

 if ~isempty(vars{ii})
   % DimensionVariable requested

   val = vars{ii}(:);
   
   did = dim_nr(ii);
   vid = var_nr(ii);  
   len = dim{did+1,2};

    vv = ncmex('vargetg',fid,var_nr(ii), 0 , len , 1 )';

  vv01 = [ min(vv) max(vv) ]; 
   nvv = size(vv,1);

  val0 = val(:);  nv = size(val,1);

  if nv ~= 1
    is_flip(ii) =  any( diff(val0) < 0 );
  
    if is_flip(ii)
     val0 = flipud(val0);  
    end

    if any( diff(val0) < 0 ) 
     msg = [' LOAD_CDF : DimensionVariables must be monotonic. '];
     ncmex('close',fid);
     return
    end
  end

  val01(ii,:) = [ val0(1)  val0(nv) ];

  vars{ii} = val0;

  is_transf(ii) = ( ii == 1 ); 

%**********************************************************************
% Transformation
%
% SPECIAL PROCEDURE FOR X (LONGITUDE)
%  ii == 1
%                                                    floor(min(vv)/180)
%  Possibilities: [ -540      -180                             -3
%                        -360         0                        -2
%                             -180        180                  -1
%                                     0        360              0
%                                         180       540   )     1
% 
%
%

  is_transf(ii) = ( ii == 1 ); 

   i2 = 180;  % half Intervall, OverLap


if is_transf(ii)
%-------------------------------------------------

   is_transf(ii) = ( ( vv01(2)- vv01(1)) < 2*i2  & ...
                     (val01(ii,2)-val01(ii,1)) <= 2*i2         );



    % Only, if DimensionVariable in Range 


   % Begin and End Connected ?

   is_conect = ( abs(diff(vv01)-2*i2) <= 2*mean(diff(vv)) );

 
   ff      = floor(vv01(1) / i2 );           % Intervall   [ ff   ff+2 ) * i2 
   ff      = ff - 1 * ( vv01(2) < (ff+1)*i2 );

   intv = [ ff   ff+2 ] * i2;

   is_intv(ii) = ( ( intv(1)<= vv01(1)  &  vv01(2) <  intv(2) ) | ...
                   ( intv(1)<  vv01(1)  &  vv01(2) <= intv(2) )      );

   is_intv(ii) = is_intv(ii) + is_intv(ii)*(vv01(2) == intv(2)); 

   % 1 == i2*[ ff  ff+2)  
   % 2 == i2*( ff  ff+2]


   is_transf(ii) = ( is_transf(ii) & is_intv(ii) );

end

if ~is_transf(ii)
%-------------------------------------------------
    val_new = val0;

    % StartIndex for Segment
      ind01 = [ 1 ; size(val_new,1)+1 ];

    is_conect = 0;

else
%-------------------------------------------------
 
   sum1 = - (2*i2) * floor(val0/(2*i2)) * is_transf(ii);
  
   val = val0 + sum1;  % [  val0   ]   --> [ 0 .. 360 )


   part = ( floor(val/i2) - 0 ) * ( is_intv(ii) == 1 ) + ...
          (  ceil(val/i2) - 1 ) * ( is_intv(ii) == 2 ) ;

   sum2 = (2*i2) *  ceil( (ff-part) / 2 ) *  is_transf(ii);

   val = val  + sum2;  % [ 0 .. 360 )  -->  [ vv ]



   % Search for Segments and fill SegmentSteps with Borders
 
   ind_seg = find( abs(diff(val)-diff(val0)) > 1e-10 ) + 1;

   nseg = length(ind_seg);

   ind_new          = ones(nv,1);
   ind_new(ind_seg) = 3;

   ind_new = cumsum(ind_new);
 
   val_new = NaN*ones(nv+2*nseg,1);
   val_new(ind_new) = val;
   
   ind_seg = [ind_seg;1] + 2*(1:nseg+1)'-1;   ind_seg(nseg+1) = [];

   val_new = val_new - ff*i2;   % --> [ 0  360 ]
 
   val_new(ind_seg-1) = 2*i2;  %  End   of Segment
   val_new(ind_seg  ) =    0;  %  Start of follw. Segment

   val_new = ( val_new + ff*i2 );


   % StartIndex for Segments
   ind01 = [ 1 ; ind_seg ; size(val_new,1)+1 ];    

end
%-------------------------------------------------



   nseg = size(ind01,1)-1;


   %  [ start count stride ]
   %  Stride specified by User or Default 
   dims{ii}{1} = ones(2*nseg,1) * dims{ii}{1} ;

   str = dims{ii}{1}(1,3);   % Stride

   zz = 0;   % Counter for the multiple Coords


   for jj = 1:nseg

      fak = sign(mean(diff(vv)));

     vn01 =  val_new([ ind01(jj)   ind01(jj+1)-1 ]);

     n1 = find( ( vn01(1) < vv )   &   ( vv < vn01(2) ) );

     n = [];
     if     ~isempty(n1)
         n = [ n1(1)  n1(end) ];
         n = n  +  [ -1  1 ] .* ( ([1 -1].*n) > [1 -nvv] );    
     else

        n2 = find( fak*vv <= fak*vn01(1) ); n2 = max(n2);
        n3 = find( fak*vv >= fak*vn01(2) ); n3 = min(n3);

        ff = [ 0  0 ];
        ff = ff( 1 : (2-(~isempty(n2))-(~isempty(n3))) );

        n = [ n2 n2 n3 n3 ];
        n = n( ( 1 : ( 2 * ((~isempty(n2) | ~isempty(n3))) )) + ...
                ((~isempty(n2)) & (~isempty(n3))) );

     end
        
     if ~isempty(n)
       
       if is_transf(ii)    &    is_conect    &    ...
          (   jj == 1 )    &   ~isempty(n1)  &    ...
          ( n(1) == 1 )    &   ( vn01(1) < vv(n(1+(fak==(-1)))) )
         % First Segment found one at beginning, add the End
           nn = nvv+1 - str*(nvv+1>str) - nvv*(nvv+1<=str);
           zz = zz+1;
           dims{ii}{1}(zz,1:3) = [ nn-1  1  1 ];       
       end

       zz = zz+1;

       % Index in CDF-Dimension 

        cc0 = n(2)-n(1)+1;       % Number of Counts for stride == 1

        cc  = floor(cc0/str);    % Counts ( length((n(1):str:n(2))) )
               

       % ende = start  + stride*(count-1)
           nn =  n(1)  +  str  *( cc - 1) ;

        cc = cc + 1*( ( nn+str <= nvv )  &  ( nn < n(2) ) );

        dims{ii}{1}(zz,1:2) = [ n(1)-1  cc ];

        n(2) = nn;

       if is_transf(ii)      &   is_conect    &    ...
          (   jj == nseg )   &  ~isempty(n1)  &    ...
          ( n(2) == nvv  )   &   ( vv(n(2-(fak==(-1)))) < vn01(2) )   
          % Last Segment found one at the End, add the Begin
            nn = 0 + str*(nvv>str) + nvv*(nvv<=str);
            zz = zz+1;
            dims{ii}{1}(zz,1:3) = [ nn-1  1  1 ];       
       end

     end

   end
   % jj 

      dims{ii}{1} = dims{ii}{1}(1:zz,:);

    if zz == 0   % isempty( dims{ii}{1} )
         msg = [ msg  msg0 ' DimensionVariable '  upper([pre{ii} suf{2}]) ...
                           ' out of range.' nl ]; 
       dims{ii}{1} = [ 0  1  1 ];
       vars{ii}    = [];
    end



       dimd(dim_nr(ii)+1) = dims{ii}(1);


    clear vv val_new val val0 ind_seg ind01



 
 elseif dims{ii}{2}
   % Dimension requested
   if isnan( dims{ii}{1}(2) )
     dims{ii}{1}(2) = floor( ( (dimd{dim_nr(ii)+1}(2)-1) - dims{ii}{1}(1)) / ... 
                             dims{ii}{1}(3) ) + 1;
     % count        = floor( (  MaxCount-1 - start )  /  stride )+1;
   end

    nn = find( isnan(dims{ii}{1}) );
    dims{ii}{1}(nn) = dimd{dim_nr(ii)+1}(nn);

    end_count = dims{ii}{1}(1) + (dims{ii}{1}(2)-1) * dims{ii}{1}(3);
                %  start       + ( count        -1) * stride

    
    if end_count > dim{ dim_nr(ii)+1 , 2 }
         msg = [ msg  msg0  ' Dimension '  upper([pre{ii} suf{1}]) ...
                            ' out of range.' nl ]; 
       dims{ii}{1} = [ 0  1  1 ];
    end
    
    dimd(dim_nr(ii)+1) = dims{ii}(1);

 end
end
% ii



% Multiple Counts of X-Dimension (dim_nr(1))
c_si  = size(dimd{dim_nr(1)+1},1);      % CountSize  of X-Dimension 


dim_si = zeros(1,ndim);
for ii = 1:ndim
 dim_si(ii) = sum(dimd{ii}(:,2));
end


%****************************************************
% Read DimensionVariables

ndat  = length(dat_nr);

dat = cell(N+ndat,1);

for ii = 1:N
 if ~isnan(var_nr(ii))
 
   dat{ii} = NaN * ones(dim_si(dim_nr(ii)+1),1);

   Nint = 1+(c_si-1)*(ii==1);

    zz0 = zeros( Nint+1 , 1 );
  
   for jj = 1 : Nint
    
     % STartCountSTRide  

     stcstr = dimd{dim_nr(ii)+1}(jj,:);

     val = ncmex('vargetg',fid,var_nr(ii), ...
                  stcstr(1),stcstr(2),stcstr(3));
     val = val(:);

             ind          =  ( 1 : size(val,1) )';
     dat{ii}(ind+zz0(jj)) =  val;

     zz0(jj+1) = zz0(jj) + max(ind);

   end


   if(  ii == 1 )  &  is_transf(ii)
   % reverse LongitudeTransformation

     sum1 = - (2*i2) * floor(dat{ii}/(2*i2));
     dat1 = dat{ii} + sum1;  % [  dat   ]   --> [ 0 .. 360 )

     ffv1 = floor(val01(1,1) / i2 );

     ffv = [ ffv1 ; 
             ffv1 - 1*( val01(1,2) <  (ffv1(ii,1)+1)*i2 ) ; 
             ffv1 + 1*( val01(1,1) >= (ffv1(ii,1)+1)*i2 )    ];

    ok = 0;
    jj = 0;

    while ~ok   &  jj < length(ffv) 

     jj = jj + 1;
   
     % is_intv(ii) == 1  ==>  i2*[ ff  ff+2)  
     % is_intv(ii) == 2  ==>  i2*( ff  ff+2]

     part = ( floor(dat1/i2) - 0 ) * ( is_intv(ii) == 1 ) + ...
            (  ceil(dat1/i2) - 1 ) * ( is_intv(ii) == 2 ) ;

     sum2 =   (2*i2) * ceil( (ffv(jj)-part) / 2 );
     dat2 = dat1  + sum2;  % [ 0 .. 360 )  -->  [ origin ]

     % dat2 should be monotonicly increased !!!
     %  we could get Problems at Begin and End  

     n = length(dat2);
     
     if  n > 2
       dat2(1) = dat2(1) - 2*i2 * ( dat2(2) <= dat2(1) );
       dat2(n) = dat2(n) + 2*i2 * ( dat2(n) <= dat2(n-1) );
      
       ok = xor( any( diff(dat2) <= 0 ) , any( diff(dat2) >= 0 ) );  
     else
       ok = 1;
     end
   
    end

     if ok
      dat{ii} = dat2;
     end
     
   end 

  if length(dat{ii}) > 1
   if ~xor( any(diff(dat{ii})<=0) , any(diff(dat{ii})>=0 ) )   &  ~is_origin
     is_interp = 0;
     msg = [ msg  msg0 'DimensionVariable for ' upper(pre{ii}) ...
                       ' not monoton, No Interpolation.' nl ]; 
   end
  end

  
 end

end



% Prepare X,Y,Z,T for Interpolate 

if ~is_origin  &  is_interp

   XYZ0 = cell(ndim,2);
   for ii = 1:ndim
%    XYZ0(ii,:) = {(1:dim_si(dim_nr(ii)+1))'}; % !!!!!!!!!!!!!
     XYZ0(ii,:) = {(1:dim_si(ii))'};
   end
   for ii = 1:N
    if ~isempty(dat{ii})
     XYZ0(ii,:) = { dat{ii}  dat{ii} };
    end 
    if ~isempty(vars{ii})
     XYZ0(ii,2) = { vars{ii} } ;
     if ~is_origin
      dat{ii} = vars{ii};
     end
    end
   end

 % Permutation: [ X Y Z T ] !!!!!!!


 %  [ ... k j i ] = ncmex("get", ... , [ i j k ... ] );
 %  XYZ = XYZ([ndim:-1:1],:);  % !!!!!!!!!!!!!!!!!!!!

end


%****************************************************
% Read DataVariables


for ii = 1:ndat
 
 if ~isnan(dat_nr(ii))

   vid = dat_nr(ii);

   dim_ii = var{vid+1,4}+1;

   dim_n  = var{vid+1,3};

   dim_perm = (N+1) + 0*dim_ii;

   for dd = 1:N
    if ~isnan(dim_nr(dd))
     dim_perm( find(dim_ii==dim_nr(dd)+1) ) = dd;  % later: + 1*(dd==1) - 1*(dd == 2); 
    end
   end

   %  [ ... k j i ] = ncmex("get", ... , [ i j k ... ] );
 
    dim_perm = fliplr(dim_perm);   % !!!!!!!!!!!!!!!!!

 
   [dim_org,dim_perm] = sort(dim_perm); % [ X Y Z T ... ]


   dim_len0 = dim_si(dim_ii); 
   dim_len0 = fliplr(dim_len0);     % "fliplr" !!!!!!!!!!
   dim_len  = dim_len0(dim_perm);       % [ X Y Z T ... ]

   % Define DataVariable
   dat{N+ii} = NaN * ones( dim_len );


   % [ Start Count Stride ]  without X - Dimensions

   dimx_ii = dim_nr(1)+1;
   
   stcstr1 = dimd(1:dimx_ii-1) ;
   if isempty(stcstr1)
      stcstr1 = zeros(3,0);
   else
      stcstr1 = cat(1,stcstr1{:} )';
   end

   stcstr2 = dimd(dimx_ii+1:ndim) ;
   if isempty(stcstr2)
      stcstr2 = zeros(3,0);
   else
      stcstr2 = cat(1,stcstr2{:} )';
   end

 
   zz = 0;  % IndexCounter for X-Dimension (1.) in DataVariable

   
   for jj = 1 : 1+(c_si-1)*(any(dim_ii==1))
    
     % STartCountSTRide
     stcstr = [ stcstr1  dimd{dimx_ii}(jj,:)'  stcstr2 ];

     val = ncmex('vargetg',fid,vid, ...
                  stcstr(1,dim_ii),stcstr(2,dim_ii),stcstr(3,dim_ii));
     val = permute(val,dim_perm);

               ind       =  ( 1 : size(val,1) );
     dat{N+ii}(ind+zz,:) =  val(ind,:);      clear val

     zz = zz + max(ind);

   end



   % Change X-Y-... to Y-X-...

   dim_perm1 = [ 2  1  3:dim_n];  % [ Y  X  Z  T ... ]

   dat{N+ii} = permute(dat{N+ii},dim_perm1(1:dim_n));   



   is_int = 0;  % True for UINT8


   % Fill Missing Values
  
   if ~isempty(att{vid+1})

     is_miss = [];
     if ~isempty(miss_fill)
         is_miss = find( ( strcmp( att{vid+1}(:,1) , 'missing_value' ) | ...
                           strcmp( att{vid+1}(:,1) , '_FillValue' )  ) );
         if ~isempty(is_miss)
            if prod(size(att{vid+1}{is_miss(1),3})) == 1
               miss_val = att{vid+1}{is_miss(1),3};
               is_miss  = find( dat{N+ii} == miss_val );
            end
         end
     end

        %--------------------------------------------------------------
        % Scale and Add Offset

        jj = strcmp( att{vid+1}(:,1) , 'scale_factor' );
        scl = 1;
        if any(jj)
           jj = find(jj);
           scl = att{vid+1}{jj,3};
           if ( prod(size(scl)) == 1 ) & ~isequal(scl,1)
              dat{N+ii} = dat{N+ii} * scl;
%              if is_char(vid+1)
%                 val = char(val);
%              end
           end
        end

        jj = strcmp( att{vid+1}(:,1) , 'add_offset' );
        off = 0;
        if any(jj)
           jj = find(jj);
           off = att{vid+1}{jj,3};
           if ( prod(size(off)) == 1 ) & ~isequal(off,0)
              dat{N+ii} = dat{N+ii} + off;
%              if is_char(vid+1)
%                 val = char(val);
%              end
           end
        end

        is_int = ( strcmp(var{vid+1,2},'byte') & ...
                   isequal(scl,1) & isequal(off,128)  );

        %--------------------------------------------------------------

      if ~isempty(is_miss)
          dat{N+ii}(is_miss) = miss_fill;
      end

  end


   % Interpolate 

   if ~is_origin  &  is_interp


     dim_perm2 = dim_perm(dim_perm1);  % NETCDF->Matlab --> [ Y X Z T ... ]

     dim_len   =  dim_len(dim_perm1);     %  [ X Y Z T ... ] -->
                                        %  [ Y X Z T ... ]
       XYZ     =      XYZ0(dim_perm1,:);


      % Interpolate only over nonsingleton Dimensions

       kk1 = find(dim_len ~= 1);
       kk2 = find(dim_len == 1);  

       if ~isempty(kk1)
        XYZI = XYZ(kk1,2);
        for jj = 1:length(kk1)
         XYZI(jj) = { permute( XYZI{jj} , ( jj : -1+2*(jj==1) : 1+(jj==1) ) ) };
        end
        if length(kk1) > 1
         val = interpm(XYZ{kk1,1}, permute(dat{N+ii},[kk1 kk2]) ,XYZI{:});
        else
         val = interp1(XYZ{kk1,1}, permute(dat{N+ii},[kk1 kk2]) ,XYZI{:});
        end

        [hilf,kk] = sort([kk1 kk2]);
        dat{N+ii} = permute(val,kk);      

       end
       
   end    

   if is_int
      dat{N+ii} = uint8(dat{N+ii});
   end

 end
 % ~isnan(dat_nr(ii))

end


for ii = 1:N
 if is_flip(ii)
  
  dat{ii} = flipud(dat{ii});

  for jj = 1:ndat
   dat{jj+N} = flipdim( dat{jj+N} , ii+1*(ii==1)-1*(ii==2) );
  end

 end
end


% close NetCDF file
%

ncmex('close',fid);

%*******************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cells,var_req,zz,is_origin,miss_fill,ncfile] = checkin(v,pre,suf)

nv = prod(size(v));

msg = cell(0,1);

%----------------------------------------------------------------
% Read InputArguments (varargin)

N = size(pre,1);
M = size(suf,1);

cells = cell(N,2);

% DimensionDefault: { [ start  count  stride ]   VarArg } 
cells(:,1) = {{ [0 NaN 1]  0 }};  

   var_req = cell(nv,1);  % Requested DataVariables: 'NAME' or [Number]
        zz = 0;           % Counter for var_req;

miss_fill = NaN;
is_origin = 0;
ncfile    = '';

for ii = 1 : nv

    str = v{ii};

    if ischar(str)  &  ( size(str,1) == 1 )

       str = lower(deblank(str));

       % Look vor Dim or Var

       ok = 0;

       for jj = 1 : N
           for kk = 1 : M   
               if strcmp( str , [ pre{jj} suf{kk} ] ) & ( ii < nv )
                  vv = v{ii+1};
                  if isnumeric(vv) & ~isempty(vv)   
                     vv = vv(:)';
                     if strcmp(suf{kk},'dim')
                        % Dimension [ start  count  stride ]
                        mv = min(size(vv,2),3);        
                        vv = round( vv(1:mv) );        
                        nn = find( ~( isnan(vv)  |  ( vv <= 0 ) |  ~finite(vv) ) );   
                        % Single Value ==> Set Stride !!!
                        cells{jj,kk}{1}(nn+2*(mv==1)) = vv(nn);
                        cells{jj,kk}{2}     = 1;
                     elseif strcmp(suf{kk},'var')
                        % DimensionVariable
                        cells(jj,kk) = { vv };
                     end
                     ok = 1;
                  else
                     msg = cat(1,msg,{sprintf('Invalid Data for %s%s.',upper(pre{jj}),upper(suf{kk}))});
                  end
                  %  isnumeric  & ~isempty
               end
               % strcmp
           end
           % kk
       end
       % jj 

       if ~ok
           % An DataVariable required
           if strcmp(str,'data') 
              if ii < nv 
                 vv = v{ii+1};
                 if ischar(vv) | iscellstr(vv);     
                    zz  = zz+1; 
                    var_req{zz} = char( vv );
                 else
                    msg = cat(1,msg,{'Values for DATA must be a CharArray or CellArray of Strings.'});
                 end
                 % ischar/iscellstr   
              end
              % ii < nv
           elseif strcmp(str,'orig')
              is_origin = 1;
           elseif strcmp(str,'fill')
              if ii < nv
                 vv = v{ii+1};
                 if isempty(vv)
                    miss_fill = [];
                 elseif isnumeric(vv) & ( prod(size(vv)) == 1 )
                    miss_fill = vv;
                 else
                    msg = cat(1,msg,{'Value for FILL must be EMPTY or a single Numeric.'});
                 end
              end
           elseif strcmp(str,'file')
              if ii < nv
                 vv = v{ii+1};
                 if ~isempty(vv)
                     if ischar(vv) & ( prod(size(vv)) == size(vv,2) )
                        ncfile = vv;
                     else
                        msg = cat(1,msg,{'Value for FILE must be a String.'});
                     end
                 end
              end            
           end  % 'data' | 'orig' | 'fill'
       end % ~ok

    end
    % ischar(str)


end
% ii


