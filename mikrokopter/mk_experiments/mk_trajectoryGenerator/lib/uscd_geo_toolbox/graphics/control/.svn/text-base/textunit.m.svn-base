function org = textunit(h,uni)

% TEXTUNIT   Set Property UNITS of TextObjects
%
% TEXTUNIT preserve the Problem of PositionShift
%  while change from 'data' to normalized Units.
%
% Note: During set to Units 'pixels' a PositionShift
%       may be occur cause of round of PositionValues.
%
% OrgUnit = TEXTUNIT( TextHandle , NewUnit  )
% OrgUnit = TEXTUNIT( TextHandle , NewUnits )
%
% NewUnit(s) can be a CharacterArray or CellArray of Strings.
%
% To get the possible Values of the TextProperty Unit type:
%
%    >> display(set(0,'DefaultTextUnits'))
%
% The Number of Elements of TextHandles and NewUnits 
%  must be equal.
%

Nin  = nargin;
Nout = nargout;

org = cell(0,1);

if Nin < 1
   if Nout == 0, clear org, end
   return
end

%-----------------------------------------------------
% Check Inputs

msg = cell(0,1);

if ~isempty(h)

    h = h(:);
    n = size(h,1);

    ok = isnumeric(h);
    if ok
       ok = all(ishandle(h));
       if ok
          ok = all(strcmp(get(h,'type'),'text'));
       end
    end

    if ~ok
        msg = cat( 1 , msg , {'Invalid TextHandles.'} )
    else
        org = cell(n,1);
        for ii = 1 : n
            org{ii} = get(h(ii),'units');
        end
    end
end

if Nin == 2
   [ok,uni] = chkcstr(uni);
   if ok
      uni = uni(:);
      ok = ~any(strcmp(uni,''));
   end
   if ~ok
       msg = cat( 1 , msg , {'Units must be String(s).'} );
   elseif ~isempty(h)
       def = set(0,'defaulttextunits');
       for ii = 1 : size(uni,1)
           if isempty(strmatch(uni{ii},def))
               m = sprintf(' ''%s'' |',def{:});
               msg = cat( 1 , msg , {sprintf('Units must be any of:%s',m(1:(end-1)))} );
               break
           end
       end 
       if ~any( size(uni,1) == [ 1  n ] )
           msg = cat( 1 , msg , {'UnitNumber must be equal to HandleNumber.'} );
       elseif ( size(uni,1) == 1 ) & isempty(msg) 
              uni = uni(ones(1,n));
       end    
   end
end

if ~isempty(msg)
    error(sprintf('%s\n',msg{:}));
end

if ( Nin < 2 ) | isempty(h)
   if n == 1, org = org{1}; end
   if Nout == 0, clear org, end
   return
end

%-----------------------------------------------------
% Set Units

c = 'xyz';
m = size(c,2);

for ii = 1 : n

    
    if strcmp(org{ii},'data')

       % 'data' --> 'norm'

       a = get( h(ii) , 'parent' );

       p = get( h(ii) , 'position' );

       for jj = 1 : m

           l = get( a , [ c(jj)  'lim' ] );

           if strcmp( get( a , [ c(jj)  'dir' ] ) , 'reverse' )
              l = l([2 1]);
           end

           p(jj) = ( p(jj) - l(1) ) / ( l(2) -l(1) );

       end 

       set( h(ii) , 'units' , 'normalized' , ...
                 'position' , p                    )

    end

    set( h(ii) , 'units' , uni{ii} );

end

if n == 1, org = org{1}; end

if Nout == 0, clear org, end
