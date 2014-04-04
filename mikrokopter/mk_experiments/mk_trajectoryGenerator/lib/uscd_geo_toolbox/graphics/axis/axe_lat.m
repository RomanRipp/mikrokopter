function axe_lat(axe,xy)

% AXE_LAT   Nice Latitudal AxesLabels
%
% AXE_LAT( Axe , [ 'x' | 'y' ])


if nargin < 1
 axe = gca;
end
if nargin < 2
 xy = 'y';
end

if ischar(axe)
  if any( lower(axe) == 'xy' )
   xy = axe;
  end
  axe = gca;
end


 tick = get(axe,[ xy  'tick' ]);

 set(axe,[ xy  'tick' ],tick);


lab = 'SN';


   grad =   fix(tick);
   minu =   fix(  60*(tick-grad));
   secu = round(3600*(tick-grad-minu/60));

    minu_secu =  fix(secu/60);
   secu = secu - 60*minu_secu;
   minu = minu + minu_secu;

    grad_minu =  fix(minu/60);
   minu = minu - 60*grad_minu;
   grad = grad + grad_minu;

   dec_secu = [];
   sec_form = '';

  if any(diff(1e3*grad+minu)==0)
   dec_secu = round(1e2*( 60*(tick-grad) - minu)) ;
   sec_form = '.%2.2d';
  end
  if ~isempty(dec_secu)
   if any(diff(1e3*grad+minu+1e-3*dec_secu)==0)
    dec_secu = round(1e3*( 60*(tick-grad) - minu)) ;
    sec_form = '.%3.3d';
   end
  end

  is_dec = ~isempty(dec_secu);

  nt = prod(size(tick));

  str = cell(nt,1);
  str(:) = { '' };

  for ii = 1 : nt

     str{ii} = [sprintf('%3.0f',abs(grad(ii))) char(176) , ...
                sprintf('%2.2d',abs(minu(ii)))           , ...
                sprintf(sec_form,abs(dec_secu(ii*(1:is_dec)))) char(39)  , ...
                lab( 1 + ( tick(ii) > 0 ) )  ];
  end

  set( axe , [ xy 'ticklabel' ] , str )
