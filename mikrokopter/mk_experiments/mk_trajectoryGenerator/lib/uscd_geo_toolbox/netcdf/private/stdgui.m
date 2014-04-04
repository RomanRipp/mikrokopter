function s = stdgui(varargin);

% STDGUI returns DefaultStructure for single GUI-Item
%
%        Type: ''
%       Width: NaN
%       Color: [1 1 1]
%      String: {''}
%       Value: [1 0 1 0.0100 0.1000]
%       Style: ''
%        Text: ''
%    UserData: []
%       CBFcn: ''
%      Option: []


Nin = nargin;

%       Value  Min  Max  SliderStep
val = [  1     0    1    0.01  0.1  ];

 cl = [ 1  1  1 ]; 


s = struct( 'Type'   , {'text' } , ...
            'Width'  , { NaN } , ...  % ButtoWidth in Character
            'Color'  , { cl  } , ...  % BackGroundColor
            'String' , { {''}} , ...  
            'Value'  , { val } , ...
            'Style'  , { ''  } , ...  % Stlye used for WIdth 
            'Text'   , { ''  } , ...  % TextLabel above
          'UserData' , { []  } , ...
            'CBFcn'  , { ''  } , ...
            'Option' , { []  }       );


Nin = 2*floor(Nin/2);
ind = ( 1 : 2 : Nin );

if ~( ( Nin >= 2 ) & iscellstr( varargin(ind) ) )
   return
end

fs = fieldnames( s );

lv = lower( varargin(ind) );


for ii = 1 : length(fs)

    jj = find( strcmp( lower(fs{ii}) , lv ) );

    if ~isempty(jj)

       jj = max(jj);

       v = varargin{2*jj};

       if strcmp( fs{ii},'Value')
          v0 = v(:)';
          n  = min( length(v) , length(val) );
          v  = val;
          v(1:n) = v0(1:n);
       end

       s = setfield( s , fs{ii} , v );

       switch fs{ii}
        case 'Type'

             switch v
              case 'text'
                   s.Color  =  NaN;
             end

        case 'Width'

             % MultipleLine-Edit
             s.Value(3) = s.Value(3) + 1 * ( imag(v) > 1 ) * ...
                                            strcmp(s.Type,'edit'); 

        case 'String'

         if isnan(real(s.Width))  & ...
            ~any( strcmp( s.Type , {'tag_frame' 'tab_list' 'sel_list'} ) )

           s.Width = size(char(s.String),2) + imag(s.Width);

         end             

       end

    end

end
