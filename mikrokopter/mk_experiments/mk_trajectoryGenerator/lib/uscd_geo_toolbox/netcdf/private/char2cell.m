function  [Msg,bb,ref] = char2cell(bb,mnl,ref01);

% CHAR2CELL Converts CharArray to CellStringArray
%
% [Msg,CellString] = CHAR2CELL( String )
%
% Converts the CharacterArray String into a CellStringArry.
%
%---------------------------------------------------------
% CHAR2CELL( String , MarkerNewLine )  
%
% Replace the String's defined by MarkerNewLine with NewLine,
%   default: EMPTY
%
%---------------------------------------------------------
% [Msg,CellString,Reference] = ...
%    CHAR2CELL( String , MarkerNewLine , ReferenceCell)
%  
% Returns References in String, using GET_REF, 
%   ReferenceCell = { S1 S2 C1 C2 }, type: >> help get_ref
%


 Msg = '';
 
 Msg0 = 'CHAR2CELL: ';

 ref = cell(0,2);

 Nin  = nargin;
 Nout = nargout;
 
 nl = char(10);

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  % Marker for NewLine
  if Nin < 2
    mnl = '';
  end

%*****************************************************
% Check Inputs

  if iscellstr(bb)
     bb = cellstr(char(bb));
     bb = strhcat(bb,'',1);
  end

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  ok = ( isnumeric(bb) | ischar(bb) );
  if ok

    if ischar(bb)
       bb = double(bb);
    end

    ok = all( ( mod(bb,1) == 0 )  & ( bb >= 0 ) & isfinite(bb)  );

  end

  if ~ok

      Msg = [ Msg0 ...
              'Input String must be a CharacterArray or ASCII-Codes.'];
 
      bb = cell(0,1);

      return

  end

  %---------------------------------------------------
  % Check MarkerNewLine

  if ~( isempty(mnl) |  ...
        ( ischar(mnl) &  ( prod(size(mnl)) == size(mnl,2) ) ) )

    Msg = [ Msg0  'Input MarkerNewLine must be a String.'];

    return

  end

  %---------------------------------------------------
  % Check Reference

  if ( Nin == 3 )  &  ( Nout == 3 )

     [Msg,ref] = get_ref('',ref01{:});

    if ~isempty(Msg)
       Msg = [ Msg0  'Invalid Input for Reference.' nl Msg ];
       return
    end

  end

%*****************************************************

  if ( size(bb,1) > 1 )  &  ( size(bb,2) > 1 )
     bb = cat( 2 , bb , 10*ones(size(bb,1),1) );
     bb = permute( bb , [ 2 1 ] );
     bb = bb(:);
  end

  if ( size(bb,1) > 1 ) 
     bb = permute( bb , [ 2 1 ] );
  end

  %---------------------------------------------------
  % Check Characters

    if any( ( 126 < bb ) & ( bb < 160 ) );

       % Old DOS !!!
       old = [ 132   148   129   142   153   154 ];
       new = [ 'ä'   'ö'   'ü'   'Ä'   'Ö'   'Ü' ];

       for ii = 1 : size(old,2)
           bb(find(bb==old(ii))) = double(new(ii));
       end 

    end

  ok = all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

  if ~ok
    Msg = [Msg0 'Invalid Characters in String.' ];
    return
  end


%*****************************************************
 

  %---------------------------------------------------
  % Remove CR
  bb( find( bb == 13 ) ) = [];


  bb = char(bb);


  %---------------------------------------------------
  % TAB  --> 8 Blanks 
  bb = strrep( bb , char(9) , char(32*ones(1,8)) ); 

  %---------------------------------------------------
  % mnl --> NewLine   % !!!!!
  if ~isempty(mnl)
    bb = strrep( bb , mnl , char(10) ); 
  end

  %---------------------------------------------------
  % Reference
  if ( Nin == 3 )  &  ( Nout == 3 )

     [MsgR,ref] = get_ref(bb,ref01{:});

  end

  %---------------------------------------------------
  % Form CellString

  % 1. "'"     --> "''"
  bb = strrep( bb , char(39) , char([39  39]) ); 

  % 2. NL --> "';'"
  bb = strrep( bb , char(10) , char([39  59  39]) );
  

  bb = eval([  '{'''  bb  '''}' ]);
