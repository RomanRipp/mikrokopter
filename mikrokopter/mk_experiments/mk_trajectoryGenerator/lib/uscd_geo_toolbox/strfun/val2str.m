function [str,val,form,msg] = val2str(varargin)

% VAL2STR  Converts between Formated Strings and values
%
% [ String , Value , Format , Msg ] = VAL2STR( Argument , Format )
%
% VAL2STR is the same as STR2VAL with permuted Outputs !!!
%
% Argument can be a String or a Value, Format gives the Conversion 
%   between them:   String <--- Format ---> Value
%
% To Replace MissingValues or Dummies by NaN use a 3. and 4. Input:
%
%    VAL2STR( Argument , Format , Dummy , Accuracy )
% 
% All Values which are in Range of Accuracy around the Dummies
%  will be replaced by NaN. Dummy can be a Vector or Matrice.
%  Example:   STR2VAL( ... , [ -9999  1e32 ] , 1e-10 )
% 
% Valid Formats are:
%
%  'char'    String == Value
%
%  '%#.#*'   for using by SPRINTF
%
%  'glon#'   for Geographic Coordinates Longitude    **°**.**' <E|W>
%  'glat#'   for Geographic Coordinates  Latitude    **°**.**' <N|S>
%  'gpos#'   for Geographic Coordinates  Lat & Lon: [ "glat#"  "glon#" ] 
%
%            # gives the Number of decimal Minutes
%
%            use the first upper Character "G" in Positionsformat
%              to return a short Form:  ***<E|W>**.**  or  ***<N|S>**.**
%
%  'dgms#'   for Sexagesimal Conversion: [ Degree Minute Second ]
%                # == '-'   Degree = ( -180 .. 180 ]   ±DDD° MM' SS"
%                # == '+'   Degree = [    0 .. 360 )    DDD° MM' SS"
%
%  'time#'   for TimeConversation, first Value  Day: # == 1    DD HH:MM
%                                              Hour: # == 2       HH:MM
%                                              Hour: # == 3       HH:MM:SS
%                                            Minute: # == 4          MM:SS
%
%            increase the Value by 4, to use Duplicate TimeValues:
%
%                # == 5    DD HH:MM - DD HH:MM
%                # == 6       HH:MM - HH:MM
%                # == 7       HH:MM:SS - HH:MM:SS
%                # == 8          MM:SS - MM:SS
%
%              use 'time#<seperator>' to define the Seperator between the Values
%               default: ' - '
%
%  'date#*'  for DateConversion, # == 0       Month as Number
%                                     1       Month as Short String (3 Characters)
%                                     2       Month as long String
%                                * == e | g   english / german Date only
%                                     E | G   english / german Date and HH:MM
%
%                                     d | D   YYYY-MM-DD Notation
%
%            Note:  Year < 100  will expand with 2000 !!!
%
%            Example:
%
%            'date0E'   MM/DD/YYYY HH:MM       12/01/1995 15:45
%            'date0G'   DD/MM/YYYY hh:mm       01.12.1995 15:45
%
%            'date1E'   DD-MMM-YYYY hh:mm      01-Dec-1995 15:45
%            'date1G'   DD.MMM.YYYY hh:mm      01.Dez.1995 15:45
%
%            'date2E'   DD. Month YYYY hh:mm   01. December 1995 15:45
%            'date2G'   DD. Monat YYYY HH:MM   01. Dezember 1995 15:45
%
%            'date0D'   YYYY/MM/DD hh:mm       1995/12/01 15:45
%            'date1D'   YYYY-MM-DD hh:mm       1995-12-01 15:45
%            'date2D'   YYYY MM DD hh mm       1995 12 01 15 45
%
%
% { Seperator }      Converts between String and CellArray of Strings,
%                     splitted by the single SeperatorCharacter
%
% { Value  Format }  a 2-Column CellArray gives allowed Values 
%                      and the corresponding Format
%
%
%  '@'            String returns the Size and Class of Argument
%
%  '@ClassName'   Specify the Class of Value
%
% Special ClassNames are:
%
% 'numeric'   True if ISNUMERIC
% 'string'    True for ONE-Row CharacterArray or ONE-Element CellArray of it
% 'cellstr'   True for CellArray or CharacterArray of Strings
% 'object'    True if ISOBJECT
%
% Special additional Handlings for Matlab-defined Class's:
%
% 'double'    True if ISNUMERIC
% 'single'    True if ISNUMERIC
% 'sparse'    True if ISNUMERIC
% 'logical'   True if ISNUMERIC  with Values [ 0 | 1 ]
% 'int<B>'    True if ISNUMERIC  with Values [ -2^(B-1) .. 2^(B-1)-1 ]
% 'uint<B>'   True if ISNUMERIC  with Values [  0       .. 2^(B)-1   ]
% 'char'      True if ISCELLSTR
% 'cell'      True if CharacterArray
% 'struct'    True if Object
% ClassName   True if Object is of given Class
% @ClassName  True if Object inherits of given Class, using ISA, returns Object
% @ClassName@ True if Object inherits of given Class, using ISA, returns Parent
%
%
%-------------------------------------------------------------------------------
%
% See also: STR2VAL, STR2VEC, SPRINTF, VAR2MSTR, STRHCAT, SEPNAME, CLASSCHK
%
% requires: STR2VAL
%


[val,str,form,msg] = str2val(varargin{:});


