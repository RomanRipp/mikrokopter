function mkpcode(pfad,varargin)

% MKPCODE  Creates pre-parsed pseudo-code file
%
% mkpcode(pfad,'-r','>Zielpfad')
% mkpcode(File,'-r','>Zielpfad')
%
%  '-r'  works recursivly
%
% Output to ZielPfad, not given ==> Option INPLACE
%


VarIn = varargin;
VarIn = VarIn(:);

Nin = size(VarIn,1);


nl = char(10);
fs = filesep;

ext = '.m';

%--------------------------------------------
% Get Options

rec  = 0;
ziel = '';
lid  = [];

for ii = 1 : Nin

  %--------------------------------------------
  % Check for Recurse, Ziel
  if ischar(VarIn{ii}) & ~isempty(VarIn{ii}) & ...
     ( prod(size(VarIn{ii})) == size(VarIn{ii},2) )

    if strcmp(VarIn{ii}(1),'>');
      ziel = VarIn{ii}(2:end);
      if isempty(ziel)
         ziel = cd;
      end
    end

    rec = ( rec | strcmp(VarIn{ii},'-r') );

  %--------------------------------------------
  % Check for LogID
  elseif isnumeric(VarIn{ii}) & ~isempty(VarIn{ii}) & ...
           ( prod(size(VarIn{ii})) == 1 )

    if ~isempty(fopen(VarIn{ii}))
       lid = VarIn{ii};
    end

  end

end


%--------------------------------------------
% Check Pfad

if isempty(pfad)
 pfad = cd;
end

if any(strcmp(pfad,{ '..' }))
  return
end

if any(strcmp(pfad,{ '.'  '*' }))
  pfad = cd;
end


%--------------------------------------------
% Try to get absolute Pfad, Ziel

p00 = cd;
if isempty(p00)
   p00 = fs;
end

p00 = cat( 2, p00 , fs(1:(end*(~strcmp(p00(end),fs)))) );


for p = { 'pfad'  'ziel' }

  
  p0 = eval(p{1});

  if ~isempty(p0)

    d1 = dir(p0);

    p0 = cat(2,p00,eval(p{1}));

    d0 = dir(p0);

    if isequal( d0 , d1 )

       p0 = cat( 2 , p0 , fs(1:(end*(~strcmp(p0(end),fs)))) );
    
       eval([ p{1} ' = p0;' ])

    end

  end

end


%--------------------------------------------
% Change to Ziel

OrgPfad = cd;
new_log = 0;


if ~isempty(ziel)

   try
     cd(ziel)
   catch
     fprintf(['Error call CD('''  ziel ''').' lasterr nl ]);
     return
   end

   new_log = isempty(lid);

   if new_log
      lid = fopen([ ziel  '.pcode.log' ],'wt');
      if lid == -1
        fprintf(['Can''t write into '''  ziel '''.']);
        cd(OrgPfad);
        return
      end
   end       

end
   

%--------------------------------------------

fprintf(pfad);

d = dir(pfad);

if isempty(d)
  fprintf([' ... can''t read File or Directory ' nl ]);
  cd(OrgPfad)
  return
end

IsFile = ( (size(d,1)==1)  &  (d(1).isdir==0) );
if ~IsFile
  if ~strcmp(pfad(length(pfad)),fs)
    pfad = [ pfad fs ];
  end
end


fprintf(nl)

is_file = find(~cat(1,d.isdir));

for ii = is_file(:)'

  name = d(ii).name;

  if ~strcmp( name , 'mkpcode.m' )

    % Check for m-File
    kk = findstr(name,ext);
  
    if  length(kk) == 1  

      if  kk == length(name)-length(ext)+1

        file = [ pfad((1:(end*(~IsFile))))  name ];

        PIn = { '-inplace'};
        PIn = PIn( 1 : isempty(ziel) );


        try
          pcode( file , PIn{:} );
          if ~isempty(lid)
             fprintf(lid,[file nl]);
          end
        catch
          fprintf(['Error call PCODE('''  file ''').' nl lasterr nl ]);
        end
 
      end

    end
  end

end

if rec

  is_dir = find(cat(1,d.isdir));

  for ii = is_dir(:)'

    if ~any(strcmp(d(ii).name,{ '.'  '..' }))

      try
         mkpcode([pfad d(ii).name],'-r',['>' ziel],lid);
      catch
          fprintf(['Error call MKPCODE('''  [pfad d(ii).name] ''',''-r'',''' ...
                   ['>' ziel] ''','  sprintf('%g',lid)  ').' nl lasterr nl ]);
      end

    end

  end

end      

   
if new_log
   fclose(lid);
end

cd(OrgPfad);
