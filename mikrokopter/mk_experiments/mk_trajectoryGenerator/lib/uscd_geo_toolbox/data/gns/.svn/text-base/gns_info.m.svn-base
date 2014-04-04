function ini = gns_info

% GNS_INFO  Returns GNS-InfoData from GNS_INFO.TXT
%
% Info = GNS_INFO
%
% Info.src = { #Source }
% Info.inf = { #Contents }
%
% Info.fld = { Field  Name  Description }  #FieldDescription
%
% Info.dsg.[A:Z] = { Code  Name }          #DesignationCodes
%
% Info.cnt = { CountyCode   CountryName }  #CountryCodes
%

ini = { 'src'  'Source'
        'inf'  'Contents'
        'fld'  'FieldDescription'
        'dsg'  'DesignationCodes'
        'cnt'  'CountryCodes'     };

[m,v] = load_asc('gns_info.txt',ini(:,2),'#','%');

if ~isempty(m)
    error(m)
end

%---------------------------------------------------
% Build Structure

for ii = 1 : size(ini,1)
    ini{ii,2} = v(ii,2);
end

ini = permute(ini,[2 1]);

ini = struct(ini{:});

%---------------------------------------------------
% Source

if chkstr(ini.src,1)

    [m,v] = read_par(ini.src,{'' '@'},'','>');

    if isempty(m) & ~isempty(v)
       ini.src = v;
    end

else

  ini.src = [];

end

%---------------------------------------------------
% Content

if chkstr(ini.inf,1)

    [m,v] = read_par(ini.inf,{'' '>'},'','%');

    if isempty(m) & ~isempty(v)
       ini.inf = v;
    end

else

  ini.inf = [];

end


%---------------------------------------------------
% Fields

if chkstr(ini.fld,1)

    [m,v] = read_par(ini.fld,{'=' ''},'+','>');

    if isempty(m) & ~isempty(v)
       v = v(:,[1 2 2]);
       for ii = 1 : size(v,1)
           c = sepname(v{ii,2},NaN,char(10));
           v{ii,2} = rmblank(c{1});
           if size(c,2) == 1
              v{ii,3} = '';
           else
              v{ii,3} = sprintf('%s\n',c{2:end});
           end
       end
       ini.fld = v;
    end

else

  ini.fld = [];

end

%---------------------------------------------------
% DesignationCodes

if chkstr(ini.dsg,1)

    k = cellstr(char(abs('A'):abs('Z'))');

    [m,v] = load_asc(ini.dsg,k,'@','%');

    if isempty(m) & ~isempty(v)
       v = v(:,[1 2]);
       for ii = 1 : size(v,1)
           if chkstr(v{ii,2},1)
              [m,d] = read_par(v{ii,2},{'=' ''},'+','>');
              if isempty(m) & ~isempty(d)
                 for jj = 1 : size(d,1)
                     if any( d{jj,2} == 10 )
                        d{jj,2}(find(d{jj,2}==10)) = ' ';
                        d{jj,2} = remfollw(d{jj,2},' ');
                     end
                 end
                 v{ii,2} = {d};
              else
                 v{ii,2} = v(ii,2);
              end
           else
              v{ii,2} = {cell(0,2)};
           end
       end
       v = permute(v,[2 1]);
       ini.dsg = struct(v{:});
    end

else

  ini.dsg = [];

end

%---------------------------------------------------
% CountryCodes

if chkstr(ini.cnt,1)

    [m,v] = read_par(ini.cnt,{'=' ''},'','%');

    if isempty(m) & ~isempty(v)
       ini.cnt = v;
    end

else

  ini.cnt = [];

end
