function main_msg(par,varargin)

% MAIN_MSG  Set Message to MessageList


[Msg,fig] = recpar(par,'Root');

if ~isempty(Msg)
    return
end

fig = fig(1);


ud = get( fig               , 'userdata' );
ud = get( ud.Children.Frame , 'userdata' );

hl = ud.Children.Logo.Message;

try
  Msg = msg_list(hl,'Message',varargin{:});
catch 
  Msg = lasterr;
end


if ~isempty(Msg)

  Msg = ['Error call MSG_LIST( Message ).' char(10) Msg ];

  warndlg(Msg,'Warning','warn');

end