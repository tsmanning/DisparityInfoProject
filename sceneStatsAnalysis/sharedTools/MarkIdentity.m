function h=MarkIdentity(hndl,varargin)
% function MarkIdentity(hndl,varargin)
% Marks the identity line on a set of axes.
% Example: plot(x,y); MarkIdentity(gca,'k')
% hndl must be the handle of a set of axes (usually 'gca').
% varargin allows you to supply arguments as for plot.
% The 'hold' status of the plot is left in its original condition.
% Example: x=rand(1,100); y = x + 2*rand(1,100); plot(x,y,'k.');
% Now to draw the identity line on with a thick dotted green line:
% MarkIdentity(gca,'color','g','linewidth',2,'linestyle',':');
% Jenny Read 2/4/2003

% Remember whether axes are st to hold on or hold off
holdstatus = get(hndl,'NextPlot');
% Set holdto on so that identity does not redraw data
hold on
xlm = get(gca,'xlim');
ylm = get(gca,'ylim');
tmp = [min([xlm ylm]) max([xlm ylm])];
if nargin==1
    plot(tmp,tmp)
else % use user-specified inputs to plot
    str = [];
    for j=1:nargin-1
        if ischar(varargin{j})
            str = [str ',''' varargin{j} ''''];
        else
            str = [str ',' num2str(varargin{j}) ''];
        end
    end
    eval(['h=plot(tmp,tmp' str ');']);
end


% Restore hold status and axis limits
set(gca,'NextPlot',holdstatus,'xlim',xlm,'ylim',ylm)