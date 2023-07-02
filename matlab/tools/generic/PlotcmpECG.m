function  fig = PlotcmpECG(x, y, L,varargin)
% Plot two multi channel ECG for comparison
%
% inputs:
% x: first input data matrix
% y: second input data matrix
% L: number of panels per figure
% fs: sampling rate (1 by default)
%
% The Open Source Electrophysiological Toolbox, 
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

if(nargin>3 && ~isempty(varargin{1}))
    fs = varargin{1};
    lbl = 'time(s)';
else 
    fs = 1;
    lbl = 'index';
end

if(nargin>4  && ~isempty(varargin{2}))
    ttl = varargin{2};
else
    ttl = '';
end

L1 = size(y, 1); L2 = size(y, 2);
% h of each axis in normalized units
axish = (1 / L) * 0.65;

t = (0:L2-1)/fs;
for i = 1:L1
    
    if(mod(i,L)==1 || L==1)
        fig = figure;
        lbwh = get(1, 'position');
        figw = lbwh(3);
        figh = lbwh(4);
    end
    row = mod( i-1, L );
    axisb = (axish+0.025) * (L - row);
    h= subplot('position', [0.1, axisb, 0.85, axish] );
        
    plot(t,x(i,1:L2),'Color', [.5 .5 .5]);
    hold on;
    plot(t,y(i,1:L2),'r');
    if (i ~= L)
    end
    xlabel('Time(s)');
    fig = gcf;
    ax = fig.CurrentAxes;
    ax.FontSize = 11;
    ylabel(num2str(i));
    grid;
    if(mod(i,L)==1 || L==1)
        title(ttl);
    end
    if(mod(i,L)==0 || L==1)
        xlabel(lbl);
    end
end
