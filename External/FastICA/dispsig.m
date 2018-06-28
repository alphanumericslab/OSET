function dispsig(signalMatrix, range, titlestr);
%DISPSIG - deprecated!
%
% Please use icaplot instead.
%
%   See also ICAPLOT

% 24.8.1998
% Hugo Gävert

fprintf('\nNote: DISPSIG is now deprecated! Please use ICAPLOT.\n');

if nargin < 3, titlestr = ''; end
if nargin < 2, range = 1:size(signalMatrix, 1); end

icaplot('dispsig',signalMatrix',0,range,range,titlestr);
