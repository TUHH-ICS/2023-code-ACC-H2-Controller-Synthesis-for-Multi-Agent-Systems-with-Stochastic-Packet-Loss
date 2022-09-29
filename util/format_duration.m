%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function str = format_duration(duration)
%FORMAT_DURATION Format a duration given in seconds into a more easily
%human readible string
%   For long times measured with tic/toc, the output becomes very hard to
%   read. This function improves the formating drastically.

% Split duration into larger quantities
days     = floor(duration / (60*60*24));
duration = mod(duration, 60*60*24);
hours    = floor(duration / (60*60));
duration = mod(duration, 60*60);
minutes  = floor(duration / 60);
duration = mod(duration, 60);
seconds  = floor(duration);
duration = duration - seconds;

% Only print large quantities if they are nonzero
if days > 0
    str = sprintf('%dd ', days);
else
    str = '';
end
if hours > 0
    str = [str, sprintf('%dh ', hours)];
end
if minutes > 0
    str = [str, sprintf('%dm ', minutes)];
end
if days == 0
    if seconds > 0
        str = [str, sprintf('%ds ', seconds)];
    end
    if minutes == 0
        str = [str, sprintf('%dms', floor(duration*1e3))];
    end
end
str = strtrim(str);
