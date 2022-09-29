%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function G = line_graph(nvert, nedge, directed)
%LINE_GRAPH Generates a Matlab graph object that represents a line
%graph with a given number of vertices.
%   This function genenerates a Matlab graph object that contains a
%   line graph with the given amount of vertices. It can be configured if
%   the connections should be directed or not and to how many of the next
%   vertices the connection should be established.
%
%   Arguments:
%       nvert    -> Number of vertices
%       nedge    -> Number of forward edges per vertex
%       directed -> [optional] Directed edges or not

if nargin <= 1 || isempty(nedge)
    nedge = 1;
elseif nedge >= nvert
    error('So many neighbours are not possible!')
end
if nargin <= 2
    directed = true;
end

% Define closure that keeps the index in the ring [1, nvert]
lim = @(i) min(nvert, max(1, i));

%% Generate graph structure
A = zeros(nvert);
for i = 1:nvert-1
    A(i, (i+1):lim(i+nedge)) = 1;
end

% If not direction, also add the other direction
if ~directed   
    for i = 2:nvert
        A(i, lim(i-nedge):(i-1)) = 1;
    end
end

% G needs only to be a digraph if it is directed
if directed
    G = digraph(A);
else
    G = graph(A);
end
end
