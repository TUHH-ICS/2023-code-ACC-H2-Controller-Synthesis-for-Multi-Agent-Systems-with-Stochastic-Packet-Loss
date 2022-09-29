%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function ret = iszero(A, tol)
%ISZERO Checks whether a matrix is completely zero
%   Direct comparison with 0 will lead to spurios false negatives, since
%   floating point arithmetic is not exact. For that reason this function
%   has a tolerance build in, until which point an entry counts as zero.

if nargin <= 1
    tol = 1e-8;
end

ret = all(abs(A) < tol, 'all');
end
