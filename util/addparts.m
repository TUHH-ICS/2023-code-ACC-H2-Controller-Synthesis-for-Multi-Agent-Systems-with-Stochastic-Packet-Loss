%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function sys = addparts(sysD, sysC, lambda)
%ADDPARTS Function that linearly combines the matrices of the two systems
%   The system sysC is scaled by lambda before being added onto sysD

[Ad, Bd, Cd, Dd, Ts] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
sys = ss(Ad+lambda*Ac, Bd+lambda*Bc, Cd+lambda*Cc, Dd+lambda*Dc, Ts);     
end
