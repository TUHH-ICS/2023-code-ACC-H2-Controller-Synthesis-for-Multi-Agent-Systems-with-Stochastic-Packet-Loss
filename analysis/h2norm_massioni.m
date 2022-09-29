%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
%
% Contains a self-implemented version of the approach described in
% P. Massioni and M. Verhaegen, “Distributed control for identical dynamically coupled systems:
% A decomposition approach,” IEEE Transactions on Automatic Control, 2009
%
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, time] = h2norm_massioni(sysD, sysC, Kd, Kc, L)
%H2NORM_MASSIONI Calculates the H2-norm of a decomposable LTI system
%   This function calculates the H2-norm of a decomposable LTI systme in a
%   scalable manner by utilizing the approach proposed by Massioni and
%   Verhaegen. The system is decomposed into modal subsystem, which are
%   then analysed one by one, leading to linear complexity in the number of
%   agents.
%
%   Arguments:
%       sysD -> Decoupled part of the open-loop plant
%       sysC -> Coupled part of the open-loop plant
%       Kd   -> Decoupled part of the controller
%       Kc   -> Coupled part of the controller
%       L    -> Pattern matrix
%   Returns:
%       H2   -> H2-norm of the closed-loop system
%       time -> Time taken to for the total calculations

timer  = tic;

lambda = eig(L);
H2     = 0;
for i = 2:length(lambda)
    sys = addparts(sysD, sysC, lambda(i));
    K   = addparts(Kd, Kc, lambda(i));
    CL  = lft(sys, K);
    H2  = H2 + norm(CL)^2;
end

H2   = sqrt(H2);
time = toc(timer);
end
