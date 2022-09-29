%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [CL_d, CL_c, CL_p] = assemble_cl(sysD, sysC, sysP, Kd, Kc)
%ASSEMBLE_CL Combine the plant and controller into the closed-loop system
%   Because we treat decomposable system in decomposed form, assembling the
%   closed-loop cannot be done using feedback() or lft(). This function
%   assembles the closed loop under the assumption that it is affine in the
%   Laplacian and will fail otherwise

[Ad, Bd, Cd, Dd, Ts] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP);
[Ad_K, Bd_K, Cd_K, Dd_K] = ssdata(Kd);
[Ac_K, Bc_K, Cc_K, Dc_K] = ssdata(Kc);

[nu, ny] = size(Dd_K);
nw = size(Bd,2) - nu;
nz = size(Cd,1) - ny;

Bd_w  = Bd(:,1:nw);
Bc_w  = Bc(:,1:nw);
Bp_w  = Bp(:,1:nw);
Bd_u  = Bd(:,nw+1:end);
Bc_u  = Bc(:,nw+1:end);
Bp_u  = Bp(:,nw+1:end);
Cd_z  = Cd(1:nz,:);
Cc_z  = Cc(1:nz,:);
Cp_z  = Cp(1:nz,:);
Cd_y  = Cd(nz+1:end,:);
Cc_y  = Cc(nz+1:end,:);
Cp_y  = Cp(nz+1:end,:);
Dd_zw = Dd(1:nz,1:nw);
Dc_zw = Dc(1:nz,1:nw);
Dp_zw = Dp(1:nz,1:nw);
Dd_zu = Dd(1:nz,nw+1:end);
Dc_zu = Dc(1:nz,nw+1:end);
Dp_zu = Dp(1:nz,nw+1:end);
Dd_yw = Dd(nz+1:end,1:nw);
Dc_yw = Dc(nz+1:end,1:nw);
Dp_yw = Dp(nz+1:end,1:nw);
Dd_yu = Dd(nz+1:end,nw+1:end);
Dc_yu = Dc(nz+1:end,nw+1:end);
Dp_yu = Dp(nz+1:end,nw+1:end);

if ~iszero(Dd_yu) || ~iszero(Dc_yu) || ~iszero(Dp_yu)
    error('The plant may not have direct feedthrough from control input to measured output')
end

P_in  = ~iszero(Bc_u) || ~iszero(Bp_u) || ~iszero(Dc_zu) || ~iszero(Dp_zu);
P_out = ~iszero(Cc_y) || ~iszero(Cp_y) || ~iszero(Dc_yw) || ~iszero(Dp_yw);
K_in  = ~iszero(Bc_K) || ~iszero(Dc_K);
K_out = ~iszero(Cc_K) || ~iszero(Dc_K);
if (P_in && P_out) || (P_in && K_out) || (P_out && K_in)
    error('The resulting closed-loop system is not linear in the Laplacian')
end

%% Decoupled part
Ad_CL = [ Ad+Bd_u*Dd_K*Cd_y  Bd_u*Cd_K ;
          Bd_K*Cd_y          Ad_K      ];
Bd_CL = [ Bd_w+Bd_u*Dd_K*Dd_yw ;
          Bd_K*Dd_yw           ];
Cd_CL = [ Cd_z+Dd_zu*Dd_K*Cd_y  Dd_zu*Cd_K ];
Dd_CL = Dd_zw+Dd_zu*Dd_K*Dd_yw;
CL_d  = ss(Ad_CL, Bd_CL, Cd_CL, Dd_CL, Ts);

%% Stochastically connected part
Ac_CL = [ Ac+Bc_u*Dd_K*Cd_y+Bd_u*Dc_K*Cd_y+Bd_u*Dd_K*Cc_y  Bc_u*Cd_K+Bd_u*Cc_K ;
          Bc_K*Cd_y+Bd_K*Cc_y                              Ac_K                ];
Bc_CL = [ Bc_w+Bc_u*Dd_K*Dd_yw+Bd_u*Dc_K*Dd_yw+Bd_u*Dd_K*Dc_yw ;
          Bc_K*Dd_yw+Bd_K*Dc_yw                                ];
Cc_CL = [ Cc_z+Dc_zu*Dd_K*Cd_y+Dd_zu*Dc_K*Cd_y+Dd_zu*Dd_K*Cc_y  Dc_zu*Cd_K+Dd_zu*Cc_K ];
Dc_CL = Dc_zw+Dc_zu*Dd_K*Dd_yw+Dd_zu*Dc_K*Dd_yw+Dd_zu*Dd_K*Dc_yw;
CL_c  = ss(Ac_CL, Bc_CL, Cc_CL, Dc_CL, Ts);

%% Deterministically connected part
Ap_CL = [ Ap+Bp_u*Dd_K*Cd_y+Bd_u*Dd_K*Cp_y  Bp_u*Cd_K         ;
          Bd_K*Cp_y                         zeros(size(Ad_K)) ];
Bp_CL = [ Bp_w+Bp_u*Dd_K*Dd_yw+Bd_u*Dd_K*Dp_yw ;
          Bd_K*Dp_yw                           ];
Cp_CL = [ Cp_z+Dp_zu*Dd_K*Cd_y+Dd_zu*Dd_K*Cp_y  Dp_zu*Cd_K ];
Dp_CL = Dp_zw+Dp_zu*Dd_K*Dd_yw+Dd_zu*Dd_K*Dp_yw;
CL_p  = ss(Ap_CL, Bp_CL, Cp_CL, Dp_CL, Ts);
end
