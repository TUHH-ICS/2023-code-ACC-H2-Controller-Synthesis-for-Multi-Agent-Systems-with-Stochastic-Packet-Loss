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

function [Fd, Fc, H2, P, solver_stats] = h2syn_nominal_dual(sysD, sysC, ny, nu, L0, backoff)
%H2SYN_NOMINAL_DUAL Function to synthesize a dynamic output-feedback
%controller for a decomposable LTI system that is suboptimal with respect
%to the H2-norm of the closed loop.
%   This function implements the decomposable output-feedback controller
%   synthesis that was proposed by Massioni and Verhaegen in 2009. In
%   addition, a backoff scheme proposed in
%
%   C. W. Scherer and S. Weiland, “Linear matrix inequalities in control,”
%
%   is implemented to improve the numerical conditioning of the resulting
%   controller at the cost of slightly worse performance bounds.
%
%   Arguments:
%       sysD    -> Decoupled part of the plant
%       sysC    -> Coupled part of the plant
%       ny      -> Number of measured outputs
%       nu      -> Number of control inputs
%       L0      -> Laplacian of the nominal graph
%       backoff -> [optional] Percentage of backoff. Default: 0.05
%   Returns:
%       Fd           -> Decoupled part of the controller
%       Fc           -> Coupled part of the controller
%       H2           -> Upper bound on the H2-norm of the closed-loop
%       P            -> Gramian matrix of the solution
%       solver_stats -> Timing information about the algorithm

% Setup the SDP solver
% The offset is used to convert strict LMI into nonstrict ones by pushing
% the solution away from the boundary.
offset = 1e-8;
opts   = sdpsettings('verbose', 0);

% Allocate storage for time measurements
solver_stats        = struct;
solver_stats.prep   = NaN;
solver_stats.solver = NaN;
solver_stats.trans  = NaN;
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

% Initialize return values with reasonable defaults
Fd    = [];
Fc    = [];
H2    = Inf;
P     = [];
timer = tic;

% Set default value for backoff
if nargin <= 5
    backoff = 0.05;
elseif backoff < 0
    error('Backoff must be non-negative!')
end

%% Prepare the system description
[Ac, Bc, Cc, Dc, T] = ssdata(sysC);
[Ad, Bd, Cd, Dd] = ssdata(sysD);

% Read performance channel size
nx = size(Ac,1);
nw = size(Bc,2) - nu;
nz = size(Cc,1) - ny;
nL = size(L0,1);

Bd_w  = Bd(:,1:nw);
Bc_w  = Bc(:,1:nw);
Bd_u  = Bd(:,nw+1:end);
Bc_u  = Bc(:,nw+1:end);
Cd_z  = Cd(1:nz,:);
Cc_z  = Cc(1:nz,:);
Cd_y  = Cd(nz+1:end,:);
Cc_y  = Cc(nz+1:end,:);
Dd_zw = Dd(1:nz,1:nw);
Dc_zw = Dc(1:nz,1:nw);
Dd_zu = Dd(1:nz,nw+1:end);
Dc_zu = Dc(1:nz,nw+1:end);
Dd_yw = Dd(nz+1:end,1:nw);
Dc_yw = Dc(nz+1:end,1:nw);
Dd_yu = Dd(nz+1:end,nw+1:end);
Dc_yu = Dc(nz+1:end,nw+1:end);

if ~iszero(Dd_yu) || ~iszero(Dc_yu)
    error('There may not be feedthrough from input to output (D_yu)!')
end

P_in  = ~iszero(Bc_u);
P_out = ~iszero(Cc_y);
if P_in && P_out
    error('The plant may not be distributed on both control input and measured output at the same time')
end

lambda = eig(L0);

%% Solve the SDP from Theorem 15, Massioni and Verhaegen, 2009
% The LMIs utilized in this function are dualized compared to the original
% ones from the paper. Instead of the controllability Gramien, the
% observability Gramian is used. This is done for the purpose of better
% comparability with the other results from the current paper.

alpha = sdpvar(1);
beta  = sdpvar(1);
Q     = sdpvar(nx, nx, nL-1, 'symmetric');
H     = sdpvar(nx, nx, nL-1, 'symmetric');
W     = sdpvar(nw, nw, nL-1, 'symmetric');
J     = sdpvar(nx, nx, nL-1, 'full');
S     = sdpvar(nx, nx, 'full');
X     = sdpvar(nx, nx, 'full');
Y     = sdpvar(nx, nx, 'full');
Kd    = sdpvar(nx, nx, 'full');
Kc    = sdpvar(nx, nx, 'full');
Ld    = sdpvar(nx, ny, 'full');
Lc    = sdpvar(nx, ny, 'full');
Md    = sdpvar(nu, nx, 'full');
Mc    = sdpvar(nu, nx, 'full');
Nd    = sdpvar(nu, ny, 'full');
Nc    = sdpvar(nu, ny, 'full');
cost  = 0;

% Check the additional constraints that need to be placed on the
% controller. These are taken from Table 1 in the above paper, which
% however contains an error of one missing but required constraint.
if ~P_in && ~P_out
    Constraints = [];
elseif ~P_in
    Constraints = [ Nc == 0, Lc == 0 ];
else
    Constraints = [ Nc == 0, Mc == 0 ];
end

% Abbreviation of a common term
Phi = alpha*eye(nx);

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 2:nL
    li  = lambda(i);
    Qi  = Q(:,:,i-1);
    Hi  = H(:,:,i-1);
    Wi  = W(:,:,i-1);
    Ji  = J(:,:,i-1);
    
    % Construct system and controller matrices for this modal subsystem
    A   = Ad    + li*Ac;
    Bu  = Bd_u  + li*Bc_u;
    Bw  = Bd_w  + li*Bc_w;
    Cy  = Cd_y  + li*Cc_y;
    Cz  = Cd_z  + li*Cc_z;
    Dyw = Dd_yw + li*Dc_yw;
    Dzu = Dd_zu + li*Dc_zu;
    Dzw = Dd_zw + li*Dc_zw;
    K   = Kd    + li*Kc;
    L   = Ld    + li*Lc;
    M   = Md    + li*Mc;
    N   = Nd    + li*Nc;
    
    % Main LMIs
    LMI = [ Qi          Ji           Y'*A'+M'*Bu'   K'             Y'*Cz'+M'*Dzu'  ;
            Ji'         Hi           A'+Cy'*N'*Bu'  A'*X'+Cy'*L'   Cz'+Cy'*N'*Dzu' ;
            A*Y+Bu*M    A+Bu*N*Cy    Y+Y'-Qi        eye(nx)+S'-Ji  zeros(nx,nz)    ;
            K           X*A+L*Cy     eye(nx)+S-Ji'  X+X'-Hi        zeros(nx,nz)    ;
            Cz*Y+Dzu*M  Cz+Dzu*N*Cy  zeros(nz,nx)   zeros(nz,nx)   eye(nz)         ];
    TRC = [ Wi             Bw'+Dyw'*N'*Bu'  Bw'*X'+Dyw'*L'  Dzw'+Dyw'*N'*Dzu' ;
            Bw+Bu*N*Dyw    Y+Y'-Qi          eye(nx)+S'-Ji   zeros(nx,nz)      ;
            X*Bw+L*Dyw     eye(nx)+S-Ji'    X+X'-Hi         zeros(nx,nz)      ;
            Dzw+Dzu*N*Dyw  zeros(nz,nx)     zeros(nz,nx)    eye(nz)           ];
    cost = cost + trace(Wi);
    
    % Add additional conditions to improve numerical conditioning
    NUM_1 = [ K, L; M, N ];
    NUM   = [ alpha * eye(nx+nu), NUM_1; NUM_1', alpha * eye(nx+ny) ]; 
    EIG   = [ Y+Y'-Qi             beta*eye(nx)+S'-Ji ;
              beta*eye(nx)+S-Ji'  X+X'-Hi            ];

    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        Constraints = [ Constraints                                ,...
                        Qi  <= alpha  * eye(nx)                    ,...
                        Hi  <= alpha  * eye(nx)                    ,...
                        LMI >= offset * eye(size(LMI))             ,...
                        TRC >= offset * eye(size(TRC))             ,...
                        NUM >= offset * eye(size(NUM))             ,...
                        EIG >= offset * eye(size(EIG))             ,...
                        [ Phi, Ji; Ji', Phi] >= offset * eye(2*nx) ];
    catch ME
        warning(ME.message)
        solver_stats.total = toc(timer);
        return
    end
end

% Add additional conditions to improve numerical conditioning
Constraints = [ Constraints                               ,...
                beta >= 1                                 ,...
                [ Phi, X; X', Phi ] >= offset * eye(2*nx) ,...
                [ Phi, Y; Y', Phi ] >= offset * eye(2*nx) ];

%% Implement three stage procedure to improve numerical conditioning
% This section implements the three stage procedure proposed by Scherer and
% Weiland in their lecture notes to improve the numerical conditioning of
% the controller.

solver_stats.prep = toc(timer);
sol = optimize(Constraints, cost, opts);

if backoff > 0
    if sol.problem ~= 0
        warning('YALMIP return an error: %s', sol.info)
        solver_stats.total = toc(timer);
        return
    end
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
    
    % Backoff is squared since the cost is gamma^2
    c_opt = (1+backoff)^2 * value(cost);
    sol   = optimize([Constraints, cost <= c_opt], alpha, opts);
    if sol.problem ~= 0
        warning('YALMIP return an error: %s', sol.info)
        solver_stats.total = toc(timer);
        return
    end
    solver_stats.yalmip = solver_stats.yalmip + sol.yalmiptime;
    solver_stats.solver = solver_stats.solver + sol.solvertime;
    
    a_opt = (1+backoff) * value(alpha);
    sol   = optimize([Constraints, cost <= c_opt, alpha <= a_opt], -beta, opts);
else
    solver_stats.yalmip = 0;
    solver_stats.solver = 0;
end

%% Calculate the controller from the final solution variables
if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
else
    solver_stats.yalmip = solver_stats.yalmip + sol.yalmiptime;
    solver_stats.solver = solver_stats.solver + sol.solvertime;
    trans = tic;
    
    % Calculate U & V from the identity UV' + XY = S
    [U, Sigma, V] = svd(value(S - X*Y));
    U = U*sqrtm(Sigma);
    V = V*sqrtm(Sigma);
    
    % Recover the Gramian from the new set of variables
    H2  = sqrt(value(cost));
    P   = zeros(2*nx,2*nx,nL-1);
    Psi = value([ Y   eye(nx)   ;
                  V'  zeros(nx) ]);
    for i = 1:nL-1
        P(:,:,i) = value([ Q(:,:,i)   J(:,:,i) ;
                           J(:,:,i)'  H(:,:,i) ]);
        P(:,:,i) = Psi' \ P(:,:,i) / Psi;
    end
    
    % Recover the decoupled part of the controller
    Dd_K = value(Nd);
    Cd_K = value(Md - Nd*Cd_y*Y)/V';
    Bd_K = U\value(Ld - X*Bd_u*Nd);
    Ad_K = U\value(Kd - X*(Ad+Bd_u*Nd*Cd_y)*Y)/V' - Bd_K*value(Cd_y*Y)/V' - U\value(X*Bd_u)*Cd_K;
    Fd = ss(Ad_K, Bd_K, Cd_K, Dd_K, T);
    
    % Recover the coupled part of the controller. This block is only true
    % under the assumption that the closed-loop system is affine in the
    % Laplacian.
    Dc_K = value(Nc);
    Cc_K = value(Mc - (Nc*Cd_y+Nd*Cc_y)*Y)/V';
    Bc_K = U\value(Lc - X*(Bc_u*Nd+Bd_u*Nc));
    Ac_K = U\value(Kc - X*(Ac+Bc_u*Nd*Cd_y+Bd_u*Nc*Cd_y+Bd_u*Nd*Cc_y)*Y)/V' ...
           - value((Bc_K*Cd_y+Bd_K*Cc_y)*Y)/V' - U\value(X*(Bc_u*Cd_K+Bd_u*Cc_K));
    Fc = ss(Ac_K, Bc_K, Cc_K, Dc_K, T);
    
    solver_stats.trans = toc(trans);
end

solver_stats.total = toc(timer);
end
