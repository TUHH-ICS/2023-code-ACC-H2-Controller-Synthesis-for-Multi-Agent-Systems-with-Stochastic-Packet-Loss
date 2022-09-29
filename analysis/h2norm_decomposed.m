%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Q, solver_stats] = h2norm_decomposed(sysD, sysC, sysP, L0, p, robust)
%H2NORM_DECOMPOSED Calculate an upper bound on the H2-norm of a
%decomposable jump system in a scalable manner
%
%   Arguments:
%       sysD   -> Decoupled part of the system
%       sysC   -> Stochastically coupled part of the system
%       sysP   -> Deterministically coupled part of the system
%       L0     -> Laplacian describing the nominal graph of the system
%       p      -> Probability of a successful package transmission
%       robust -> [optional] Use single Z matrix for robust analyis. Default: false
%   Returns:
%       H2           -> Upper bound on the H2-norm of the system
%       Q            -> Storage function matrix of the solution
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
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

% Initialize return values with reasonable defaults
H2    = inf;
timer = tic;

% Set default value for robust synthesis
if nargin <= 5
    robust = false;
end

%% Prepare the system description
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP);

nx = size(Ac,1);
nu = size(Bc,2);
N  = size(L0,1);

lambda = eig(L0);

%% Solve the SDP from Theorem 2
Q    = sdpvar(nx, nx);
cost = 0;
Constraints = Q >= offset * eye(nx);

% If robust flag is set, a single Z variable will be used. Since the LMIs
% are concex in lambda(i), its enough to check the boundary in that case.
if robust
    Z   = repmat(sdpvar(nu, nu), 1, 1, N-1);
    set = [2, N];
else
    Z   = sdpvar(nu, nu, N-1);
    set = 2:N;
end

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = set
    li  = lambda(i);
    lit = p*(1-p)*li;
    Zi  = Z(:,:,i-1);
    
    Acl = Ad + p*li*Ac + li*Ap;
    Bcl = Bd + p*li*Bc + li*Bp;
    Ccl = Cd + p*li*Cc + li*Cp;
    Dcl = Dd + p*li*Dc + li*Dp;
    
    LMI = Acl'*Q*Acl + Ccl'*Ccl - Q  + 2*lit * (Ac'*Q*Ac + Cc'*Cc);
    TRC = Bcl'*Q*Bcl + Dcl'*Dcl - Zi + 2*lit * (Bc'*Q*Bc + Dc'*Dc);
    cost = cost + trace(Zi);
    
    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        Constraints = [ Constraints                     ,...
                        LMI <= -offset * eye(size(LMI)) ,...
                        TRC <= -offset * eye(size(TRC)) ];
    catch ME
        warning(ME.message)
        Q = [];
        solver_stats.total = toc(timer);
        return
    end
end

% Since we only checked the boundary, we need to change the cost to reflect
% that
if robust
    cost = (N-1)*trace(Z(:,:,1));
end

solver_stats.prep = toc(timer);
sol = optimize(Constraints, cost, opts);

if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
    Q = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt(value(cost));
    
    Q  = value(Q);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
end

solver_stats.total = toc(timer);
end
