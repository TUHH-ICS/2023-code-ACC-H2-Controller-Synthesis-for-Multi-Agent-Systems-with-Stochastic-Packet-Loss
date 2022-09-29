%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Q, solver_stats] = h2norm_enumerated(sysD, sysC, sysP, Kd, Kc, L0, p, symmetric)
%H2NORM_ENUMERATED Calculates the H2-norm of the given decomposable jump
%system with Bernoulli packet loss.
%   This function calculates the H2-norm of the given decomposable jump
%   system by application of Theorem 1 from the accompanying paper. To do
%   so, it needs to enumerate all possible modes, the number of which grows
%   exponentially with the number of agents in the system.
%   In contrast to the other performance functions, we must distinguish
%   between symmetric and asymmetric packet loss, because this will change
%   the performance measure in this case.
%
%   Arguments:
%       sysD      -> Decoupled part of the plant
%       sysC      -> Coupled part of the plant
%       sysP      -> Deterministic interconnected part of the plant
%       Kd        -> Decoupled part of the controller
%       Kc        -> Coupled part of the controller
%       L0        -> Laplacian describing the nominal graph of the system
%       p         -> Probability of a successful package transmission
%       symmetric -> [Optional] Handle symmetric or asymmetric loss
%   Returns:
%       H2           -> H2-norm of the system
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

% Set defaults
if nargin <= 7
    symmetric = false;
end

%% Prepare the system description
[AKc, BKc, CKc, DKc] = ssdata(Kc);
[AKd, BKd, CKd, DKd] = ssdata(Kd);

nx = order(sysC) + size(AKc,1);
ny = size(BKc,2);
nu = size(CKc,1);
nw = size(sysC,2) - nu;
nL = size(L0,1);

hat = @(L, K, M) kron(eye(nL-1), K) + kron(L, M);
assemble_controller = @(L) ss(hat(L, AKd, AKc), hat(L, BKd, BKc),...
                              hat(L, CKd, CKc), hat(L, DKd, DKc), 1);

% Calculate adjacency matrix and truncated diagonalizing transformation
A0 = diag(diag(L0)) - L0;
[T, ~] = eig(L0);
T  = T(:,2:end);
Lr = T'*L0*T;

%% Split network into separete channels
if symmetric
    A0 = tril(A0);
end

channels = find(A0);
m = length(channels);

%% Solve the SDP from Theorem 1
% To generate the required SDP, we need to enumerate all modes of the jump
% system. This is done by splitting the communication into channels and
% disableing them individually.

Q   = sdpvar(nx*(nL-1));
Z   = sdpvar(nw*(nL-1));
LMI = -Q;
TRC = -Z;

for i = 1:2^m
    % Calculate new adjacency matrix with certain channels disabled
    A       = A0;
    mask    = channels(dec2bin(i-1, m) ~= '0');
    A(mask) = 0;
    if symmetric
        A = A + A';
    end
    L = T'*(diag(sum(A, 2)) - A)*T;
    
    % Assemble the closed loop with the centre of gravity removed
    sys = assemble_mas(sysD, sysC, sysP, ny, nu, L, Lr);
    K   = assemble_controller(L);
    CL  = lft(sys, K);
    [Acl, Bcl, Ccl, Dcl] = ssdata(CL);
    
    loss   = length(mask);
    chance = p^(m-loss)*(1-p)^loss;
    LMI    = LMI + chance * (Acl'*Q*Acl + Ccl'*Ccl);
    TRC    = TRC + chance * (Bcl'*Q*Bcl + Dcl'*Dcl);
end

% For certain LMIs that are obviously infeasible, Yalmip will refuse to
% construct the SDP and issue an error. This try catch will convert that 
% error into a warning so that we can successfully finish running the
% remainder of the script.
try
    Constraints = [ Q   >=  offset * eye(size(Q))   ;
                    LMI <= -offset * eye(size(LMI)) ;
                    TRC <= -offset * eye(size(TRC)) ];
catch ME
    warning(ME.message)
    Q = [];
    solver_stats.total = toc(timer);
    return
end
            
solver_stats.prep = toc(timer);
sol = optimize(Constraints, trace(Z), opts);

if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
    Q = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt(value(trace(Z)));
    
    Q  = value(Q);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
end

solver_stats.total = toc(timer);
end
