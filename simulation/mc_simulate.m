%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [h2, converged] = mc_simulate(graph, simconf, netconf)
%MC_SIMULATE Perform a Monte-Carlo simulation of the H2-norm of the
%network of connected masses
%   This function repeatedly runs simulations of networks of masses trying
%   to maintain a formation under disturbances. By varying the properties
%   of the network, this can be exploited to test controllers under varying
%   circumstances.
%
%   Arguments:
%       graph   -> Nominal communication graph
%       simconf -> Configuration for the simulation
%       netconf -> Might be struct array, configuration for network. Either
%                  SinrConfiguration or plain struct for Bernoulli
%   Returns:
%       h2        -> Mean H2-norm of the system for each type of network
%       converged -> True, if all simulations converged

% Extract problem data
N   = height(graph.Nodes);
dim = size(simconf.positions, 1);

% Calculate all required graph matrices
L0   = full(laplace_matrix(graph));
A0   = full(adjacency(graph));
Ld   = kron(L0, eye(dim));
Rhat = kron(eye(N), sqrtm(simconf.R));
Pi   = kron(eye(N) - ones(N)/N, eye(dim));

% Allocate storage for the simulation parameters and results
n = 1:simconf.samples;
m = 1:N*dim;
[Net, M, ~] = meshgrid(netconf, m, n);
perf = zeros(size(M));
conv = zeros(size(M), 'logical'); % Convergence flag

parfor i = 1:numel(M)
    % Initialize the network
    if isa(Net, 'SinrConfiguration')
        Network = SinrNetwork(Net(i));
    else
        Network = BernoulliNetwork(N, simconf.dT, dim,...
                                    Net(i).range, Net(i).p, Net(i).sym);
    end

    % Initialize the agents
    Agents = cell(N, 1);
    for j = 1:length(Agents)
        pos    = simconf.positions(:,j);
        nghbrs = find(A0(j,:)); % All nonzero entries numbered
        dist   = (mod(M(i)-1, dim)+1) * (mod(M(i)-1, N) == j-1);
        Agents{j} = DisturbedAgent(Network.getId(), pos, pos, nghbrs, dist);
    end
    Agents = [Agents{:}];
   
    % Initialize simulation
    sim   = SimulationManager(Network, Agents);
    steps = sim.estimateSteps(simconf.Tf);
    leech = DataLeech(Agents, steps, 'position', 'ref', 'u');

    % Simulate!
    t = sim.step();
    leech.save(t)
    while formation_error(Agents) > simconf.tol && t < simconf.Tf
        t = sim.step();
        leech.save(t)
    end
    
    % Deallocate network object. Required for SinrNetwork
    delete(Network)
    
    % Check whether the simulation converged prematurely
    conv(i) = t < simconf.Tf;
   
    % Evaluation
    pos = leech.data.position;
    ref = leech.data.ref;
    u   = leech.data.u;
    err = zeros(steps, dim*N);
    ctr = zeros(steps, dim*N);
    for j = 1:steps
        err(j,:) = Ld*(pos(j,:) - ref(j,:))';
        ctr(j,:) = Pi*Rhat*u(j,:)';
    end
    
    % Exclude first step, since that will only contain the disturbance
    perf(i) = sum(err.^2, 'all') + sum(ctr(2:end,:).^2, 'all');
end

%% Evaluate
% Check whether the simulations for each netconf converged
converged = all(conv, [1,3]);

% Calculate mean performance for each netconf
h2 = sqrt(sum(perf, [1,3]) / simconf.samples);
h2 = reshape(h2, size(netconf));
end

%% Helper functions
function err = formation_error(Agents)
    pos = [Agents.position];
    ref = [Agents.ref];
    err = norm(pos - mean(pos - ref,2) - ref, 'fro');
end
