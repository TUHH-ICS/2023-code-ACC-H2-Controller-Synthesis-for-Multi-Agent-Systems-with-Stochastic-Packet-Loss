%---------------------------------------------------------------------------------------------------
% For Paper
% "Distributed H2 Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner and H. Frey
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

clear

addpath('analysis', 'synthesis')
addpath('graphs', 'simulation', 'util')
addpath(genpath('mas-simulation/lib'))

%% Settings
% The script will design a lossy controller with this fixed transmission
% probability in addition to the nominal controller
p_design = 0.5;

% For both the theoretical and the empirical analysis, this is the number
% of points the probability axis will be divided into
points = 41;

% Settings for the Monte-Carlo simulation
simconf = struct;
simconf.Tf      = 1e3;  % Maximum simulation time [s]
simconf.tol     = 1e-5; % Tolerance before stopping the simulation
simconf.samples = 100;  % Number of samples per probability

%% Problem definition
m   = 1;  % Model mass [kg]
b   = 10; % Coefficient of friction [kg/s]
dT  = 1;  % Sampling time [s]
dim = 1;  % Problem dimension

% Discretized state-space model of a mass with friction
A = [ 0   1   ;
      0  -b/m ];
B = [ 0 ; 1/m ];
C = [ 1  0 ];
G = c2d(ss(kron(eye(dim), A), kron(eye(dim), B), kron(eye(dim), C), 0), dT);

% Controller tuning
R = kron(eye(dim), 1);

% Communicaton structure
graph = line_graph(6, 2, false);

% Additional properties of the communication links
range = Inf;  % Maximum range of communication
sym   = false; % If set to true, failure is always symmetric

%% Prepare synthesis and analysis steps
% Generate dependend settings
N = height(graph.Nodes);
simconf.R         = R;   % Penalty on control signal
simconf.dT        = dT;  % Sampling time during simulation
simconf.positions = [1:N; zeros(dim-1,N)]; % Initial positions of the agents

% Fill probability grid
p_swp  = linspace(0, 1, points)';

% Assemble the generalized plant
[sysD, sysC, sysP, ny, nu] = prepare_generalized_plant(G, R);

% Calculate graph Laplacian
L0 = full(laplace_matrix(graph));

figure()
for j = 1:2
    lossy = j == 2;
    
    %% Controller synthesis
    tic
    if lossy
        disp('Phase 2 of 2, lossy design')
        [Kd, Kc] = h2syn_lossy(sysD, sysC, sysP, ny, nu, L0, p_design);
    else
        disp('Phase 1 of 2, nominal design')
        [Kd, Kc] = h2syn_nominal_dual(sysD, addparts(sysC, sysP, 1), ny, nu, L0);
    end
    disp(['Controller synthesis completed in ' format_duration(toc)])

    % Save controller for the Monte-Carlo simulation
    save('controller.mat', 'dT', 'm', 'b', 'Kd', 'Kc')
    
    %% System analysis
    % Calculate closed-loop system
    [clD, clC, clP] = assemble_cl(sysD, sysC, sysP, Kd, Kc);
    
    % Allocate storage
    H2_swp = zeros(size(p_swp));
    H2_dec = zeros(size(p_swp));
    
    % Sweep over the probability range to analyse the system performance
    tic
    parfor i = 1:length(p_swp)
        H2_swp(i) = h2norm_enumerated(sysD, sysC, sysP, Kd, Kc, L0, p_swp(i), sym);
        H2_dec(i) = h2norm_decomposed(clD, clC, clP, L0, p_swp(i));
    end
    disp(['System analysis completed in ' format_duration(toc)])
    
    %% Monte-Carlo validation
    netconf = struct('p', num2cell(p_swp), 'range', range, 'sym', sym);
        
    tic
    [H2_mtc, converged] = mc_simulate(graph, simconf, netconf);
    disp(['Monte-Carlo simulation completed in ' format_duration(toc)])
    
    % Remove values that have not converged yet
    H2_mtc(~converged) = NaN;
    
    %% Visualize results
    plot(p_swp, H2_swp, p_swp, H2_dec, '-.', p_swp, H2_mtc, '--')
    hold on
    
    %% Export results
    if lossy
        name = sprintf('bernoulli_lossy_%g_%d.csv', p_design, uint32(posixtime(datetime())));
    else
        name = sprintf('bernoulli_nominal_%d.csv', uint32(posixtime(datetime())));
    end
    tbl = table(p_swp, H2_swp, H2_dec, H2_mtc);
    writetable(tbl, name)
end

% Add formating to the figure
hold off
xlim([0,1])
ylim padded

xlabel('Probability p')
ylabel('H_2-Performance')
title('Transmission Probability Sweep')
legend('Nominal Analysis', 'Nominal Decoupled', 'Nominal Empirical',...
       sprintf('Lossy (p=%g) Analysis', p_design),...
       sprintf('Lossy (p=%g) Decoupled', p_design),...
       sprintf('Lossy (p=%g) Empirical', p_design))
