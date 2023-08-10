% Script running example implementation of control methods in
% "Foresight and relaxation enable efficient control of nonlinear complex systems", 
% by Xin Zhang, Xiaozhu Zhang, Gang Yan and Jack Murdoch Moore,
% currently under review at Physical Review Research
%
% Jack Moore and Xin Zhang, August 2023
%
close all;
clear;

fprintf('Finding fixed points ');

ticTotal = tic;

% Generate network defining coupling
randSeed = 0;%Seed for repeatability
rng(randSeed);%Seed random number generator for repeatability
m = 5; p = 0.5; % Number of neurons; Edge density
L = round(m*(m-1)*p/2);%Number of links
fUnd = 0.4; %Fraction of links which are undirected
J = gen_adj_dir_undir(m, L, fUnd);%Adjacency matrix

sysStr = 'bist';%Bistable neuronal system: "Dynamics of random neural networks with bistable units" by Stern, Sompolinsky and Abbott (2014)
s = 1.5; g = 2.0;
% Governing differential equation:
evol_dyn = @(x) - x + s*tanh(x) + g*J*tanh(x);
% Jacobian matrix:
jac_mat = @(x) (s - 1)*eye(m) + g*J - (s*eye(m) + g*J)*diag(tanh(x).^2);

inputStructure.sysStr = sysStr;%Type of dynamical system we are investigating

tol1 = 10^-3; inputStructure.tol1 = tol1;%Tolerance to consider two fixed points the same; minimum error for successful control
tol2 = 10^-6; inputStructure.tol2 = tol2;%Tolerance for return point
tol3 = 10^-6; inputStructure.tol3 = tol3;%Used when defining vectors with the colon operator (to ensure number of elements of resulting vector is predictable)
tol4 = tol1; inputStructure.tol4 = tol4;%Tolerance to terminate control
inputStructure.randSeed = randSeed;

% Obtain the stable fixed point allowing system to evolve from randomly
% chosen points:

%First estimates of stable fixed points:
Xs0 = [];%List of stable fixed points
opt = odeset('RelTol', 10^-10, 'AbsTol', 10^-12);

inputStructure.opt = opt;%ODE solver options

for rnd = 1:1000
    rng(rnd);
    x_init = 10*(rand(m,1) - 0.5);
    [~, xs] = ode45(@(t,x) evol_dyn(x), 0:(10^3), x_init, opt);
    xs = xs(end, :)';
    Xs0 = [Xs0, xs];
end
[Xs0, ~, ~] = uniquetol(Xs0', 10^-2, 'ByRows', true);
Xs0 = Xs0';
%Check that Xs0 does not simply comprise instances of a fixed point at zero,
%with small numerical differences:
if (max(abs(Xs0(:))) < 10^-6) && all(evol_dyn(zeros(m, 1)) < 10^-6); Xs0 = zeros(m, 1); end

% %Refine estimates of stable fixed points:
Xs = Xs0; iterList = ones(1, size(Xs0, 2));
for ii = 1:size(Xs0, 2)
    [~, xs] = ode45(@(t,x) evol_dyn(x), 0:(10^3), Xs0(:, ii), opt);
    Xs(:, ii) = xs(end, :)';
    iterList(ii) = iterList(ii) + 1;
    while (norm(Xs(:, ii) - Xs0(:, ii)) >= 10^-9)
        Xs0(:, ii) = Xs(:, ii);
        [~, xs] = ode45(@(t,x) evol_dyn(x), 0:(10^3), Xs0(:, ii), opt);
        Xs(:, ii) = xs(end, :)';
        iterList(ii) = iterList(ii) + 1;
        if iterList(ii) >= 100
            Xs(:, ii) = NaN*Xs(:, ii);
            continue;
        end
    end
end
Xs = Xs(~isnan(sum(Xs, 2)), :);
[Xs, ~, ~] = uniquetol(Xs', tol1, 'ByRows', true);
Xs = Xs';
Xs = Xs(:,~all(isnan(Xs),1)); % remove Nan in Xs

%Check that Xs does not simply comprise instances of a fixed point at zero,
%with small numerical differences:
if (max(abs(Xs(:))) < 10^-6) && all(evol_dyn(zeros(m, 1)) < 10^-6); Xs = zeros(m, 1); end

% %Check the identified fixed points correspond to small dx/dt
iterList = ones(1, size(Xs, 2));
for ii = 1:size(Xs, 2)
    xs = Xs(:, ii);
    if (norm(evol_dyn(xs)) > 10^-6)
        Xs(:, ii) = NaN*Xs(:, ii);
    end
end
Xs = Xs(~isnan(sum(Xs, 2)), :);
[Xs, ~, ~] = uniquetol(Xs', tol1, 'ByRows', true);
Xs = Xs';
Xs = Xs(:,~all(isnan(Xs),1)); % remove Nan in Xs
%Check that Xs does not simply comprise instances of a fixed point at zero,
%with small numerical differences:
if (max(abs(Xs(:))) < 10^-6) && all(evol_dyn(zeros(m, 1)) < 10^-6); Xs = zeros(m, 1); end

num_fp = size(Xs, 2);%Number of fixed points

timeFixedPoints = toc(ticTotal);
disp(['took ', num2str(timeFixedPoints), ' s.']);

dtfMin = 0;
dx = 10^-1;
dt = 10^-2;
dt00 = 10^-2*dt;
DtDistExp = 0.5;%Exponent with which Dt depends on distance from final target (strategy 43)
pertDistExp = 1;%Exponent with which perturbation depends on distance from final target

% T0 = 10;%Fixed T0
T0 = NaN;%T0 = sqrt(dist0/dx)*dt;
% T0 = Inf;%T0 = (dist0/dx)*dt;

maxNumIter = 10000;

rng(randSeed); %Reset seed for convenience
ii0_ii_f_pairs = nchoosek(1:num_fp, 2); %Control between which pairs of fixed points?
ii0_ii_f_pairs = [ii0_ii_f_pairs; fliplr(ii0_ii_f_pairs)];
ii0_ii_f_pairs = ii0_ii_f_pairs(randperm(size(ii0_ii_f_pairs, 1)), :); %Randomly reorder pairs of fixed points
% If there are too many pairs of fixed points then we can just choose the first few
maxNumPairs = 1;
numPairs = size(ii0_ii_f_pairs, 1);
if (numPairs > maxNumPairs)
    numPairs = maxNumPairs;
    ii0_ii_f_pairs = ii0_ii_f_pairs(1:numPairs, :);
end

numPairs = size(ii0_ii_f_pairs, 1);
    
ticControlStrategy = tic;

for comparisonType = [1, 2]
    
    switch comparisonType
        case 1
            fprintf('Comparing EILOCS, DDLOCS and ALITE control strategies ');
            strategyList = [45, 46, 43];%Control strategies to consider: 43 - ALITE; 45 - EILOCS; 46 - DDLOCS
            r0 = 10^-4; rList = r0*[1, 1, 1];%Relaxation parameter r
            deList = [10^5, 1, 1];%Energy increment Delta E
            controlNodes = [3, 4];%Driver nodes/control set
        case 2
            fprintf('Comparing relaxation values for ALITE control strategy ');
            strategyList = [43, 43, 43, 43];%Control strategies to consider: 43 - ALITE; 45 - EILOCS; 46 - DDLOCS
            rList = [10^-6, 0.001, 0.1, 1];%Relaxation parameter r
            deList = [1, 1, 1, 1];%Energy increment Delta E
            controlNodes = [3, 5];%Driver nodes/control set
    end
    
    ticControl = tic;
    
    inputStructure.T0 = T0;
    inputStructure.maxNumIter = maxNumIter;%Maximum allowed number of iterations
    inputStructure.dtfMin = dtfMin;%Minimum allowed time increment
    inputStructure.rList = rList;%Perturbation added to Gramian matrix when calculating v_1 (modes 17 and 18 only)
    inputStructure.dx = dx;%Each step, aim to move dx closer.
    inputStructure.deList = deList;%Each step, aim to use energy de.
    inputStructure.dt = dt;
    inputStructure.dt00 = dt00;
    inputStructure.DtDistExp = DtDistExp;%Exponent by which time increment Dt scales with distance
    inputStructure.pertDistExp = pertDistExp;%Exponent by which perturbation scales with distance
    
    for jjPair = 1:numPairs % parfor jj = 1:numPairs
        line_para = cell(1,10);
        ii0_i_f = ii0_ii_f_pairs(jjPair, :);
        ii0 = ii0_i_f(1); ii_f = ii0_i_f(2); % Initial, final stable fixed point
        
        control_func(ii0, ii_f, Xs, strategyList, evol_dyn, jac_mat, controlNodes, inputStructure);
        
    end
    
    timeControl = toc(ticControl);
    disp(['took ', num2str(timeControl), ' s.']);
    
end

timeTotal = toc(ticTotal);
disp(['Total time: ', num2str(timeTotal), ' s.']);


%% Functions
%%%%%%%%%%%%%%%%%
function A = gen_adj_dir_undir(N, L, fUnd)

%Generate the adjacency matrix of a network with N nodes and L links. A
%fraction fU are undirected and a fraction (1 - fU) are directed.
%
%Example. Generate a network with 50 nodes and 200 links, 40% of which are
%undirected:
%N = 50; L = 200; fUnd = 0.4; A = gen_adj_dir_undir(N, L, fUnd);
%
% Jack Murdoch Moore, November 2020
assert(L <= N*(N - 1)/2, 'The number of links L must be no more than the number of node pairs N*(N-1)/2.');
A = zeros(N);%Adjacency matrix
uppTriInds = NaN(1, N);%To store indices of all node pairs
for ii = 2:N%Find indices of all node pairs in part of adjacency matrix above main diagonal
    uppTriInds((ii - 2)*(ii - 1)/2 + (1:(ii - 1))) = N*(ii - 1) + (1:(ii - 1));
end
%%A more intuitive, but possibly slower, way to find indices of all node pairs in part of adjacency matrix above main diagonal
% uppTriInds = find(triu(ones(N), 1));%Find indices of all node pairs
linkInds = datasample(uppTriInds, L, 'Replace', false);%Choose all node pairs involved in links
numUnd = round(fUnd*L);%Number of undirected links
undLinkInds = datasample(linkInds, numUnd, 'Replace', false);%Find indices of node pairs involved in undirected links
A(undLinkInds) = 1; A = A + A';%Place undirected links in matrix
dirLinkInds = setdiff(linkInds, undLinkInds);%Find indices of node pairs involved in directed links
numDir = L - numUnd;%Number of directed links
numUppTriDirLink = randi(numDir + 1) - 1;%Number of directed links i -> j for which j > i
uppTriDirLinkInds = datasample(dirLinkInds, numUppTriDirLink, 'Replace', false);%Choose directed links in the upper triangle of the adjacency matrix; node pairs i, j with j > i for which there is a directed link i -> j and for which
A(uppTriDirLinkInds) = 1;%Add directed links i -> j for which j > i
lowTriDirLinkInds = setdiff(dirLinkInds, uppTriDirLinkInds);%Choose directed links in the lower triangle of the adjacency matrix;
[lowTriDirLinkRowNum, lowTriDirLinkColNum] = ind2sub([N, N], lowTriDirLinkInds);
lowTriDirLinkInds = sub2ind([N, N], lowTriDirLinkColNum, lowTriDirLinkRowNum);
A(lowTriDirLinkInds) = 1;%Add %Choose directed links in the lower triangle of the adjacency matrix; directed links j -> i for which j > i
end

%%%%%%%%%%%%%%%%%