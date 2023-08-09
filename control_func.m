% Function containing example implementation of control methods in
% "Foresight and relaxation enable efficient control of nonlinear complex systems", 
% by Xin Zhang, Xiaozhu Zhang, Gang Yan and Jack Murdoch Moore,
% currently under review at Physical Review Reesearch
% 
% This function is used by the script "gen_data_comparison"
%
% Jack Moore and Xin Zhang, August 2023
% 

function control_func(ii0, ii_f, Xs, strategyList, evol_dyn, jac_mat, controlNodes, inputStructure)

sysStr = inputStructure.sysStr;%Type of dynamical system we are investigating

opt = inputStructure.opt;%ODE solver options
tol2 = inputStructure.tol2;%Tolerance for return point
tol3 = inputStructure.tol3;%Used when defining vectors with the colon operator (to ensure number of elements of resulting vector is predictable)
tol4 = inputStructure.tol4;%Tolerance for detecting whether system is stuck
randSeed = inputStructure.randSeed;
dx = inputStructure.dx;%Each step, aim to move dx closer (strategy 46) or move distance dx (strategy 43).
deList = inputStructure.deList;%Each step, aim to use energy de (strategy 45).
dt = inputStructure.dt;%Each iteration lasts time dt (strategies 45, 46)
T0 = inputStructure.T0;%Estimated time to reach final target (strategy 43)
dt00 = inputStructure.dt00;%Sampling time (strategy 43)
dtfMin = inputStructure.dtfMin;%Minimum allowed time increment
maxNumIter = inputStructure.maxNumIter;%Maximum allowed number of iterations
rList = inputStructure.rList;%Perturbation added to Gramian matrix when calculating v_1
DtDistExp = inputStructure.DtDistExp;%Exponent with which Dt depends on distance from final target (strategy 43)
pertDistExp = inputStructure.pertDistExp;%Exponent with which perturbation depends on distance from final target

k = 10;%Number of samples per time step (probably no longer relevant)

rng(randSeed);%Seed for repeatibility

m = size(Xs,1);
x0 = Xs(:,ii0);
xf = Xs(:,ii_f);
t0 = 0;
dist0 = norm(xf - x0);%Distance between initial point and final target

if isnan(T0)
    T0 = round((dist0/dx)*dt/dt00)*dt00;%Choose T0 according to distance between initial point and final target; make T0 divisible by dt00
elseif (T0 == Inf)
    T0 = round(sqrt(dist0/dx)*dt/dt00)*dt00;%Choose T0 according to distance between initial point and final target; make T0 divisible by dt00
end

controlDimension = numel(controlNodes);%Number of nodes controlled

%Make control matrix B
B = zeros(m, controlDimension);
for j = 1:controlDimension
    B(controlNodes(j), j) = 1;
end

num_r = numel(rList);

x_totalCell = cell(num_r, 1);
t_totalCell = cell(num_r, 1);
d_totalCell = cell(num_r, 1);
s_totalCell = cell(num_r, 1);
e_totalCell = cell(num_r, 1);

for ii_r = 1:num_r

r = rList(ii_r);

strategy = strategyList(ii_r);

de = deList(ii_r);

% % Attempt control with the chosen control set
% Initial time and initial conditions
t_1 = t0;
x_t = x0;%Set current point to initial point
iter = 0;
s_1 = 0;
e_1 = 0;

t_total = [t0];
tList = [];
dist_t = norm(xf - x_t);%Distance between current point and final target
x_total = [x0];%States
s_total = [0];%Arc lengths
e_total = [e_1];%Control energies
d_total = [dist_t];%Distances from final target

% Mode switching (switching control strategy)
prevStateList = x0;%The list of previous states

closestReturn = Inf;%Closest return to a previous state (for detecting recurrence)

while ((dist_t > tol4) && (~isnan(dist_t)) && (iter < maxNumIter)) && (closestReturn >= tol2)%Keep incrementing time while the number of iterations is less than the maximum allowed and the distance between current point and final target is greater than threshold
    iter = iter + 1;
    
    tList = [tList, t_1];
    
    A = jac_mat(x_t); % Matrix A in linearisation dx/dt = Ax + f + Bu
    f = evol_dyn(x_t) - A*x_t; % Constant vector f in linearisation dx/dt = Ax + f + Bu

    dist2 = norm(x_t - xf);%Error (distance between current point x_1 and intermediate target xf)
    
    %Work out appropriate initial conditions for DE which implements control (strategies 43, 45, 46) and time increment dtf (strategy 43): 
    
    if ismember(strategy, [45, 46])%EILOCS, DDLOCS
        dtf = dt;%Time increment

        W_t = gram_list(A, B, dtf, dtf + tol3); % Gramian at times 0:dt:(t_f-t_0)
        W_tf = W_t(:, :, end);%Gramian at final time
        W_tf = 0.5*(W_tf + W_tf');
        try
            [V, D] = eig(W_tf);
        catch%Numerical problem
            V = NaN(size(W_tf)); D = diag(NaN(m, 1));
        end
        eigList = diag(D);
        eig2List = eigList + r*2*(dist2/dist0)^pertDistExp*mean(eigList);%Perturb Gramian (i.e., use relaxation)
        D2 = diag(eig2List);
        W_tf2 = V*D2*V'; W_tf2 = 0.5*(W_tf2 + W_tf2');
        g_t = free_state(A, x_t, f, 0, dtf, dtf + tol3);% Free state (without control) under linear approximation at times 0:dt:(t_f-t_0)
        g_tf = g_t(:, end);% Free state (without control) under linear approximation at final time

        tSpan = linspace(0, dtf, k+1)';
    end
    
    if (strategy == 43)%ALITE
        %For arc length dx, follow minimum-energy trajectory
        %to final target according to local linearisation
        %
        %If there is a problem calculating minimum energy control
        %input at current location then stop
        
        %If (pertDistExp == 0.5), then target remaining time is proportional to square root of
        %remaining distance 
        
        DT = (dist2/dist0)^DtDistExp*T0;
        
        W_t = gram_list(A, B, DT, 1.5*DT);
        W_T = W_t(:, :, end); %Gramian after time DT
        W_T = 0.5*(W_T + W_T');
        try
            [V, D] = eig(W_T);
        catch%Numerical problem
            V = NaN(size(W_T)); D = diag(NaN(m, 1));
        end
        eigList = diag(D);
        eig2List = eigList + r*2*(dist2/dist0)^pertDistExp*mean(eigList);
        D2 = diag(eig2List);
        W_T2 = V*D2*V'; W_T2 = 0.5*(W_T2 + W_T2');
        g_t = free_state(A, x_t, f, 0, DT, 1.5*DT);
        g_T = g_t(:,end); %Free state after time DT
        dtList = (0:dt00:max(1.5*dt00, (DT + tol3)))';
        ode_fun_lin = @(t, xvse) lin_dyn_with_cont_fun_2(xvse(1:m), xvse((m + 1):(2*m)), f, A, B);
        v_1 = expm(A'*DT)/(W_T2)*(xf - g_T);%Calculate initial condition for ODE which implements control
        if ~any(isnan(v_1))
            [dtList, XVSELin] = ode45(ode_fun_lin, dtList, [x_t; v_1; 0; 0], opt);
            ssL = XVSELin(:, 2*m + 1)';%Distance travelled after each time for linearised system
            [ssL2, ind] = unique(ssL);
            keepIndBool = (isfinite(ssL2) & ~isnan(ssL2));
            ind = ind(keepIndBool); ssL2 = ssL2(keepIndBool);
            dtList2 = dtList(ind);
            if (ssL2(end) < dx)
                dtf = dtList2(end);
            else
                dtf = interp1(ssL2, dtList2, dx);%Time required to use energy for linear system
            end
            if ~isempty(dtf) && (dtf >= dtfMin)%Time increments which are too small could cause numerical problems
                tSpan = linspace(0, dtf, k + 1)';
            else
                v_1 = NaN*v_1;
            end

        end

        if any(isnan(v_1))
            v_1 = NaN(m, 1);
            dtf = 0;
            tSpan = linspace(0, dtf, k+1)';
        end

    elseif (strategy == 45)%EILOCS: Use (at most) energy de in moving, according to local linearisation, directly towards final target xf from destination g_tf under free evolution

        def = (xf - g_tf)'/(W_tf2)*(xf - g_tf);%Energy required to reach final target xf in time dt
        if (def <= de)
            x_tf = xf;%Immediate target
        else
            x_tf = g_tf + sqrt(de/def)*(xf - g_tf);%Immediate target
        end
        v_1 = expm(A'*dtf)/(W_tf2)*(x_tf - g_tf);%Calculate initial condition for ODE which implements control

    elseif (strategy == 46)%DDLOCS: Move distance dx directly towards final towards final target xf according to local linearisation

        dxf = norm(xf - x_t);%Distance from final target xf
        if (dxf <= dx)
            x_tf = xf;%Immediate target
        else
            x_tf = x_t + dx*(xf - x_t)/norm(xf - x_t);%Immediate target
        end
        v_1 = expm(A'*dtf)/(W_tf2)*(x_tf - g_tf);%Calculate initial condition for ODE which implements control

    end

    dtf = max(tSpan) - min(tSpan);%Time increment for this iteration

    %Evolution of controlled nonlinear system with selected v_1:            
    if (dtf > 0)
        ode_fun = @(t, xvse) dyn_with_cont_fun_2(xvse(1:m), xvse((m + 1):(2*m)), evol_dyn, A, B);%Nonlinear system evolution under optimal linear control solution
        [~, XVSE] = ode45(ode_fun, tSpan, [x_t; v_1; 0; 0], opt);
        if size(XVSE, 1) ~= numel(tSpan); XVSE = nan(numel(tSpan), 2*m + 2); end
    else
        XVSE = repmat([x_t', v_1', 0, 0], [k+1, 1]);
    end

    xx = XVSE(:, 1:m)';%States
    ss = XVSE(:, 2*m + 1)';%Arc lengths
    ee = XVSE(:, 2*m + 2)';%Control energies
    dd = sqrt(sum((xx - xf).^2, 1));%istances from final target

    s_total = [s_total, ss(:,2:end) + s_1]; s_1 = s_1 + ss(end); 
    e_total = [e_total, ee(:,2:end) + e_1]; e_1 = e_1 + ee(end); 
    d_total = [d_total, dd(:,2:end)];
    x_total = [x_total, xx(:, 2:end)];
    t_total = [t_total, t_1 + tSpan(2:end)'];%Small tSpan can be numerically equivalent to zero
    
    % Update data for next iteration:
    t_1 = t_total(end);
    x_t = xx(:, end);
    
    % Calculate variables to check terrmination conditions:
    dist_t = norm(x_t - xf);%Error (distance between current point x_1 and final target xf)
    closestReturn = min(sqrt(sum((prevStateList - x_t).^2, 1)))/(10*dx*dtf);%How close is the current state to the closest earlier point in the current control mode?
    prevStateList = [prevStateList, x_t];

end

x_totalCell{ii_r} = x_total;%Store record of states for control strategies currently considered
t_totalCell{ii_r} = t_total;%Store record of times for control strategies currently considered
d_totalCell{ii_r} = d_total;%Store record of distances from final target for control strategies currently considered
s_totalCell{ii_r} = s_total;%Store record of arc length for control strategies currently considered
e_totalCell{ii_r} = e_total;%Store record of control energy for control strategies currently considered

end

%Save results
saveCellPatience = {...
    'sysStr', 'controlNodes', 'strategyList',...
    'rList', 'x0', 'xf', 'x_totalCell', 't_totalCell', 'd_totalCell', 's_totalCell', 'e_totalCell', 'saveCellPatience'};
save(['comparison-', vector_to_string(rList), '_', sysStr, '_cont-', vector_to_string(controlNodes), '_strat-', vector_to_string(strategyList), '_energy-', vector_to_string(deList), '.mat'],...
    saveCellPatience{:});

end


%% Functions
%
% Calculate the trajectory of the linear system in the absence of control.
% The necessary integral
% \int_0^t e^(A (t - \tau)) d \tau
% for each t in dt:dt:T is calculated based on
% Loan (1978) Computing integrals involving the matrix exponential.
%
% Jack Murdoch Moore, December 2020
%
function II = free_state(A, x0, f, T0, dt, T)
tt = T0:dt:T;
num_t = numel(tt);
n = size(A, 1);
II = NaN(n,  num_t);
for ii_t = 1:num_t
    Dt = tt(ii_t);
    FG = expm([A, eye(n); zeros(n), zeros(n)]*Dt);
    I = FG(1:n, (n + 1):end);
    I_e = expm(A*Dt);
    II(:,ii_t) = I*f + I_e*x0;
end

end
%
% Calculate the controllability Grammian
% \int_0^t e^(A (t - \tau)) B B^T e^(A^T (t - \tau)) d \tau
% for each t in dt:dt:T .
%
% The method is based on Loan (1978) Computing integrals involving the
% matrix exponential. The code is based on the following answer of Sheng
% Cheng (16, 20 Feb 2017) in the MATLAB Answers forum:
% https://www.mathworks.com/matlabcentral/answers/286516-integral-of-controllability-gramian
%
% Jack Murdoch Moore, December 2020
%
function WW = gram_list(A, B, dt, T)
tt = 0:dt:T;
num_t = numel(tt);
n = size(A, 1);
expAdt = expm(A*dt);
expADt = eye(n);
AA = -A';%So that AA' = -A
WW = NaN(n, n, num_t);
for ii_t = 1:numel(tt)
    BB = expADt*B;
    Dt = tt(ii_t);
    expADt = expADt*expAdt;
    Q_c = BB*BB';
    FG = expm([-AA', Q_c; zeros(n), AA]*Dt); % Coming from the first equation below Eq. (1.4) of Loan (1978)
    W = FG((n + 1):end, (n + 1):end)'*FG(1:n,(n + 1):end); % From the second equation below (1.4) of Loan (1978)
    WW(:, :, ii_t) = W;
end

end

% ODE for system under linear optimal control.
%
% Output:
% dxdvds = [dx/dt; dv/dt; ds/dt] is column vector concatenating derivatives
% of state (column) vector x, column vector 
% v = e^(A' (tf - t)) (W(t_f)^-1*(x_f - g_f)
% such that B'*v is control (column) vector, and derivative ds/dt of length
% of trajectory s such that
% (ds/dt)^2 = \sum_i (dx_i/dt)^2
%
% Arguments:
% x is state (column)
% B'*v is control (column) vector
% dyn_fun is ODE defining self-dynamics
%
% Jack Murdoch Moore, June 2022
%
function dxdvdsde = dyn_with_cont_fun_2(x, v, self_dyn_fun, A, B)

dx = self_dyn_fun(x) + B*B'*v;
dv = -A'*v;
ds = sqrt(sum(dx.^2));
u = B'*v;
de = sum(u.^2);
dxdvdsde = [dx; dv; ds; de];

end

% ODE for linear system under linear optimal control.
%
% Output:
% dxdv = [dx/dt; dv/dt; ds/dt] is column vector concatenating derivatives
% of state (column) vector x, column vector 
% v = e^(A' (tf - t)) (W(t_f)^-1*(x_f - g_f)
% such that B'*v is control (column) vector, and derivative ds/dt of length
% of trajectory s such that
% (ds/dt)^2 = \sum_i (dx_i/dt)^2.
%
% Arguments:
% x is state (column)
% B'*v is control (column) vector
% dyn_fun is ODE defining self-dynamics
%
% Jack Murdoch Moore, June 2022
%
function dxdvdsde = lin_dyn_with_cont_fun_2(x, v, f, A, B)

dx = A*x + f + B*B'*v;
dv = -A'*v;
ds = sqrt(sum(dx.^2));
u = B'*v;
de = sum(u.^2);
dxdvdsde = [dx; dv; ds; de];

end

% A function to turn a vector into a string format useful for file names
% 
% Jack Murdoch Moore, June 2022
% 
function str = vector_to_string(vec)
str = num2str(vec(1));
for ii = 2:numel(vec)
    str = [str, ',', num2str(vec(ii))];
end
end