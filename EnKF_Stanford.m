
%   ===================================================================   %
%   Created by Joseph N. Squeo, University of Connecticut
%   September 2019
%   ===================================================================   %
%	
%   Testing Stanford EnKF formulation for 1D heat conduction case
%
%   ===================================================================   %


function [x_EnKF, x_tr] = EnKF_Stanford(x0,N,C,H,w,v,q,y_meas,sigma,tPlot,caseFolder_OF,solverName,tEnd,dt,solverRuns,T_FD) 
% OUTPUTS
%   x_EnKF: Predicted state forward in time (k+1) (with noise)
%
%   x_tr: Truth state (with noise)
%
%   MSE: mean-squred error between x_tr and x_EnKF for error quantification
%
% INPUTS
%   T0:     Initial temperature at t=0
%
%    N:     Number of cells
%
%    C:     Observation operator matrix maps the prior state vector into the 
%           vector of measurements (observables)---> y_k = C_k * x_k (nxn
%           matrix). If mapping the
%           state, then C = eye(nxn).
%    w:     Standard deviation of model noise with a Gaussian distribution
%
%    v:     Standard deviation of measurement noise with a Gaussian
%           distribution
%
%    q:     Number of ensemble members (greater q, less noise in x_EnKF
%           output but greater computational requirements)
%
%  y_meas:  Measurement matrix (n,k) ---> rows are number of states (or
%           cells), columns are number of time iterations (k). If no
%           measurement is available, use "nan"). Initiate as y_meas = nan
%           * ones(N,tEnd) then fill in availble measurements at
%           corresponding times.
%
%   sigma:  standard deviation of random sampling error noise to initialize
%           the ensemble
%
%   tplot:  time interval "k" to be plotted to show x_EnKF, x_tr and meas
%
%   caseFolder_OF: directory path to the OpenFOAM case folder
%                   i.e.) caseFolder_OF = '~/OpenFOAM/1D_Heat_Conduction/'
%
%   solverName: name of the solver to be called in OpenFOAM
%               i.e.) solverName = 'myLaplacianFoam'
%
%   tEnd:   end time
%
%   dt:     time step 
%
%   solverRuns: number of times OF solver should run per EnKF time
%               advancement
%   ===================================================================   %

% ======================== INITIALIZATION =============================   %
num_iterations = tEnd/(dt*solverRuns);
[y_rows,y_col] = size(y_meas);     % # rows, # col (nxk)
[C_rows,C_col] = size(C);          % # rows (m), # col (nx1)
x_est = x0 .* ones(N,q);       % generate ensemble for estimate state
if y_col ~= num_iterations         % y_meas must have a measurement or "nan" for each state (cell) and time iteration
    error('Must have a measurement for each individual time iteration. y_meas must have a column for each time iteration, k. Use "nan" if no measurement available.')
end

%   ==================== 1.ITERATE THROUGH TIME =======================   %

for k = 2:num_iterations + 1             % time loop (k=1 --> t=0, k=2 --> t=1*dt, k=3 --> t=2*dt)
   
    
   %  =============== 2.ITERATE THROUGH ENSEMBLE MEMBERS ==============   %
   for j = 1:q                              % loop through ensemble members one at a time
     W(:,j) = w .* x_est(:,j) .* randn(N,1);         % random Gaussian distributed noise with standard deviation input "w" (array size nx1)
     V(:,j) = v .* randn(C_rows,1);         % random Gaussian distributed noise with standard deviation input "z" (array size nx1)
     V(:,j) = v .* randn(C_rows,1) .* y_meas(:,k-1);
     samp_err(:,j) = sigma .* y_meas(:,k-1) .* randn(N,1);
     
     % ======== INITIALIZE THE ENSEMBLE WITH RANDOM SAMPLE ERROR ======== %
     if k == 2          % sample error only added during first time iteraiton
        x_est(:,j) = x_est(:,j) + samp_err(:,j);    % add sample error to each ensemble member
     end
     
     y_for(:,j) = C * x_est(:,j);           % forecast measurement (nxn * nxq = nxq matrix) 
     y(:,j) = y_meas(:,k-1) + V(:,j);         % add noise to measurements y_k,i = y_k + v_k,i (equations 3.6) (nxq matrix)
   end
   
   % Remove noise from B.C. cells (1st and last row in A matrix are B.C.s)
   W(1,:) = 0;
   W(end,:) = 0;
   V(1,:) = 0;
   V(end,:) = 0;
   
   
   %  =============== 3.AVERAGE ALL ENSEMBLE MEMBERS ==================   %
   x_estbar = mean(x_est,2);                % mean of ensemble of forecasted state (nx1)  
   y_bar = mean(y,2);                % mean of ensemble of forecasted measurement (nx1)  
   Vmean = mean(V,2);
   
   %  ================= 4.ERROR COVARIANCE MATRICES ===================   %
   P = zeros(N);
   R = zeros(C_rows);
   
   for j = 1:q
     Ex(:,j) = [x_est(:,j) - x_estbar(j)];  % ensemble error matrix --> (n x q) matrix
     P = P + Ex(:,j)*Ex(:,j)';
     
     Ey(:,j) = [y(:,j) - y_bar(j)];  % ensemble of output error matrix --> (p x q) matrix
%      Ey(:,j) = [V(:,j) - Vmean(j)];  % ensemble of output error matrix --> (p x q) matrix
     R = R + Ey(:,j)*Ey(:,j)';
   end 
   P = P./(q-1);
   R = R/(q-1);
   
   
   %  ==================== 6.KALMAN GAIN MATRIX =======================   %
%    K = P*H'/(H*P*H' + R);
   K = P*H'*pinv(H*P*H' + R);
   
   %  ======================= 7.ANALYSIS STEP =========================   %
   x_est = x_est + K * (y - y_for);         % new state estimate (nxp * pxq = nxq)
   
   % Send analsis step output to OpenFOAM and run the solver for one dt
   tFolder = write_controlDict(k,caseFolder_OF,solverName,dt,tEnd,solverRuns); % write new controlDict file to start at the current time t
   varname = 'T';
   T = zeros(N,q);    % rows are each second of time, columns are cell centers
   
   for j = 1:q   % iterate through for each ensemble member
       changeFolderOF(varname,tFolder,caseFolder_OF);      % change directories to current OF time folder to create 'T' file from analysis step
       Matlab2OF(x_est,j,varname,tFolder,x0);   % generate the T file for OpenFOAM from analysis step
       Tsolver = Matlab_callOF_myLaplacianFoam(k,varname,caseFolder_OF,solverName,dt,solverRuns);  % change to t+dt time folder, run the solver for one time step
       T(:,j) = Tsolver;    % store values from OF solver output to a matrix (each row is a new cell, each column a new ensemble)

   
%       x_est(:,j) = T(:,j) + W(:,j);  % forecast state forward in time (nxn * nxq = nxq)                      
   end
   
%  ======================= 8. FORECAST STEP ===========================   %
   x_est = T + W;   %forecast state forward in time (nxn * nxq = nxq)     
   
%  ======================= OUTPUT VARIABLES ===========================   %
    x_EnKF(:,k) = mean(x_est,2);             % EnKF state output (each row is a state (cell), each column is an iteraiton in time)
    x_tr(:,k) = mean(T,2);     % true state = OpenFOAM output without noise
    x_EnKF(:,1) = x0;                 % set the initial condition for k=1
    x_tr(:,1) = x0;
end

% ======================== POST-PROCESSING ============================   %
ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt,solverRuns,T_FD)
end
