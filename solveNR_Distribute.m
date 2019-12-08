function [V, success, et, niter] = solveNR_Distribute(mpc,Method)
%% Newton-Raphson Power Flow with distribution
% ECE 6320 Fall 2019 project 1
% Method 1 basing on the cost function.
% Method 2 basing on the generation capacity of the buses.
tic;

%% Specify parameters
tol = 1e-6; % Convergence tolerance
max_iter = 30;

%% Load the test case
mpc = loadcase(mpc);

%% Store placeholders for various columns in the Matpower data format.
% See Appendix B of the Matpower manual.
define_constants; 

Sbase = mpc.baseMVA; % base power

nbus = size(mpc.bus,1);

%% Compute the admittance matrix
Y = makeYbus(mpc);

%% Get bus indices
pv = find(mpc.bus(:,BUS_TYPE) == PV);
pq = find(mpc.bus(:,BUS_TYPE) == PQ);
ref = find(mpc.bus(:,BUS_TYPE) == REF);

%% Newton-Raphson iteration

% Store voltage initialization
V = mpc.bus(:,VM).*exp(1i*mpc.bus(:,VA)*pi/180); % Complex voltage phasor

% Store the power injection vector
S0 = -(mpc.bus(:,PD) + 1i*mpc.bus(:,QD))/Sbase; % Store the load demands
for i=1:size(mpc.gen,1)
    S0(mpc.gen(i,GEN_BUS)) = S0(mpc.gen(i,GEN_BUS)) + mpc.gen(i,PG)/Sbase; % Store the power injections from the generator
end

delta = sum(mpc.gen(:,PG)) - sum(mpc.bus(:,PD));
% Compute the initial participant
if Method == 1
    c_2 = mpc.gencost(:,COST);
    for i = 1:length(c_2)
        part(i,1) = sum(c_2)/c_2(i);
    end
elseif Method == 2
    PG = mpc.gen(:,PG);
    part = PG/(sum(PG));
end


% Compute the initial mismatches
S = V.*conj(Y*V);
mis = [real(S([ref; pv]) + delta * part - S0([ref; pv]));
       imag(S(pq) - S0(pq));
       real(S(pq) - S0(pq))];

   
niter = 0;
while any(abs(mis) > tol) && niter < max_iter

    % Update the iteration counter
    niter = niter + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the Jacobian at this iteration:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Constructing the matrices directly using the derivative function dSbus_dV:
    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Y, V);
    j11 = real(dSbus_dVa([ref; pv], (pv)));
    j12 = real(dSbus_dVa([ref; pv], (pq)));
    j13 = real(dSbus_dVm([ref; pv], (pq)));
    j14 = part;  
    j21 = imag(dSbus_dVa((pq), (pv)));
    j22 = imag(dSbus_dVa((pq), (pq)));
    j23 = imag(dSbus_dVm((pq), (pq)));
    j24 = zeros(length(pq),1);  
    j31 = real(dSbus_dVa((pq), (pv)));
    j32 = real(dSbus_dVa((pq), (pq)));
    j33 = real(dSbus_dVm((pq), (pq)));
    j34 = zeros(length(pq),1);    
    J = [   j11 j12 j13 j14;
            j21 j22 j23 j24;
            j31 j32 j33 j34;   ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the Newton-Raphson update step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default Matlab linear algebra 
%     dx = -J \ mis;
    [L,U,P,Q] = lu(J);         
    b_hat = P*mis;
    y = L \ b_hat;
    dx_hat = U \ y;
    dx = -Q * dx_hat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the voltage phasors and loss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vmag = zeros(nbus,1);
    Vmag([ref; pv]) = abs(V([ref; pv])); % Store the voltage magnitudes that are fixed (ref and pv)
    Vmag(pq) = abs(V(pq)) + dx(length(pv)+length(pq)+1:length(pv)+2*length(pq)); % Store the pq bus voltage magnitudes
    
    Vang = zeros(nbus,1);
    Vang([pv; pq]) = angle(V([pv; pq])) + dx(1:length(pv)+length(pq));
    
    delta = delta + dx(length(pv)+2*length(pq)+1:length(pv)+2*length(pq)+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruct the voltage phasors and store in the Matpower data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V = Vmag.*exp(1i*Vang); 
    mpc.bus(:,VM) = abs(V);
    mpc.bus(:,VA) = angle(V)*180/pi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the mismatch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = V.*conj(Y*V);
    dPG = PG + part * delta - real(S([ref; pv]));
    PG = PG + dPG;
    if Method == 2
        part = PG/(sum(PG)); %(Method 3)
    end

    mis = [real(S([ref; pv]) + delta * part - S0([ref; pv]));
           imag(S(pq) - S0(pq));
           real(S(pq) - S0(pq))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
end


%% Update the success flag
if all(abs(mis) <= tol)
    success = true;
else
    success = false;
end

%% Store the elapsed time
et = toc;