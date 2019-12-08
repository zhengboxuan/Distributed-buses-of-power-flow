function [V, success, et, niter] = solveFDPF_Distribute(mpc,Method)
%% Newton-Raphson Power Flow
% ECE 6320 Fall 2019 HW 3 Problem 1
%
% Name: SOLUTION

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
Pref = real(S(ref)) * part(1,1); 
Pmis = [real(S(pv) .* part(2:length(part),1) - S0(pv));
        real(S(pq) - S0(pq))];
Qmis = imag(S(pq) - S0(pq));

% Construct the Bp and Bpp matrices
[Bp,Bpp] = makeB(mpc,'FDXB');

% Restrict to rows and columns used in the Jacobian
Bp = Bp([pv;pq],[pv;pq]);
Bpp = Bpp(pq,pq);

% Compute LU factors
[Lp,Up,Pp,Qp] = lu(Bp);
[Lpp,Upp,Ppp,Qpp] = lu(Bpp);

niter = 0;
while any(abs([Pmis; Qmis]) > tol) && niter < max_iter

    % Update the iteration counter
    niter = niter + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the dVa updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve Bp*dVa = -Pmis ./ V using the precomputed LU factorization
    % dVa = -Bp \ (Pmis./abs(V([pv;pq])));
    dVa = -Qp*(Up\(Lp\(Pp*(Pmis./abs(V([pv;pq]))))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the voltage phasors with the new angles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Va = angle(V);
    Va([pv;pq]) = Va([pv;pq]) + dVa;
    
    V = abs(V).*exp(1i*Va);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the mismatch phasors with the new angles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = V.*conj(Y*V);
    Qmis = imag(S(pq) - S0(pq));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the dVm updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve Bpp*dVa = -Qmis ./ V using the precomputed LU factorization
    % dVm = -Bpp \ (Qmis./abs(V(pq)));
    dVm = -Qpp*(Upp\(Lpp\(Ppp*(Qmis./abs(V(pq))))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the voltage phasors with the new magnitudes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vm = abs(V);
    Vm(pq) = Vm(pq) + dVm;
    
    V = Vm.*exp(1i*angle(V));
    if Method == 2
        part = PG/(sum(PG)); %(Method 3)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the mismatch with the new voltage magnitudes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = V.*conj(Y*V);
    Pmis = real(S([pv; pq]) - S0([pv; pq]));
    Qmis = imag(S(pq) - S0(pq));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% Update the success flag
if all(abs([Pmis; Qmis]) < tol)
    success = true;
else
    success = false;
end

%% Store the elapsed time
et = toc;