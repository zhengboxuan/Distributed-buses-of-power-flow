
function res = solveGS_Distribute(mpc,Method)
define_constants;
mpopts = mpoption;
mpopts.out.all = 0;
mpopts.verbose = 0;

mpopts.pf.gs.max_it = 30;
mpopts.pf.tol = 1e-6;
mpopts.pf.alg = 'GS';
mpopts.pf.v_cartesian = 0;
mpc = loadcase(mpc);
delta = sum(mpc.gen(:,2)) - sum(mpc.bus(:,3));
% Compute the initial participant
if Method == 1
    c_2 = mpc.gencost(:,COST);
    for k = 1:length(c_2)
        part(k,1) = c_2(k)/sum(c_2);
    end
elseif Method == 2
    PG = mpc.gen(:,PG);
    part = PG/(sum(PG));
end
mpc.gen(:,2) = mpc.gen(:,2) + delta * part;
res = runpf(mpc,mpopts);
