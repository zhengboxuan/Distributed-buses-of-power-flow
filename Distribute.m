function res = Distribute(mpc, method, part)
if method == 1
    [V, success, et, nither] = solveFDPF_Distribute(mpc, part);
    res.V = V;
    res.success = success;
    res.et = et;
    res.nither = nither;
elseif method == 2
    [V, success, et, nither] = solveNR_Distribute(mpc, part);
    res.V = V;
    res.success = success;
    res.et = et;
    res.nither = nither;
else
    fprintf('fail to find such a method')
end
