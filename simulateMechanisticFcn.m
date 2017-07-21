function p_hat = simulateMechanisticFcn(visContrast,params,areaIdx,bilateral)
CL = visContrast(1); CR = visContrast(2);
r = [CR CL 0 0 CR CL]' + params.baselineActivity;

if areaIdx == 0
    zeta = params.weights*r;
    pL = exp(zeta(1)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    pR = exp(zeta(2)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    pNG = 1 / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    p_hat = [pL pR pNG];
else
    %Inactivation trial
    if bilateral == 1
        areaIdx = [areaIdx-1 areaIdx];
    end
    
    r_inact = r ;
    r_inact(areaIdx) = r_inact(areaIdx) * params.InactivationScalingFactor;
    zeta = params.weights*r_inact;
    
    pL = exp(zeta(1)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    pR = exp(zeta(2)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    pNG = 1 / ( 1 + exp(zeta(1)) + exp(zeta(2)));
    p_hat = [pL pR pNG];
end
end

