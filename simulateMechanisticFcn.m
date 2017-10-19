function p_hat = simulateMechanisticFcn(visContrast,p,areaIdx,bilateral)
CL = visContrast(1); CR = visContrast(2);
r = [CR CL 0 0 CR CL]' + p.baseline;

weights = [p.V1i p.V1c p.S1i p.S1c p.M2i p.M2c;
           p.V1c p.V1i p.S1c p.S1i p.M2c p.M2i];

if areaIdx == 0
    zeta = weights*r;
else    %Inactivation trial
    if bilateral == 1
        areaIdx = [areaIdx-1 areaIdx];
    end
    
    r_inact = r ;
    r_inact(areaIdx) = r_inact(areaIdx) * p.inact;
    zeta = weights*r_inact;
end

pL = exp(zeta(1)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
pR = exp(zeta(2)) / ( 1 + exp(zeta(1)) + exp(zeta(2)));
pNG = 1 / ( 1 + exp(zeta(1)) + exp(zeta(2)));
p_hat = [pL pR pNG];
end

