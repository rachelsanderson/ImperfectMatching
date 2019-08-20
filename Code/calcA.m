%given a1 and believes pi, returns a2 and a3 that give unbiased muHat

function a = calcA(p, varX1, varX2)
    A = fmincon(@(X) varX1*(X(1)^2) + varX2*(X(2)^2), [1;1], [], [], [p 1-p], [1]);
    a = [A; p*A(2)+(1-p)*A(1)];
end
