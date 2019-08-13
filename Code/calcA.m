%given a1 and believes pi, returns a2 and a3 that give unbiased muHat

function [a_2,a_3] = calcA(a1,pi)
    syms a2 a3;
    eqns = [pi*a1 + (1-pi)*a2 - 1, pi*a2+(1-pi)*a1 - a3];
    [a_2, a_3] = solve(eqns,a2,a3);
end
