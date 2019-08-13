function var_muHat = calc_var(a1, a2, p, sig, omeg, kappa, mu, n)
    varX1 = p*(sig^2-mu^2) + (1-p)*(omeg^2 - kappa^2) - (p*mu + (1-p)*kappa)^2;
    varX2 = (1-p)*(sig^2-mu^2) + p*(omeg^2 - kappa^2) - ((1-p)*mu + p*kappa)^2;
    var_muHat = ((a1^2)*varX1 + (a2^2)*varX2)/n;
end



