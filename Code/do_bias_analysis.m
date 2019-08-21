function f = do_bias_analysis(mu, sig, kappa, omeg)

% range of P to try, numP set  for graphing purposes
numP = 9;
trueP = linspace(0, 1, numP);

% variance of x1,x2 for different values of p
varX1 = trueP*sig^2 + (1.-trueP)*omeg^2 + trueP.*(1.-trueP)*(kappa-mu)^2;
varX2 = (1.-trueP)*sig^2 + trueP*omeg^2 + trueP.*(1.-trueP)*(kappa-mu)^2;

% define bias function with parameter values
calc_bias = @(a, p1) a(1)*(p1*mu + (1-p1)*kappa) + a(2)*((1-p1)*mu + p1*kappa) - a(3)*kappa;

% define variables for plotting
formatSpec = 'true p = %f';
fig1 = figure;
fig2 = figure;

for t=1:numP            % t indicates true prob
    bias = zeros(numP,1);
    vars = zeros(numP,1);
    
    for i=1:numP        % i indicates believed prob
        a = calcA(trueP(i),varX1(i), varX2(i));
        bias(i) = calc_bias(a, trueP(t));
        vars(i) = calc_var(a, varX1(t), varX2(t));
    end
    
    % benchmark variance and MSE for a = (1,1,1)
    varComp = calc_var([1,1,1], varX1(t), varX2(t));
    MSE = bias.^2 + vars;
    
    % plot bias/variance tradeoff
    figure(fig1);
    subplot(3,3,t);
    hold on
    title(sprintf(formatSpec, round(trueP(t),1)))
    plot(trueP, bias)
    plot(trueP, vars)
    plot([trueP(1), trueP(numP)], [varComp, varComp], '--k')
    plot([0.5, 0.5], [min(bias), max(vars)], '--k')
    hold off
    
    % plot MSE
    figure(fig2);
    subplot(3,3,t);
    hold on
    title(sprintf(formatSpec, round(trueP(t),1)))
    plot(trueP, MSE)
    plot([trueP(1), trueP(numP)], [varComp, varComp], '--k')
    hold off
end
paramSpec = '%d_%d_%d_%d';
str = sprintf(paramSpec,mu, sig, kappa, omeg);

figure(fig1);
legend('Location','southoutside', 'bias','variance', 'benchmark a = (1,1,1)')
saveas(fig1,strcat('../Figures/bias_var_tradeoff_',str,'.png'));
figure(fig2);
legend('Location','southoutside', 'implied MSE', 'benchmark MSE for a = (1,1,1)')
saveas(fig2,strcat('../Figures/mse_',str,'.png'));
end

