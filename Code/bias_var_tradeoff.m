clear;

% Bias-Variance Tradeoff
% Imperfect Matching Paper w/Bo Honoré
% Rachel Anderson, 2019-07-30

%% Settings

n = 100;          % sample size
L = 2;          % number of matches

% true X_il ~ N(mu, sig^2)
mu = 0;
sig = 1;

% false X_il ~ N(kappa, omeg^2)
kappa = 1;
omeg = 2;

% range of p1 values to try -- true and beliefs
minP = 0.1;
maxP = 0.9;
numP = 9;

% range of P and a1 to try
trueP = linspace(minP, maxP, numP);
a1List = linspace(0, 1, 10);

%% Calculate variance of X1 and X2

varX1 = trueP*sig^2 + (1.-trueP)*omeg^2 + trueP.*(1.-trueP)*(kappa-mu)^2;
varX2 = (1.-trueP)*sig^2 + trueP*omeg^2 + trueP.*(1.-trueP)*(kappa-mu)^2;

%% Functions (so far for L=2, but to be made for arbitrary L)

% calculates bias of muHat given choices of a1, a2
calc_bias = @(a, p1) a(1)*(p1*mu + (1-p1)*kappa) + a(2)*((1-p1)*mu + p1*kappa) - a(3)*kappa;

% calcA in other file -- calculates optimal A given a1, believs piHat
% calc_var in other file -- calculates variance of muHat given a2, a3, p

%% Plot results

formatSpec = 'true p = %f';
fig1 = figure;
MSE = zeros(numP, numP);
varComp = zeros(numP);
for t=1:numP            % t indicates true prob
    subplot(3,3,t);
    bias = zeros(numP,1);
    vars = zeros(numP,1);
    varComp(t) = calc_var([1,1,1], varX1(t), varX2(t)); % variance from setting a = [1 1 1]
    for i=1:numP        % i indicates believed prob
        a = calcA(trueP(i),varX1(i), varX2(i));
        bias(i) = calc_bias(a, trueP(t));
        vars(i) = calc_var(a, varX1(t), varX2(t));
    end
    MSE(t,:) = bias.^2 + vars;
    hold on
    title(sprintf(formatSpec, round(trueP(t),1)))
    plot(trueP, bias)
    plot(trueP, vars)
    plot([trueP(1), trueP(numP)], [varComp(t), varComp(t)], '--k') 
    plot([0.5, 0.5], [min(bias), max(vars)], '--k')
    hold off
end
legend('Location','southoutside', 'bias','variance', 'benchmark a = (1,1,1)')
saveas(fig1,'bias_var_tradeoff.png');

%% DO MSE 
fig2 = figure
for t=1:numP            % t indicates true prob
    subplot(3,3,t);
    hold on
    title(sprintf(formatSpec, round(trueP(t),1)))
    plot(trueP, MSE(t,:))
    plot([trueP(1), trueP(numP)], [varComp(t), varComp(t)], '--k') 
    hold off
end
legend('Location','southoutside', 'implied MSE', 'benchmark MSE for a = (1,1,1)')
saveas(fig2,'bias_var_tradeoff.png');