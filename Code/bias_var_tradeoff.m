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

%% Gen data -- not needed for now

% trueMatch = mu + sig*randn(n,1);
% falseMatch = kappa + omeg*randn(n,1);
% trueMatchID = (rand(n,1) < p1);  % entry = 1 if X1 is true match
%
% X1 = trueMatchID.*trueMatch + (1.-trueMatchID).*falseMatch;
% X2 = (1-trueMatchID).*trueMatch + (trueMatchID).*falseMatch


%% Functions (so far for L=2, but to be made for arbitrary L)

% calculates muHat = a1/n sum X_1i + a2/n sum X_2i - a3*kappa
calc_muHat = @(a1, a2, a3, X1, X2) a1*mean(X1) + a2*mean(X2) - a3*kappa;

% calculates bias of muHat given choices of a1, a2
calc_bias = @(a1, a2, a3, p1) a1*(p1*mu + (1-p1)*kappa) + a2*((1-p1)*mu + p1*kappa) - a3*kappa;

% calcA in other file -- calculates optimal A given a1, believs piHat
% calc_var in other file -- calculates variance of muHat given a2, a3, p

%% Plot results


for a = 1:10
    a1 = a1List(a);
    formatSpec = 'true p = %f';
    fig = figure;
    for t=1:numP
        subplot(3,3,t);
        bias = zeros(numP,1);
        vars = zeros(numP,1);
        for i=1:numP
            [a2,a3]= calcA(a1, trueP(i));
            bias(i) = vpa(calc_bias(a1,a2,a3,trueP(t)));
            vars(i) = vpa(calc_var(a1, a2, trueP(t), sig, omeg, kappa, mu, n));
        end
        hold on
        title(sprintf(formatSpec, round(trueP(t),1)))
        plot(trueP, bias)
        plot(trueP, vars)
        hold off
    end
    legend('Location','bestoutside','orientation','horizontal','bias','variance')
    saveas(fig,sprintf('a1_%d.png',a));
end

%%
testP = 0.8;
title(sprintf('trueP = %d', testP))
for j = 2:length(a1List)
    subplot(3,3,j-1);
    hold on 
    for i=1:numP
        [a2,a3]= calcA(a1List(j), trueP(i));
        bias(i) = vpa(calc_bias(a1List(j),a2,a3,testP));
        vars(i) = vpa(calc_var(a1List(j), a2, testP, sig, omeg, kappa, mu, n));
    end
    plot(trueP, bias)
    plot(trueP, vars)
    hold off
end
legend('Location','bestoutside','orientation','horizontal','bias','variance')

