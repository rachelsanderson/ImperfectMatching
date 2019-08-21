function f = do_bias_analysis(correct_mean, correct_var, incorrect_mean, incorrect_var)

% range of P to try, numP set  for graphing purposes
numP = 5;
numBeliefs = 20;
trueP = linspace(0, 0.5, numP);
beliefs = linspace(0,1, numBeliefs);

% variance of x1,x2 for different values of p
varX1 = beliefs*correct_var+ (1.-beliefs)*incorrect_var+ beliefs.*(1.-beliefs)*(incorrect_mean - correct_mean)^2;
varX2 = (1.-beliefs)*correct_var + beliefs*incorrect_var+ beliefs.*(1.-beliefs)*(incorrect_mean - correct_mean)^2;

% define bias function with parameter values
calc_bias = @(a, p1) a(1)*(p1*correct_mean + (1-p1)*incorrect_mean) + a(2)*((1-p1)*correct_mean + p1*incorrect_mean) - a(3)*incorrect_mean;

% define variables for plotting
formatSpec = 'true p = %f';
fig1 = figure;

for t=1:numP            % t indicates true prob
    bias = zeros(numP,1);
    vars = zeros(numP,1);
    
    for i=1:numBeliefs        % i indicates believed prob
        a = calcA(beliefs(i),varX1(i), varX2(i));
        bias(i) = abs(calc_bias(a, trueP(t)));
        vars(i) = calc_var(a, varX1(t), varX2(t));
    end
    
    % benchmark variance and MSE for a = (1,1,1)
    varComp = calc_var([1,1,1], varX1(t), varX2(t));
    MSE = bias.^2 + vars;
    
    % make last subplot span two columns
    if t == numP
        subplot(3,2,[numP:(numP+1)])
    else
        subplot(3,2,t);
    end
    
    hold on
    title(sprintf(formatSpec, round(trueP(t),1)))
    h(1) = plot(beliefs, bias, 'r', 'DisplayName','Bias (Abs Value)');
    h(2) = plot(beliefs, vars, 'b', 'DisplayName','Variance');
    h(3) = plot(beliefs, MSE, 'g', 'DisplayName', 'MSE');
    h(4) = plot([beliefs(1), beliefs(numBeliefs)], [varComp, varComp], '--k', 'DisplayName', 'Equal weights');
    h(5) = plot([0.5, 0.5], [min(bias), max(MSE)], '--k');
    hold off
    
end

legend(h(1:4), 'Location', 'southoutside','orientation','horizontal')

titleSpec = 'Correct Distribution: (%d,%d)  Incorrect Distribution: (%d,%d)';
titleName = sprintf(titleSpec, correct_mean, correct_var, incorrect_mean, incorrect_var);
fig1.Name = titleName;

% uncomment this with Matlab 2019a
sgtitle(sprintf(titleSpec, correct_mean, correct_var, incorrect_mean, incorrect_var))

paramSpec = '%d_%d_%d_%d';
str = sprintf(paramSpec,correct_mean, correct_var, incorrect_mean, incorrect_var);
saveas(fig1,strcat('../Figures/bias_var_tradeoff_',str,'.png'));

end

