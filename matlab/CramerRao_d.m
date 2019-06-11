function d_CRstd = CramerRao_d(correlation, defaultBarrier, nObligors,T)
% expectedLL = CRAMERRAO_D(correlation, defaultBarrier, nObligors,T)
%
% Computes the Cramer-Rao standard deviation limit for the d default
% barrier threshold in the points specified by defaultBarrier, correlation
%
% @inputs:       - correlation: scalar or Ax1x1 matrix
%                - defaultBarrier: scalar or 1xBx1 matrix
%                - nObligors: scalar
%                - T: time parameter
% @outputs       - d_CRstd: scalar or AxB matrix, each row corresponds to a
%                   correlation, each column to a value of defaultBarrier

epsilon = 1e-6;                        % increment of independent variable

X(1,1,1:nObligors+1) = (0:nObligors);  % matrix 1x1xC number of defaults
binomialcoef = factorial(nObligors)/...% matrix 1x1xC binomial coefficients
    (factorial(X).*factorial(nObligors-X));

p = @(s,y) normcdf((defaultBarrier+s*epsilon-sqrt(correlation)*y)./...
            (sqrt(1-correlation)));    % AxBx1 matrix of handle functions
                                       % Probabilities of 1 default given y
integrand = @(s,y) (p(s,y).^X).*(1-p(s,y)).^(nObligors-X).*binomialcoef;
                                       % AxBxC matrix of handle functions
                                       % Probabilities X defaults given y
P = @(s) integral(@(y) normpdf(y).*integrand(s,y),-30,30,...
            'ArrayValued',true);       % AxBxC matrix of handle functions
                                       % Probabilities of X defaults

% AxBxC matrices of probabilities, where each dimension corresponds to
% A-> correlation, B-> defaultBarriers, C-> number of defaults X
P_central = P(0);    % without change on defaultBarrier
P_plus = P(1);       % with indrement on defaultBarrier of epsilon
P_minus = P(-1);     % with decrement on defaultBarrier of -epsilon

% Logarithm of the above matrices -> LogLikelihood
LL_central = log(P_central);
LL_plus = log(P_plus);
LL_minus = log(P_minus);

secondDerivativeLL = (LL_plus-2.*LL_central+LL_minus)./epsilon^2;

if sum(sum(sum(isnan(secondDerivativeLL))))>0
    error("loglikelihood value goes to infinity, choose other input value")
end

% AxBx1 matrix of expected second derivative

expected = sum(secondDerivativeLL.*P_central,3);
% Cramer-Rao limit, standard deviation

d_CRstd = sqrt(-expected.^-1./T);

end
