function correlation = correlationFromBasel2(defaultRate)
% correlation = correlationFromBasel2(defaultRate)
% Computes the correlation as a function of the defualt Rate as in Basel II
% requiremnts.
%
% @inputs: defaultRate.
%
% @outputs: correlation.
%

correlation = 0.12*(1-exp(-50*defaultRate))/(1-exp(-50))+...
    +0.24*(1-(1-exp(-50*defaultRate))/(1-exp(-50)));

end

