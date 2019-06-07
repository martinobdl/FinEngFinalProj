function correlation = correlationFromBasel2(defaultRate)
% correlation = CORRELATIONFROMBASEL"(defaultRate)
%
% Computes the correlation as a function of the defualt rate as in Basel II
% requiremnts.
%
% @inputs:          - defaultRate: scalar
%
% @outputs:         - correlation: scalar


correlation = 0.12*(1-exp(-50*defaultRate))/(1-exp(-50))+...
            + 0.24*(1-(1-exp(-50*defaultRate))/(1-exp(-50)));

end

