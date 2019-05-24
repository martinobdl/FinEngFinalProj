function correlation = correlationFromBasel2(defaultRate)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

correlation = 0.12*(1-exp(-50*defaultRate))/(1-exp(-50))+...
    +0.24*(1-(1-exp(-50*defaultRate))/(1-exp(-50)));

end

