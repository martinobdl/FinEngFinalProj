function [CapitalRequrement] = CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

CapitalRequrement = ...
    (1-recoveryRate)*normcdf((norminv(defaultRate)-sqrt(correlation)*norminv(1-confidenceLevel))/sqrt(1-correlation))+...
    -(1-recoveryRate)*defaultRate;

end

