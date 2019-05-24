function data = readData()
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

    fileName = '../data/Data.csv';
    matrixData = dlmread(fileName,',',1,0);
    
    data.year   = matrixData(:,1);
    data.DG_SG  = matrixData(:,2);
    data.DG_All = matrixData(:,3);
    data.RR     = matrixData(:,4);
    
end

