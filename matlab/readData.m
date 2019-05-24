function data = readData(fileName)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: fileName: file name of the file to read
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

    matrixData = dlmread(fileName,',',1,0);
    
    data.years   = matrixData(:,1);
    data.DG_SG  = matrixData(:,2);
    data.DG_All = matrixData(:,3);
    data.RR     = matrixData(:,4);
    
end

