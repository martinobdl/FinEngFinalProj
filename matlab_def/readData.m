function data = readData(fileName)
% data = READDATA()
%
% Read the data file into a struct with year, specultaive grade default
% rate, all rates defalut rate and recovery rate.
%
% @inputs:            - fileName: file name of the file to read
% @outputs:           - data: struct with year, DG_SG, DG_AR, RR


    matrixData = dlmread(fileName,',',1,0);
    
    data.years  = matrixData(:,1);
    data.DR_SG  = matrixData(:,2);
    data.DR_AR  = matrixData(:,3);
    data.RR     = matrixData(:,4);
    
end

