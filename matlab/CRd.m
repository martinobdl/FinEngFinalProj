function d_CRstd = CramerRao_d(defaultBarrier_vect,correlation_vect, nObligors)
% d_CRstd = CRAMERRAO_D(defaultBarrier_vect,correlation_vect, nObligors)
%
% Computes the Cramer-Rao standard deviation limit for the d default
% barrier threshold in the points specified by defaultBarrier_vect,
% correlation_vect 
n = length(correlation_vect);
m = length(defaultBarrier_vect);
CR = zeros(n,m);
f = waitbar(0,'');
for i = 1:n
    for j= 1:m
      
        CR(i,j)=ExpectedSecondDerivativeLL(correlation_vect(i), defaultBarrier_vect(j), nObligors);
        waitbar((i*j)/(n*m),f,sprintf('Please wait %2.2f%%',(i*j*100)/(n*m)))

    end
end
delete(f)
d_CRstd = sqrt(-CR.^(-1));
end

    