function d_CRstd = CRd(defaultBarrier_vect,correlation_vect, nObligors)

n = length(correlation_vect);
m = length(defaultBarrier_vect);
CR = zeros(n,m);
f = waitbar(0,'');
for i = 1:n
    for j= 1:m
      
        CR(i,j)=ExpectedSecondDerivativeLL(correlation_vect(i), defaultBarrier_vect(j), nObligors);
    end
    waitbar(i/n,f,sprintf('Please wait %2.2f%%',i*100/n))
end

d_CRstd = sqrt(-CR.^(-1)/20);
delete(f)
end

    