function x_distr = bayesianPrediction(x_vect,param_vect,param_distr,handlef)

n = length(x_vect);
x_dd = zeros(n,1);
for i = 1:n
    x_dd(i) = trapz(param_vect,handlef(x_vect(i)).*param_distr);
end

x_distr = x_dd./trapz(x_vect,x_dd);

end

    