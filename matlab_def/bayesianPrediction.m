function x_distr = bayesianPrediction(x_vect,p_vect,p_distr,handlef)

n = length(x_vect);
x_dd = zeros(n,1);
for i = 1:n
    x_dd(i) = trapz(p_vect,handlef(x_vect(i)).*p_distr);
end

x_distr = x_dd./trapz(x_vect,x_dd);

end

    