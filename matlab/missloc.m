N=60; 
syms d m
f=m*log(normcdf(d))+(N-m)*log(normcdf(-d))
dLL=diff(f,d,2)
p=nchoosek(N,m)*normcdf(d)^m*(1-normcdf(d))^(N-m)
ELL = subs(dLL*p, m, 0:60)
ELL_sum = sum(ELL)
figure 
fplot(ELL_sum, [-7 7])
% CR_ao=-1/ELL_sum
% fplot(CR_ao, [-4 4])