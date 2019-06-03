rho_vect = linspace(0.05,0.995,1980);

rho_sim = samplingFromPosterior(1e4,rho_vect,H_rho);
r=zeros(1e4*1e4,1);
d_sim=zeros(1e4*1e4,1);
for i = 1:1e4
    r(1+(i-1)*1e4:1e4*i)=rho_sim(i);
    CR_tmp = interp2(CRsurface.d,CRsurface.rho,CRsurface.surface,d_vect,rho_sim(i),'linear','extrap');
    tmp = normpdf(d_vect,0,1).*prod(normpdf(d,d_vect,CR_tmp));
    d_H = tmp./trapz(d_vect,tmp);
    d_sim(1+(i-1)*1e4:1e4*i)=samplingFromPosterior(1e4,d_vect,d_H);
end

 