function X = samplingFromPosterior(nSim,x,y)

    cdf = cumtrapz(x,y);
    [cdf_unique, idx] = unique(cdf);
    u = rand(nSim,1);
    
    X = interp1(cdf_unique,x(idx),u);
    
end

