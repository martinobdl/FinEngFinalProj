function ccd_normalized = bayesianPrediction(x_vect,param_vect,...
                                             posterior,model_pdf)
% X = BAYESIANPREDICTION(x_vect,param_vect,posterior,model_pdf)
%
% Computes numerically the class conditional denisty (ccd) given a
% posterior and a model pdf.
%
% @inputs:        - x_vect:     vector of where the ccd will be defined
%                 - param_vect: vector of where the posterior is defined
%                 - posterior:  values of the posterior (same dim as
%                               param_vect)
%                 - model_pdf:  handle function of the pdf
% @outputs:       - ccd_normalized: vector class conditional density
%                               defined on x_vect
% 

n   = length(x_vect);
ccd = zeros(n,1);

for i = 1:n
    ccd(i) = trapz(param_vect,model_pdf(x_vect(i)).*posterior);
end

ccd_normalized = ccd./trapz(x_vect,ccd);

end

    