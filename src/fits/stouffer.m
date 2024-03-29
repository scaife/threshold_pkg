function pcomb = stouffer(p)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
    if length(p)==0
        error('pfast was passed an empty array of p-values')
        pcomb=1;
    else
        pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
    end