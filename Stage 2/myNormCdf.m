function p = myNormCdf(x)
    %MYNORMCDF Compute the standard normal CDF using the error function.
    p = 0.5*(1 + erf(x./sqrt(2)));
end
