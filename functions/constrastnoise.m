function CNR = constrastnoise(sROI,bROI)
% CONTRASTNOISE  Computes contrast to noise ratio
%
% CNR = constrastnoise(sROI,bROI) computes the contrast-to-noise
% ratio as defined in [1].
% 
% References: 
% 
% [1] Warr, Ametova,  Cernik, Fardell, Handschuh, Jørgensen, Papoutsellis,
%     Pasca and Withers: "Enhanced hyperspectral tomography for 
%     bioimaging by spatiospectral reconstruction", Scientifc Reports 11, 
%     20818 (2021). 
%     DOI: https://doi.org/10.1038/s41598-021-00146-4

N1 = length(sROI(:));
N2 = length(bROI(:));

fs = 1/N1 * sum (sROI(:));
fb = 1/N2 * sum (bROI(:));
sigmas = 1/(N1-1) * sum((sROI(:)-fs).^2);
sigmab = 1/(N2-1) * sum((bROI(:)-fb).^2);

CNR = 2 * abs(fs-fb)/(sigmas+sigmab);

end