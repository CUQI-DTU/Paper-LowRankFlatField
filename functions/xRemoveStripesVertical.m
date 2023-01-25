function [nima]=xRemoveStripesVertical(ima,decNum,wname,sigma)

% Removes vertical stripes in image, for sinogram ring removal.
%
% This file is implemented based on code provided in 
% https://doi.org/10.1364/OE.17.008567
% by Beat Münch, Pavel Trtik, Federica Marone, and Marco Stampanoni
% 
% The code is redistributed here by permission of authors. 
% License: Apache 2.0

% wavelet decomposition
for ii=1:decNum
[ima,Ch{ii},Cv{ii},Cd{ii}]=dwt2(ima,wname);
end

% FFT transform of horizontal frequency bands
for ii=1:decNum
% FFT
fCv=fftshift(fft(Cv{ii}));
[my,mx]=size(fCv);
% damping of vertical stripe information
damp=1-exp(-[-floor(my/2):-floor(my/2)+my-1].^2/(2*sigma^2));
fCv=fCv.*repmat(damp',1,mx);

% inverse FFT
Cv{ii}=ifft(ifftshift(fCv));
end

% wavelet reconstruction
nima=ima;
for ii=decNum:-1:1
nima=nima(1:size(Ch{ii},1),1:size(Ch{ii},2));
nima=idwt2(nima,Ch{ii},Cv{ii},Cd{ii},wname);
end
return