function tsino_nlm = tsino_nonlocal_medfilt(tsino, sz)
% TSINO_NONLOCAL_MEDFILT  Applies non-local median filter to transmission
% sinogram.
%
% tsino_nlm = TSINO_NONLOCAL_MEDFILT(tsino, sz) applies algorithm 3 from
% [1] to a transmission sinogram. The input tsino must be a parallel beam
% transmission sinogram; the rows of tsino are the projections.
%
% References: 
%
% [1] Vo, Atwood, and Drakopoulos, "Superior techniques for eliminating 
%     ring artifacts in X-ray micro-tomography," Optics Express, 2018.
%     DOI: 10.1364/OE.26.028396

assert(isnumeric(sz) && sz > 1 && floor(sz) == sz, 'sz must be a positive integer')
tsino = tsino';
[p,r] = size(tsino);

tsino_nlm = tsino;
tsino_idx = zeros(p,r);

% Sort columns
for k = 1:r
    [tmp,idx] = sort(tsino(:,k));
    tsino_idx(:,k) = idx;
    tsino_nlm(:,k) = tmp;
end

% Apply median filter to each row of sorted transmission sinogram
for k = 1:p
    tsino_nlm(k,:) = medfilt1(tsino_nlm(k,:), sz);
end

% Reverse sort
for k = 1:r
    tsino_nlm(tsino_idx(:,k),k) = tsino_nlm(:,k);
end
tsino_nlm = tsino_nlm';
end
