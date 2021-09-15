function [firstobsfeatures, LnewfeatWithNoise, feathist] = ...
    updatelandmarks(Zout, feathist, L, tt, n, featNoiseSigma)
%LnewfeatWithNoise Update landmark history with new features and outputs
%noisy landmark measurement at newely observed features.
newobsfeat = unique(nonzeros(Zout(:,tt,1:n)));
firstobsfeatures = setdiff(newobsfeat,feathist);
feathist = [feathist ; firstobsfeatures];

Lnewfeat = L(:,firstobsfeatures);
LnewfeatWithNoise = Lnewfeat + featNoiseSigma^2*randn(size(Lnewfeat));

end

