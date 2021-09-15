function [Lgrad, Posegrad, costVal, PoseGrad_ft, posegrad_ht,ftval, grdft, hmtnrm2val, gradhArray_Tt] = symmultiValComp(X,tTldVec,tf, zid, uvi)

% X=Xc;
% tf=2;
% tTldVec = Ttild;
% zid = Z{1};
% uvi=u{1};
%%
XT = X.T;
XL = X.L;
Ind16 = reshape(1:16,4,4);
Ind12 = Ind16(1:3,1:4);
gradfInd = reshape(Ind12(:) +[0,16],[],1);
dhhdTtInd = Ind12(:);
grdft=[];
%Evaluate fts
for tt=2:tf
    Ttm1 = reshape(XT((1:16)+(tt-2)*16),4,4);
    Tt   = reshape(XT((1:16)+(tt-1)*16),4,4);
    tTtldtm1 = reshape(tTldVec((1:16)+(tt-2)*16),4,4);    
    [ftval(tt), grdft(gradfInd+(tt-2)*16,tt)] = symslam_ft(Ttm1, Tt, tTtldtm1);
end


%Evaluate hts
Lgrad = zeros(3,199);
for tt=2:tf
    Tt   = reshape(XT((1:16)+(tt-1)*16),4,4);
    % since the landmarks numbering starts at 0, we add 1 to all of them
    zmtID = zid{tt};
    %%
    % note that the data here is not indexed by m but just listed in
    % the order the landmarks id appear in zid_t
    ut = uvi{tt};
    dhdT = [];
%     dhdL = zeros(3,199);
    grdhmt = zeros(15,size(ut,1));
    for m =1:numel(zmtID)
        mid = zmtID(m)+1;
        umt = ut(m,:).';
        lmbar = [XL(:,mid);1];
        [hmtnrm2val(m,tt), grdhmt(:,m)] = symslam_ht(Tt,umt,lmbar);
    end
    dhdT(dhhdTtInd+(tt-1)*16,:) = grdhmt(4:end,:);
    gradhArray_Tt{tt} = grdhmt(4:end,:); %dhdT;
%     dhdL(1:3,zmtID+1) =   grdhmt(1:3,:);
    Lgrad(1:3,zmtID+1) = Lgrad(1:3,zmtID+1)+ grdhmt(1:3,:);
%     gradhArray_Lt{tt} = dhdL;
end


%% Compute totatl derivatives
PoseGrad_ft = sum(grdft,2);
posegrad_ht=zeros(16*tf,1);
for tt=2:tf
    posegrad_ht(dhhdTtInd +(tt-1)*16,:) = posegrad_ht(dhhdTtInd +(tt-1)*16,:) + sum(gradhArray_Tt{tt},2);
end
Posegrad = [PoseGrad_ft;0] + posegrad_ht;
costVal = sum(ftval(:)) + sum(hmtnrm2val(:));
end