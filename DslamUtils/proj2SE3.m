function [Xout] = proj2SE3(X,n ,Tt, posLndmrkSepInd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T_L_sep = 16*n*Tt;
Xposes = X(1:posLndmrkSepInd);
Xlm = X(posLndmrkSepInd+1:end);
Pvecmat = reshape(Xposes(1:T_L_sep),[],3);
%% Computing Poses loss
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
% TBD - can speedup by direct assignment and avoid copy of array
PX = zeros(size(Xposes));
for i=1:n
    Ti  = Pvecmat(:,i);
    %     tildTi = tildPveccmat(:,i);
    Tit = posvecret(Ti,1);
    Ki = numel(Ti)/16;
    for t=1:(Ki)
        % squeeze(posesarray{1}(t,:,:))  return Tit :4x4 pose at time t
        [Uit,Dit,Vit] = svd(Tit(1:3,1:3));
        if det(Dit)>0
            P_O3_Tit = Uit*eye(3)*Vit.';
        else
            [~,minind]=min(Dit(:));
            D = eye(3);
            D(minind)= -1 ;
            P_O3_Tit = Uit*D*Vit.';
        end
        
        Pit = Tit; %copy original matrix (for translation values)
        Pit(1:3,1:3) = P_O3_Tit; %replace w/ projected rotation values
        PX((1:16)+(i-1)*16*Ki+(t-1)*16)=Pit(:); % copy to output Vecpose
        % squeeze(tTarray{1}(t,:,:))return tildeTit :4x4 IMU from Tit to Titp1
        %         tildTit = squeeze(tTarray{i}(t,:,:));
        %         G = Titp1\(tildTit*Tit);
        %         Test_ind=(1:16)+(i-1)*16*Ki+(t-1)*16
        
    end
end

Xout = [PX(:); Xlm];
end




