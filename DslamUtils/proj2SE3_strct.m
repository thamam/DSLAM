function [Xout] = proj2SE3_strct(XT,n,tf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Computing Poses loss
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
% TBD - can speedup by direct assignment and avoid copy of array
PX = zeros(size(XT));
for i=1:n
    Ti  = XT(:,i);
    %     tildTi = tildPveccmat(:,i);
    for t=1:tf
        Tit = posvecret(Ti,t);
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
        PX((1:16)+(t-1)*16,i)=Pit(:); % copy to output Vecpose
        % squeeze(tTarray{1}(t,:,:))return tildeTit :4x4 IMU from Tit to Titp1
        %         tildTit = squeeze(tTarray{i}(t,:,:));
        %         G = Titp1\(tildTit*Tit);
        %         Test_ind=(1:16)+(i-1)*16*Ki+(t-1)*16
        
    end
end

Xout = PX;
end




