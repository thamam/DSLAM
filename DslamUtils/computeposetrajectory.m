function [posetrajec] = computeposetrajectory(T0, Ttild, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);

%% Initializing

[tf] = size(Ttild,1)/16+1;
posetrajec = zeros(tf*16,n);
posetrajec(1:16,1:n) = T0;
%% Computing poses

for i=1:n
    Ttilde_i = Ttild(:,i);
    Titm1  = reshape(T0(:,i),4,4);
    temp=[];
    for t=2:tf
        Ttilde_itm1 = posvecret(Ttilde_i,t-1);
        Tit = Titm1/ Ttilde_itm1 ;
        temp = [temp;Tit(:)];
        Titm1=Tit;
    end
    posetrajec(17:end,i)=temp;
end
end

