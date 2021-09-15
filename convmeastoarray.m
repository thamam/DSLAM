function [Lout,Rout] = convmeastoarray(cin,SZ)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Test code
% cin = Zarray;
% SZ = [M,397];

%%

[ii,jj] = size(cin);

Lout = zeros([SZ, jj]);
Rout = zeros([SZ, jj]);


assert(ii==1);
for j=1:jj
    cij = cin{ii,j};
    [rr,cc]=size(cij);
    assert(cc==1);
    for r=1:rr
        cijrc = cij{r,cc};
        Lout(1:size(cijrc,1),r,j)=(cijrc(:,1));
        Rout(1:size(cijrc,1),r,j)=(cijrc(:,2));
    end
end

end
    
