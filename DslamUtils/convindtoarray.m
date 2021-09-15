function [aout] = convindtoarray(cin,SZ)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Test code
% cin = Zarray;
% SZ = [M,397];

%%

[ii,jj] = size(cin);

aout = zeros([SZ, jj]);


assert(ii==1);
for j=1:jj
    cij = cin{ii,j};
    [rr,cc]=size(cij);
    assert(cc==1);
    for r=1:rr
        cijrc = cij{r,cc};
        aout(1:numel(cijrc),r,j)=(cijrc(:))+1; %plus 1:landmarks number start from 0
    end
end

end
    
