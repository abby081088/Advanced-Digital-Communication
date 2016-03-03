function [out] = decoder_2belief(degree,in)

%This is the check node

out=ones(degree,1);

for ii = 1:degree
    for jj=1:degree
        if(ii~=jj)
            out(ii) = out(ii)*in(jj); %product for the indexes different of ii
        end
    end
end
