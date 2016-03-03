function [out] = decoder_2(degree,in)
  
%This is the check node

tmp = prod(in);%% BPSK maps the sum into the product
for ii = 1:degree    
    out(ii) = tmp*in(ii); %in(ii)*in(ii)==1 so it is the product for everything but the one at index ii
end
  
