function [out] = decoder_2belief2(degree,in)
  
%This is the check node

for ii = 1:degree    
    ind = find([1:degree]~=ii);
    out(ii) = prod(in(ind));
end
  
