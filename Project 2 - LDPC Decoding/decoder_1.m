function [out, b_hat] = decoder_1(degree,y,in)
  
% This implements the variable node
   
for ii = 1:degree %For all edges of variable node
    ind = find([1:degree]~=ii); %indices different from ii (extrinsic infos)
    
    % flip input if all other agree its wrong
    if all(in(ind)==-y)
      out(ii) = -y;
    else
      out(ii) = y;
      
    end
 end
  
% decoder output after this iteration
if all(in == -y)    
    b_hat = -y;
else    
    b_hat = y;
end
  
  
  
