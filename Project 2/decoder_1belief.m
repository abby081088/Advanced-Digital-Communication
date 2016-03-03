function [out, b_hat] = decoder_1belief(degree,y,in)

% This implements the variable node


for ii = 1:degree
    out(ii) = y;
    b_hat = y;
    if(y==0)
        for jj = 1:degree
            if(jj~=ii && in(jj)~=0)
                out(ii) =in(jj);
                b_hat = in(jj);
            end
        end
    end
end
