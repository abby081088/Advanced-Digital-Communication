function [matrix] = GenerateMatrix(pulse,L_o,k0,k1,k2,m)
  
  L_f = length(pulse);
 
  matrix = zeros(L_o,k1+k2+1);
    
  tmp = pulse(1:min(L_f,L_o-(k0-1)));
  matrix(k0+[0:length(tmp)-1],k1+1) = tmp;

    
  for ii=1:k2

    tmp = pulse(1:min(L_f,L_o-(k0-1)-ii*m));
    matrix(k0+ii*m+[0:length(tmp)-1],k1+1+ii) = tmp;
 
  end
   
  
  for ii=1:k1
    
    if k0-ii*m>0
      tmp = pulse(1:min(L_f,L_o-(k0-1)+ii*m));
      
    else
      tmp = pulse( -(k0-ii*m)+2:min(L_f,L_o-(k0-ii*m)+1));
    end
    
    matrix(max(1,k0-ii*m)+[0:length(tmp)-1],k1+1-ii) = tmp;
  
  end
  
