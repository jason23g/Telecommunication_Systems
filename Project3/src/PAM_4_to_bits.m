function est_X = PAM_4_to_bits(X,A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

temp = zeros(2*length(X),1);
for i=1:1:length(X)
   if(X(i) == -A)
       temp(2*i-1) = 1;
       temp(2*i) = 1;
   elseif(X(i) == -3*A)
       temp(2*i-1) = 1;
       temp(2*i) = 0;
   elseif(X(i) == 3*A)
       temp(2*i-1) = 0;
       temp(2*i) = 0;
   elseif(X(i) == A)
       temp(2*i-1)= 0;
       temp(2*i) = 1;
   end
    
end

est_X = temp;
end

