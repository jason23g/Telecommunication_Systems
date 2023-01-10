function est_X = detect_4_PAM(Y,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X = [-3*A,-A,A,3*A];


 
 
for i = 1 : length(Y)
 temp(1) = abs(Y(i) - X(1));
 temp(2) = abs(Y(i) - X(2));
 temp(3) = abs(Y(i) - X(3));
 temp(4) = abs(Y(i) - X(4));
 
    tmp = temp(1);
    est_X(i) = -3*A;
    for j = 2:4
        if(tmp > temp(j))
   tmp = temp(j);
   if(j == 2)
      est_X(i) = -A;
   end
   if(j == 3)
      est_X(i) = A;

   end
   if(j == 4)
      est_X(i) = 3*A;
   end
        end
    end
    
    
    
end
end

