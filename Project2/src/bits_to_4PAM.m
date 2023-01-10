function X =bits_to_4PAM(b)

	for k = 1:2:length(b)
		if(b(k) == 0 && b(k+1) == 0)
			X(k) = 3;
		elseif(b(k) == 0 && b(k+1) == 1)
			X(k) = 1;
		elseif(b(k) == 1 && b(k+1) == 1)
			X(k) = -1;
		else
			X(k) = -3;
		end
	end
end