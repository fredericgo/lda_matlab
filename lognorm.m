function x = lognorm(v_log)
	x = v_log(1);
	for i = 2:length(v_log)
		x = logsum1d(x, v_log(i));
	end
end

function x = logsum1d(log_a, log_b) 
	if (log_a < log_b)
		x = log_b + log(1 + exp(log_a - log_b));
	else
		x = log_a + log(1 + exp(log_b - log_a));
	end
end
