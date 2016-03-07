function opt = opt_alpha(ss, D, K)
	NEWTON_THRESH = 1e-5;
	MAX_ALPHA_ITER = 1000;

	init_a = 100;
    iter = 0;
    log_a = log(init_a);

    do
        iter++;
        a = exp(log_a);
        if (isnan(a))
            init_a = init_a * 10;
            printf("warning : alpha is nan; new init = %5.5f\n", init_a);
            a = init_a;
            log_a = log(a);
        end

        f = alhood(a, ss, D, K);
        df = d_alhood(a, ss, D, K);
        d2f = d2_alhood(a, D, K);
        log_a = log_a - df/(d2f * a + df);
        printf("alpha maximization : %5.5f   %5.5f\n", f, df);
    until ((abs(df) < NEWTON_THRESH) || (iter > MAX_ALPHA_ITER));
    opt = exp(log_a);
end


function value = alhood(a, ss, D, K)
	value = D * (log_gamma(K * a) - K * log_gamma(a)) + (a - 1) * ss; 
end

function value = d_alhood(a, ss, D, K)
	value = D * (K * digamma(K * a) - K * digamma(a)) + ss; 
end

function value = d2_alhood(a, D, K)
	value = D * (K * K * trigamma(K * a) - K * trigamma(a)); 
end

