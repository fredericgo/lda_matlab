% conversion from c to matlab
function [likelihood, var_gamma, phi] = lda_inference(document, model, var_gamma, phi)
    %
    % phi(N, K)

    global VAR_CONVERGENCE;
    global MAX_ITER;

    converged = 1;
    likelihood = 0;
    likelihood_old = 0;    

    % compute posterior dirichlet
    var_gamma = ones(1, model.num_topics) * model.alpha + document.total / model.num_topics;
    digamma_gam = psi(var_gamma);

    phi(:) = 1 / model.num_topics;

    var_iter = 0;
    while (converged > VAR_CONVERGENCE) && (var_iter < MAX_ITER),

        var_iter += 1;
        for n=1:document.length

            oldphi = phi(n, :);
            phi(n, :) = digamma_gam + model.log_beta(:, document.words(n))';

            % normalize phis
            % phisum = log_sum(phisum, phi(n, :));
            phisum = lognorm(phi(n,:));

            phi(n, :) = exp(phi(n, :) - phisum);

            % var_gamma[k] =
            %        var_gamma[k] + doc->counts[n]*(phi[n][k] - oldphi[k]);

            var_gamma += document.counts(n) * (phi(n, :) - oldphi);
            digamma_gam = psi(var_gamma);
        end

        likelihood = compute_likelihood(document, model, phi, var_gamma);
        %assert(!isnan(likelihood));
        if likelihood_old == 0,
            converged = 1;
        else 
            converged = (likelihood_old - likelihood) / likelihood_old;
        end
        % printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
        likelihood_old = likelihood;

    end
end

function likelihood = compute_likelihood(document, model, phi, var_gamma);

    likelihood = 0;
    digsum = 0;
    var_gamma_sum = 0;

    dig = psi(var_gamma);
    var_gamma_sum = sum(var_gamma);
    digsum = psi(var_gamma_sum);


    likelihood = log_gamma(model.alpha * model.num_topics) - model.num_topics * log_gamma(model.alpha) ...
                 - log_gamma(var_gamma_sum);

    for k=1:model.num_topics
        likelihood += (model.alpha - 1)*(dig(k) - digsum) + log_gamma(var_gamma(k)) ...
                     - (var_gamma(k) - 1)*(dig(k) - digsum);

        for n=1:document.length
        
            if (phi(n, k) > 0)
                likelihood += document.counts(n) * phi(n,k)*((dig(k)-digsum)-log(phi(n,k)) + model.log_beta(k, document.words(n)));
            end
        end
    end

end
