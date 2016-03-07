function run_em(corpus, directory)
 	global EM_CONVERGENCE;
 	global EM_MAX_ITER;
 	global VAR_MAX_ITER;
 	global INITIAL_ALPHA;
 	global LAG;
 	global NTOPICS;

 	% allocate variational parameters
 	phi = zeros(corpus.max_doc_length, NTOPICS);
	var_gamma = zeros(corpus.num_docs, NTOPICS);

    % initialize model
	model = struct('num_topics', 0, 'num_terms', 0, 'alpha', 1, 'log_beta', []);
	model.num_topics = NTOPICS;
	model.num_terms = corpus.num_terms;
	model.log_beta = zeros(model.num_topics, model.num_terms);
	model.alpha = INITIAL_ALPHA;

	% sufficient stats
	ss = struct('class_word', [], 'class_total', [], 'num_docs', 0, 'alpha_suffstats', 0);
	ss.class_word = (1.0 / model.num_terms) + rand(model.num_topics, model.num_terms);
	ss.class_total = sum(ss.class_word, 2);
	ss.alpha_suffstats = 0;

	for i=1:model.num_topics
		for j=1:model.num_terms
		    if (ss.class_word(i, j) > 0)
				model.log_beta(i, j) = ss.class_word(i, j);
			end
		end
		model.log_beta(i, :)  = log(model.log_beta(i, :)) - log(sum(model.log_beta(i, :)));
	end

	model = lda_mle(model, ss, 0);

    filename = sprintf("%s/000",directory);
    save_lda_model(model, filename);


    % run expectation maximization

    likelihood = 0, likelihood_old = 0, converged = 1;

    filename = sprintf("%s/likelihood.dat", directory);
    likelihood_file = fopen(filename, "w");

    iter = 0;
    while (((converged < 0) || (converged > EM_CONVERGENCE) || (i <= 2)) && (iter <= EM_MAX_ITER))
        iter += 1; 
        printf("**** em iteration %d ****\n", iter);
        likelihood = 0;
        ss = zero_initialize_ss(ss, model);

        % e-step

		for d=1:corpus.num_docs 
            if (rem(d, 10) == 0) 
            	printf("document %d\n",d);
            end
			[l, ss] = doc_e_step(corpus.docs{d}, var_gamma(d, :), ...
                                     phi, model, ss);
			likelihood += l;
        end

        % m-step

        model = lda_mle(model, ss, 1);

        % check for convergence

        converged = (likelihood_old - likelihood) / (likelihood_old);
        if (converged < 0) 
        	VAR_MAX_ITER = VAR_MAX_ITER * 2;
        end
        likelihood_old = likelihood;

        % output model and likelihood

        fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
        fflush(likelihood_file);
        if (rem(iter, LAG) == 0)
            filename = sprintf("%s/%03d",directory, iter);
            save_lda_model(model, filename);
            filename = sprintf("%s/%03d.gamma",directory, iter);
            save_gamma(filename, var_gamma, corpus.num_docs, model.num_topics);
        end
    end

    % output the final model

    filename = sprintf("%s/final",directory);
    save_lda_model(model, filename);
    filename = sprintf("%s/final.gamma",directory);
    save_gamma(filename, var_gamma, corpus.num_docs, model.num_topics);

    % output the word assignments (for visualization)

    filename = sprintf("%s/word-assignments.dat", directory);
    w_asgn_file = fopen(filename, "w");
    for d=1:corpus.num_docs
        if (rem(d, 100) == 0) 
        	printf("final e step document %d\n",d);
        end
        likelihood += lda_inference(corpus.docs{d}, model, var_gamma(d), phi);
        % write_word_assignment(w_asgn_file, corpus.docs(d), phi, model);
    end
    fclose(w_asgn_file);
    fclose(likelihood_file);

end

function model = lda_mle(model, ss, estimate_alpha)

    for k=1:model.num_topics
    	for w=1:model.num_terms

    		if (ss.class_word(k,w) > 0)
    			model.log_beta(k, w) = log(ss.class_word(k, w)) - log(ss.class_total(k));
    		else
    			model.log_beta(k, w) = -100;
    		end
    	end
    end

    if (estimate_alpha)
    	model.alpha = opt_alpha(ss.alpha_suffstats, ss.num_docs, model.num_topics);
        printf("new alpha = %5.5f\n", model.alpha);
    end

end


function ss = zero_initialize_ss(ss, model)
    for k=1:model.num_topics
		ss.class_total(k) = 0;
        for w=1:model.num_terms
			ss.class_word(k, w) = 0;
        end
    end
    ss.num_docs = 0;
    ss.alpha_suffstats = 0;
end

function [likelihood, var_gamma, phi, ss] = doc_e_step(document, var_gamma, phi, model, ss)

    likelihood = 0;

    % posterior inference

    [likelihood, var_gamma, phi] = lda_inference(document, model, var_gamma, phi);

    % update sufficient statistics

    gamma_sum = 0;
    for k=1:model.num_topics
        gamma_sum += var_gamma(k);
        ss.alpha_suffstats += psi(var_gamma(k));
    end
    ss.alpha_suffstats -= model.num_topics * psi(gamma_sum);

    for n=1:document.length
        for k=1:model.num_topics
            ss.class_word(k, document.words(n)) += document.counts(n) * phi(n, k);
            ss.class_total(k) += document.counts(n)*phi(n, k);
        end
    end

    ss.num_docs = ss.num_docs + 1;
end

function save_gamma(filename, var_gamma, num_docs, num_topics)
    fileptr = fopen(filename, "w");

	for d=1:num_docs
    	% fprintf(fileptr, "%5.10f", var_gamma(d,0));
		for k=1:num_topics
		    fprintf(fileptr, " %5.10f", var_gamma(d,k));
		end
		fprintf(fileptr, "\n");
    end
    fclose(fileptr);
end

