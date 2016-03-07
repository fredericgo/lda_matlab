function save_lda_model(model, model_root)
	filename = sprintf("%s.beta", model_root);
	file = fopen(filename, "w");
	disp(filename)

	for i=1:model.num_topics
		for j=1:model.num_terms
			fprintf(file, '%5.10f ', model.log_beta(i, j));
		end
	end

	fprintf(file, '\n');
	fclose(file);

	filename = sprintf("%s.ohter", model_root);
	file = fopen(filename, "w");
	fprintf(file, "num_topics %d\n", model.num_topics);
    fprintf(file, "num_terms %d\n", model.num_terms);
    fprintf(file, "alpha %5.10f\n", model.alpha);
    fclose(file);

end

