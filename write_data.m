function write_data()
	data = load('nips12raw_str602.mat');
	filename = 'nips5';
	out = fopen(filename, 'w');

	for d=1:5,
		[words, docs, counts] = find(data.counts(:,d));
		fprintf(out, '%d ', length(words));
		for i=1:(length(words)-1)
			fprintf(out, '%d:%d ', words(i), counts(i));
		end
		fprintf(out, '%d:%d', words(length(words)), counts(length(words)));

		fprintf(out, '\n');

	end

	fclose(out);
end