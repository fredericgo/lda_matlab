% global variables
global MAX_ITER = 20;
global VAR_CONVERGENCE = 1e-6;
global EM_MAX_ITER = 100;
global EM_CONVERGENCE = 1e-4;
global INITIAL_ALPHA = 1;
global LAG = 5;
global NTOPICS = 3;

% set documents
data = load('nips12raw_str602.mat');
directory = "nips";
mkdir(directory);

corpus = struct('docs', [], 'max_doc_length', 0, 'num_docs', 0, 'num_terms', 0);
corpus.docs = {};
corpus.num_docs = 10; % data.Np;
corpus.num_terms = data.Nw;

for d=1:corpus.num_docs,
	doc = struct('words', []);
	[doc.words, _, doc.counts] = find(data.counts(:,d));
	doc.total = data.plengths(d);
	doc.length = length(doc.words);
	if (doc.length > corpus.max_doc_length) 
		corpus.max_doc_length = doc.length;
	end
	corpus.docs{d} = doc;
end


run_em(corpus, directory);