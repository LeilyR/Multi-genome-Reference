#include "dynamic_mc.hpp"

template<typename reader> 
alignment_contexts<reader>::alignment_contexts(reader & r, size_t level):rd(r), level(level)  {
}

template<typename reader> 
alignment_contexts<reader>::~alignment_contexts() {}



/*
	This function computes the next modification (strand 1 to 2) encoded in a 
	pairwise alignment starting from column 'at'.

	Possible modifications are
	Keep (several characters, number is a power of two for binary encoding), on_char is the first character kept
	Modify (single character, on char is the changed character)
	Insert (we insert a single character before the current position, on_char is a character on strand 1 after the insertion)
	Delete (we delete several characters, number is a power of two for binary encoding, on_char is the first character that is deleted)

	in all cases num_cols returns the number of alignment columns used by this modification, the next call to this function should be at 'at+num_cols'

*/
template<typename reader>
void alignment_contexts<reader>::next_modification_1_2(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols) {
	char ch1, ch2;
	al.alignment_col(at, ch1, ch2);


	if(ch1=='-') { // Insertion
// TODO look at insertions at alignment end again: on_char could be after the current alignment, but in decoding we do not know that the alignment has ended. Maybe this problem is not so bad because we usually cut end gaps
		on_char = (size_t) dnastring::base_to_index('A');
		for(size_t next_r1base = at+1; next_r1base < al.alignment_length(); ++next_r1base) {
			char n_ch1, n_ch2;
			al.alignment_col(next_r1base, n_ch1, n_ch2);
			if(n_ch1!='-') {
				on_char = (size_t) dnastring::base_to_index(n_ch1);
				break;
			}
		}
		modification = dynamic_mc_model::modification_character(dnastring::base_to_index(ch2), -1, -1, -1);
		num_cols = 1;
		return;
	}

	on_char = (size_t) dnastring::base_to_index(ch1); // in all other cases there is a current char base for the modification

	if(ch2=='-') { // Deletion
		size_t max = (1<<(NUM_DELETE_DYN-1));
		size_t extra_del = 0; // number of deleted chars after the first one
		for(; at + extra_del + 1 < al.alignment_length(); ++extra_del) {
			char n_ch1, n_ch2;
			al.alignment_col(at + extra_del+1, n_ch1, n_ch2);
			if(n_ch2!='-') {
				break; // break before increase of extra_del
			}
			// no break, now we know we can add one to extra_del
			if(extra_del + 2 == max) {
				extra_del++;
				break;
			}
		}
		size_t num_deleted = extra_del + 1;
		size_t log_del = 1;
		for(; log_del < NUM_DELETE_DYN; ++log_del) {
			if(num_deleted < (1<<log_del)) {
				break;
			}
		}
		num_cols = (1<<(log_del-1));
		modification = dynamic_mc_model::modification_character(-1, log_del -1, -1, -1);

		return;
	}

	if(ch1==ch2) { // Keep
		size_t max = (1<<(NUM_KEEP_DYN-1));
		size_t extra_keep = 0;
		for(; at + extra_keep + 1 < al.alignment_length(); ++extra_keep) {
			char n_ch1, n_ch2;
			al.alignment_col(at + extra_keep+1, n_ch1, n_ch2);
			if(n_ch1!=n_ch2) {
				break;
			}
			if(extra_keep + 2 == max) {
				extra_keep++;
				break;
			}
		}
		size_t num_kept = extra_keep + 1;
		size_t log_kep = 1;
		for(; log_kep < NUM_KEEP_DYN; ++log_kep) {
			if(num_kept < (1<<log_kep)) {
				break;
			}
		}
		num_cols = (1<<(log_kep-1));
		modification = dynamic_mc_model::modification_character(-1, -1, -1, log_kep -1);
		
	} else { // Modify
		num_cols = 1;
		modification = dynamic_mc_model::modification_character(dnastring::base_to_index(ch2), -1 , -1, -1);
	}

}

template<typename reader>
void alignment_contexts<reader>::next_modification_2_1(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols) {
	char ch1, ch2;
	al.alignment_col(at, ch1, ch2);


	if(ch2=='-') { // Insertion
// TODO look at insertions at alignment end again: on_char could be after the current alignment, but in decoding we do not know that the alignment has ended. Maybe this problem is not so bad because we usually cut end gaps
		on_char = (char) dnastring::base_to_index('A');
		for(size_t next_r2base = at+1; next_r2base < al.alignment_length(); ++next_r2base) {
			char n_ch1, n_ch2;
			al.alignment_col(next_r2base, n_ch1, n_ch2);
			if(n_ch2!='-') {
				on_char = (char) dnastring::base_to_index(n_ch2);
				break;
			}
		}
		modification = dynamic_mc_model::modification_character(dnastring::base_to_index(ch1), -1, -1, -1);
		num_cols = 1;
		return;
	}

	on_char = (char) dnastring::base_to_index(ch2); // in all other cases there is a current char base for the modification

	if(ch1=='-') { // Deletion
		size_t max = (1<<(NUM_DELETE_DYN-1));
		size_t extra_del = 0; // number of deleted chars after the first one
		for(; at + extra_del + 1 < al.alignment_length(); ++extra_del) {
			char n_ch1, n_ch2;
			al.alignment_col(at + extra_del+1, n_ch1, n_ch2);
			if(n_ch1!='-') {
				break; // break before increase of extra_del
			}
			// no break, now we know we can add one to extra_del
			if(extra_del + 2 == max) {
				extra_del++;
				break;
			}
		}
		size_t num_deleted = extra_del + 1;
		size_t log_del = 1;
		for(; log_del < NUM_DELETE_DYN; ++log_del) {
			if(num_deleted < (1<<log_del)) {
				break;
			}
		}
		num_cols = (1<<(log_del-1));
		modification = dynamic_mc_model::modification_character(-1, log_del -1, -1, -1);

		return;
	}

	if(ch1==ch2) { // Keep
		size_t extra_keep = 0;
		size_t max = (1<<(NUM_KEEP_DYN-1));
		for(; at + extra_keep + 1 < al.alignment_length(); ++extra_keep) {
			char n_ch1, n_ch2;
			al.alignment_col(at + extra_keep+1, n_ch1, n_ch2);
			if(n_ch1!=n_ch2) {
				break;
			}
			if(extra_keep + 2 == max) {
				extra_keep++;
				break;
			}
		}
		size_t num_kept = extra_keep + 1;
		size_t log_kep = 1;
		for(; log_kep < NUM_KEEP_DYN; ++log_kep) {
			if(num_kept < (1<<log_kep)) {
				break;
			}
		}
		num_cols = (1<<(log_kep-1));
		modification = dynamic_mc_model::modification_character(-1, -1, -1, log_kep -1);
		
	} else { // Modify
		num_cols = 1;
		modification = dynamic_mc_model::modification_character(dnastring::base_to_index(ch1), -1 , -1, -1);
	}

}




template<typename reader>
void alignment_contexts<reader>::read_alignment_1_2(const pw_alignment & pw) {

	std::string context;
	for(size_t i=level; i>0; --i) {
		size_t numk = i-1;
		if(numk > NUM_KEEP_DYN-1) numk = NUM_KEEP_DYN;
		char keepc = dynamic_mc_model::modification_character(-1, -1, -1, numk);
		context.append(1, keepc);
	}
	size_t at_al_col = 0;
	while(at_al_col < pw.alignment_length()) {

		char modc;
		char onc;
		size_t len;
		next_modification_1_2(pw, at_al_col, modc, onc, len);
		
		// add current known base to context
		std::string newcontext = context;
		newcontext.append(1, onc);

		rd.see_context(newcontext, (size_t)modc);
		
		// update mod context for next round
		context = context.substr(1);
		context.append(1, modc);
			
		at_al_col += len;
	}

	assert(at_al_col == pw.alignment_length());
}

template<typename reader>
void alignment_contexts<reader>::read_alignment_2_1(const pw_alignment & pw) {

	std::string context;
	for(size_t i=level; i>0; --i) {
		size_t numk = i-1;
		if(numk > NUM_KEEP_DYN-1) numk = NUM_KEEP_DYN;
		char keepc = dynamic_mc_model::modification_character(-1, -1, -1, numk);
		context.append(1, keepc);
	}
	size_t at_al_col = 0;
	while(at_al_col < pw.alignment_length()) {

		char modc;
		char onc;
		size_t len;
		next_modification_2_1(pw, at_al_col, modc, onc, len);
		
		// add current known base to context
		std::string newcontext = context;
		newcontext.append(1, onc);

		rd.see_context(newcontext, (size_t)modc);
		
		// update mod context for next round
		context = context.substr(1);
		context.append(1, modc);
			
		at_al_col += len;
	}

	assert(at_al_col == pw.alignment_length());
}

void reader_counter::see_context(const std::string & context, size_t modification) {
	
	std::map<std::string, std::vector<size_t> >::iterator findc = countsmap.find(context);
	if(findc == countsmap.end()) {
		countsmap.insert(std::make_pair(context, vector<size_t>(NUM_MODIFICATIONS, 0)));
		findc = countsmap.find(context);
	}
	findc->second.at(modification)++;
}

dynamic_mc_model::dynamic_mc_model(all_data & d, size_t nt): data(d), num_threads(nt) {
	// use (1<<i) should be faster	
/*
	size_t numberOfPowers = 32;
	powersOfTwo = std::vector<size_t>(numberOfPowers, 1);
	for(size_t i=1; i< numberOfPowers; ++i) {
		powersOfTwo.at(i) = powersOfTwo.at(i-1)*2;
	}
*/
}
dynamic_mc_model::~dynamic_mc_model(){}

char dynamic_mc_model::modification_character(int modify_base, int num_delete, int insert_base, int num_keep) {
		//return the enc
	if(num_delete != -1){
	//	std::cout<< "there is a delete of length "<< num_delete <<std::endl; 
	//	return NUM_DELETE+num_delete;
		return 5 + num_delete;
	}
	if(num_keep != -1){
	//	std::cout<< "there is a keep of length" << num_keep << std::endl;
	//	return NUM_KEEP+num_keep;
		return 5 + NUM_DELETE + num_keep;
	}
	if(modify_base != -1) {
	//	std::cout<< "there is a modification to " << dnastring::index_to_base(modify_base) << std::endl;
		return modify_base;
	}
	if(insert_base != -1){
	//	std::cout<< " there is a insertion of " << dnastring::index_to_base(insert_base) <<std::endl;
		return insert_base + NUM_KEEP + NUM_DELETE + 5;
	}
	assert (false);
	return -1;
}

std::string dynamic_mc_model::tranlsate_modification_character(char enc) {
	int modify_base = -1;
	int num_delete =-1;
	int insert_base = -1;
	int num_keep = -1;
	modification(enc, modify_base, num_delete, insert_base, num_keep);
	std::stringstream s;
	if(num_delete != -1){
		s<< " a delete of length " <<  num_delete; 
	}
	if(num_keep != -1){
		s<< " a keep of length" << num_keep;
	}
	if(modify_base != -1) {
		s<< " a modification to " << dnastring::index_to_base(modify_base);
	}
	if(insert_base != -1){
		s<< " an insertion of " << dnastring::index_to_base(insert_base);
	}
	return s.str();
}

size_t dynamic_mc_model::modification_length(char mod) {
	int modify_base = -1;
	int num_delete =-1;
	int insert_base = -1;
	int num_keep = -1;
	size_t length = 0;
	modification(mod, modify_base, num_delete, insert_base, num_keep);
	if(modify_base != -1){
		length = 1;
	}
	if(insert_base != -1){
		length = 0;
	}
	if(num_delete != -1){
		length = num_delete;
	}
	if(num_keep != -1){
		length = num_keep;
	}
	return length;
}

void dynamic_mc_model::modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep) {
	modify_base = -1;
	num_delete =-1;
	insert_base = -1;
	num_keep = -1;

	if(enc < 5) {
		modify_base = enc;
		return;	
	} 
	if(enc < 5 + NUM_DELETE) {
		num_delete = (1<<(enc - 5));
		return;
	}
	if (enc < 5 + NUM_DELETE + NUM_KEEP){
		num_keep = (1<<(enc-NUM_DELETE-5));
		return;
	}
	if(enc< 5+ 5 + NUM_KEEP + NUM_DELETE){
		insert_base = enc- 5 - NUM_KEEP - NUM_DELETE;
		return;
	}
	assert(0);
}


/*
	computes total cost in that map of counts
*/
double dynamic_mc_model::total_sequence_cost(const std::map<std::string, std::vector<size_t> > & all_context_counts) const {

	double information_cost = 0;
	
	for(std::map<std::string, std::vector<size_t> >::const_iterator it = all_context_counts.begin(); it!=all_context_counts.end(); ++it) {
		size_t all = 1;
		const std::vector<size_t> & counts = it->second;
		for(size_t i=0; i<5; ++i) {
			all+=counts.at(i);
		}
		for(size_t i=0; i<5; ++i) {
			information_cost -= counts.at(i) * log2((double)counts.at(i) / (double) all);
		}
	}
	return information_cost;
}

/*
	add a count of base num times after context into the map

*/

void dynamic_mc_model::context_count(std::map<std::string, std::vector<size_t> > & countmap, const std::string & context, size_t base, size_t num) const {
	std::map<std::string, std::vector<size_t> >::iterator findc = countmap.find(context);
	if(findc == countmap.end()) {
		countmap.insert(std::make_pair(context, std::vector<size_t>(5,1)));
		findc = countmap.find(context);
	}
	findc->second.at(base)+=num;
}

void dynamic_mc_model::compute_sequence_model() {
	sequence_model_widths.resize(data.numAcc());
	vector<unsigned char> alphabet;
	get_sequence_alphabet(alphabet);


	for(size_t acc = 0; acc<data.numAcc(); ++acc) {
		model_bits_to_full_high_values(sequence_models.at(acc), acc_sequence_all_bases_total.at(acc), MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, alphabet, sequence_model_widths.at(acc));
	}

}
void dynamic_mc_model::compute_alignment_model() {
	vector<unsigned char> alignment_alphabet;
	get_alignment_alphabet(alignment_alphabet);
	alignment_model_widths = std::vector<std::vector<std::map<std::string, std::vector<uint32_t> > > > (data.numAcc(), std::vector<std::map<std::string, std::vector<uint32_t> > >(data.numAcc()));
	for(size_t a1 = 0; a1<data.numAcc(); ++a1) {
		for(size_t a2 = 0; a2<data.numAcc(); ++a2) {
			model_bits_to_full_high_values(alignment_models.at(a1).at(a2), alignments_all_modification_total.at(a1).at(a2), MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL + 1, alignment_alphabet, alignment_model_widths.at(a1).at(a2));
		}
	}

}


void dynamic_mc_model::train_sequence_model() {

	vector<size_t> alignments_on_acc(data.numAcc(), 0);
#pragma omp parallel for schedule(static) num_threads(num_threads)
	for(size_t a=0; a<data.numAlignments(); ++a) {
		const pw_alignment & al = data.getAlignments().at(a);
		size_t r1 = al.getreference1();
		size_t r2 = al.getreference2();
#pragma omp critical(count)
{
		alignments_on_acc.at(data.accNumber(r1))++;
		alignments_on_acc.at(data.accNumber(r2))++;
}
	}



	std::vector<unsigned char> alphabet;
	get_sequence_alphabet(alphabet);

	std::string astring("");
	for(size_t i=0; i<MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; ++i) {
		astring.append("A");
	}

	sequence_models.resize(data.numAcc());
	acc_sequence_all_bases_total.resize(data.numAcc());
	acc_sequence_alignment_flag.resize(data.numAcc());
#pragma omp paralell for schedule(dynamic) num_threads(num_threads)
	for(size_t acc=0; acc<data.numAcc(); ++acc) {
		std::map< std::string, std::vector<size_t> > context_counts;
		size_t acc_num_bases = 0;
		for(size_t k = 0; k <data.getAcc(acc).size(); k++) {
			std::string sequence = data.getSequence(data.getAcc(acc).at(k)).str();
			acc_num_bases += sequence.length();



			for(size_t i=0; i<MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; ++i) {
				std::string cont = astring.substr(i);
				cont.append(sequence.substr(0, i));
				context_count(context_counts, cont, dnastring::base_to_index(sequence.at(i)), 1);
				assert(cont.length() == MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
			}
			for(size_t i=MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; i<sequence.length(); ++i) {
				std::string cont = sequence.substr(i - MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
				context_count(context_counts, cont, dnastring::base_to_index(sequence.at(i)), 1);
			}
		} // for sequences in accession

		// we estimate the probability of having an alignment on this accession by the total number of half alignments on that accession
		// in the end we could have less (we select only some alignments) or more (we split alignments to remove overlap)
		size_t al_on_this_acc = alignments_on_acc.at(acc);
		size_t acc_sequence_ends = data.getAcc(acc).size()+1;
#pragma omp critical(count_results)
{
		std::cout << "on acc " << acc << " there are " << al_on_this_acc << " half alignments " << acc_num_bases << " bases and " << acc_sequence_ends << " sequence ends " << std::endl;
		if(al_on_this_acc < 1 ) al_on_this_acc = 1;
		assert(acc_num_bases>0);
		size_t total = acc_sequence_ends + al_on_this_acc + acc_num_bases;
		double acc_al_prob = (double) al_on_this_acc / (double) total;
		double acc_sequence_end_prob = (double) acc_sequence_ends / (double) total;
		double acc_base_prob = (double) acc_num_bases/ (double) total;
	//	std::cout << "approximate encoding cost on acc: any base: " << -log2( acc_base_prob) << "bit half alignment begin: " << -log2(acc_al_prob) << "bit sequence end: " << -log2(acc_sequence_end_prob) << "bit"<< std::endl; 
		// compute flags weight on this accession, this means the total value available for encoding bases decreases
		vector<double> distri(3);
		distri.at(0) = acc_base_prob;
		distri.at(1) = acc_al_prob;
		distri.at(2) = acc_sequence_end_prob;
		// we use total for max width because event 0 is the sum of all sequence events
		vector<uint32_t> flag_bits;
		distribution_to_bits(distri, MAX_ENCODING_WIDTH_BITS, flag_bits);
		vector<uint32_t> flag_widths;
		bits_to_width_encoding(flag_bits, TOTAL_FOR_ENCODING, flag_widths);
		std::cout <<"acc " << acc << ": encoding widhts: any base: "<< flag_widths.at(0) << " half alignment begin: " << flag_widths.at(1) << " sequence end: " << flag_widths.at(2) << std::endl;
		std::cout <<"acc " << acc << ": encoding costs: any base: " << -log2((double)flag_widths.at(0)/(double)TOTAL_FOR_ENCODING) << 
			 "bit half alignment begin: " << -log2((double)flag_widths.at(1)/(double)TOTAL_FOR_ENCODING) << "bit sequence end: " << -log2((double)flag_widths.at(2)/(double)TOTAL_FOR_ENCODING) <<"bit "<< std::endl;
		
		uint32_t new_total = flag_widths.at(0);
		std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > acc_seq_model;

		std::cout << " run context selector on " << context_counts.size() << " contexts"<< std::endl;
		double enc_cost = context_selector(context_counts, alphabet, new_total, acc_seq_model, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL );
		std::cout << " select " << acc_seq_model.size() << " contexts, total sequence encoding needs " << enc_cost << "bit (" << enc_cost/(double)acc_num_bases << " bit/base)" <<  std::endl;

		sequence_models.at(acc) = acc_seq_model;
		
		acc_sequence_all_bases_total.at(acc) = flag_widths.at(0);
		acc_sequence_alignment_flag.at(acc) = flag_widths.at(0) + flag_widths.at(1);


		for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::iterator it = acc_seq_model.begin(); it!=acc_seq_model.end(); ++it) {
			bool use_sqroot;
			size_t num_bits;
			get_model_type(it->second.first, use_sqroot, num_bits);

//			std::cout << "res " << it->first << " sqroot " << use_sqroot << " bits " << num_bits << std::endl;
		}
}

	}

	
	
		
	compute_sequence_model();

}





unsigned char dynamic_mc_model::get_model_index(bool use_sqroot, size_t num_bits) {
	unsigned char res = 0;
	if (use_sqroot) {
		res+= MAX_ENCODING_WIDTH_BITS - MIN_ENCODING_WIDTH_BITS + 1;
	}
	res+=num_bits - MIN_ENCODING_WIDTH_BITS;
	return res;
}

void dynamic_mc_model::get_model_type(unsigned char model_index, bool & use_sqroot, size_t & num_bits) {
	if(model_index > MAX_ENCODING_WIDTH_BITS - MIN_ENCODING_WIDTH_BITS) {
		num_bits = model_index - (MAX_ENCODING_WIDTH_BITS - MIN_ENCODING_WIDTH_BITS + 1) + MIN_ENCODING_WIDTH_BITS;
		use_sqroot = 1;
	} else {
		num_bits = model_index + MIN_ENCODING_WIDTH_BITS;
		use_sqroot = 0;
	}
}



void dynamic_mc_model::counts_to_distribution(const std::vector<size_t> & counts, std::vector<double> & distri) {
	size_t sum = 0;
	distri.resize(counts.size());
	for(size_t i=0; i<counts.size(); ++i) {
		sum+=counts.at(i);
	}
	if(sum==0) sum++;
	for(size_t i=0; i<counts.size(); ++i) {
		distri.at(i) = (double) counts.at(i) / (double) sum;
	}
}


double dynamic_mc_model::counts_to_bits_eval(const std::vector<size_t> & counts, uint32_t target_total, unsigned char & model_index, std::vector<uint32_t> & bits) {
	
	double best_cost = std::numeric_limits<double>::max();
	std::vector<uint32_t> best_bits;
	unsigned char best_index;

	vector<double> distri;
	vector<uint32_t> widths;
	counts_to_distribution(counts, distri);
	
	for(model_index = 0; model_index < NUM_ENCODING_TYPES; ++model_index) {
		size_t num_bits;
		bool use_sqroot;
		get_model_type(model_index, use_sqroot, num_bits);
		if(use_sqroot) {
			distribution_to_sqroot_bits(distri, num_bits, bits);
			sqroot_bits_to_width_encoding(bits, target_total, widths);
		} else {
			distribution_to_bits(distri, num_bits, bits);
			bits_to_width_encoding(bits, target_total, widths);
		}

//		std::cout << " mi " << (size_t)model_index;
//		for(size_t i=0; i<bits.size(); ++i) {
//			std::cout << " " << distri.at(i) << " " << bits.at(i) <<std::endl;
//		}


		// encoding cost
		double cost =  count_cost_eval(counts, widths);
		// model parameter cost

//		std::cout << "c " << cost;
			cost += widths.size() * num_bits;
//		std::cout << " cp " << cost;
		if(cost < best_cost) {
//			std::cout << " * ";
			best_cost = cost;
			best_bits = bits;
			best_index = model_index;
		}
//		std::cout << std::endl;
	}
	model_index = best_index;
	bits = best_bits;
	return best_cost;
}


double dynamic_mc_model::count_cost_eval(const std::vector<size_t> & counts, const std::vector<uint32_t> & enc_widths) {
	assert(counts.size() == enc_widths.size());
	size_t total = 0;
	for(size_t i=0; i<counts.size(); ++i) {
		total+=enc_widths.at(i);
	}
	if(total==0) total++;
	double total_cost = 0;
	for(size_t i=0; i<counts.size(); ++i) {
//		std::cout << " " << enc_widths.at(i);
		if(counts.at(i) > 0) {
			double icost = -log2((double) enc_widths.at(i) / (double) total  );
			total_cost += counts.at(i) * icost;
		}
	}
//	std::cout << std::endl;
	return total_cost;
}


/*
	preorder recursion function to sum up all context counts

*/
void dynamic_mc_model::context_subtree_cost(const std::string & cur_context, const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, size_t max_length, 
	std::map< std::string, std::vector<size_t> > & sum_counts) {

/*
	std::cout << " callxc " << cur_context.length();
	for(size_t i=0; i<cur_context.length(); ++i) {
		std::cout << ":" << (size_t) cur_context.at(i);
	}
	std::cout << std::endl;

*/

	if(cur_context.length() < max_length) {

		std::vector<size_t> cur_sums(alphabet.size(), 0);
		// add subtrees
		for(size_t i=0; i<alphabet.size(); ++i) {
			std::string lcont;
			lcont.append(1, alphabet.at(i));
			lcont.append(cur_context);
			context_subtree_cost(lcont, context_counts, alphabet, max_length, sum_counts);
			std::map<std::string, std::vector<size_t> >::iterator findc = sum_counts.find(lcont);
			assert(findc != sum_counts.end());
			std::vector<size_t> & lcounts = findc->second;
			for(std::size_t j=0; j<lcounts.size(); ++j) {
				cur_sums.at(j)+=lcounts.at(j);
			}
 
		}
		// add current counts
		std::map<std::string, std::vector<size_t> >::const_iterator fcur = context_counts.find(cur_context);
		if(fcur != context_counts.end()) {
			const std::vector<size_t> & cc = fcur->second;
			for(size_t i=0; i<cc.size(); ++i) {
				cur_sums.at(i)+=cc.at(i);
			}
	//		std::cout << " * ";
		}

		sum_counts.insert( std::make_pair( cur_context, cur_sums) );

	
		
		// recursion end, sum of length 1
	} else if(cur_context.length() == max_length) {
		std::map<std::string, std::vector<size_t> >::const_iterator findc = context_counts.find(cur_context);
		if(findc == context_counts.end()) {
	//		std::cout << " -" << context_counts.size() << "- "; 
			sum_counts.insert( std::make_pair( cur_context, std::vector<size_t>(alphabet.size(), 0)) );
		} else {
	//		std::cout << " . "; 
			sum_counts.insert( std::make_pair( cur_context, findc->second) );
			
			std::vector<size_t> c = findc->second;
/*
		std::cout << " csc " << cur_context;
		for(size_t i=0; i<c.size(); ++i) {
			std::cout << ": " << c.at(i);
		}
		std::cout << std::endl;
*/



		}
	} else {
		assert(0);
	}
	


}



double dynamic_mc_model::context_selector(const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, uint32_t target_total, 
	std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & result_bits, size_t max_length ) {
	double sum_cost = 0;	
	// context -> ((cost, counts),(type_index, bits))
	typedef std::pair< std::pair<double, std::vector<size_t> >, std::pair<unsigned char, std::vector<uint32_t> > > pairstype;
	std::map<std::string, pairstype > all_evals;

	std::map< std::string, std::vector<size_t> > subtree_context_counts;
	std::string empty;
	
	context_subtree_cost(empty, context_counts, alphabet, max_length, subtree_context_counts);

	std::set<std::string> cur_set;
	cur_set.insert(empty);
// empty eval
	std::map<std::string, std::vector<size_t> >::iterator finde = subtree_context_counts.find(empty);
	assert(finde!=subtree_context_counts.end());
	std::vector<size_t> ecounts = finde->second;
	unsigned char empty_model_index;
	std::vector<uint32_t> ebits;
	double ecost = counts_to_bits_eval(ecounts, target_total, empty_model_index, ebits);
	all_evals.insert(std::make_pair(empty, std::make_pair(std::make_pair(ecost, ecounts), std::make_pair(empty_model_index, ebits))));

/*
	std::cout << " ecounts ";
	for(size_t i=0; i<ecounts.size(); ++i) {
		std::cout << " " << ecounts.at(i);
	}
	std::cout << " cost " << ecost << std::endl;
*/

	
	for(size_t cur_length=0; cur_length<=max_length; ++cur_length) {

		std::set<std::string> next_set;
		for(std::set<std::string>::iterator it = cur_set.begin(); it!=cur_set.end(); ++it) {
			std::string cont = *it;
			std::map<std::string, pairstype>::iterator findc = all_evals.find(cont);
			pairstype cont_eval = findc->second;
			double cont_cost = cont_eval.first.first;
			double cur_children_cost = 0;
			std::set<std::string > cur_children;
			if(cur_length < max_length) {
			// all child nodes:
			for(size_t i=0; i<alphabet.size(); ++i) {
				std::string lcont;
				lcont.append(1, alphabet.at(i));
				lcont.append(cont);
				std::map< std::string, std::vector<size_t> >::iterator findcounts = subtree_context_counts.find(lcont);
				if(findcounts == subtree_context_counts.end()) assert(0);
				vector<size_t> lcounts = findcounts->second;
				unsigned char model_index;
				std::vector<uint32_t> lbits;
				double lcost = counts_to_bits_eval(lcounts, target_total, model_index, lbits);
				cur_children_cost+=lcost;
				cur_children.insert(lcont);
				all_evals.insert(std::make_pair(lcont, std::make_pair(std::make_pair(lcost, lcounts), std::make_pair(model_index, lbits))));
			}
			}


			if(cont_cost <= cur_children_cost || cur_length == max_length) {
//				std::cout << "SELECT " << cont << " costs " << cont_cost << " " << cur_children_cost<<  std::endl;
//				std::vector<size_t> ccount = cont_eval.first.second;
//				for(size_t i=0; i<ccount.size(); ++i) {
//					std::cout << " " << ccount.at(i);
//				}
//				std::cout << std::endl;

				sum_cost += cont_cost;

				result_bits.insert(std::make_pair(cont, std::make_pair(cont_eval.second.first, cont_eval.second.second)));

			} else {
				for(std::set<std::string>::iterator cit = cur_children.begin(); cit!=cur_children.end(); ++cit) {
					next_set.insert(*cit);
				}

			} 

		}
	
		cur_set = next_set;
	}

	

	return sum_cost;


}

/*
double dynamic_mc_model::context_selector(const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<char> & alphabet, uint32_t target_total, 
	std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & cpref_bits ) {
	double sum_cost = 0;	
	// context -> ((cost, counts),(type_index, bits))
	typedef std::pair< std::pair<double, std::vector<size_t> >, std::pair<unsigned char, std::vector<uint32_t> > > pairstype;
	typedef std::map< std::string, pairstype >  ce_type;
	std::vector<ce_type> context_evals; // context length -> context -> ((cost, counts),(type_index, bits))
	for(std::map< std::string, std::vector<size_t> >::const_iterator it = context_counts.begin(); it!=context_counts.end(); ++it) {
		std::string cont = it->first;
		std::vector<size_t> counts = it->second;
		
		unsigned char model_index;
		std::vector<uint32_t> bits;
		double cost = counts_to_bits_eval(counts, target_total, model_index, bits);
		if(context_evals.size() < cont.length()+1) context_evals.resize(cont.length()+1);
		context_evals.at(cont.length()).insert(std::make_pair(cont, std::make_pair(std::make_pair(cost, counts), std::make_pair(model_index, bits)))); 
			
		bool use_sqroot;
		size_t num_bits;
		get_model_type(model_index, use_sqroot, num_bits);
	//	std::cout << " init eval of " << cont << " bits " << num_bits << " sqroot " << use_sqroot << " cost " << cost << std::endl;



	}


	std::cout << "Cont SEL lens " << std::endl;
	for(size_t i=0; i<context_evals.size(); ++i) {
		std::cout << " l " << i<< " size " << context_evals.at(i).size() << std::endl;
	}

	
	bool go = true;

	while(go) {
		// get a set of sibling leaves (char 0 variable, 1..n fixed) with the longest possible context still available
		std::string cur_patt("");
		go = false;
		for(size_t c = context_evals.size(); c>0; --c) {
			size_t cl = c-1;
			if(!context_evals.at(cl).empty()) {
				cur_patt = context_evals.at(cl).begin()->first;
				go = true;
				break;
			}


		}
		if(0!=cur_patt.compare("")) {	 
			std::vector<std::string> sibling_leaves;
			std::vector< pairstype > sibling_evals;
			double leaves_cost = 0;
			for(size_t i=0; i<alphabet.size(); ++i) {
				std::string searchv = cur_patt;
				searchv.at(0) = alphabet.at(i);
				std::map< std::string, pairstype >::iterator findcur = context_evals.at(searchv.length()).find(searchv);
				if(findcur != context_evals.at(searchv.length()).end()) {
					sibling_leaves.push_back(searchv);
					sibling_evals.push_back(findcur->second);
					leaves_cost += sibling_evals.at(sibling_evals.size()-1).first.first;
				}	
			}
			// TODO alle suffixe rein, wo kommen kürzere contexte her?

			// build a common model (reducing context length by one)
			vector<size_t> sums(sibling_evals.at(0).first.second.size(), 0);
			for(size_t i=0; i< sums.size() ; ++i) {
				for(size_t j=0; j<sibling_evals.size(); ++j) {
					sums.at(i)+=sibling_evals.at(j).first.second.at(i);
				}
			}
			unsigned char common_index;
			std::vector<uint32_t> common_bits;
			double common_cost = counts_to_bits_eval(sums, target_total, common_index, common_bits);

			// remove all the leaves
			for(size_t i=0; i<sibling_leaves.size(); ++i) {
				size_t num_erased = context_evals.at(sibling_leaves.at(i).length()).erase(sibling_leaves.at(i));
				assert(num_erased==1);
			}	

			if(common_cost < leaves_cost) {
				// remove first char of context string
				std::string common_context("");
				for(size_t i=1; i<cur_patt.length(); ++i) {
					common_context.append(1, cur_patt.at(i));
				}
				bool use_sqroot;
				size_t num_bits;
				get_model_type(common_index, use_sqroot, num_bits);
			//	std::cout << "join contexts to " << common_context << " cost " << common_cost << " less than leaves costs " << leaves_cost <<  " bits " << num_bits << " sqroot: " << use_sqroot << std::endl;
					// add joint leaves as a shorter context node
				if(common_context.length()>0) {
					context_evals.at(common_context.length()).insert(std::make_pair(common_context, std::make_pair( std::make_pair(common_cost , sums), std::make_pair( common_index, common_bits) )));
				} else {

			
					cpref_bits.insert(std::make_pair(common_context, std::make_pair( common_index, common_bits ) ) );
					sum_cost += common_cost;
				}
			} else {
			//	std::cout << " keep all sibling leaves of " << cur_patt << " cost " << leaves_cost <<  " number " << sibling_leaves.size();
				// add the leaves into result set
				for(size_t i=0; i<sibling_leaves.size(); ++i) {
	
			//		std::cout << " " << sibling_leaves.at(i) << ": " << sibling_evals.at(i).first.first;
					cpref_bits.insert(std::make_pair(sibling_leaves.at(i), std::make_pair(sibling_evals.at(i).second.first, sibling_evals.at(i).second.second ) ) );
				}
			//	std::cout << std::endl;
				sum_cost+=leaves_cost;
			}
		} else {
			// we should only get here if the function gets called with only empty context
			std::map<std::string, pairstype>::iterator finde = context_evals.at(0).find(cur_patt);
			assert(finde!=context_evals.at(0).end() || !go);
			if(go) {
				pairstype evals = finde->second;
				unsigned char common_index = evals.second.first;
				std::vector<uint32_t> common_bits = evals.second.second;
				cpref_bits.insert(std::make_pair(cur_patt, std::make_pair(common_index, common_bits ) ) );
				sum_cost += evals.first.first;
			}
	

		}

	
	} // while true
	return sum_cost;
}
*/


/*
	distribution to storable parameters
	bits bit per parameter (encoding width), zero width is ok for very small width

	bit encoding is used for efficient storage of parameters

	we store square roots of widths with bits precision	


*/


void dynamic_mc_model::distribution_to_sqroot_bits(const std::vector< double > & distri, size_t bits, std::vector<uint32_t> & sqbits) {
	size_t max_value = (1<<bits) - 1; // no encoded sqroot width larger than this

	double maxprob = 0;
	for(size_t i=0; i<distri.size(); ++i) {
		if(maxprob <distri.at(i)) maxprob = distri.at(i);
	}

	double scale = (double) max_value / sqrt(maxprob);
	sqbits.resize(distri.size());
	for(size_t i=0; i<distri.size(); ++i) {
		double nv = sqrt(distri.at(i)) * scale;
		uint32_t nvi = (uint32_t) (nv+ 0.5);
		if(nvi > max_value) nvi = max_value;
		sqbits.at(i) = nvi;
	}
} 


void dynamic_mc_model::distribution_to_bits(const std::vector<double> & distri, size_t bits, std::vector<uint32_t> & vbits) {
	size_t max_value = (1<<bits) -1; // no encoded width larger than this
	double maxprob = 0;
	for(size_t i=0; i<distri.size(); ++i) {
		if(maxprob <distri.at(i)) maxprob = distri.at(i);
	}

	double scale = (double) max_value / maxprob;
	vbits.resize(distri.size());
	for(size_t i=0; i<distri.size(); ++i) {
		double nv = distri.at(i) * scale;
		uint32_t nvi = (uint32_t) (nv+ 0.5);
		if(nvi > max_value) nvi = max_value;
		vbits.at(i) = nvi;
	}
}




void dynamic_mc_model::sqroot_bits_to_width_encoding(const std::vector<uint32_t> & sqbits, uint32_t target_total, std::vector<uint32_t> & widths) {
	std::vector<size_t> vbits(sqbits.size());
	widths.resize(sqbits.size());
	size_t total = 0;
	for(size_t i=0; i<sqbits.size(); ++i) {
		vbits.at(i) = sqbits.at(i) * sqbits.at(i);
		total += vbits.at(i);
	}

	double scale = (double)  target_total / (double) total;

	size_t real_total = 0;
	for(size_t i=0; i<vbits.size(); ++i) {
		widths.at(i) = (uint32_t)(scale * (double)(vbits.at(i)));
		if(widths.at(i) < ENCODING_MIN_WIDTH) widths.at(i) = ENCODING_MIN_WIDTH;
		real_total+=widths.at(i);
	}

	int correction = target_total - real_total; // we have to add that number to real_total
	int correction_done = 0;
	double correction_factor = (double) correction / (double) real_total;
	int extra_correction = 1; // correct some more to be sure to get it done
	if(correction < 0) extra_correction = -1;
	for(size_t i=0; i<widths.size(); ++i) {
		int this_correct = (int) widths.at(i) * correction_factor + extra_correction;
		if(this_correct < 0 && widths.at(i) + this_correct < ENCODING_MIN_WIDTH) {
			this_correct = ENCODING_MIN_WIDTH - widths.at(i);
		}
		if(correction < 0) {
			if(correction_done + this_correct <= correction) {
				this_correct = correction - correction_done;
			}
		} else {
			if(correction_done + this_correct >= correction) {
				this_correct = correction - correction_done;
			}
		}
		

		correction_done += this_correct;
		widths.at(i) += this_correct;
		if(correction_done == correction) break;
	}
}


void dynamic_mc_model::bits_to_width_encoding(const std::vector<uint32_t> & bits, uint32_t target_total, std::vector<uint32_t> & widths) {
	widths.resize(bits.size());
	size_t total = 0;
	for(size_t i=0; i<bits.size(); ++i) {
		total += bits.at(i);
	}

	double scale = (double) target_total / (double) total;
		
//	std::cout << " sc " << scale << " total " << total;

	size_t real_total = 0;
	for(size_t i=0; i<bits.size(); ++i) {
		widths.at(i) = (uint32_t)(scale * (double)(bits.at(i)));
		if(widths.at(i) < ENCODING_MIN_WIDTH) widths.at(i) = ENCODING_MIN_WIDTH;
		real_total+=widths.at(i);
//		std::cout << " b " << bits.at(i) << " w " << widths.at(i) << std::endl;
	}

	int correction = target_total - real_total; // we have to add that number to real_total

//	std::cout << " real total " << real_total << " corr " << correction << std::endl;

	int correction_done = 0;
	double correction_factor = (double) correction / (double) real_total;
	int extra_correction = 1; // correct some more to be sure to get it done
	if(correction < 0) extra_correction = -1;
	for(size_t i=0; i<widths.size(); ++i) {
		int this_correct = (int) widths.at(i) * correction_factor + extra_correction;
		if(this_correct < 0 && widths.at(i) + this_correct < ENCODING_MIN_WIDTH) {
			this_correct = ENCODING_MIN_WIDTH - widths.at(i);
		}
		if(correction < 0) {
			if(correction_done + this_correct <= correction) {
				this_correct = correction - correction_done;
			}
		} else {
			if(correction_done + this_correct >= correction) {
				this_correct = correction - correction_done;
			}
		}
		

		correction_done += this_correct;
		widths.at(i) += this_correct;
		if(correction_done == correction) break;
	}
}







/*
void dynamic_mc_model::markov_chain(size_t & accession, size_t  level){
	assert(level != 0);
	for(size_t k = 0; k <data.getAcc(accession).size(); k++){
		const dnastring & sequence = data.getSequence(data.getAcc(accession).at(k));
		std::cout << "seq length "<<sequence.length() << " level "<< level << std::endl;
		for(size_t i = 0 ; i < sequence.length(); i++ ){
			char schr = sequence.at(i);
			size_t s = dnastring::base_to_index(schr);
			std::stringstream context;		
			for(size_t j = level; j>0; j--){
				char rchr;
				signed long temp = i - j;
				if(temp >= 0){
					rchr = sequence.at(i-j);
				}else{
					rchr = 'A';

				}
				context << rchr;
			}
			std::string seq;
			context>>seq;
			std::map <std::string, std::vector<double> >::iterator it= sequence_successive_bases.at(accession).find(seq);
			if(it==sequence_successive_bases.at(accession).end()){
				sequence_successive_bases.at(accession).insert(std::make_pair(seq, std::vector<double>(5,1)));
				it= sequence_successive_bases.at(accession).find(seq);
			}
			it->second.at(s)++;	
		}
	}
	for(std::map <std::string, std::vector<double> >::iterator it= sequence_successive_bases.at(accession).begin();it!=sequence_successive_bases.at(accession).end();it++){
		std::string seq = it->first;
		int total = 0;
		std::vector<double> & base = sequence_successive_bases.at(accession).at(seq);
		for(size_t j = 0; j<5;j++){
			total += base.at(j);
		}
		for(size_t k=0; k<5;k++){
			base.at(k) = -log2(base.at(k)/total);
			create_cost.at(accession).at(k) += base.at(k);
		}
	}
}
*/
void dynamic_mc_model::calculate_sequence_high(size_t & accession, const std::map<std::string, std::vector<size_t> > & counts){//TODO level0
/*	size_t level = accLevel.at(accession);
	if(level > 0){
		//make_all_patterns(level);
		for(std::map<std::string, std::vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){
			std::string pattern = it ->first;
			std::map<std::string,std::vector<double> >::iterator it1=sequence_successive_bases.at(accession).find(pattern);
			if(it1 != sequence_successive_bases.at(accession).end()){
				for(size_t n = 0; n <5; n ++){
					it->second.at(n) = it1 ->second.at(n);
				}
			}else{
				for(size_t n =0; n < 5 ; n++){
					it -> second.at(n) = -log2(0.2);
				}
			}
		}
		std::map <std::string, std::vector<unsigned int> >lower_bound;
		for(std::map<std::string, std::vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){	
			std::vector<double> num(5,0);
			std::vector<unsigned int> low(5,0);
			std::vector<unsigned int> high_value(5,0);
			unsigned int l = 0;
			size_t bit = 12;
			std::string current_pattern = it->first;
                        high.at(accession).insert(std::make_pair(current_pattern,std::vector<unsigned int>(5,0)));
			std::map<std::string, std::vector<unsigned int> >::iterator it1=high.at(accession).find(current_pattern);
			assert(it1 != high.at(accession).end());
			for (size_t j=0; j < 5;j++){
				low.at(j) = l;
				num.at(j)=it->second.at(j);
				int power_of_two = exp((-num.at(j))*log(2))*powersOfTwo.at(bit);
				high_value.at(j) = l + power_of_two;
				l = high_value.at(j);
			}
			for(size_t j = 0; j < 5 ; j++){
				if(high_value.at(j)==low.at(j)){
					for(size_t m = 0; m < 4 ; m++){
						int power_of_two = exp((-num.at(m))*log(2))*(powersOfTwo.at(bit)-5);
						high_value.at(m)= low.at(m)+power_of_two+1;
						low.at(m+1)=  high_value.at(m);
					}
					high_value.at(4) = low.at(4)+exp((-num.at(4))*log(2))*(powersOfTwo.at(bit)-5)+1;
					break;
				}
			}
			for(size_t j = 0; j < 5; j++){
				it1->second.at(j)=high_value.at(j);
			}
		}
	}else{// TODO if level 0 needs to be added
		

	}

*/
}



void dynamic_mc_model::model_bits_to_full_high_values(const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & bitmodel, uint32_t target_total, size_t target_pattern_length, const std::vector<unsigned char> & alphabet,
	std::map< std::string, std::vector<uint32_t> > & hv_results) {
	std::string seq = "";
	std::set<std::string> pattern_set;

	std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > lsorted_inp(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL + 1);
	for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = bitmodel.begin(); it!=bitmodel.end(); ++it) {
		const std::string & cont = it->first;
		size_t len = cont.length();
		lsorted_inp.at(len).insert(make_pair(it->first, it->second));
	}

	for(size_t l=0; l < MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; l++) {
		// for each context of length l, copy data to all derived contexts of length l+1 (any character added to the front)
		for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::iterator it = lsorted_inp.at(l).begin(); it!=lsorted_inp.at(l).end(); ++it) {
			const std::string & cont = it->first;
			std::pair<unsigned char, std::vector<uint32_t> > pdata = it->second;

			for(size_t i=0; i<alphabet.size(); ++i) {
				std::string newcont;
				newcont.append(1, alphabet.at(i));
				newcont.append(cont);
				lsorted_inp.at(newcont.length()).insert(std::make_pair(newcont, pdata)); 
			}
		}
	}

	// copy all to result map
	for(size_t l=0; l < MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; l++) {
		for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::iterator it = lsorted_inp.at(l).begin(); it!=lsorted_inp.at(l).end(); ++it) {
			const std::string & cont = it->first;
			const std::vector<uint32_t> & pdata = it->second.second;

			hv_results.insert(std::make_pair(cont, pdata));
		}
	}
} 




/*
void dynamic_mc_model::make_all_patterns(const size_t& level, ){
	std::string seq = "";
	std::set<std::string> pattern;
	assert(level !=0 );
	for(size_t i = 0 ; i < level ; i++){
		seq += dnastring::index_to_base(0);
	}
	pattern.insert(seq);
	for(size_t i = 0; i < level;i++){	
		std::set<std::string> intermediate_pattern;
		for(size_t j = 0; j <5; j++){
			for(std::set<std::string>::iterator pat = pattern.begin();pat != pattern.end();++pat){
				std::string seq1 = *pat;
				seq1.at(level-1-i)=dnastring::index_to_base(j);	
				intermediate_pattern.insert(seq1);
			}
		}
		for(std::set<std::string>::iterator pat1 = intermediate_pattern.begin(); pat1 != intermediate_pattern.end();++pat1){
			std::string seq2 = *pat1;
			pattern.insert(seq2);
		}
	}	
	for(std::set<std::string>::iterator it = pattern.begin();it !=pattern.end();it++){
		std::string seq3 = *it;
		all_the_patterns.insert(std::make_pair(seq3,std::vector<double>(5,0)));
	}	
}*/


/*
	Write all bits to out

	If bits.size() is not divisible by 8 we add some zero bits

*/

void dynamic_mc_model::bits_write(const std::vector<bool> & bits, std::ofstream & out) {
	for(size_t n=0; n < bits.size(); n+=8) {
		unsigned char outc = 0;
		// this loop computes a byte from 8 bits
		for(size_t m=n; m<n+8; m++) {
			unsigned char nextbit = 0;
			if(m<bits.size()) nextbit = bits.at(m);
			outc += (1<<(m-n))*nextbit;
		}
		out.put(outc);
	}




}

/*
	read num_bits_to_read from in and put them into bits

	if num_bits_to_read is not divisible by 8 we will still read whole bytes and discard a few bits

*/
void dynamic_mc_model::bits_read(std::ifstream & in, std::vector<bool> & bits, size_t num_bits_to_read) {
	bits.resize(num_bits_to_read + 7);
	for(size_t num_read = 0; num_read < num_bits_to_read; num_read+=8) {
		unsigned char inc;
		binary_read(in, inc);
		for(size_t m=0; m<8; ++m) {
			if(((1<<m) & inc) == (1<<m)) {
				bits.at(num_read + m) = 1;
			} else {
				bits.at(num_read + m) = 0;
			}
		}
	}



	bits.resize(num_bits_to_read);
}

/*
	convert a vector of integers into a vector of bits
	using bits_per_int lowest significance bits per input integer

*/
void dynamic_mc_model::ints_to_bits(const std::vector<uint32_t> & ints, std::vector<bool> & bits, size_t bits_per_int) {
	bits.resize(ints.size() * bits_per_int);
	for(size_t i=0; i<ints.size(); ++i) {
		for(size_t m=0; m<bits_per_int; ++m) {
			if(((1<<m) & ints.at(i)) == (1<<m)) {

				bits.at(i*bits_per_int+m) = 1;
			} else {
				bits.at(i*bits_per_int+m) = 0;
			}
		}
#ifndef NDEBUG
		for(size_t m=bits_per_int; m<32; ++m) {
			assert(0== ( ints.at(i) & (1<<m) ) );
		}
#endif
	}


}

/*
	convert vector of bits into integers
	we create num_ints integers using bits_per_int bits for each, we use those bits as lowest significance bits
		
*/
void dynamic_mc_model::bits_to_ints(const std::vector<bool> & bits, std::vector<uint32_t> & ints, size_t num_ints, size_t bits_per_int) {
	assert(bits.size()>=num_ints*bits_per_int);
	ints.resize(num_ints);
	for(size_t i=0; i<num_ints; ++i) {
		uint32_t next_int = 0;
		for(size_t m=0; m<bits_per_int; ++m) {
			if(bits.at(i*bits_per_int + m)) {
				next_int += (1<<m);
			}
		}
		ints.at(i) = next_int;

	}
}



/*
	TODO these binary read/write functions have to be change if we want to transfer files to big endian systems
*/
void dynamic_mc_model::binary_write(std::ofstream & out, size_t v) {
	out.write(reinterpret_cast<char*>(&v), sizeof(size_t));
}

void dynamic_mc_model::binary_write(std::ofstream & out, uint32_t v) {
	out.write(reinterpret_cast<char*>(&v), sizeof(uint32_t));
}


void dynamic_mc_model::binary_write(std::ofstream & out, const std::string & str) {
	out.write(str.c_str(), str.length());
	out.put((char) 0);

}


void dynamic_mc_model::binary_read(std::ifstream & in, size_t & v) {
	char buf[sizeof(size_t)];
	in.read(buf, sizeof(size_t));
	check_read_error(in);
	size_t * vp = reinterpret_cast<size_t*>(buf);
	v = *vp;	
}


void dynamic_mc_model::binary_read(std::ifstream & in, uint32_t & v) {
	char buf[sizeof(uint32_t)];
	in.read(buf, sizeof(uint32_t));
	check_read_error(in);
	uint32_t * vp = reinterpret_cast<uint32_t*>(buf);
	v = *vp;	
}


void dynamic_mc_model::binary_read(std::ifstream & in, std::string & str) {
	std::getline(in, str, (char)0);
	check_read_error(in);
}



void dynamic_mc_model::binary_read(std::ifstream & in, unsigned char & c) {
	char buf[sizeof(unsigned char)];
	in.read(buf, sizeof(unsigned char));
	check_read_error(in);
	unsigned char * cp = reinterpret_cast<unsigned char *>(buf);
	c= *cp;
}



void dynamic_mc_model::mapmodel_write(std::ofstream & outs, const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level) {
	// sequence contexts:
	std::set< std::string> contexts; //  (this means we have to look at all extensions in the next level
		
	std::string empty_context;
	contexts.insert(empty_context);


	size_t nomodel_counter = 0;
	// model_index - NUM_ENCODING_TYPES + 1 is number of nomodel events
	size_t MAX_NOMODEL_NUM = 256 - NUM_ENCODING_TYPES; // nomodel counter fits into same byte as all encoding types

	for(size_t cur_len = 0; cur_len <= max_level; ++cur_len) {
		std::set<std::string> next_contexts;
		for(std::set<std::string>::iterator it = contexts.begin(); it!=contexts.end(); ++it) {
			std::string cont = *it;
			assert(cont.length() == cur_len);
			std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator findcont = mapmodel.find(cont);
			if(findcont == mapmodel.end()) { // no model for this context, look for larger contexts
					for(size_t i=0; i<alphabet.size(); ++i) {
						std::string lcont;
						lcont.append(1, alphabet.at(i));
						lcont.append(cont);
						next_contexts.insert(lcont);

					}
					

					if(nomodel_counter < MAX_NOMODEL_NUM - 1) {
						nomodel_counter++;
					} else {
						outs.put((unsigned char) 255);
						nomodel_counter = 0;
					}
					// write no model symbol
				//	outs.put((unsigned char) NUM_ENCODING_TYPES);
				} else {

					if(nomodel_counter > 0) {
						size_t nmchar = nomodel_counter -1 + NUM_ENCODING_TYPES;
						outs.put((unsigned char) nmchar);
						nomodel_counter = 0;
					}


					unsigned char model_index = findcont->second.first;
					std::vector<uint32_t> model_values = findcont->second.second;
					// write model type (how many bits per value)
					outs.put(model_index);
					bool use_sqroot;
					size_t num_bits;
					get_model_type(model_index, use_sqroot, num_bits);
					// write the values
					vector<bool> bits;
					ints_to_bits(model_values, bits, num_bits);
					bits_write(bits, outs);




				}
			}
			contexts = next_contexts;
		}

}


void dynamic_mc_model::mapmodel_read(std::ifstream & in,  std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level)  {
	// sequence contexts:
	std::set< std::string> contexts; //  (this means we have to look at all extensions in the next level
		
	std::string empty_context;
	contexts.insert(empty_context);


	size_t nomodel_counter = 0;
	// model_index - NUM_ENCODING_TYPES + 1 is number of nomodel events
	size_t MAX_NOMODEL_NUM = 256 - NUM_ENCODING_TYPES; // nomodel counter fits into same byte as all encoding types




	for(size_t cur_len = 0; cur_len <= max_level; ++cur_len) {
		std::set<std::string> next_contexts;
		for(std::set<std::string>::iterator it = contexts.begin(); it!=contexts.end(); ++it) {
			std::string cont = *it;
			assert(cont.length() == cur_len);


			if(nomodel_counter > 0) {
				nomodel_counter--;
	// insert all child nodes so that we will exspect them later in correct order
				for(size_t i=0; i<alphabet.size(); ++i) {
					std::string lcont;
					lcont.append(1, alphabet.at(i));
					lcont.append(cont);
					next_contexts.insert(lcont);
				}



			} else {




			// now we have to test if there is a model for cont or if we have to proceed with its child nodes later
			unsigned char model_index;
			binary_read(in, model_index);

			if(model_index >= NUM_ENCODING_TYPES) {
				nomodel_counter = model_index - NUM_ENCODING_TYPES + 1;
				nomodel_counter--;




				// insert all child nodes so that we will exspect them later in correct order
				for(size_t i=0; i<alphabet.size(); ++i) {
					std::string lcont;
					lcont.append(1, alphabet.at(i));
					lcont.append(cont);
					next_contexts.insert(lcont);
				}

			} else {
					
					

				bool use_sqroot;
				size_t num_bits;
				get_model_type(model_index, use_sqroot, num_bits);
				size_t read_bits = num_bits * alphabet.size();
				std::vector<bool> bits;
				bits_read(in, bits, read_bits);
				std::vector<uint32_t> values;
				bits_to_ints(bits, values, alphabet.size(), num_bits);
				mapmodel.insert(std::make_pair(cont, std::make_pair(model_index, values ) ) );
			}

			}

		}
			contexts = next_contexts;
	}







}

/*
	this function writes all parameters necessary to recompute an identical model to a stream

	Structure

	* number of accessions (size_t)
	For each accession (
		* accession name, zero terminated
		* 2 sequence flag widths (uint32_t)
		* number of acc alignment flag widths

	)


*/
void dynamic_mc_model::write_parameters(std::ofstream & outs) const {//write accessions, their levels and high values to a file
	binary_write(outs, data.numAcc());
	
		 

	std::vector<unsigned char> salphabet;
	get_sequence_alphabet(salphabet);
	std::vector<unsigned char> aalphabet;
	get_alignment_alphabet(aalphabet);


	for(size_t acc = 0; acc < data.numAcc(); ++acc) {
//	for(size_t acc = 0; acc < 1; ++acc) {
		binary_write(outs, data.get_acc(acc));
		binary_write(outs, acc_sequence_all_bases_total.at(acc));
		binary_write(outs, acc_sequence_alignment_flag.at(acc));
		mapmodel_write(outs, sequence_models.at(acc),  salphabet, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
		
		for(size_t a2 = 0; a2<data.numAcc(); ++a2) {
			binary_write(outs, alignments_all_modification_total.at(acc).at(a2));
			mapmodel_write(outs, alignment_models.at(acc).at(a2), aalphabet, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);


		}	
	}
	


}

/*
	this is a test function for correctnes of paramter read/write

	it reads parameters from a file and verifies that those values are identical to
	those store in the class members
*/
void dynamic_mc_model::check_parameters_in_file(const std::string & file) const {
	std::ifstream in(file.c_str(), std::ios::binary);
	std::vector<std::string> acc_names;
	std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > >  sequence_models;
	std::vector<uint32_t> acc_sequence_all_bases_total;
	std::vector<uint32_t> acc_sequence_alignment_flag;
	std::vector<std::vector<uint32_t> > alignments_all_modification_total;
	std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > alignment_models;

	read_paramters(in, acc_names, 
			sequence_models,
			acc_sequence_all_bases_total,
			acc_sequence_alignment_flag,
			alignments_all_modification_total,
			alignment_models);

// check all flags and names:
	assert(data.numAcc() == acc_names.size());
	size_t na = data.numAcc();
	assert(na == acc_sequence_all_bases_total.size());
	assert(na == acc_sequence_alignment_flag.size());
	assert(na == alignments_all_modification_total.size());
	for(size_t acc=0; acc<acc_names.size(); ++acc) {
		assert( 0 == data.get_acc(acc).compare(acc_names.at(acc)) );
		assert(na == alignments_all_modification_total.at(acc).size());
		for(size_t a2 = 0; a2 < na; ++a2) {
			assert( this->alignments_all_modification_total.at(acc).at(a2) == 
                                     alignments_all_modification_total.at(acc).at(a2) );
		}
	}

// check dynamic context length sequence model
	assert(this->sequence_models.size() == sequence_models.size());
	for(size_t acc = 0; acc < na; ++acc) {
		assert(this->sequence_models.at(acc).size() == sequence_models.at(acc).size());
		for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = this->sequence_models.at(acc).begin();
			it!=this->sequence_models.at(acc).end(); ++it) {
			std::string cont = it->first;
			unsigned char model_index = it->second.first;
			std::vector<uint32_t> model_values = it->second.second;
			std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator t_it = sequence_models.at(acc).find(cont);
			assert(t_it!=sequence_models.at(acc).end());
			unsigned char t_model_index = t_it->second.first;
			std::vector<uint32_t> t_model_values = t_it->second.second;
			assert(model_index == t_model_index);
			assert(model_values.size() == t_model_values.size());
			
			for(size_t i=0; i<model_values.size(); ++i) {
				assert(t_model_values.at(i) == model_values.at(i));
			}
		}
	}

// check alignment models
	assert(this->alignment_models.size() == alignment_models.size());
	for(size_t acc = 0; acc < na; ++acc) {
		assert(this->alignment_models.at(acc).size() == alignment_models.at(acc).size());
		for(size_t a2=0; a2 < na; ++a2) {
			assert(this->alignment_models.at(acc).at(a2).size() == alignment_models.at(acc).at(a2).size());
			for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = this->alignment_models.at(acc).at(a2).begin();
			it!=this->alignment_models.at(acc).at(a2).end(); ++it) {
				std::string cont = it->first;

				unsigned char model_index = it->second.first;
				std::vector<uint32_t> model_values = it->second.second;
				std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator t_it = alignment_models.at(acc).at(a2).find(cont);
				assert(t_it != alignment_models.at(acc).at(a2).end());
				unsigned char t_model_index = t_it->second.first;
				std::vector<uint32_t> t_model_values = t_it->second.second;
				assert(model_index == t_model_index);
				assert(model_values.size() == t_model_values.size());
			
				for(size_t i=0; i<model_values.size(); ++i) {
					assert(t_model_values.at(i) == model_values.at(i));
				}
			}
		}
	}




	in.close();

}



void dynamic_mc_model::check_read_error(std::ifstream & in) {
	if(in.bad()) {
		std::cerr << "Error while reading from file." << std::endl;
		std::exit(1);
	}
}
	

void dynamic_mc_model::get_sequence_alphabet(std::vector<unsigned char> & alphabet) {
	alphabet.resize(5);
	alphabet.at(dnastring::base_to_index('A')) = 'A';
	alphabet.at(dnastring::base_to_index('T')) = 'T';
	alphabet.at(dnastring::base_to_index('C')) = 'C';
	alphabet.at(dnastring::base_to_index('G')) = 'G';
	alphabet.at(dnastring::base_to_index('N')) = 'N';
}


void dynamic_mc_model::get_alignment_alphabet(std::vector<unsigned char> & alphabet) {
	alphabet.resize(NUM_MODIFICATIONS);
	for(size_t i=0; i< NUM_MODIFICATIONS; ++i) {
		alphabet.at(i) = i;
	}

}






/*
	read parameters and store them into temporary variables

	written this way so that we can have a test function which verifies that writing and reading again results in identical data
*/
void dynamic_mc_model::read_paramters(std::ifstream & in, 
	std::vector<std::string> & acc_names, 
	std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > & sequence_models,
	std::vector<uint32_t> & acc_sequence_all_bases_total,
	std::vector<uint32_t> & acc_sequence_alignment_flag,
	std::vector<std::vector<uint32_t> > & alignments_all_modification_total,
	std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > & alignment_models
) const {  // TODO set accession names in data
	size_t num_acc;
	binary_read(in, num_acc);
	acc_names.resize(num_acc);
	alignments_all_modification_total.resize(num_acc);
	acc_sequence_all_bases_total.resize(num_acc);
	acc_sequence_alignment_flag.resize(num_acc);

	std::vector<unsigned char> sequence_alphabet;
	get_sequence_alphabet(sequence_alphabet);
	std::vector<unsigned char> alignment_alphabet;
	get_alignment_alphabet(alignment_alphabet);
	sequence_models.resize(num_acc);
	alignment_models.resize(num_acc);

	for(size_t acc=0; acc < num_acc; ++acc) {

		alignment_models.at(acc).resize(num_acc);
		// read accession name
		binary_read(in, acc_names.at(acc));
		
		// flag values
		binary_read(in, acc_sequence_all_bases_total.at(acc));
		binary_read(in, acc_sequence_alignment_flag.at(acc));
		alignments_all_modification_total.at(acc).resize(num_acc);

		mapmodel_read(in, sequence_models.at(acc), sequence_alphabet, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);

		for(size_t a2 = 0; a2<data.numAcc(); ++a2) {
			binary_read(in, alignments_all_modification_total.at(acc).at(a2));
			mapmodel_read(in, alignment_models.at(acc).at(a2), alignment_alphabet, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);


		}	

	}
}



/*
 Read parameters and store them into this object
*/
void dynamic_mc_model::set_patterns(std::ifstream & in){
	std::vector<std::string> acc_names; // TODO set in data
	read_paramters(in, acc_names, 
			this->sequence_models,
			this->acc_sequence_all_bases_total,
			this->acc_sequence_alignment_flag,
			this->alignments_all_modification_total,
			this->alignment_models);













/*	size_t bit = 12;
	char c;
	char h;
	c= in.get();	
	while(c != 8){
		size_t accession = 0;
		std::string acc;
		std::stringstream s;
		while(c != 0){
			s << c;
			c = in.get();
		}
		s >> acc;
		size_t retrieved_level;
		in>>retrieved_level;
		data.set_accession(acc);//During the decoding we have no access to the fasta file, thus we need to set the accession names and numbers in the data class
		accession = data.get_acc_id(acc);
		set_acc_level(retrieved_level);
		size_t level = get_acc_level(accession);//just to test and see if returns the same thing.
		std::cout<< "ret_level: " << retrieved_level << " level "<< level <<std::endl;
		//make_all_patterns(level);// TODO later on change it to retrieved_level	
		for(std::map<std::string,std::vector<double> >::const_iterator it= all_the_patterns.begin(); it!= all_the_patterns.end();it++){
			std::string pattern = it ->first;
			high.at(accession).insert(std::make_pair(pattern, std::vector<unsigned int>(5,0)));
		}
		for(std::map<std::string,std::vector<unsigned int> >::iterator it= high.at(accession).begin(); it!= high.at(accession).end();it++){
			std::vector<bool> binary_high_value(0);
			size_t bound = (5*bit)/8;
			for(size_t j = 0 ; j < bound ; j++){ 
				h=in.get();
				size_t H = size_t(h); 
				for(size_t k = 0; k < 8 ; k++){
					binary_high_value.push_back(H%2);
					H = H/2;
				}
			}
			size_t counter = 0;
			for(size_t i = 0; i < binary_high_value.size()-bit;i++){
				unsigned int high_value = 0;					
				for(size_t j =i; j < i+bit; j++){
					high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
				}
				i=i+bit-1;
				it -> second.at(counter)=high_value;
				counter = counter + 1;
			}
		}
		c= in.get();	
	}

*/
}


void dynamic_mc_model::alignment_counts_add(std::map<std::string, std::vector<size_t> > & countsmap, const std::map<std::string, std::vector<size_t> > & newcounts) {
	for(std::map<std::string, std::vector<size_t> >::const_iterator it = newcounts.begin(); it!=newcounts.end(); ++it) {
		std::map<std::string, std::vector<size_t> >::iterator findit = countsmap.find(it->first);
		if(findit==countsmap.end()) {
			std::string context = it->first; 
			countsmap.insert(std::make_pair(context, it->second));
		} else {
			std::vector<size_t> & c = findit->second;
			const std::vector<size_t> & d = it->second;
			for(size_t i=0; i<c.size(); ++i) {
				c.at(i)+=d.at(i);
			}
		}
	}
}


/*
	train alignment/modification with a dynamic model
*/


void dynamic_mc_model::train_alignment_model() {

// for better multithreading we sort alignments using accession their accession indices
	size_t numa = data.numAcc();
	std::vector<std::vector<std::vector<const pw_alignment *> > > sorted_als(numa, std::vector<std::vector<const pw_alignment * > >(numa));
	vector<unsigned char> alphabet;
	get_alignment_alphabet(alphabet);

	alignments_all_modification_total = std::vector<std::vector<uint32_t> > (numa, std::vector<uint32_t>(numa));


	for(size_t i=0; i<data.numAlignments(); ++i) {
		const pw_alignment & al = data.getAlignment(i);
		size_t r1 = al.getreference1();
		size_t r2 = al.getreference2();
		size_t a1 = data.accNumber(r1);
		size_t a2 = data.accNumber(r2);
		sorted_als.at(a1).at(a2).push_back(&al); // al from a1 to a2
	}

	// all the (directed) pairs of accessions (including within each) where alignments can be found
	std::vector<std::pair<size_t, size_t> > acc_compare_order;
	for(size_t i=0; i<numa; ++i) {
		for(size_t j=0; j<numa; ++j)  acc_compare_order.push_back(std::make_pair(i,j));
	}

	
	size_t max_batch_size = 10000;
	std::vector< std::vector< std::map<std::string, std::vector<size_t> > > > result_counts(
		numa,std::vector< std::map<std::string, std::vector<size_t> > > (numa, 
				  std::map<std::string, std::vector<size_t> >()));

	alignment_models = std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > >(numa, 
			   std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > >(numa));

	std::vector<std::vector<size_t> > alignment_length_sums(numa, std::vector<size_t>(numa, 0));

	for(size_t comp=0; comp < acc_compare_order.size(); ++comp) {
		size_t a1 = acc_compare_order.at(comp).first; // we look at modification from a1 to a2
		size_t a2 = acc_compare_order.at(comp).second;
		// modification from a1 to a2 consist either of a 1-2 reading of a1,a2 or a 2-1 reading of a2,a1 alignments


		size_t lensum1 = 0;	
		size_t lensum2 = 0;

		size_t num = sorted_als.at(a1).at(a2).size();
		size_t batch_size = num / num_threads + 1;
		if(batch_size > max_batch_size) batch_size = max_batch_size;
		size_t num_batchs = num / batch_size + 1;

		size_t alcounter = 0;
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
		for(size_t b = 0; b < num_batchs; ++b) {
			size_t from = b*batch_size;
			size_t to = from + batch_size;
			if(to > num) to = num;
			
			reader_counter local_counts;
			context_counter local_counter(local_counts,  MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			



			for(size_t i=from; i<to; ++i) {
				const pw_alignment * al = sorted_als.at(a1).at(a2).at(i);

				lensum1+= al->alignment_length();
				alcounter++;
				local_counter.read_alignment_1_2(*al);
			}
#pragma omp critical(counts)
{
			alignment_counts_add(result_counts.at(a1).at(a2), local_counts.countsmap);
			alignment_length_sums.at(a1).at(a2)+=lensum1;
}

		}
		assert(alcounter == sorted_als.at(a1).at(a2).size());


		num = sorted_als.at(a2).at(a1).size();
		batch_size = num / num_threads + 1;
		if(batch_size > max_batch_size) batch_size = max_batch_size;
		num_batchs = num / batch_size + 1;
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
		for(size_t b = 0; b < num_batchs; ++b) {
			size_t from = b*batch_size;
			size_t to = from + batch_size;
			if(to > num) to = num;
			reader_counter local_counts;
			context_counter local_counter(local_counts,  MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			



			for(size_t i=from; i<to; ++i) {
				const pw_alignment * al = sorted_als.at(a2).at(a1).at(i);
				local_counter.read_alignment_2_1(*al);

				lensum2+=al->alignment_length();
			}
#pragma omp critical(counts)
{
			alignment_counts_add(result_counts.at(a2).at(a1), local_counts.countsmap);
			alignment_length_sums.at(a2).at(a1)+=lensum2;
}


		}

	} // for accession comparison

// #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
	for(size_t comp=0; comp < acc_compare_order.size(); ++comp) {
		size_t a1 = acc_compare_order.at(comp).first; // we look at modification from a1 to a2
		size_t a2 = acc_compare_order.at(comp).second;
// TODO XXX it would make sense to recompute the width of alignment end flags per accession pair after clustering
/*
	the ideal width of alignment end flag is (average alignment length after cutting)^-1
	at the moment we assume after cutting alignments we have lengths between 25 and 50 length percentiles

*/

		std::multiset<size_t>  al_len;
		for(size_t i=0; i<sorted_als.at(a1).at(a2).size(); ++i) {
			al_len.insert(sorted_als.at(a1).at(a2).at(i)->alignment_length());
		}
		size_t num = 1; // safeguard values to always have sensible results
		size_t len = 100; 
		size_t pos = 0;
		for(std::multiset<size_t>::iterator it = al_len.begin(); it!=al_len.end(); ++it) {
			double percentile = (double) pos / (double) al_len.size();
			if(percentile > 0.25 && percentile < 0.5) {
				num++;
				len+=*it;
			}			
			pos++;
		}
		double average_length = (double) len / (double) num;
		double end_prob = 1 / average_length;
		vector<double> distri(2);
		distri.at(0) = 1 - end_prob;
		distri.at(1) = end_prob;
		vector<uint32_t> flag_bits;
		distribution_to_bits(distri, MAX_ENCODING_WIDTH_BITS, flag_bits);
		vector<uint32_t> flag_widths;
		bits_to_width_encoding(flag_bits, TOTAL_FOR_ENCODING, flag_widths);
		std::cout <<"align acc " << a1 << " to " << a2 << ": encoding widhts: any mod: "<< flag_widths.at(0) << " (" << distri.at(0)<<") " <<   " alignment end: " << flag_widths.at(1) << " ("<<distri.at(1)<<") "<< std::endl;
		

		uint32_t new_total = flag_widths.at(0);
		alignments_all_modification_total.at(a1).at(a2) = new_total;

		std::map<std::string, std::vector<size_t> > & counts = result_counts.at(a1).at(a2);

/*
		size_t checks = 0;
		for(std::map<std::string, std::vector<size_t> >::iterator it = counts.begin(); it!=counts.end(); ++it) {
			std::vector<size_t> c = it->second;
			for(size_t i=0; i<c.size(); ++i) checks+=c.at(i);
			std::string cont = it->first;
			std::cout << "CO " << cont.size();
			for(size_t i=0; i<cont.length(); ++i) std::cout << ": " << (size_t) cont.at(i);
			std::cout << " === ";
			for(size_t i=0; i<c.size(); ++i) std::cout << ": " << c.at(i);
			std::cout << std::endl;
		}
*/		
		std::map<std::string, std::pair<unsigned char, std::vector<uint32_t> > > mod_model;

		std::cout << " run context selector on " << counts.size() << " contexts " << std::endl;
		double enc_cost = context_selector(counts, alphabet, new_total, mod_model, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL + 1);
		std::cout << " select " << mod_model.size() << " contexts, cost: " << enc_cost << " ( " << (enc_cost/(double) alignment_length_sums.at(a1).at(a2))<< " bits/alignment column)" << std::endl;
		alignment_models.at(a1).at(a2) = mod_model;
	}
	compute_alignment_model();

}

/*
	new functors - separated per acc

	new computing_modif function 
		- why is it slow

	go over al - shift pattern seq1 (now copy from seq)

	find s1n, s2n -> skip gaps on ref1, why?  for seq1 context? why at level?


	keep: find length, set pattern and n
*/


void dynamic_mc_model::recursive_al_model(){//TODO again the same question about level 0
/*	counting_functor functor(data);
	std::vector<std::vector<std::vector<std::vector<double> > > >modification(data.numAcc(), std::vector<std::vector<std::vector<double> > >(data.numAcc(),std::vector<std::vector<double> > (6, std::vector<double>(6,1) )));
	//set an initial mod_cost for the alignments and go through the while loop as long as it decreases!
	for(size_t i = 0 ; i< data.numAlignments() ; i++){
		const pw_alignment & p = data.getAlignment(i);
		size_t acc1 = data.accNumber(p.getreference1());	
		size_t acc2 = data.accNumber(p.getreference2());
		for(size_t j = 0 ; j< p.alignment_length() ; j++){
			char s1ch;
			char s2ch;
			p.alignment_col(j, s1ch, s2ch);
			size_t s1 = dnastring::base_to_index(s1ch);
			size_t s2 = dnastring::base_to_index(s2ch);
			modification.at(acc1).at(acc2).at(s1).at(s2)+=1.0;
			modification.at(acc2).at(acc1).at(s2).at(s1)+=1.0;
		}
	}
	for(size_t acc1=0; acc1 < data.numAcc(); ++acc1) {
		for(size_t acc2=0; acc2 < data.numAcc(); ++acc2) {
			double old_cost = 0.0;
			for(size_t j=0; j<6; ++j) {
				double sum = 0;
				for(size_t k=0; k<6; ++k) {
					sum+=modification.at(acc1).at(acc2).at(j).at(k);
				}	
				for(size_t k=0; k<6; ++k) {
					modification.at(acc1).at(acc2).at(j).at(k) = -log2(modification.at(acc1).at(acc2).at(j).at(k)/sum);	
				}
			}
			for(size_t j=0; j<6; ++j) {
				for(size_t k=0; k<6; ++k) {
				old_cost +=modification.at(acc1).at(acc2).at(j).at(k);
				}
			}		
			double new_cost =0.0;
			size_t level = 1;
			while(new_cost<old_cost){//what is the reasonable way of defining new_cost? Decreasing mod_cost or increasing gain? Increasing gain is more complicated though
				if(level != 1){//update the mod_cost
					old_cost = new_cost;
				}
				make_all_alignment_pattern(level);			
				for(std::set<std::string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){
					std::string pattern = *it;	
					functor.create_context(acc1,acc2,pattern);			
				}
				// count all contexts over all the alignments
				for(size_t k = 0; k < data.numAlignments(); k++){
					const pw_alignment & p = data.getAlignment(k);
					if( data.accNumber(p.getreference1())==acc1 &&  data.accNumber(p.getreference2())==acc2){
						count_context_OneToTwo(p,level,functor);
						count_context_TwoToOne(p,level,functor); 
					}
				}
				functor.total_context(acc1,acc2);
				for(std::map <std::string, std::vector<size_t> >::const_iterator it= functor.get_context(acc1,acc2).begin();it!=functor.get_context(acc1,acc2).end();it++){
					std::string seq1 = it->first;
					const std::vector<size_t> & base = functor.get_context(acc1,acc2).at(seq1);
					for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++) {
						std::map <std::string, std::vector<double> >::iterator it1= mod_cost.at(acc1).at(acc2).find(seq1);
						if(it1==mod_cost.at(acc1).at(acc2).end()) {
							mod_cost.at(acc1).at(acc2).insert(std::make_pair(seq1, std::vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
							it1= mod_cost.at(acc1).at(acc2).find(seq1);
						}
						it1->second.at(k)=-log2(base.at(k)/functor.get_total(acc1,acc2,seq1));
					}
				}
				level = level + 1;
			}
			al_level.at(acc1).at(acc2)=level-2;
			std::cout<< "level between acc1 " << acc1 << " and acc2 "<< acc2 <<" is "<< level -2 <<std::endl;
			if(level -2 > 0){
				make_all_alignment_pattern(al_level.at(acc1).at(acc2));
			}
			calculate_alignment_high(acc1,acc2);
		}
	}
*/
}
void dynamic_mc_model::make_all_alignment_pattern(const size_t & level){ // it is only used when level > 0
/*	assert(level != 0);
	std::string context; 
	std::set<std::string>pattern;
	// Alignment_level is makov chain level for alignments
	for(size_t i = 0 ; i < level ; i++){
		context += (char)0;
	}
	pattern.insert(context);
	// this will create about (ND+NK+10)^Alignment_length patterns:
	for(size_t i =0; i < level; i++) {
		std::set<std::string> intermediate_pattern;
		for(size_t j = 0; j <NUM_DELETE+NUM_KEEP+10; j++){
			// For each current pattern: modify position i to character j
			for(std::set<std::string>::iterator it = pattern.begin(); it!= pattern.end();it++){
				std::string seq = *it;
				seq.at(i)=j;	
				// throw away some patterns to enforce decreasing order of binary encoding for num keep/delete
				size_t keepthispattern  = 1;
				if(i>0) {
					int mod, del, ins, keep;
					int mod_prev, del_prev, ins_prev, keep_prev;
					modification(seq.at(i), mod, del, ins, keep);
					modification(seq.at(i-1), mod_prev, del_prev, ins_prev, keep_prev);
					if(del > -1 && del_prev > -1) {
						if(del <= del_prev) { // enforce increasing order binary encoding
							keepthispattern = 0;
						}
					
					} else if (keep >-1 && keep_prev > -1) {
						if(keep <= keep_prev) {
							keepthispattern = 0;
						}
					}
				}
				if(keepthispattern)
					intermediate_pattern.insert(seq);
			//	std::cout<< "pattern: " << seq << std::endl;
			}
		}
		pattern.clear();
		for(std::set<std::string>::iterator it1 = intermediate_pattern.begin(); it1 != intermediate_pattern.end();++it1){
			std::string seq1 = *it1;
			pattern.insert(seq1);
		}
	} // for Alignment_level
	std::set<std::string> intermediate1_pattern;
	for(std::set<std::string>::iterator it = pattern.begin(); it != pattern.end();++it){
		for(size_t i = 0 ; i < 6 ; i++){
			std::string seq = *it;
			char c = i;
			std::string seq1 = seq + c;
		//	std::cout << "seq1: " << seq1 <<std::endl;
			intermediate1_pattern.insert(seq1);
		}
	}
	pattern.clear();
	for(std::set<std::string>::iterator it1 = intermediate1_pattern.begin(); it1 != intermediate1_pattern.end();++it1){
			std::string seq1 = *it1;
			pattern.insert(seq1);
	}
	all_alignment_patterns.clear();
	for(std::set<std::string>::iterator it = pattern.begin();it !=pattern.end();it++){
		std::string seq = *it;
		all_alignment_patterns.insert(seq);
	}
*/
}
void dynamic_mc_model::calculate_alignment_high(size_t & acc1, size_t & acc2){//TODO level0?
/*
	counting_functor functor(data);
	size_t level = al_level.at(acc1).at(acc2);
	if(level != 0){
		// Set all counts to 1 for all context/accession pairs
		make_all_alignment_pattern(level);
		for(std::set<std::string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){
			std::string pattern = *it;	
			functor.create_context(acc1, acc2, pattern);			
		}	
		// count all contexts over all the alignments
		for(size_t k = 0; k < data.numAlignments(); k++){
			const pw_alignment & p = data.getAlignment(k);
			if( data.accNumber(p.getreference1())==acc1 &&  data.accNumber(p.getreference2())==acc2){
				count_context_OneToTwo(p,level,functor);
				count_context_TwoToOne(p,level,functor); 
			}
		}
		functor.total_context(acc1,acc2);
		for(std::set<std::string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){
			std::vector<double> num(NUM_DELETE+NUM_KEEP+10,0);
			std::vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
			std::vector<unsigned int> high_value(NUM_DELETE+NUM_KEEP+10,0);
			unsigned int l = 0;
			size_t bit = 12; // number of bits to use for encoding event width todo shall i use a smaller number?
			std::string current_pattern= *it;
			// it1: high values for current pattern/accession pair
	                highValue.at(acc1).at(acc2).insert(std::make_pair(current_pattern,std::vector<unsigned int>(NUM_DELETE+NUM_KEEP+10,0)));
			std::map<std::string, std::vector<unsigned int> >::iterator it1=highValue.at(acc1).at(acc2).find(current_pattern);
			assert(it1 != highValue.at(acc1).at(acc2).end());
			// it3: get counts for current pattern/accession pair
			std::map <std::string, std::vector<size_t> >::const_iterator it3= functor.get_context(acc1,acc2).find(current_pattern);
			assert(it3 != functor.get_context(acc1,acc2).end());
			double total =  functor.get_total(acc1,acc2,current_pattern);
			for (size_t f=0; f < NUM_DELETE+NUM_KEEP+10;f++){
				low.at(f) = l;
				num.at(f)=it3->second.at(f);
				size_t rescaledNum = (num.at(f)/total)*(powersOfTwo.at(bit) - NUM_DELETE - NUM_KEEP - 11) + 1;
				assert(rescaledNum >= 1);
				assert(rescaledNum < powersOfTwo.at(bit));
				high_value.at(f) = l + rescaledNum;
				l = high_value.at(f);
			}
			// store high values
			for(size_t f = 0; f < NUM_DELETE+NUM_KEEP+10; f++){
				it1->second.at(f)=high_value.at(f);
			}
		}
	}else{ //In case level 0

	}
*/
}
void dynamic_mc_model::write_alignment_parameters(std::ofstream & outs){//since we alreay wrote the accession names on the file and can retrieve them from set parameters here we can only save acc_number not the name anymore.
/*	size_t bit = 12;
	for(size_t i = 0 ; i < data.numAcc(); i++){
		for(size_t j =0; j < data.numAcc(); j++){
			outs<<char(0);
			outs << i;
			outs<< j;
			size_t level = get_al_level(i,j);
			outs<< level;
			for(std::map<std::string, std::vector<unsigned int> >::iterator it= highValue.at(i).at(j).begin(); it != highValue.at(i).at(j).end(); it++){	
				std::vector<bool> bit_to_byte(0);
				for(size_t j = 0; j < NUM_DELETE+NUM_KEEP+10; j++){
					int h =it->second.at(j);
					for(size_t m = 0; m < bit; m++){
						bit_to_byte.push_back(h%2);
						h = h/2;
					}
				}
				for(size_t n =0; n < bit_to_byte.size()-8; n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= powersOfTwo.at(m-n)* bit_to_byte.at(m);
					}
					n= n+7;
					outs<< a;
				}
			}
		}
	}	
	outs<<(char) 8;
*/
}

void dynamic_mc_model::set_alignment_pattern(std::ifstream & in){
/*	size_t bit = 12;
	char c;
	char h;
	c = in.get();
	while(c != 8){
		size_t acc1;
		size_t acc2;
		size_t level;
		in >> acc1;
		in >> acc2;
		in >> level;
		set_al_level(acc1,acc2,level);// because in decoding we don't have access to level vector values
		make_all_alignment_pattern(level);
		for(std::set<std::string>::const_iterator it= all_alignment_patterns.begin(); it!= all_alignment_patterns.end();it++){
			std::string pattern = *it;
			highValue.at(acc1).at(acc2).insert(std::make_pair(pattern, std::vector<unsigned int>(NUM_KEEP+NUM_DELETE+10,0)));
		}
		for(std::map<std::string,std::vector<unsigned int> >::iterator it= highValue.at(acc1).at(acc2).begin(); it!= highValue.at(acc1).at(acc2).end();it++){
			std::vector<bool> binary_high_value(0);
			size_t bound = ((NUM_KEEP+NUM_DELETE+10)*bit)/8;
			for(size_t j = 0 ; j < bound ; j++){ 
				h=in.get();
				size_t H = size_t(h);
				for(size_t k = 0; k < 8 ; k++){
					binary_high_value.push_back(H%2);
					H = H/2;	
				}
			}
			size_t counter =0;
			for(size_t i = 0; i < (binary_high_value.size())-bit;i++){
				unsigned int high_value = 0;					
				for(size_t j =i; j < i+bit; j++){
					high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
				}
				i=i+bit-1;
				it -> second.at(counter)= high_value;
				counter = counter +1;
			}
			it -> second.at(NUM_KEEP+NUM_DELETE+9)=powersOfTwo.at(bit); 
		}
		c=in.get();
	}
*/
}
/*
void dynamic_mc_model::count_context_OneToTwo(const pw_alignment & p ,size_t level, abstract_context_functor & functor)const{
	std::string seq = "";
	size_t left_1; 
	size_t right_1;
	p.get_lr1(left_1,right_1);
	size_t acc1 = data.accNumber(p.getreference1());
	size_t acc2 = data.accNumber(p.getreference2());
	size_t first_patterns = level;
	size_t power = powersOfTwo.at(NUM_KEEP-1);//we break keeps longer than 'power'
	for(size_t j = 0; j < level; j++){//Ex. Two keeps of length 2^1 and 2^0 is created at the begining of each alignment if level == 2.
		first_patterns--;
		seq+=modification_character(-1,-1,-1,first_patterns);//returns "keep"
	}
	for (size_t i = 0; i< p.alignment_length(); i++){
		size_t n = 0;
		int modify_base =-1;
		int num_delete=-1;
		int insert_base=-1;
		int num_keep=-1;
		std::string seq1(" ",level+1);
		char seq2;
		for(size_t w = level; w>0 ;w--){
			seq1.at(level-w)=seq.at(seq.size()-w);
		}
		char s1chr;
		char s2chr;
		size_t s1;
		size_t s2;
		char s1nchr;
		char s2nchr;
		size_t s1n;
		size_t s2n;
		p.alignment_col(i, s1chr, s2chr);				
		s1 = dnastring::base_to_index(s1chr);
		s2 = dnastring::base_to_index(s2chr);
		for(size_t j =i; j < p.alignment_length(); j++){
			p.alignment_col(j, s1nchr, s2nchr);				
			s1n = dnastring::base_to_index(s1nchr);
			s2n = dnastring::base_to_index(s2nchr);
			if(s1n == 5){
				
			}else break;	
		}
		seq1.at(level)=s1n;
		if(s1 == s2){
			size_t klength = 0;
			for(size_t j = i; j<p.alignment_length(); j++){
				char q1chr;
				char q2chr;
				p.alignment_col(j, q1chr, q2chr);				
				size_t q1 = dnastring::base_to_index(q1chr);
				size_t q2 = dnastring::base_to_index(q2chr);
				if(q1 == q2){
					klength +=1;
				}else{	
					break;
				}
			}
			if(klength > power){
				num_keep = NUM_KEEP-1;
				n=powersOfTwo.at(num_keep)-1;
				seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
			}else{
				for (size_t m = NUM_KEEP; m > 0; m--){
					if((klength & powersOfTwo.at(m-1)) != 0){
						num_keep=m-1;
						n= powersOfTwo.at(num_keep)-1;
						seq += modification_character(modify_base,num_delete,insert_base,num_keep);
						break;
					}
				}
			}
			seq2 = seq.at(seq.size()-1);
			functor. see_context(acc1,acc2,p,i,seq1,seq2);
		}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s2;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
				}
				if(s1 == 5){
					insert_base = s2;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);						
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
				}
				if(s2 == 5){
					size_t dlength = 0;
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
						size_t q2 = dnastring::base_to_index(q2chr);
						if(q2 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					if(dlength > powersOfTwo.at(NUM_DELETE-1)){
						num_delete = NUM_DELETE-1;
						n= powersOfTwo.at(num_delete)-1;
						seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
					}else{
						for(size_t m = NUM_DELETE; m > 0; m--){
							if((dlength & powersOfTwo.at(m-1)) != 0){
								num_delete = m-1;
								n= powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);		
								break;
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
				}

			}
			i=i+n;
		}

}
void dynamic_mc_model::count_context_TwoToOne(const pw_alignment & p ,size_t level, abstract_context_functor & functor)const{
		std::string seq = "";
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		size_t first_patterns = level;
		for(size_t j = 0; j < level; j++){
			first_patterns --;
			seq+=modification_character(-1,-1,-1,first_patterns);
		}
		for (size_t i = 0; i< p.alignment_length(); i++){
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			size_t n = 0;			
			std::string seq1(" ",level+1);
			char seq2;
			for(size_t w = level; w>0 ;w--){
				seq1.at(level-w)=seq.at(seq.size()-w);
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			char s1nchr;
			char s2nchr;
			size_t s2n;
			if(s2 == 5){
				for(size_t j =i; j < p.alignment_length(); j++){
					p.alignment_col(j, s1nchr, s2nchr);				
					s2n = dnastring::base_to_index(s2nchr);
					if(s2n != 5){
						seq1.at(level)=s2n;
						break;
					}else continue;
				}
			}else seq1.at(level)=s2;
			if(s1 == s2){
				size_t klength = 0;
				for(size_t j = i; j<p.alignment_length(); j++){
					char q1chr;
					char q2chr;
					p.alignment_col(j, q1chr, q2chr);				
					size_t q1 = dnastring::base_to_index(q1chr);
					size_t q2 = dnastring::base_to_index(q2chr);
					if(q1 == q2){
						klength +=1;
					}else {	
						break;
					}
				}
				if(klength > powersOfTwo.at(NUM_KEEP-1)){
					num_keep = NUM_KEEP-1;
					n=powersOfTwo.at(num_keep)-1;
					seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
				}else{
					for(size_t m = NUM_KEEP; m > 0; m--){
						if((klength & powersOfTwo.at(m-1)) != 0){
							num_keep=m-1;
							n=powersOfTwo.at(num_keep)-1;						
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);
							break;
						}
					}
				}
				seq2 = seq.at(seq.size()-1);
				functor. see_context(acc2,acc1,p,i,seq1,seq2);
			}else{

				if((s1!=5) & (s2!=5)){
					modify_base = s1;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);
				}
				if(s2 == 5){
					insert_base = s1;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);					
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);						
				}
				if(s1 == 5){
					size_t dlength = 0;
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
						size_t q1 = dnastring::base_to_index(q1chr);
						if(q1 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					if(dlength > powersOfTwo.at(NUM_DELETE-1)){
						num_delete = NUM_DELETE-1;
						n = powersOfTwo.at(num_delete)-1;
						seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
					}else{
						for(size_t m = NUM_DELETE ; m > 0; m--){
							if((dlength & powersOfTwo.at(m-1)) != 0){
								num_delete = m-1;
								n = powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);
								break;			
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);
				}
			}
			i=i+n;
		}		
}
*/

/*
void dynamic_mc_model::modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep)const {
	modify_base = -1;
	num_delete =-1;
	insert_base = -1;
	num_keep = -1;
	if(enc < 5) {
		modify_base = enc;
		return;	
	} 
	if(enc < 5 + NUM_DELETE) {
		num_delete = powersOfTwo.at(enc - 5);
		return;
	}
	if (enc < 5 + NUM_DELETE + NUM_KEEP){
		num_keep = powersOfTwo.at(enc-NUM_DELETE-5);
		return;
	}
	if(enc< 5+ 5 + NUM_KEEP + NUM_DELETE){
		insert_base = enc- 5 - NUM_KEEP - NUM_DELETE;
		return;
	}
	assert(0);
}
char dynamic_mc_model::modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const {
	//return the enc
	if(num_delete != -1){
		return 5 + num_delete;
	}
	if(num_keep != -1){
		return 5 + NUM_DELETE + num_keep;
	}
	if(modify_base != -1) {
		return modify_base;
	}
	if(insert_base != -1){
		return insert_base + NUM_KEEP + NUM_DELETE + 5;
	}
	assert (false);
	return -1;
}
*/


void dynamic_mc_model::train(std::ofstream & outs){//call it before training mc model to find the best level.
	train_sequence_model();	
	train_alignment_model();


// after training we write training results to a file, read it again and verify that this results in identical class members
#ifndef NDEBUG
	std::string tfile("graph_parameters_test_file");
	::pid_t pid = ::getpid();
	size_t spid = pid;
	char buf[200];
	sprintf(buf, "%u.dat", spid);
	tfile+=(buf);
	std::ofstream outf(tfile.c_str(), std::ios::binary | std::ios::out);
	if(outf) {
		write_parameters(outf);
		outf.close();
		check_parameters_in_file(tfile);
	} else {
		std::cerr << "Warning: could not write to test file: "<< tfile << std::endl;
	}

#endif


exit(0);
	recursive_al_model();
	write_parameters(outs);
	write_alignment_parameters(outs);
	std::cout<< "dynamic train is done! "<<std::endl;
}

void dynamic_mc_model::cost_function(const pw_alignment& p , double & c1 , double & c2 , double & m1 , double & m2)const{//TODO add the possibilty of level = 0
/*	if(p.is_cost_cached()) {
		c1 = p.get_create1();
		c2 = p.get_create2();
		m1 = p.get_modify1();
		m2 = p.get_modify2();
		return; 
	}
	cost_functor f(data,mod_cost);
	std::vector<double> cost_on_sample(2,0);
	std::vector<double> modify_cost(2,0);
	size_t acc1 = data.accNumber(p.getreference1());
	size_t acc2 = data.accNumber(p.getreference2());
	size_t level = get_al_level(acc1,acc2);
	count_context_OneToTwo(p,level,f);
	count_context_TwoToOne(p,level,f);
	char s1chr;
	char s2chr;
	size_t left1;
	size_t right1;	
	size_t left2;
	size_t right2;
	p.get_lr1(left1,right1);
	p.get_lr2(left2,right2);
	for(size_t i = left1; i< right1; i++){
		s1chr = data.getSequence(p.getreference1()).at(i);
		size_t s1 = dnastring::base_to_index(s1chr);
		std::stringstream context1;
		for (size_t j = level; j>0; j--){
			
			if(i<j){
				char r1chr = 'A';
				context1 << r1chr;					
			}else{
				char r1chr = data.getSequence(p.getreference1()).at(i-j);
				context1 << r1chr;
			}
		}
		std::string seq1;
		context1>>seq1;
		std::map <std::string, std::vector<double> >::const_iterator it= sequence_successive_bases.at(acc1).find(seq1);
		assert(it != sequence_successive_bases.at(acc1).end());
		cost_on_sample.at(0) += it->second.at(s1);
	}	
	for(size_t i = left2; i<right2; i++){
		s2chr = data.getSequence(p.getreference2()).at(i);
		size_t s2 = dnastring::base_to_index(s2chr);
		std::stringstream context2;
		for(size_t j = level; j>0; j--){
			if(i<j){
				char r2chr = 'A';
				context2 << r2chr;	
			}else{
				char r2chr = data.getSequence(p.getreference2()).at(i-j);
				context2 << r2chr;
			}
		}
		std::string seq2;
		context2>>seq2;
		std::map <std::string, std::vector<double> >::const_iterator it1= sequence_successive_bases.at(acc2).find(seq2);
		assert(it1 != sequence_successive_bases.at(acc2).end());
		cost_on_sample.at(1) += it1->second.at(s2);
	}
	c1 = cost_on_sample.at(0);
	c2 = cost_on_sample.at(1);
	m1 = f.get_modify(p,acc1,acc2);
	m2 = f.get_modify(p,acc2,acc1);
	modify_cost.at(0) = m1;
	modify_cost.at(1) = m2;
	p.set_cost(cost_on_sample, modify_cost);
	assert(p.is_cost_cached());
	assert(m1==p.get_modify1());
	assert(m2==p.get_modify2());
	assert(c1==p.get_create1());
	assert(c2==p.get_create2());
*/
}

void dynamic_mc_model::gain_function(const pw_alignment& p , double & g1 , double & g2)const{
	double c1;
	double c2;
	double m1;
	double m2;
	cost_function(p, c1, c2, m1,m2);

	
	g1 = c2 - m1;
	g2 = c1 - m2;
}

// TODO we should not need that
size_t dynamic_mc_model::get_acc_level(size_t & acc)const{
//	return accLevel.at(acc);
	return -1;
}
void dynamic_mc_model::set_acc_level(size_t & level){
//	accLevel.push_back(level);

}
size_t dynamic_mc_model::get_al_level(size_t& acc1 , size_t & acc2)const{
//	return al_level.at(acc1).at(acc2);
	return -1;

}
void dynamic_mc_model::set_al_level(size_t& acc1 , size_t & acc2, size_t & level){
//	al_level.at(acc1).at(acc2) = level;
}

const std::map<std::string, std::vector<unsigned int> >&  dynamic_mc_model::get_high(size_t acc)const{
//	return high.at(acc);
}
std::string dynamic_mc_model::get_firstPattern(size_t& accession)const{
/*	std::string first_pattern;
	size_t level = accLevel.at(accession);
	for(size_t i = 0; i< level; i++){
		first_pattern += "A";
	}
	return first_pattern;
*/
}
std::string dynamic_mc_model::get_firstAlignmentPattern(size_t& acc1, size_t& acc2)const{
/*		std::string first_pattern;
		size_t level = al_level.at(acc1).at(acc2);
		size_t firstPatterns = level;
		for(size_t j = 0; j < level; j++){
			firstPatterns --;
			first_pattern += modification_character(-1,-1,-1,firstPatterns);
		}
		return first_pattern;

*/
	}	
void dynamic_mc_model::get_encoded_member(pw_alignment & al, size_t center_ref, size_t center_left, encoding_functor & functor,std::ofstream& outs)const{
/*	size_t acc1 = data.accNumber(al.getreference1());
	size_t acc2 = data.accNumber(al.getreference2());
	size_t accession = data.accNumber(center_ref);
	size_t left_1; 
	size_t left_2;
	size_t right_1;
	size_t right_2;
	size_t level = al_level.at(acc1).at(acc2);
	al.get_lr1(left_1,right_1);
	al.get_lr2(left_2,right_2);
	if(al.getreference1()==center_ref && left_1 == center_left){
		count_context_OneToTwo(al, level, functor);	
	}else{
		count_context_TwoToOne(al, level, functor);
	}
*/
}

std::vector<unsigned int> dynamic_mc_model::get_high_at_position(size_t seq_index, size_t position)const{
/*	const dnastring & sequence = data.getSequence(seq_index);
	size_t accession = data.accNumber(seq_index);
	size_t level = accLevel.at(accession);
	std::stringstream context;
	for(size_t j = level; j>0; j--){
		if(position < j){
			char chr = 'A';
			context<<chr;
		}else{
			char chr = sequence.at(position-j);
			context<<chr;
		}
	}
	std::string current_pattern;
	context >> current_pattern;
	std::map<std::string, std::vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
	assert(it!=high.at(accession).end());
	return it->second;
*/
}
const std::map<std::string, std::vector<unsigned int> > & dynamic_mc_model::get_highValue(size_t acc1, size_t acc2)const{
//	return highValue.at(acc1).at(acc2);
}


