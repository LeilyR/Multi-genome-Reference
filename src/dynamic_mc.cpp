#include "dynamic_mc.hpp"


template<typename reader>
sequence_contexts<reader>::sequence_contexts(reader & r, size_t level): rd(r), level(level) {
	astring.append(level, 'A');

}

template<typename reader>
sequence_contexts<reader>::~sequence_contexts() {}

template<typename reader>
void sequence_contexts<reader>::read_sequence(const std::string & sequence) {
	for(size_t i=0; i<level && i<sequence.length(); ++i) {
		std::string cont = astring.substr(i);
		cont.append(sequence.substr(0, i));
		rd.see_context(cont, dnastring::base_to_index(sequence.at(i)));
		assert(cont.length() == level);
	}
	for(size_t i=level; i<sequence.length(); ++i) {
		std::string cont = sequence.substr(i - level, level);
		rd.see_context(cont, dnastring::base_to_index(sequence.at(i)));
	}
}


template<typename reader>
void sequence_contexts<reader>::read_sequence(const dnastring & sequence, const size_t & from, const size_t & to) {
	assert(from != to); // should only be called for sequences longer than 1 so that we know which strand we are on
	if(from < to) {
		size_t len = to - from + 1;
		std::string seq(len, ';');
		for(size_t i=0; i<len; ++i) {
			seq.at(i) = sequence.at(from + i);
		}
		read_sequence(seq);
	} else {
		size_t len = from - to + 1;
		std::string seq(len, ';');
		for(size_t i=0; i<len; ++i) {
			seq.at(i) = dnastring::complement(sequence.at(from - i));
		}
		read_sequence(seq);
	}
}

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
	//	on_char = (size_t) dnastring::base_to_index('A');
		for(size_t next_r1base = at+1; next_r1base < al.alignment_length(); ++next_r1base) {
			char n_ch1, n_ch2;
			al.alignment_col(next_r1base, n_ch1, n_ch2);
			if(n_ch1!='-') {
				on_char = (size_t) dnastring::base_to_index(n_ch1);
				break;
			}
		}
	/*	for(size_t last_r1base = at-1; last_r1base > 0; last_r1base--){
			char n_ch1, n_ch2;
			al.alignment_col(last_r1base, n_ch1, n_ch2);
			if(n_ch2!='-') {
				on_char = (char) dnastring::base_to_index(n_ch2);
				break;
			}
		}*/
		modification = dynamic_mc_model::modification_character(-1, -1, dnastring::base_to_index(ch2),-1);
	//	std::cout << "on base " << size_t(on_char) <<  " modification is " << size_t(modification) << std::endl;
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
//		std::cout << "log del " << log_del << std::endl;
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
	//	on_char = (char) dnastring::base_to_index('A');
		for(size_t next_r2base = at+1; next_r2base < al.alignment_length(); ++next_r2base) {
			char n_ch1, n_ch2;
			al.alignment_col(next_r2base, n_ch1, n_ch2);
			if(n_ch2!='-') {
				on_char = (char) dnastring::base_to_index(n_ch2);
				break;
			}
		}
	/*	for(size_t last_r2base = at-1; last_r2base > 0; last_r2base--){
			char n_ch1, n_ch2;
			al.alignment_col(last_r2base, n_ch1, n_ch2);
			if(n_ch2!='-') {
				on_char = (char) dnastring::base_to_index(n_ch2);
				break;
			}
		}*/
		modification = dynamic_mc_model::modification_character(-1, -1, dnastring::base_to_index(ch1),-1);
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
//	pw.print();
	std::string context;
	for(size_t i=level; i>0; --i) {
		size_t numk = i-1;
		if(numk > NUM_KEEP_DYN-1) numk = NUM_KEEP_DYN;
		char keepc = dynamic_mc_model::modification_character(-1, -1, -1, numk);
		context.append(1, keepc);
	}
//	std::cout << "first al context is "<<int(context.at(0)) << " " << int(context.at(1)) << std::endl;
	size_t at_al_col = 0;
	while(at_al_col < pw.alignment_length()) {

		char modc;
		char onc;
		size_t len;
		next_modification_1_2(pw, at_al_col, modc, onc, len);
		
		// add current known base to context
		std::string newcontext = context;
		newcontext.append(1, onc);


	//	std::cout << " found con at " << at_al_col << " is " ;
	//	for(size_t i=0; i<newcontext.length(); ++i) {
	//		std::cout << " " << (size_t)newcontext.at(i);
	//	}
	//	std::cout << std::endl;

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

void a_reader_counter::see_context(const std::string & context, size_t modification) {
	
	std::map<std::string, std::vector<size_t> >::iterator findc = countsmap.find(context);
	if(findc == countsmap.end()) {
		countsmap.insert(std::make_pair(context, vector<size_t>(NUM_MODIFICATIONS, 0)));
		findc = countsmap.find(context);
	}
	findc->second.at(modification)++;
}
s_reader_counter::s_reader_counter():numchars(5) {}

s_reader_counter::~s_reader_counter() {}

void s_reader_counter::see_context(const std::string & context, size_t modification) {
	
	std::map<std::string, std::vector<size_t> >::iterator findc = countsmap.find(context);
	if(findc == countsmap.end()) {
		countsmap.insert(std::make_pair(context, vector<size_t>(numchars, 0)));
		findc = countsmap.find(context);
	}
	findc->second.at(modification)++;
}



a_reader_costs::a_reader_costs(const std::map<std::string, std::vector<double> > & mcost):cost_sum(0), model_costs(mcost) {

}

a_reader_costs::~a_reader_costs() {}

void a_reader_costs::see_context(const std::string & context, size_t modification) {
//	for(std::map<std::string, std::vector<double> >::const_iterator it = model_costs.begin(); it != model_costs.end(); it++){
//		std::cout << it->first.size() << std::endl;
//	}
	std::map<std::string, std::vector<double> >::const_iterator it = model_costs.find(context);
//	std::cout << "context size " <<context.size() <<std::endl;
	assert(it!=model_costs.end());
	const std::vector<double> & mc = it->second;
	cost_sum+=mc.at(modification);
}



s_reader_costs::s_reader_costs(const std::map<std::string, std::vector<double> > & mcost):cost_sum(0), model_costs(mcost) {

}

s_reader_costs::~s_reader_costs() {}


void s_reader_costs::see_context(const std::string & context, size_t modification) {
	std::map<std::string, std::vector<double> >::const_iterator it = model_costs.find(context);
	assert(it!=model_costs.end());
	const std::vector<double> & mc = it->second;
	cost_sum+=mc.at(modification);
}


a_reader_encode::a_reader_encode(const std::map<std::string, std::vector<uint32_t> > & mhigh, wrapper & wrap):model_high(mhigh),wrappers(wrap){}

a_reader_encode::~a_reader_encode(){}

void a_reader_encode::see_context(const std::string & context,size_t modification){
//	std::cout << "context size "<< context.size()<<std::endl;
//	std::cout << "context: "; 
	for(size_t i = 0; i < context.size();i++){
//		std::cout << int(context.at(i)) << " ";
		int con = int(context.at(i));
		wrappers.context(con);
	}
//	std::cout << " " <<std::endl;
//	std::cout << "modification "<< modification << std::endl;
	std::map<std::string, std::vector<uint32_t> >::const_iterator it = model_high.find(context);
	assert(it!=model_high.end());
	std::vector<unsigned int> low_high;
	const std::vector<uint32_t> & high = it->second;
	unsigned int l = 0;
	unsigned int h = high.at(modification);
	if(modification != 0){
		l = high.at(modification-1);
	}
//	std::cout<< "l_model "<< l << "h_model: "<< h <<std::endl;
	low_high.push_back(l);
	low_high.push_back(h);
	low_high_values.push_back(low_high);
}
dynamic_mc_model::dynamic_mc_model(all_data & d, wrapper & wr, size_t nt): data(d),wrappers(wr), num_threads(nt) {
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
		return 5 + NUM_DELETE_DYN + num_keep;
	}
	if(modify_base != -1) {
	//	std::cout<< "there is a modification to " << dnastring::index_to_base(modify_base) << std::endl;
		return modify_base;
	}
	if(insert_base != -1){
	//	std::cout<< " there is a insertion of " << dnastring::index_to_base(insert_base) <<std::endl;
		return insert_base + NUM_KEEP_DYN + NUM_DELETE_DYN + 5;
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

size_t dynamic_mc_model::modification_length(char mod){
	int modify_base = -1;
	int num_delete =-1;
	int insert_base = -1;
	int num_keep = -1;
	size_t length = 0;
	modification(mod, modify_base, num_delete, insert_base, num_keep);
	if(modify_base != -1){
		length = 1;
	}
	if(insert_base != -1){//on the reference that '-' happens.
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
	if(enc < 5 + NUM_DELETE_DYN) {
		num_delete = (1<<(enc - 5));
		return;
	}
	if (enc < 5 + NUM_DELETE_DYN + NUM_KEEP_DYN){
		num_keep = (1<<(enc-NUM_DELETE_DYN-5));
		return;
	}
	if(enc< 5+ 5 + NUM_KEEP_DYN + NUM_DELETE_DYN){
		insert_base = enc- 5 - NUM_KEEP_DYN - NUM_DELETE_DYN;
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
	sequence_model_highs.resize(data.numAcc());
	sequence_model_costs.resize(data.numAcc());
	vector<unsigned char> alphabet;
	get_sequence_alphabet(alphabet);

	std::cout << "number of acc " << data.numAcc() << std::endl;
	for(size_t acc = 0; acc<data.numAcc(); ++acc) {
// storage widhts (few bits) to encoding widths (32 bits)
		model_bits_to_full_high_values(sequence_models.at(acc), acc_sequence_all_bases_total.at(acc), MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, alphabet, sequence_model_highs.at(acc), sequence_model_costs.at(acc));
	//	std::cout << sequence_model_highs.size() << std::endl;
	}

}


void dynamic_mc_model::compute_alignment_model() {
	vector<unsigned char> alignment_alphabet;
	get_alignment_alphabet(alignment_alphabet);
	alignment_model_highs = std::vector<std::vector<std::map<std::string, std::vector<uint32_t> > > > (data.numAcc(), std::vector<std::map<std::string, std::vector<uint32_t> > >(data.numAcc()));
	alignment_model_costs = std::vector<std::vector<std::map<std::string, std::vector<double> > > > (data.numAcc(), std::vector<std::map<std::string, std::vector<double> > >(data.numAcc()));
	for(size_t a1 = 0; a1<data.numAcc(); ++a1) {
		for(size_t a2 = 0; a2<data.numAcc(); ++a2) {
			// TODO can we make this faster?
			std::cout << "al end in model " << a1 << " " << a2 << " " <<  alignments_all_modification_total.at(a1).at(a2) <<std::endl;			
			model_bits_to_full_high_values(alignment_models.at(a1).at(a2), alignments_all_modification_total.at(a1).at(a2), MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+1, alignment_alphabet, alignment_model_highs.at(a1).at(a2), alignment_model_costs.at(a1).at(a2));
		}
	}

// compute base cost
/*
	Using an alignment in encoding has these costs:
	
	* alignment begin flag on 2 sequences
	* 2 * cluster center adress 
	* the modifications are contained in using modify - create cost (0 base cost)
	* alignment end flag on 2 sequences

*/

	alignment_adress_cost = log2(data.numAlignments());
	alignment_begin_cost.resize(data.numAcc());
	alignment_end_cost.resize(data.numAcc());
	alignment_base_cost.resize(data.numAcc());
	std::cout << "alignment base cost estimate: " << alignment_adress_cost << std::endl;
	for(size_t a = 0; a<data.numAcc(); ++a) {
		alignment_end_cost.at(a).resize(data.numAcc());
		alignment_base_cost.at(a).resize(data.numAcc());
		size_t bases = acc_sequence_all_bases_total.at(a);
		size_t alflag = acc_sequence_alignment_flag.at(a);
		size_t al_begin_width = alflag - bases;
		double al_begin_cost = -log2((double) al_begin_width / (double) TOTAL_FOR_ENCODING);
		alignment_begin_cost.at(a) = al_begin_cost;
		std::cout << "acc "<< a << " alignment begin cost estimate " << alignment_begin_cost.at(a) << std::endl;
		for(size_t a2 =0; a2<data.numOfAcc(); ++a2) {
			size_t al_mods = alignments_all_modification_total.at(a).at(a2);
			size_t al_end_width = TOTAL_FOR_ENCODING - al_mods;
			double al_end_cost = -log2((double) al_end_width / (double) TOTAL_FOR_ENCODING);
			alignment_end_cost.at(a).at(a2) = al_end_cost;
			alignment_base_cost.at(a).at(a2) = alignment_adress_cost + alignment_begin_cost.at(a) + alignment_end_cost.at(a).at(a2);
			std::cout << " acc " << a << " to acc " << a2 << " alignment end cost estimate " << alignment_end_cost.at(a).at(a2) <<  " total base cost " << alignment_base_cost.at(a).at(a2) << std::endl;
		}
	}

}
void dynamic_mc_model::compute_sequence_model_decoding(){
	sequence_model_highs.resize(data.numOfAcc());
	sequence_model_costs.resize(data.numOfAcc());
	vector<unsigned char> alphabet;
	get_sequence_alphabet(alphabet);

	std::cout << "number of acc " << data.numOfAcc() << std::endl;
	for(size_t acc = 0; acc<data.numOfAcc(); ++acc) {
// storage widhts (few bits) to encoding widths (32 bits)
		model_bits_to_full_high_values(sequence_models.at(acc), acc_sequence_all_bases_total.at(acc), MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, alphabet, sequence_model_highs.at(acc), sequence_model_costs.at(acc));
		std::cout << sequence_model_highs.size() << std::endl;
	}
}


void dynamic_mc_model::compute_alignment_model_decoding(){
	vector<unsigned char> alignment_alphabet;
	get_alignment_alphabet(alignment_alphabet);
	alignment_model_highs = std::vector<std::vector<std::map<std::string, std::vector<uint32_t> > > > (data.numOfAcc(), std::vector<std::map<std::string, std::vector<uint32_t> > >(data.numOfAcc()));
	alignment_model_costs = std::vector<std::vector<std::map<std::string, std::vector<double> > > > (data.numOfAcc(), std::vector<std::map<std::string, std::vector<double> > >(data.numOfAcc()));
	for(size_t a1 = 0; a1<data.numOfAcc(); ++a1) {
		for(size_t a2 = 0; a2<data.numOfAcc(); ++a2) {
			std::cout << "al end in model " << a1 << " " << a2 << " " <<  alignments_all_modification_total.at(a1).at(a2) <<std::endl;			
			model_bits_to_full_high_values(alignment_models.at(a1).at(a2), alignments_all_modification_total.at(a1).at(a2), MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+1, alignment_alphabet, alignment_model_highs.at(a1).at(a2), alignment_model_costs.at(a1).at(a2));
			for(std::map<std::string, std::vector<uint32_t> >::iterator it = alignment_model_highs.at(a1).at(a2).begin();it != alignment_model_highs.at(a1).at(a2).end();it++){
				assert(it->second.at(it->second.size()-1) <=  alignments_all_modification_total.at(a1).at(a2));
			}

		}
	}



}
double dynamic_mc_model::get_alignment_base_cost(const size_t & acc1, const size_t & acc2) const {
	return alignment_base_cost.at(acc1).at(acc2);
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
//		std::map< std::string, std::vector<size_t> > context_counts;
		s_reader_counter ccounts;
		scontext_counter ccount(ccounts, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
		size_t acc_num_bases = 0;
		for(size_t k = 0; k <data.getAcc(acc).size(); k++) {
			std::string sequence = data.getSequence(data.getAcc(acc).at(k)).str();
			acc_num_bases += sequence.length();
			ccount.read_sequence(sequence);



		} // for sequences in accession

		// we estimate the probability of having an alignment on this accession by the total number of half alignments on that accession
		// in the end we could have less (we select only some alignments) or more (we split alignments to remove overlap)
		size_t al_on_this_acc = alignments_on_acc.at(acc);
		size_t acc_sequence_ends = data.getAcc(acc).size()+1;// +1 is added to be used as the end of accession flag
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
		std::cout <<"acc " << acc << ": encoding widths: any base: "<< flag_widths.at(0) << " half alignment begin: " << flag_widths.at(1) << " sequence end: " << flag_widths.at(2) << std::endl;
		std::cout <<"acc " << acc << ": encoding costs: any base: " << -log2((double)flag_widths.at(0)/(double)TOTAL_FOR_ENCODING) << 
			 "bit half alignment begin: " << -log2((double)flag_widths.at(1)/(double)TOTAL_FOR_ENCODING) << "bit sequence end: " << -log2((double)flag_widths.at(2)/(double)TOTAL_FOR_ENCODING) <<"bit "<< std::endl;
		
		uint32_t new_total = flag_widths.at(0);//This is the total value for the sequence patterns
		std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > acc_seq_model;

		std::cout << " run context selector on " << ccounts.countsmap.size() << " contexts"<< std::endl;
		double enc_cost = context_selector(ccounts.countsmap, alphabet, new_total, acc_seq_model, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL );
		enc_cost += acc_sequence_ends * -log2((double)flag_widths.at(2)/(double)TOTAL_FOR_ENCODING);
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
	distri.resize(counts.size());
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


void dynamic_mc_model::model_to_widths(const std::vector<uint32_t> bits, unsigned char model_index, uint32_t target_total, std::vector<uint32_t> & widths) {
	widths.resize(bits.size());
	size_t num_bits;
	bool use_sqroot;
	get_model_type(model_index, use_sqroot, num_bits);
	if(use_sqroot) {
		sqroot_bits_to_width_encoding(bits, target_total, widths);
	} else {
		bits_to_width_encoding(bits, target_total, widths);
	}


}

double dynamic_mc_model::model_to_costs(const std::vector<uint32_t> bits, unsigned char model_index, uint32_t target_total, std::vector<double> & costs) {
	costs.resize(bits.size());
	std::vector<uint32_t> widths(bits.size());
	model_to_widths(bits, model_index, target_total, widths);
	
	size_t total = 0;
	for(size_t i=0; i<widths.size(); ++i) {
		if(widths.at(i) < 1) widths.at(i) = 1; 
		total+=widths.at(i);
	}
	for(size_t i=0; i<widths.size(); ++i) {
		costs.at(i) = -log2((double) widths.at(i) / (double) total);
	}
}



double dynamic_mc_model::distri_to_bits(const std::vector<double> & distri, unsigned char model_index, std::vector<uint32_t> & bits) {
	bits.resize(distri.size());
	size_t num_bits;
	bool use_sqroot;
	get_model_type(model_index, use_sqroot, num_bits);
	if(use_sqroot) {
		distribution_to_sqroot_bits(distri, num_bits, bits);
	} else {
		distribution_to_bits(distri, num_bits, bits);
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
		distri_to_bits(distri, model_index, bits);
		model_to_widths(bits, model_index, target_total, widths);
		size_t num_bits;
		bool use_sqroot;
		get_model_type(model_index, use_sqroot, num_bits);

			
		// encoding cost
		double cost =  count_cost_eval(counts, widths);

//		parameter (model size) cost
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




void dynamic_mc_model::width_correct(std::vector<uint32_t> & widths, uint32_t target_total) {
	// we correct the widths several times until the sum matches target_total
	size_t it;
	for(it=0; it<widths.size(); ++it) {
		size_t total = 1;  // avoid total 0
		for(size_t i=0; i<widths.size(); ++i) {
			total += widths.at(i);
	//		std::cout << " it " << it << " w " << widths.at(i) << std::endl;
		} 
		double scale = (double)  target_total / (double) total;
	//	std::cout << " scale " << scale << std::endl;
		size_t scaled_total = 0;
		for(size_t i=0; i<widths.size(); ++i) {
			size_t nw =  (size_t)(scale * ((double)(widths.at(i) + 0.5)));
			if(nw < ENCODING_MIN_WIDTH) nw = ENCODING_MIN_WIDTH;
			if(nw > target_total) nw = target_total;
			widths.at(i) = (uint32_t) nw;
			scaled_total+=widths.at(i);
	//		std::cout << " scaled " << widths.at(i) << std::endl;
		}
	//	std::cout << " scaled total " << scaled_total << std::endl;

		int64_t correction = (int64_t)target_total - (int64_t)scaled_total;
		if(correction == 0) break;
		int64_t correction_done = 0;
		double correction_factor = (double) correction / (double) scaled_total;
	//	std::cout << " correction " << correction << " factor " << correction_factor << std::endl;
		int extra_correction = 1; // correct some more to be sure to get it done
		if(correction < 0) extra_correction = -1;
		for(size_t i=0; i<widths.size(); ++i) {
			int64_t this_correct = (int64_t) ((double)widths.at(i) * correction_factor) + extra_correction;
			int64_t prediction = (int64_t)(widths.at(i)) + this_correct;
	//		std::cout << " first correct " << this_correct <<  " on " << widths.at(i) << " predicted to " << prediction << std::endl;
			if( ( prediction < ENCODING_MIN_WIDTH ) ) {
	//			std::cout << " second correct " << this_correct << " numbers " << widths.at(i) << " " << prediction << " " <<  ENCODING_MIN_WIDTH << std::endl;
				this_correct = ENCODING_MIN_WIDTH - (int64_t)(widths.at(i));
			}
			if (prediction > target_total) {
				this_correct = target_total - (int64_t)widths.at(i);
			}

			if(correction < 0) {
				if(correction_done + this_correct <= correction) {
					this_correct = correction - correction_done;
	//		std::cout << " third correct " << this_correct << std::endl;
				}
			} else {
				if(correction_done + this_correct >= correction) {
					this_correct = correction - correction_done;
				}
			}
			correction_done += this_correct;
			widths.at(i) += this_correct;
	//		std::cout << " real correct " << this_correct << " correction done " << correction_done << std::endl;
			if(correction_done == correction) break;
		}
	} // redo algorithm iterations
#ifndef NDEBUG
	size_t test_sum = 0;
	for(size_t i=0; i<widths.size(); ++i) {
		test_sum+=widths.at(i);
	}
	assert(test_sum == target_total);
#endif


}

void dynamic_mc_model::sqroot_bits_to_width_encoding(const std::vector<uint32_t> & sqbits, uint32_t target_total, std::vector<uint32_t> & widths) {
	widths.resize(sqbits.size());
	for(size_t i=0; i<sqbits.size(); ++i) {
		widths.at(i) = sqbits.at(i) * sqbits.at(i);
	}
	width_correct(widths, target_total);

}


void dynamic_mc_model::bits_to_width_encoding(const std::vector<uint32_t> & bits, uint32_t target_total, std::vector<uint32_t> & widths) {
	widths.resize(bits.size());
	for(size_t i=0; i<bits.size(); ++i) {
		widths.at(i) = bits.at(i);
	}
	width_correct(widths, target_total);
}







/* 
 TODO we can avoid copies made by this function

one data structure like the bitmodel map to keep all different length contexts and their high values/costs
then the current one for fixed length contexts (more of them) can keep references instead of copies


*/

void dynamic_mc_model::model_bits_to_full_high_values(const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & bitmodel, uint32_t target_total, size_t target_pattern_length, const std::vector<unsigned char> & alphabet,
	std::map< std::string, std::vector<uint32_t> > & hv_results, std::map<std::string, std::vector<double> > & cost_results) {
	std::string seq = "";
	std::set<std::string> pattern_set;

	std::vector< std::map< std::string, std::pair<std::vector<uint32_t>, std::vector<double> > > > lsorted_res(target_pattern_length+1);
	
	clock_t ctime = clock();	
	for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = bitmodel.begin(); it!=bitmodel.end(); ++it) {
		const std::string & cont = it->first;
		unsigned char model_index = it->second.first;
		const std::vector<uint32_t> & bits = it->second.second;
		size_t len = cont.length();

		std::vector<uint32_t> widths;
		std::vector<uint32_t> highs;
		std::vector<double> costs;
		model_to_widths(bits, model_index, target_total, widths);
		model_widths_to_highs(widths, target_total, highs);
		model_to_costs(bits, model_index, target_total, costs);

	//	std::cout << " high and costs for " << cont << std::endl;


		lsorted_res.at(len).insert(std::make_pair(cont, std::make_pair(highs, costs)));
	}
	ctime = clock() - ctime;
	clock_t itime = clock();
	for(size_t l=0; l < target_pattern_length; l++) {
		// for each context of length l, copy data to all derived contexts of length l+1 (any character added to the front)
		for(std::map< std::string, std::pair<std::vector<uint32_t>, std::vector<double> > >::iterator it = lsorted_res.at(l).begin(); it!=lsorted_res.at(l).end(); ++it) {
			const std::string & cont = it->first;
			const std::pair<std::vector<uint32_t>, std::vector<double> > & data = it->second;

			for(size_t i=0; i<alphabet.size(); ++i) {
				std::string newcont;
				newcont.append(1, alphabet.at(i));
				newcont.append(cont);

				std::map< std::string, std::pair<std::vector<uint32_t>, std::vector<double> > >::iterator findnc = lsorted_res.at(newcont.length()).find(newcont);
				assert(findnc == lsorted_res.at(newcont.length()).end());

				lsorted_res.at(newcont.length()).insert(std::make_pair(newcont, data)); 
			}
		}
	//	std::cout << " l " << l+1 <<  " s " << lsorted_res.at(l+1).size() << std::endl;
	}

	for(std::map< std::string, std::pair<std::vector<uint32_t>, std::vector<double> > >::const_iterator it = lsorted_res.at(target_pattern_length).begin(); it!=lsorted_res.at(target_pattern_length).end(); ++it) {
		const std::string & cont = it->first;
		const std::vector<uint32_t> & widths = it->second.first;
		const std::vector<double> & costs = it->second.second;

		hv_results.insert(std::make_pair(cont, widths));
		cost_results.insert(std::make_pair(cont, costs));
	}
	itime = clock() - itime;

	std::cout << " compute model times " << (double)ctime/CLOCKS_PER_SEC << " " << (double)itime/CLOCKS_PER_SEC << std::endl;
	std::cout << " sizes " << bitmodel.size() << " " << hv_results.size() << " " << cost_results.size() << std::endl;
} 



void dynamic_mc_model::model_widths_to_highs(const std::vector<uint32_t> & widths, uint32_t target_total, std::vector<uint32_t> & high_values) {
	uint32_t sum = 0;
	high_values.resize(widths.size());
	for(size_t i=0; i<widths.size(); ++i) {
		sum+=widths.at(i);
		high_values.at(i) = sum;
	}
	assert(sum == target_total);
}


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



#define NOMODEL_COMPRESS 1

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
					
#if NOMODEL_COMPRESS
					if(nomodel_counter < MAX_NOMODEL_NUM - 1) {
						nomodel_counter++;
					} else {
						outs.put((unsigned char) 255);
						nomodel_counter = 0;
					}
#else 
					// write no model symbol
					outs.put((unsigned char) NUM_ENCODING_TYPES);
#endif
				} else {

#if NOMODEL_COMPRESS
					if(nomodel_counter > 0) {
						size_t nmchar = nomodel_counter -1 + NUM_ENCODING_TYPES;
						outs.put((unsigned char) nmchar);
						nomodel_counter = 0;
					}
#endif


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
#if NOMODEL_COMPRESS
				nomodel_counter--;
#endif
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
#if NOMODEL_COMPRESS
				nomodel_counter = model_index - NUM_ENCODING_TYPES + 1;
				nomodel_counter--;


#endif


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
			mapmodel_write(outs, alignment_models.at(acc).at(a2), aalphabet, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+1);


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

/*

			for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = this->alignment_models.at(acc).at(a2).begin();
			it!=this->alignment_models.at(acc).at(a2).end(); ++it) {
				std::string cont = it->first;
				std::cout << "t ";
				for(size_t j=0; j<cont.length(); ++j) {
					std::cout << " " << (size_t)cont.at(j);
				}
				std::cout << std::endl;
			}
			for(std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > >::const_iterator it = alignment_models.at(acc).at(a2).begin();
			it!=alignment_models.at(acc).at(a2).end(); ++it) {
				std::string cont = it->first;
				std::cout << "l ";
				for(size_t j=0; j<cont.length(); ++j) {
					std::cout << " " << (size_t)cont.at(j);
				}
				std::cout << std::endl;
			}



*/
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
) {  // TODO set accession names in data
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

		for(size_t a2 = 0; a2<num_acc; ++a2) {
			binary_read(in, alignments_all_modification_total.at(acc).at(a2));
			mapmodel_read(in, alignment_models.at(acc).at(a2), alignment_alphabet, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+1);


		}	

	}
}



/*
 Read parameters and store them into this object
*/
void dynamic_mc_model::set_patterns(std::ifstream & in){
	std::vector<std::string> acc_names; 
	read_paramters(in, acc_names, 
			this->sequence_models,
			this->acc_sequence_all_bases_total,
			this->acc_sequence_alignment_flag,
			this->alignments_all_modification_total,
			this->alignment_models);
	for(size_t i=0; i<acc_names.size(); ++i) {
		data.add_accession(acc_names.at(i));
	}
	compute_sequence_model_decoding();
	compute_alignment_model_decoding();

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

#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
		for(size_t b = 0; b < num_batchs; ++b) {
			size_t from = b*batch_size;
			size_t to = from + batch_size;
			if(to > num) to = num;
			
			a_reader_counter local_counts;
			acontext_counter local_counter(local_counts,  MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			



			for(size_t i=from; i<to; ++i) {
				const pw_alignment * al = sorted_als.at(a1).at(a2).at(i);

				lensum1+= al->alignment_length();
				local_counter.read_alignment_1_2(*al);
			}
#pragma omp critical(counts)
{
			alignment_counts_add(result_counts.at(a1).at(a2), local_counts.countsmap);
			alignment_length_sums.at(a1).at(a2)+=lensum1;
}

		}


		num = sorted_als.at(a2).at(a1).size();
		batch_size = num / num_threads + 1;
		if(batch_size > max_batch_size) batch_size = max_batch_size;
		num_batchs = num / batch_size + 1;
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
		for(size_t b = 0; b < num_batchs; ++b) {
			size_t from = b*batch_size;
			size_t to = from + batch_size;
			if(to > num) to = num;
			a_reader_counter local_counts;
			acontext_counter local_counter(local_counts,  MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			



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

		std::cout << " run context selector on " << counts.size() << " contexts " << " from " <<  alignment_length_sums.at(a1).at(a2) << " total alignment length" << std::endl;
		double enc_cost = context_selector(counts, alphabet, new_total, mod_model, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL +1);
		std::cout << " select " << mod_model.size() << " contexts, cost: " << enc_cost << " ( " << (enc_cost/(double) alignment_length_sums.at(a1).at(a2))<< " bits/alignment column)" << std::endl;
		alignment_models.at(a1).at(a2) = mod_model;
	}
	compute_alignment_model();

}


void dynamic_mc_model::train(std::ofstream & outs){	
	train_sequence_model();	
	train_alignment_model();
	write_parameters(outs);


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

/*	for(size_t i=0; i<data.numAlignments(); ++i) {
		const pw_alignment & p = data.getAlignment(i);
	//	double c1, c2, m1, m2;
		double g1, g2;
		gain_function(p,g1,g2);
		double average_gain = g1+g2/2;
		if(average_gain <= 0){
			
		}
	//	cost_function(p, c1, c2, m1, m2);
	}*/

}

void dynamic_mc_model::cost_function(const pw_alignment& p , double & c1 , double & c2 , double & m1 , double & m2)const{

	if(p.is_cost_cached()) {
		c1 = p.get_create1();
		c2 = p.get_create2();
		m1 = p.get_modify1();
		m2 = p.get_modify2();
		return;
	}



	size_t a1 = data.accNumber(p.getreference1());
	size_t a2 = data.accNumber(p.getreference2());
	vector<double> ccost(2,0);
	vector<double> mcost(2,0);
// compute alignment cost
	const std::map< std::string, std::vector<double> > & model12 = alignment_model_costs.at(a1).at(a2);
	const std::map< std::string, std::vector<double> > & model21 = alignment_model_costs.at(a2).at(a1);

	a_reader_costs costs12(model12);
	a_reader_costs costs21(model21);

	acontext_cost ccost12(costs12, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
	acontext_cost ccost21(costs21, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
	
	ccost12.read_alignment_1_2(p);
	ccost21.read_alignment_2_1(p);

	m1=costs12.cost_sum + alignment_base_cost.at(a1).at(a2);
	m2=costs21.cost_sum + alignment_base_cost.at(a2).at(a1);

	mcost.at(0) = m1;
	mcost.at(1) = m2;

// compute sequence cost
	const std::map< std::string, std::vector<double> > & model1 = sequence_model_costs.at(a1);
	const std::map< std::string, std::vector<double> > & model2 = sequence_model_costs.at(a2);

//	std::cout << " acc 1 " << p.getreference1() << " costs size " << sequence_model_costs.at(a1).size() << " hv size " << sequence_model_highs.at(a1).size() << std::endl;
	s_reader_costs costs1(model1);
	s_reader_costs costs2(model2);

	scontext_cost scost1(costs1, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
	scontext_cost scost2(costs2, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
	
	scost1.read_sequence(data.getSequence(p.getreference1()), p.getbegin1(), p.getend1());
	scost2.read_sequence(data.getSequence(p.getreference2()), p.getbegin2(), p.getend2());

	c1 = costs1.cost_sum;
	c2 = costs2.cost_sum;

	ccost.at(0) = c1;
	ccost.at(1) = c2;

	p.set_cost(ccost, mcost);

}

/*inline void dynamic_mc_model::gain_function(const pw_alignment& p , double & g1 , double & g2)const{
	double c1;
	double c2;
	double m1;
	double m2;
	cost_function(p, c1, c2, m1,m2);

	
	g1 = c2 - m1;
	g2 = c1 - m2;
}*/
double dynamic_mc_model::get_the_gain(const pw_alignment & p, std::string  & center)const{
	double g1;
	double g2;
	gain_function(p,g1,g2);
	std::vector<std::string> center_parts;
	strsep(center, ":" , center_parts);
	unsigned int dir = atoi(center_parts.at(0).c_str());
	unsigned int ref = atoi(center_parts.at(1).c_str());
	unsigned int left = atoi(center_parts.at(2).c_str());
	if(p.getreference1()== ref){
		return g1; //if center is on ref1
	}else{
		return g2; //if center is on ref1
	}
}
const std::vector<uint32_t> & dynamic_mc_model::get_seq_high_at_position(size_t & ref, size_t & position)const{
	const dnastring & sequence = data.getSequence(ref);
	size_t accession = data.accNumber(ref);
	std::stringstream context;
	for(size_t j = MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; j>0; j--){
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
/*	if(position == 0 || position < 6){
		std::cout << "pattern "<< current_pattern << std::endl;
	}*/
//	std::cout << current_pattern <<std::endl;
	std::map<std::string, std::vector<uint32_t> >::const_iterator it=sequence_model_highs.at(accession).find(current_pattern);
	assert(it->second.at(it->second.size()-1)  == acc_sequence_all_bases_total.at(accession));
	assert(it!=sequence_model_highs.at(accession).end());
	return it->second;
}

const std::vector<uint32_t> & dynamic_mc_model::get_center_high_at_position(size_t & ref, size_t & left, size_t & position)const{
	const dnastring & sequence = data.getSequence(ref);
	size_t accession = data.accNumber(ref);
	std::stringstream context;
	for(size_t j = MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; j>0; j--){
		if(position < left+j){
			char chr = 'A';
			context<<chr;
		}else{
			char chr = sequence.at(position-j);
			context<<chr;
		}
	}
	std::string current_pattern;
	context >> current_pattern;
//	std::cout << "pattern "<< current_pattern << std::endl;
	std::map<std::string, std::vector<uint32_t> >::const_iterator it=center_model_highs.at(accession).find(current_pattern);
//	std::cout << it->second.at(4) <<std::endl;
	assert(it!=center_model_highs.at(accession).end());
	return it->second;
}
void dynamic_mc_model::calculate_center_flags(vector<size_t> & counts , uint32_t & target_total , unsigned char & model_index, std::vector<uint32_t> & bits, std::vector<uint32_t> & widths){
	counts_to_bits_eval(counts,target_total, model_index , bits);
	model_to_widths(bits,model_index,target_total,widths);
	target_total = widths.at(0);
}

void dynamic_mc_model::calculate_center_high_values(uint32_t & target_total){ //calculates the high values using the given total
	center_model_highs.resize(data.numOfAcc());
	vector<unsigned char> alphabet;
	get_sequence_alphabet(alphabet);	
	std::vector<std::map< std::string, std::vector<double>  > > center_model_costs(data.numOfAcc());
	for(size_t acc =0; acc < data.numOfAcc(); acc++){
		model_bits_to_full_high_values(sequence_models.at(acc), target_total, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, alphabet, center_model_highs.at(acc), center_model_costs.at(acc));
	//	for(std::map<std::string, std::vector<uint32_t> >::iterator it = center_model_highs.at(acc).begin(); it!= center_model_highs.at(acc).end();it++){
	//		std::cout << it->second.at(4) << std::endl;
	//	}
	}
	
}

void dynamic_mc_model::get_seq_flag(size_t & acc, std::vector<uint32_t> & al_begin, std::vector<uint32_t> & seq_acc_end)const{
	al_begin.resize(2);
	seq_acc_end.resize(2);
	al_begin.at(0) =  acc_sequence_all_bases_total.at(acc);
	al_begin.at(1) = acc_sequence_alignment_flag.at(acc);
	seq_acc_end.at(0) = acc_sequence_alignment_flag.at(acc);
	seq_acc_end.at(1) =  TOTAL_FOR_ENCODING;
}

void dynamic_mc_model::get_end_al_flag(size_t & acc1 , size_t & acc2, std::vector<uint32_t> & al_end){
	al_end.resize(2);
	al_end.at(0) = alignments_all_modification_total.at(acc1).at(acc2);
	al_end.at(1) = TOTAL_FOR_ENCODING;
}

void dynamic_mc_model::get_center_flag(std::vector<uint32_t> & bits, unsigned char model_index, uint32_t & target_total, std::vector<uint32_t> & widths)const{
	model_to_widths(bits,model_index,target_total,widths);
}

const std::map< std::string, std::vector<uint32_t>  > & dynamic_mc_model::get_sequence_model_highs(size_t & acc)const{
//	std::cout << "its size is "<<std::endl;
//	std::cout << sequence_model_highs.at(acc).size() << std::endl;
	return sequence_model_highs.at(acc);
}
const std::map< std::string, std::vector<uint32_t>  > & dynamic_mc_model::get_center_model_highs(size_t & acc)const{
//	std::cout << center_model_highs.size() <<std::endl;
	return center_model_highs.at(acc);
}
const std::map< std::string, std::vector<uint32_t>  > & dynamic_mc_model::get_alignment_model_highs(size_t & acc1, size_t & acc2)const{
//	std::cout << "here! "<< std::endl;
//	std::cout << alignment_model_highs.at(acc1).at(acc2).size() <<std::endl;
/*	if(acc1 == 1 && acc2 == 0){
		std::cout << "high values from 1 to 0" <<std::endl;
		const std::map< std::string, std::vector<uint32_t> > & hi12 = alignment_model_highs.at(acc1).at(acc2);
		for(std::map<std::string, std::vector<uint32_t> >::const_iterator it = hi12.begin(); it!= hi12.end();it++){
			std::cout << it->first <<std::endl;
			for(size_t i =0 ; i < it->second.size();i++){
				std::cout << it->second.at(i)<<std::endl;
			}
		}
	}*/
	return alignment_model_highs.at(acc1).at(acc2);
}
void dynamic_mc_model::arith_encode_al(const pw_alignment & p, unsigned int & cent_ref, unsigned int & cent_left, std::vector<std::vector< unsigned int> > & low_high){
//	p.print();
	size_t l1,l2,r1,r2;
	p.get_lr1(l1,r1);
	p.get_lr2(l2,r2);
	size_t acc1 = data.accNumber(p.getreference1());
	size_t acc2 = data.accNumber(p.getreference2());
/*	if(acc1 == 1 && acc2 == 0){
		std::cout << "high values from 1 to 0" <<std::endl;
		const std::map< std::string, std::vector<uint32_t> > & hi12 = alignment_model_highs.at(acc1).at(acc2);
		for(std::map<std::string, std::vector<uint32_t> >::const_iterator it = hi12.begin(); it!= hi12.end();it++){
			std::cout << it->first <<std::endl;
			for(size_t i =0 ; i < it->second.size();i++){
				std::cout << it->second.at(i)<<std::endl;
			}
		}
	}*/
	if(p.getreference1()==cent_ref && l1 == cent_left){
	//	std::cout << "last high value "<< std::endl;
		const std::map< std::string, std::vector<uint32_t> > & high12 = alignment_model_highs.at(acc1).at(acc2);
	//	for(std::map<std::string, std::vector<uint32_t> >::const_iterator it = high12.begin(); it!= high12.end();it++){
		//	std::cout << it->second.at(it->second.size()-1) << std::endl;
		//	std::cout << "length of context " << it->first.length() << std::endl;
	//	}
		a_reader_encode encode12(high12, wrappers);
		acontext_encode enc12(encode12, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
		enc12.read_alignment_1_2(p);
		low_high = encode12.low_high_values;
	}else{
		const std::map< std::string, std::vector<uint32_t> > & high21 = alignment_model_highs.at(acc2).at(acc1);
		a_reader_encode encode21(high21,wrappers);
		acontext_encode enc21(encode21, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
		enc21.read_alignment_2_1(p);
		low_high = encode21.low_high_values;
	}
}
void dynamic_mc_model::arith_encode_long_al(const pw_alignment & p, size_t & acc1, size_t & acc2, unsigned int & cent_ref, unsigned int & cent_left, std::vector<std::vector< unsigned int> > & low_high){
		const std::map< std::string, std::vector<uint32_t> > & high21 = alignment_model_highs.at(acc2).at(acc1);
		if(p.getreference1()==0 && p.getreference2()==984){
			for(std::map<std::string, std::vector<uint32_t> >::const_iterator it = high21.begin(); it!= high21.end();it++){
				std::string patterns = it->first;
				std::cout << patterns.size() <<std::endl;
				for(size_t i =0; i < patterns.size();i++){
					std::cout << int(patterns.at(i)) << " ";
				}
				std::cout << " " << std::endl;
				std::cout << "length of context " << it->first.length() << std::endl;
			}
		}
		std::cout << high21.size() << std::endl;
		a_reader_encode encode21(high21,wrappers);
		acontext_encode enc21(encode21, MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
		enc21.read_alignment_2_1(p);
		low_high = encode21.low_high_values;
		std::cout << "size of low_high is "<< low_high.size()<<std::endl;
}

void dynamic_mc_model::write_al_high_onstream(std::ofstream & al_high_out){
	for(size_t i =0; i < data.numOfAcc();i++){
		for(size_t j =0; j < data.numOfAcc();j++){
			for(std::map< std::string, std::vector<uint32_t> >::iterator it = alignment_model_highs.at(i).at(j).begin(); it != alignment_model_highs.at(i).at(j).end();it++){
				for(size_t m =0 ; m < it->second.size();m++){
					al_high_out<<char(0);
					al_high_out << it->second.at(m);
				//	al_high_out.write(reinterpret_cast<char*>(&highValue), sizeof(uint32_t));					
				}
			}
		}
	}
}
