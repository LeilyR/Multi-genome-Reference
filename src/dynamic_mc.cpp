#include "dynamic_mc.hpp"

#ifndef DYNAMIC_MC_CPP
#define DYNAMIC_MC_CPP


dynamic_mc_model::dynamic_mc_model(all_data & d): data(d),create_cost(data.numOfAcc(),vector<double>(5,0.0)),sequence_successive_bases(data.numOfAcc()),mod_cost(data.numOfAcc(),std::vector<std::map <std::string, std::vector<double> > > (data.numOfAcc())),high(data.numOfAcc()),highValue(data.numOfAcc(),std::vector<std::map<std::string , std::vector<unsigned int> > >(data.numOfAcc())),al_level(data.numOfAcc(),std::vector<size_t>(data.numOfAcc())){
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

void dynamic_mc_model::train_sequence_model(){//TODO What should i do in case of level 0? should it even be taken into account? If yes, then i have to think about filling out the map of highi values. I may need a new map just for l = 0. 

	std::string astring("");
	astring.append(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, 'A');

	for(size_t acc=0; acc<data.numAcc(); ++acc) {
		std::map<std::string, std::vector<size_t> > context_counter;
		std::map<std::string, std::vector<size_t> > best_counter;
		size_t best_level = 0;
		double best_cost = std::numeric_limits<double>::max();
		// for the highest markov chain level we fill the map directly from the sequences
		for(size_t k = 0; k <data.getAcc(acc).size(); k++){
			std::string sequence = data.getSequence(data.getAcc(acc).at(k)).str();

			// we look for the base at i, given its sequence context. At the beginning of each sequence we fill contexts from the left with A
			// first position where we dont need filling is at context length + 1

			for(size_t p=0; p< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL && p < sequence.length() ; ++p) {
				std::string context(astring, 0, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL - p); // left side A's
				context.append(sequence.substr(0, p));
				context_count(context_counter, context, dnastring::base_to_index(sequence.at(p)), 1);
			}
			for(size_t p = MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; p<sequence.length(); ++p) {
				std::string context = sequence.substr(p-MAX_SEQUENCE_MARKOV_CHAIN_LEVEL, MAX_SEQUENCE_MARKOV_CHAIN_LEVEL);
				context_count(context_counter, context, dnastring::base_to_index(sequence.at(p)), 1);
			}
		}

		for(size_t level = MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; level>=0; --level) {
			size_t max_contexts = 1;
			for(size_t i=0; i<level; ++i) {
				max_contexts *= 5;
			}
			double parameter_cost = (double)max_contexts * 5*  BITS_PER_SEQUENCE_HIGH_VALUE;
			double encoding_cost = total_sequence_cost(context_counter);
			double all_cost = parameter_cost + encoding_cost;
			std::cout << " on accession " << acc << " level " << level<< " param " << parameter_cost << " enc " << encoding_cost << " total " << all_cost << std::endl;

			if(all_cost < best_cost) {
				std::cout << " new best!" << std::endl;
				best_counter = context_counter;
				best_cost = all_cost;
				best_level = level;
			
			}

			// context_conter one level down:
			if(level > 0) {	
				std::map<std::string, std::vector<size_t> > new_counter;							
				for(std::map<std::string, std::vector<size_t> >::iterator it = context_counter.begin(); it!=context_counter.end(); ++it) {
					std::string oldcontext = it->first;
					std::string newcontext = oldcontext.substr(1);
					for(size_t i=0; i<5; ++i) {
						context_count(new_counter, newcontext, i, it->second.at(i) );
					}
				}
				context_counter = new_counter;
			}
		}
		accLevel.push_back(best_level);
		calculate_sequence_high(acc, context_counter);	
	}
/*
   TODO go on to compute high values
	for(size_t i =0; i < data.numAcc(); i++){
		double total_cost = 0;
		//An initial cost for each accession(mc of level 0)
		std::vector<double> cost_on_acc (5,0);
		size_t total_base_number = 0;
		for(size_t k = 0; k <data.getAcc(i).size(); k++){
			const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
			std::vector<size_t> number(5,0);
			total_base_number += sequence.length();
			for(size_t j = 0 ; j< sequence.length(); j++ ){
				size_t base = dnastring::base_to_index(sequence.at(j));
				number.at(base) ++ ;
			}
			for(size_t j = 0; j< 5; j++){	
				cost_on_acc.at(j) += number.at(j);
			}
		}
		for(size_t j=0; j<5;j++){
			cost_on_acc.at(j)= -log2(cost_on_acc.at(j)/total_base_number);
		}
		for(size_t j = 0; j < 5 ; j++){
			create_cost.at(i).at(j) = cost_on_acc.at(j);
		}
		for(size_t j =0; j < 5 ; j++){
			total_cost += cost_on_acc.at(j);
		}
		//Till here
		double new_cost = 0;
		size_t level = 1;
		double old_cost = total_cost;
		while(new_cost < old_cost){
			if(level != 1){
				old_cost = new_cost;
			}
			markov_chain(i,level);
			for(size_t k=0; k<create_cost.at(i).size();k++){
				new_cost += create_cost.at(i).at(k);
			}
			level = level+1;
		}
		accLevel.push_back(level-2);
		sequence_successive_bases.clear();// Only 'seq's that are related to the final level should be kept and the rest should be deleted.
		markov_chain(i,level-2);
		std::cout<< "level is "<< level-2 <<std::endl;
		calculate_sequence_high(i);	}
*/	
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
	size_t level = accLevel.at(accession);
	if(level > 0){
		make_all_patterns(level);
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
}
void dynamic_mc_model::make_all_patterns(const size_t& level){
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
}
void dynamic_mc_model::write_parameters(std::ofstream & outs){//write accessions, their levels and high values to a file
	size_t bit = 12;
	for(size_t i = 0 ; i < data.numAcc(); i++){
		outs << data.get_acc(i);
		outs<< (char)0;
		size_t level = accLevel.at(i);
		outs << level;
		for(std::map<std::string, std::vector<unsigned int> >::iterator it= high.at(i).begin(); it != high.at(i).end() ; it++){	
			std::vector<bool> bit_to_byte(0);
			for(size_t j = 0; j < 5; j++){
				int h = it->second.at(j);
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
	outs<<(char)8;
}

void dynamic_mc_model::set_patterns(std::ifstream & in){
	size_t bit = 12;
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
		make_all_patterns(level);// TODO later on change it to retrieved_level	
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
}
/*
	train alignment/modification with a dynamic model
*/

void dynamic_mc_model::recursive_al_model(){//TODO again the same question about level 0
	counting_functor functor(data);
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
}
void dynamic_mc_model::make_all_alignment_pattern(const size_t & level){ // it is only used when level > 0
	assert(level != 0);
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

}
void dynamic_mc_model::calculate_alignment_high(size_t & acc1, size_t & acc2){//TODO level0?
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
}
void dynamic_mc_model::write_alignment_parameters(std::ofstream & outs){//since we alreay wrote the accession names on the file and can retrieve them from set parameters here we can only save acc_number not the name anymore.
	size_t bit = 12;
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
}

void dynamic_mc_model::set_alignment_pattern(std::ifstream & in){
	size_t bit = 12;
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
}
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

void dynamic_mc_model::train(std::ofstream & outs){//call it before training mc model to find the best level.
	train_sequence_model();	
	recursive_al_model();
	write_parameters(outs);
	write_alignment_parameters(outs);
	std::cout<< "dynamic train is done! "<<std::endl;
}

void dynamic_mc_model::cost_function(const pw_alignment& p , double & c1 , double & c2 , double & m1 , double & m2)const{//TODO add the possibilty of level = 0
	if(p.is_cost_cached()) {
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
size_t dynamic_mc_model::get_acc_level(size_t & acc)const{
	return accLevel.at(acc);
}
void dynamic_mc_model::set_acc_level(size_t & level){
	accLevel.push_back(level);
}
size_t dynamic_mc_model::get_al_level(size_t& acc1 , size_t & acc2)const{
	return al_level.at(acc1).at(acc2);

}
void dynamic_mc_model::set_al_level(size_t& acc1 , size_t & acc2, size_t & level){
	al_level.at(acc1).at(acc2) = level;
}
//std::vector<size_t> dynamic_mc_model::get_powerOfTwo()const{
//	return powersOfTwo;
//}	
const std::map<std::string, std::vector<unsigned int> >&  dynamic_mc_model::get_high(size_t acc)const{
	return high.at(acc);
}
std::string dynamic_mc_model::get_firstPattern(size_t& accession)const{
	std::string first_pattern;
	size_t level = accLevel.at(accession);
	for(size_t i = 0; i< level; i++){
		first_pattern += "A";
	}
	return first_pattern;
}
std::string dynamic_mc_model::get_firstAlignmentPattern(size_t& acc1, size_t& acc2)const{
		std::string first_pattern;
		size_t level = al_level.at(acc1).at(acc2);
		size_t firstPatterns = level;
		for(size_t j = 0; j < level; j++){
			firstPatterns --;
			first_pattern += modification_character(-1,-1,-1,firstPatterns);
		}
		return first_pattern;
	}	
void dynamic_mc_model::get_encoded_member(pw_alignment & al, size_t center_ref, size_t center_left, encoding_functor & functor,std::ofstream& outs)const{
	size_t acc1 = data.accNumber(al.getreference1());
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
}

std::vector<unsigned int> dynamic_mc_model::get_high_at_position(size_t seq_index, size_t position)const{
	const dnastring & sequence = data.getSequence(seq_index);
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
}
const std::map<std::string, std::vector<unsigned int> > & dynamic_mc_model::get_highValue(size_t acc1, size_t acc2)const{
	return highValue.at(acc1).at(acc2);
}
#endif
