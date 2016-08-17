#include "dynamic_encoder.hpp"

#ifndef DYNAMIC_ENCODER_CPP
#define DYNAMIC_ENCODER_CPP

	template<typename T>
	dynamic_decoder<T>::dynamic_decoder(T & m , decoding_wrapper & wrapper):model(m),wrappers(wrapper),center_flags(3,0){}
	template<typename T>
	dynamic_decoder<T>::~dynamic_decoder(){}
	template<typename T>
	void dynamic_decoder<T>::set_center_from_stream(std::ifstream & in){
		size_t bit = 32;
		size_t partition_size = 0;
		char c ; 
		char c1[sizeof(size_t)];
		in.read(c1,sizeof(size_t));
		size_t *temp = reinterpret_cast<size_t*> (c1);
		partition_size = *temp;
		std::vector<std::pair<size_t, vector<unsigned int> > >value;
		for(size_t i = 0; i < partition_size; i++){
			cluster_high_partition.push_back(value);
		}
		unsigned char h;
		c = in.get();
		size_t i = 0;
		size_t index = 0;
		while(c != 8){
			unsigned int low = 0;
			while(c != 7){
				size_t acc;
				in.read(c1,sizeof(size_t));
				size_t *temp = reinterpret_cast<size_t*> (c1);
				acc = *temp;
				//accession of a center
				std::map<size_t , std::vector<size_t> >::iterator it = acc_of_center.find(acc);
				if(it == acc_of_center.end()){
					acc_of_center.insert(std::make_pair(acc, std::vector<size_t>()));
				}
				std::map<size_t , std::vector<size_t> >::iterator it1 = acc_of_center.find(acc);
				it1 ->second.push_back(index);
				std::vector<unsigned int> temp1(2,0);
				std::pair<size_t,std::vector<unsigned int> >Centerhigh;
				Centerhigh.second = temp1;
				Centerhigh.first = index;				
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t j = 0; j < binary_high_value.size();j++){
					high_value += binary_high_value.at(j)*(1<<j);
				}
				Centerhigh.second.at(0) = low;
				Centerhigh.second.at(1) = high_value;
				low = Centerhigh.second.at(1);
				c=in.get();
				cluster_high_partition.at(i).push_back(Centerhigh);
				std::cout << "acc " << acc<< " center " << index << " low " << Centerhigh.second.at(0) << " high " << Centerhigh.second.at(1) << std::endl;
				index = index +1;
			}
			c = in.get();
			i = i + 1;
			std::cout<< "i: "<< i << std::endl;
		}
	}
	
	template<typename T>	
	void dynamic_decoder<T>::set_partition_high_from_stream(std::ifstream & in){
		size_t bit = 32;
		unsigned char h;
		char c;
		c = in.get();
		unsigned int low = 0;
		while(c != 8){
			while(c!=7){
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t i = 0; i < binary_high_value.size();i++){
					high_value += binary_high_value.at(i)*(1<<i);
				}
				partitionHigh.push_back(high_value);
				std::cout<< "high_read: "<< high_value<<std::endl;
				c=in.get();
			}
			c=in.get();
		}
	}

	template<typename T>
	void dynamic_decoder<T>::set_center_flags_from_the_stream(std::ifstream & in){
		std::vector<uint32_t> bits(2,0);
		unsigned char index;
		index=in.get();
		char c[sizeof(uint32_t)];
		for(size_t i =0; i < bits.size();i++){
			in.read(c,sizeof(uint32_t));
			uint32_t *temp = reinterpret_cast<uint32_t*> (c);
			bits.at(i) = *temp;
		}
		unsigned int total = TOTAL_FOR_ENCODING;
		model.get_center_flag(bits,index, total, center_flags);
		std::cout << "first flag " << std::endl;
		std::cout << center_flags.at(0)  << " " << center_flags.at(1) <<std::endl;
		std::cout << "center flag size " << center_flags.size() << std::endl;
	}

	template<typename T>
	std::string dynamic_decoder<T>::first_sequence_pattern()const{
		std::string first_pattern;
		for(size_t i = 0; i< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; i++){
			first_pattern += "A";
		}
		return first_pattern;
	}
	template<typename T>
	std::string dynamic_decoder<T>::first_alignment_pattern()const{
		std::string first_context;
	//	size_t firstPatterns = MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL;
	//	for(size_t j = 0; j < MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; j++){
	//		firstPatterns --;
	//		first_pattern += model.modification_character(-1,-1,-1,firstPatterns);
	//	}
	//	return first_pattern;
		for(size_t i= MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; i>0; --i) {
			size_t numk = i-1;
			if(numk > NUM_KEEP_DYN-1) numk = NUM_KEEP_DYN;
			char keepc = dynamic_mc_model::modification_character(-1, -1, -1, numk);
			first_context.append(1, keepc);
		}
	//	std::cout << "first al context is "<<int(first_context.at(0)) << " " << int(first_context.at(1)) << std::endl;
		return first_context;
	}
	template<typename T>
	void dynamic_decoder<T>::arithmetic_decoding_centers(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec){
		size_t total_length = 0;
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
		set_center_flags_from_the_stream(in);
		model.calculate_center_high_values(center_flags.at(0));
		std::cout << "all the parameters are set" <<std::endl;
		size_t acc = 0;
		dec.set_stream(in);	
		std::string first_pattern = first_sequence_pattern();
		unsigned int target;
		size_t lengthOfSeq = 0;
		unsigned int total = TOTAL_FOR_ENCODING;
		size_t numCenter = 0;
		size_t counter =0;
		for(size_t i = 0; i < cluster_high_partition.size(); i++){
			for(size_t j = 0; j < cluster_high_partition.at(i).size();j++){
				numCenter = numCenter +1 ;
			}
		}
		std::cout<< "partition size: "<<cluster_high_partition.size()<<std::endl;
		std::cout << " numbe of centers: "<< numCenter <<std::endl;
		while(counter < cluster_high_partition.size()){
			std::vector<size_t> centerOfAPartition;
			for(size_t i =0; i < cluster_high_partition.at(counter).size();i++){
				centerOfAPartition.push_back(cluster_high_partition.at(counter).at(i).first);
			}
			std::cout<< "size of centerOfAPartition: "<<centerOfAPartition.size()<<std::endl;
			std::map<size_t, std::string> intermediate;
			for(size_t i = 0; i < centerOfAPartition.size() ; i++){//All the centers of a partition
				size_t c_id = centerOfAPartition.at(i);
				std::cout << "center id is " << c_id <<std::endl;
				for(std::map<size_t, std::vector<size_t> >::iterator it = acc_of_center.begin();it != acc_of_center.end(); it ++){	
					for(size_t k = 0; k < it->second.size();k++){		
						if(it ->second.at(k) == c_id){
							acc = it->first;
							break;
						}
					}
				}
				std::string center_seq;
				std::vector<uint32_t>high(6,0);
				std::vector<uint32_t>low(6,0);
				lengthOfSeq = 1;
				std::string pattern = first_pattern;
				std::cout << "first pattern " << first_pattern << " acc "<< acc <<std::endl;
				std::map<std::string, std::vector<uint32_t> >::const_iterator it=model.get_center_model_highs(acc).find(first_pattern);
				std::cout << it ->second.at(4) <<std::endl;
				std::cout << "here! "<<std::endl;
				assert(it != model.get_center_model_highs(acc).end());
			//	std::cout << it->second.size()<<std::endl;
				for(size_t j=0; j<5; j++){
					high.at(j) = it->second.at(j);
					if(j!= 0){
						low.at(j)=high.at(j-1);
					}else low.at(j)= 0;
				}
				low.at(5)= high.at(4);
				std::cout<< "high.at(4)" << high.at(4)<<std::endl;
				high.at(5)= center_flags.at(0)+center_flags.at(1);
				std::cout << "high.at(5) " << high.at(5) <<std::endl;
				size_t base = 10;
				target = dec.get_target(total);
				std::cout << "target "<< target << std::endl;
				for(size_t n = 0; n < 5 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				assert(base < 5);
				std::cout << "first base " << base << std::endl;
				center_seq +=dnastring::index_to_base(base);
				dec.decode(low.at(base),high.at(base));
				wrappers.decode(low.at(base),high.at(base),total);
				target = dec.get_target(total);
				std::cout << " center_flags.at(0) " << center_flags.at(0) <<std::endl;
				while(target < center_flags.at(0)){//End of a center	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					std::string current_pattern;
					std::stringstream s;
					if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL > 1){
						for(size_t M=1; M< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; M++){
							s<<pattern.at(M);
						}
						s<<p;
						s>>current_pattern;
					}
					if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL == 1){
						s<<p;
						s>>current_pattern;
					}
				//	std::cout<<"cur pattern: "<< current_pattern<<std::endl;
				//	std::cout << "base " << base << std::endl;
					std::map<std::string, std::vector<uint32_t> >::const_iterator it1=model.get_center_model_highs(acc).find(current_pattern);
					if(it1 == model.get_center_model_highs(acc).end()){
						std::cout<<"current pattern: "<< current_pattern<<std::endl; 		
					}
					assert(it1 != model.get_center_model_highs(acc).end());
					for(size_t j=0; j<5; j++){
						high.at(j) = it1->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}	
					base = 10;
					for(size_t n = 0; n < 5 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					center_seq +=dnastring::index_to_base(base);
					dec.decode(low.at(base),high.at(base));
					wrappers.decode(low.at(base),high.at(base),total);
				//	std::cout << "base: "<< base << std::endl; 
					pattern = current_pattern;
					target = dec.get_target(total);
				//	std::cout << "target "<< target << std::endl;
				}
				std::cout << "length of center " << lengthOfSeq << std::endl;
				total_length +=lengthOfSeq;
				std::cout << "total length " << total_length << std::endl;
			//	save1 << "length "<< lengthOfSeq <<std::endl;
				base = 10;
				for(size_t n = 0; n < 6 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				assert(base == 5);
				std::cout<<"flag base : "<< base << std::endl;
				dec.decode(low.at(base),high.at(base));
				wrappers.decode(low.at(base),high.at(base),total);
				size_t center_id = centerOfAPartition.at(i);
				intermediate.insert(std::make_pair(center_id,center_seq));
			}
			decoded_center_in_partition.push_back(intermediate);
			std::cout << " counter "<< counter <<std::endl;
			counter = counter + 1 ;
		}
	}

	template<typename T>
	void dynamic_decoder<T>::decoding_directions(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, bool& reverse_center, bool& reverse_member, bool& reverse_both, bool& forward_both, uint32_t & target){
                unsigned int total = TOTAL_FOR_ENCODING;
		unsigned int lowest = 0;
		std::vector<unsigned int> high(4);
		high.at(0) = TOTAL_FOR_ENCODING/4;
		high.at(1) = TOTAL_FOR_ENCODING/2;
		high.at(2) = 3*TOTAL_FOR_ENCODING/4;
		high.at(3) = TOTAL_FOR_ENCODING;
		//Reverse center
		if(target>= lowest && target < high.at(0)){
			std::cout<<"rev_center"<<std::endl;
			dec.decode(lowest, high.at(0));
			wrappers.decode(lowest,high.at(0),total);
			reverse_center = true;
			target = dec.get_target(total);
		}
		//Reverse member
		else if(target>= high.at(0) && target< high.at(1)){
			std::cout<<"rev_member"<<std::endl;
			dec.decode(high.at(0), high.at(1));
			wrappers.decode(high.at(0), high.at(1), total);
			reverse_member = true;
			target = dec.get_target(total);
		}
		//Both reverse
		else if(target>= high.at(1) && target< high.at(2)){
			std::cout<<"rev_both"<<std::endl;
			dec.decode(high.at(1), high.at(2));
			wrappers.decode(high.at(1) ,high.at(2),total);
			reverse_both = true;
			target = dec.get_target(total);
		}else{
			assert(target>= high.at(2) && target< high.at(3));
			std::cout<<"forward_both"<<std::endl;
			dec.decode(high.at(2), high.at(3));
			wrappers.decode(high.at(2), high.at(3),total);
			forward_both = true;
			target = dec.get_target(total);
			std::cout << "target when both are forward "<< target << std::endl;


		}
	}
	template<typename T>
	std::string dynamic_decoder<T>::associatedMember(std::string & center, size_t & modificationPattern, size_t & position){
		std::string member = "";
		if(modificationPattern<5){
			char base = dnastring::index_to_base(modificationPattern);
			member += base;	
			std::cout<< "modification"<<std::endl;
		//	std::cout<< "member is "<<member<<std::endl;
			return member;

		}
		if(modificationPattern>=5 && modificationPattern<5+NUM_DELETE_DYN){
			std::cout<<"deletion happened!"<<std::endl;
			return member;
		}
		if(modificationPattern>=5+NUM_DELETE_DYN && modificationPattern<5+NUM_KEEP_DYN+NUM_DELETE_DYN){
			size_t length = model.modification_length(modificationPattern);
			std::cout << "center length "<< center.size() << " pattern length: "<< length <<" pattern index : "<< modificationPattern << std::endl;
			for(size_t i = 0 ; i < length;i++){
				assert(position+ i <center.length());
				member += center.at(position+i);
			//	std::cout<< "member is "<<member<<std::endl;
			}
			return member;
		}
		if(modificationPattern>=5 + NUM_KEEP_DYN + NUM_DELETE_DYN){
			char base = dnastring::index_to_base(modificationPattern-20);
			member += base;
			std::cout<< "insertion"<<std::endl;
		//	std::cout<< "member is "<<member<<std::endl;

			return member;

		}
		assert(0);
		return "";
	}

	template<typename T>
	void dynamic_decoder<T>::al_decoding(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out){
		model.set_patterns(in);//Read the model parameters from the stream and put them in the corresponding maps
		std::cout<< "alignment patterns are retrieved! "<<std::endl;
		arithmetic_decoding_centers(in,dec);
	/*	std::cout<<"decoding the center has been done!"<<std::endl; //Just to test the retrieved high values.
		std::ofstream al_high_decode("al_high_decode.txt");
		model.write_al_high_onstream(al_high_decode);
		al_high_decode.close();*/
		size_t base = 0;
		size_t sequence_counter = 0;
		out.write(reinterpret_cast<char*>(&sequence_counter),sizeof(uint32_t));
		unsigned int total = TOTAL_FOR_ENCODING;
		std::string first_al_pattern = first_alignment_pattern();
		std::string first_pattern = first_sequence_pattern();
		std::string al_pattern;
		unsigned int target;
		size_t acc = 0;
	//	std::ofstream save("decoded_sequences.txt");
	//	save << 0;
		while(acc<acc_of_center.size()){
			std::cout<< "accession "<< acc <<std::endl;
			std::vector<uint32_t> al_begin(2);
			std::vector<uint32_t> seq_acc_end(2);
			model.get_seq_flag(acc,al_begin,seq_acc_end);
			std::cout << "al begin low " << al_begin.at(0) << " al begin high " << al_begin.at(1) <<std::endl;
			std::cout << " seq end low " << seq_acc_end.at(0) << " seq end high " << seq_acc_end.at(1) << std::endl;
			target = dec.get_target(total); 
			while(target < seq_acc_end.at(0)){//end of an accession
				std::string decodedCenter;
				std::string pattern;
		//		low.at(5)= al_begin.at(0);//al_begin flag
		//		high.at(5)= al_begin.at(1);
		//		low.at(6)= high.at(5);//seq and acc end flag
		//		high.at(6) = seq_acc_end.at(1);
				if(target<al_begin.at(0)){//If there is no alignment on the first position
					std::cout<< "Sequence doesn't start with an alignment!" <<std::endl;
					pattern = first_pattern;
					std::cout<< "pattern size: "<<pattern.size()<<std::endl;
					std::vector<unsigned int>high(5,0);
					std::vector<unsigned int>low(5,0);
					std::map<std::string, std::vector<unsigned int> >::const_iterator it=model.get_sequence_model_highs(acc).find(first_pattern);
					for(size_t j=0; j<5; j++){
						high.at(j) = it->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}
					assert(high.at(4) == al_begin.at(0));
					base = 12;
					for(size_t n = 0; n < 5 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					if(base < 5){
						std::cout << " base < 5 " << std::endl;
						char p = dnastring::index_to_base(base);
						std::string test;
						test +=p;
						std::cout << sequence_counter<< " "<< test << std::endl;
						out<<p;	
					//	save << p;
					}
					dec.decode(low.at(base),high.at(base));	
					wrappers.decode(low.at(base),high.at(base),total);
					assert(base < 5);
//					std::cout << "base " << base << endl;
					int base_int = base;
					wrappers.decodeContext(base_int);
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern.
					target = dec.get_target(total);	
				}
				else {
					std::cout<<"Sequence starts with an alignment!" << target <<std::endl;
				}
				
				while(target< seq_acc_end.at(0)){//end of each seq
					size_t length_between_two_centers = 1;
					while(target< al_begin.at(0)){//If there is no alignment on that position (see if i can move it to the else of alignment if)
						char p = dnastring::index_to_base(base);
						std::stringstream s;
						std::string current_pattern;
						if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL > 1){// TODO: Make a function which returns the current pattern!
						//	std::cout << "pattern " << pattern << std::endl;
							for(size_t M=1; M< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; M++){
								s<<pattern.at(M);
							}
						//	std::cout << " " <<std::endl;
							s<<p;
							current_pattern = s.str();
						}
						if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL ==1){	
							s<<p;
							s>>current_pattern;
						}
						std::cout<< "current pattern: " << current_pattern<<std::endl;
						std::vector<unsigned int>high(5,0);
						std::vector<unsigned int>low(5,0);
						std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_sequence_model_highs(acc).find(current_pattern);
						assert(it1 != model.get_sequence_model_highs(acc).end());
						for(size_t j=0; j<5; j++){
							high.at(j) = it1->second.at(j);
							if(j!= 0){
								low.at(j)=high.at(j-1);
							}else low.at(j)= 0;
						}	
						base = 12;
						for(size_t n = 0; n < 5 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						assert(base<5);
						if(base < 5){
							char p1 = dnastring::index_to_base(base);
							out<<p1;
						//	save << p1;
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						int base_int = base;
						wrappers.decodeContext(base_int);
						std::cout<<"base: "<<base << " low " << low.at(base) << " high " << high.at(base) <<std::endl;
						pattern = current_pattern;
						std::cout << "pattern "<< pattern <<std::endl;
						target = dec.get_target(total);	
					//	std::cout << "target2: "<<target<<std::endl;
						length_between_two_centers ++;
						if(target>= al_begin.at(0)){
							std::cout<< "starting of an alignemnt's flag!"<<std::endl;
							std::cout << "length between two centers " << length_between_two_centers << std::endl;
						}
					}
					std::cout << "target "<<target <<std::endl;
					if(al_begin.at(0)<=target && target < seq_acc_end.at(0)){//Beginning of an al but not end of the seq
						if( al_begin.at(0)<=target && target < al_begin.at(1)){
							dec.decode( al_begin.at(0), al_begin.at(1));
							wrappers.decode( al_begin.at(0), al_begin.at(1),total);
						}else{
							std::cout<< "there is a problem in bal flag"<<std::endl;
						}
						size_t center;
						size_t number_of_par = 0;
						target=dec.get_target(total);
						std::cout<< "par target: "<< target <<std::endl;
						unsigned int low_par = 0;
						unsigned int high_par = 0;
						for(size_t j = 0; j < partitionHigh.size(); j++){
							low_par = 0;
							high_par = partitionHigh.at(j);
							if(j != 0){
								low_par = partitionHigh.at(j-1);
							}
							if(target >= low_par && target < high_par){
								number_of_par = j;
								std::cout<< "number of par: "<< number_of_par << std::endl;
								std::cout<< "low par: "<< low_par << " high par"<< high_par <<std::endl;
								break;
							}
						}
						dec.decode(low_par,high_par);//Dont need to be added to the 'out'
						wrappers.decode(low_par,high_par,total);
						target = dec.get_target(total);
						std::cout << "center target " <<target <<std::endl;
						for(size_t j =0; j < cluster_high_partition.at(number_of_par).size();j++){
							if(cluster_high_partition.at(number_of_par).at(j).second.at(0) <= target && cluster_high_partition.at(number_of_par).at(j).second.at(1)> target){
								center = cluster_high_partition.at(number_of_par).at(j).first;
								dec.decode(cluster_high_partition.at(number_of_par).at(j).second.at(0),cluster_high_partition.at(number_of_par).at(j).second.at(1));
								wrappers.decode(cluster_high_partition.at(number_of_par).at(j).second.at(0),cluster_high_partition.at(number_of_par).at(j).second.at(1),total);
								std::cout << cluster_high_partition.at(number_of_par).at(j).second.at(0) << " " << cluster_high_partition.at(number_of_par).at(j).second.at(1) << std::endl;
								break;
							}
						}
						std::cout <<"center: "<< center<< std::endl;	
						std::map<size_t, std::string>::iterator seq = decoded_center_in_partition.at(number_of_par).find(center);
						assert(seq != decoded_center_in_partition.at(number_of_par).end());
						decodedCenter =seq->second;
						std::cout<< "decoded center: "<< decodedCenter <<std::endl;
						size_t cent_acc;
						for(std::map<size_t, std::vector<size_t> >::iterator it_acc = acc_of_center.begin(); it_acc != acc_of_center.end(); it_acc++){
							for(size_t g= 0; g< it_acc->second.size();g++){
								if(it_acc->second.at(g)== center){
									cent_acc = it_acc->first;
									std::cout<< "acc of cent: "<<cent_acc<<std::endl;
									break;
								}else continue;
							}
						}
						std::vector<uint32_t> al_end(2);
					
						model.get_end_al_flag(cent_acc, acc ,al_end);
						std::cout << "al end: " << al_end.at(0) << " "<<al_end.at(1) <<std::endl;
						al_pattern = first_al_pattern;
						size_t pos = 0;
						target = dec.get_target(total);	
						std::cout << "target before direction detection " << target <<std::endl;
						bool reverse_center = false;
						bool reverse_member = false;
						bool reverse_both = false;
						bool forward_both = false;
						decoding_directions(in, dec, reverse_center, reverse_member, reverse_both, forward_both, target);
						size_t modify = 0;
						std::string member;
						std::cout << "target "<< target << " al_end.at(0) "<< al_end.at(0) << " " << al_end.at(1) <<std::endl;
						size_t modify_count = 0;
						while(target < al_end.at(0) && pos < decodedCenter.size()){//while we are on an alignment
							std::vector<unsigned int>al_high(NUM_MODIFICATIONS,0);
							std::vector<unsigned int>al_low(NUM_MODIFICATIONS,0);
							size_t last_base;
							if(reverse_center == true){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							if(reverse_both == true){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							if(forward_both == true){//both references are forward
								std::cout << "both are forward" <<std::endl;
								last_base = dnastring::base_to_index(decodedCenter.at(pos));
							}
							if(reverse_member == true){
								last_base = dnastring::base_to_index(decodedCenter.at(pos));
								std::cout << "member is reversed! " << std::endl;
							}
							std::string al_context = al_pattern;
							al_context += last_base;
						//	std::cout<< "size of context: "<< al_context.size() << std::endl;
							std::cout<< "decoding_context: ";
							for(size_t k =0; k < al_context.size(); k++){
								std::cout<< int(al_context.at(k)) << " ";
								int con =  int(al_context.at(k));
								wrappers.decodeContext(con);
							}
							std::cout<< " " <<std::endl;
							std::map<std::string, std::vector<unsigned int> >::const_iterator mod = model.get_alignment_model_highs(cent_acc,acc).find(al_context);//modifiactions are considered from center to the other reference
							assert(mod != model.get_alignment_model_highs(cent_acc,acc).end());
						//	std::cout << "cent_acc " << cent_acc << " other_acc " << i <<std::endl;
							for(size_t k = 0 ; k < mod->second.size();k ++){
								al_high.at(k) = mod->second.at(k);
								if(k > 0){
									al_low.at(k) = al_high.at(k-1);
								}else{
									al_low.at(k) = 0;
								}				
							}
							std::cout << "mod->second size: " << mod->second.size() << " target: "<< target<<std::endl;
						//	al_high.at(mod->second.size()-1)= TOTAL_FOR_ENCODING;
							assert(al_high.at(mod->second.size() -1) == al_end.at(0));
							size_t modification = NUM_MODIFICATIONS+20;
							for(size_t n = 0; n < NUM_MODIFICATIONS; n++){
								if(al_low.at(n)<=target && al_high.at(n)> target){
									modification = n;
									break;
								}
							}
							assert(modification < NUM_MODIFICATIONS);
							dec.decode(al_low.at(modification),al_high.at(modification)); 
							wrappers.decode(al_low.at(modification),al_high.at(modification),total); 
							modify = modification;
							if(center == 6){
								modify_count++;
							}
							std::cout<< " modify: " << modify << std::endl;
							stringstream s;
							if(MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL > 1){
								for(size_t M=1; M< al_context.size()-1; M++){
									std::cout << int(al_context.at(M)) <<std::endl;
									s << al_context.at(M);
								}
								al_pattern=s.str();
								al_pattern += modify;
							}else{
								al_pattern = modify;
							}
							std::cout<<"al_pattern: ";
							for(size_t al = 0; al < al_pattern.size();al++){
								std::cout<< int(al_pattern.at(al))<< " ";
							}
							std::cout<<" " << std::endl;
							if(reverse_center==true || reverse_both == true){
								std::string rev_decodedCenter;
								for(size_t rev = 0; rev< decodedCenter.size(); rev++){
									char chr= dnastring::complement(decodedCenter.at(decodedCenter.size()-1-rev));	
									rev_decodedCenter +=chr;
								}
								std::cout << "associatedMem_both reverse "<<std::endl;

								member += associatedMember(rev_decodedCenter,modification,pos);							
							}else{
								std::cout << "associatedMem_else "<<std::endl;

								member += associatedMember(decodedCenter,modification,pos);
						//		std::cout<< "member "<< member << " modification "<<modification<<std::endl;

							}
							assert(pos < decodedCenter.size());
							size_t mod_length =  model.modification_length(modification);
							std::cout << "mod length "<< mod_length << std::endl;
							pos += model.modification_length(modification);
							std::cout<< "position on center: " << pos << std::endl;
							std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
							target = dec.get_target(total);	
						}
						if(al_end.at(0) <= target && target < al_end.at(1)){
							dec.decode(al_end.at(0), al_end.at(1));
							wrappers.decode(al_end.at(0), al_end.at(1),total);
							std::cout << " End of an al! "<<std::endl;
							std::cout << modify_count << std::endl;
							std::cout << "t " << target << " " <<al_end.at(0)<< " "<< al_end.at(1)<<std::endl;
						}
						else{
							std::cout<<"there is something wrong!"<<std::endl;
						}
						std::string temp;
						std::cout<< "length of the member: "<<member.length()<<std::endl;
						if(member.length()> 0){
							if(reverse_member == true || reverse_both == true){
								std::string rev_member;
								for(size_t rev = 0; rev< member.size(); rev++){
									char chr= dnastring::complement(member.at(member.size()-1-rev));	
									rev_member +=chr;
								}
								out<<rev_member;
							//	save << rev_member;
								std::cout<< "add to the 'out' 1"<<endl;								
								for(size_t mem =0; mem<rev_member.length();mem++){
									std::cout<<rev_member.at(mem);
								}
								std::cout << " " << std::endl;
								for(size_t sl =MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl>0; sl--){
									temp += rev_member.at(rev_member.length()-1-sl);
								}
								pattern = temp; 
								std::cout << temp << std::endl;
								base = dnastring::base_to_index(rev_member.at(rev_member.length()-1));
							}else{
								std::cout<< "add to the 'out' 2"<<endl;
								for(size_t mem =0; mem<member.length();mem++){
									std::cout<< member.at(mem);
								}
								out<<member;
							//	save << member;
								std::cout << " " << std::endl;
								for(size_t sl =MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl>0; sl--){
									temp += member.at(member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(member.at(member.length()-1));
							}
						}else{
							std::cout<< "add to the 'out' 3"<<endl;
							out << decodedCenter;
						//	save << decodedCenter;
							for(size_t sl=MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl > 0; sl--){
								char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
								temp +=center_base;
							}
							pattern = temp; 
							base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
						}

					//	target = dec.get_target(total);	
					}
					target = dec.get_target(total);	
					if(target<al_begin.at(0)){
						std::cout<< "back to sequence from an alignment"<<std::endl;
						std::cout << base <<std::endl;
						std::cout << "pattern size " << pattern.size() << std::endl;
					} 
				}
				if(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1)){
					dec.decode(seq_acc_end.at(0),seq_acc_end.at(1));
					wrappers.decode(seq_acc_end.at(0),seq_acc_end.at(1),total);
					std::cout << " End of a sequence! "<<std::endl;
		//		//	out << std::endl;
				//	save << sequence_counter+ 1;
					sequence_counter = sequence_counter + 1;
					std::cout << "seq counter "<< sequence_counter <<std::endl;
					out.write(reinterpret_cast<char*>(&sequence_counter),sizeof(uint32_t));

				}else{
					std::cout<<"there is something wrong here!"<<std::endl;
				}
				target = dec.get_target(total);	
			}	
			if(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1)){
				dec.decode(seq_acc_end.at(0), seq_acc_end.at(1));
				wrappers.decode(seq_acc_end.at(0), seq_acc_end.at(1),total);
				std::cout << " End of an accession! "<<std::endl;
			}else{
				std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
			}
			acc = acc + 1;
		}
	//	save.close();
	}


	template<typename T>
	void dynamic_decoder<T>::set_long_centers_from_stream(std::ifstream & in){
		size_t bit = 32;
		int index = 1;
		size_t partition_size = 0;
		char c; 
		char c1[sizeof(size_t)];
		in.read(c1,sizeof(size_t));
		size_t *temp = reinterpret_cast<size_t*> (c1);
		partition_size = *temp;
		std::cout << "par size: " << partition_size << std::endl;
		std::vector<std::pair<std::vector<int>, std::vector<unsigned int> > >value;
		for(size_t i = 0; i < partition_size; i++){
			long_centers_high_value.push_back(value);
		}
		std::cout << "long center partition size "<<long_centers_high_value.size()<<std::endl;
		unsigned char h;
		c = in.get();
		size_t i = 0;
		bool long_center_is_reached = false;
		while(c != 8){
			size_t low = 0;
			std::vector<std::pair<std::vector<int> , std::vector<unsigned int> > >longCenterHigh; //first: indices of a center second: its low and high value
			while(c != 5 && c != 7 && long_center_is_reached == false){//As logn as there is a short center
				std::vector<int> center;
				size_t acc;
				in.read(c1,sizeof(size_t));
				size_t *temp = reinterpret_cast<size_t*> (c1);
				acc = *temp;
				std::map<size_t , std::vector<int> >::iterator it = acc_of_long_center.find(acc);
				if(it == acc_of_long_center.end()){
					acc_of_long_center.insert(std::make_pair(acc, std::vector<int>()));
					it = acc_of_long_center.find(acc);
				}
				std::cout << "acc "<< acc <<std::endl;
				it->second.push_back(index);
				center.push_back(index);
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t i = 0; i < binary_high_value.size();i++){
					high_value += binary_high_value.at(i)*(1<<i);
				}
				std::vector<unsigned int> boundry(2,0);
				boundry.at(1)=high_value;
				std::cout << "high "<< high_value << std::endl;
				boundry.at(0)=low;
				std::cout << "low "<< low << std::endl;
				low = boundry.at(1);
				longCenterHigh.push_back(std::make_pair(center, boundry));
				index ++;
				c = in.get();
			}//Now there is no more a short center
			if(c == 5){
				long_center_is_reached =true;
				std::cout << " c is 5 " << index <<std::endl;
				c=in.get();
			}
			while(c != 7 && long_center_is_reached == true){
				//Add index
				char h1;
				h1 = in.get();
			//	size_t *temp = reinterpret_cast<size_t*> (h1);
			//	size_t center_length = *temp;
				std::cout << "size 0f the center "<< int(h1)<<std::endl;
				vector<int>center;
				for(size_t i =0; i < h1 ; i++){
					char c1[sizeof(int)];
					in.read(c1,sizeof(int));
					int *temp = reinterpret_cast<int*> (c1);	
					center.push_back(*temp);
				}
				for(size_t i =0; i < center.size();i++){
					std::cout << " index "<< center.at(i)<<std::endl;
				}
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t i = 0; i < binary_high_value.size();i++){
					high_value += binary_high_value.at(i)*(1<<i);
				}
				std::vector<unsigned int> boundry(2,0);
				boundry.at(1)=high_value;
				std::cout << "high "<< high_value << std::endl;
				boundry.at(0)=low;
				std::cout << "low "<< low << std::endl;
				low = boundry.at(1);
				longCenterHigh.push_back(std::make_pair(center, boundry));
				c = in.get();
				std::cout<<"c " << int(c) <<std::endl;
			}
			long_centers_high_value.at(i) = longCenterHigh;
			i = i + 1;
			std::cout <<  "i " << i <<std::endl;
			c = in.get();
			std::cout << "c " << int(c) <<std::endl;
		}
		std::cout << "c " << int(c) <<std::endl;
		for(size_t i =0; i < long_centers_high_value.size();i++){
			std::vector<bool> binary_high_value(0);
			size_t bound = bit/8;
			for(size_t j = 0 ; j < bound ; j++){ 
				h=in.get();
				for(size_t k = 0; k < 8 ; k++){
					binary_high_value.push_back(h%2);
					h = h/2;	
				}
			}
			unsigned int high_value = 0;	
			std::cout << binary_high_value.size() <<std::endl;
			for(size_t i = 0; i < binary_high_value.size();i++){
				high_value += binary_high_value.at(i)*(1<<i);
			}
			std::cout << high_value <<std::endl;
			long_center_partitions_high_value.push_back(high_value);
		}
	}
	template<typename T>
	void dynamic_decoder<T>::set_long_center_partition_high_from_stream(std::ifstream & in){
/*		size_t bit = 32;
		int index = 1;
		size_t partition_size = 0;
		char c; 
		char c1[sizeof(size_t)];
		in.read(c1,sizeof(size_t));
		size_t *temp = reinterpret_cast<size_t*> (c1);
		partition_size = *temp;
		std::cout << "par size: " << partition_size << std::endl;
		std::vector<std::pair<std::vector<int>, std::vector<unsigned int> > >value;
		for(size_t i = 0; i < partition_size; i++){
			long_centers_high_value.push_back(value);
		}
		std::cout << "long center partition size "<<long_centers_high_value.size()<<std::endl;
		unsigned char h;
		c = in.get();
		size_t i = 0;
		bool long_center_is_reached = false;
		while(c != 8){
			size_t low = 0;
			std::vector<std::pair<std::vector<int> , std::vector<unsigned int> > >longCenterHigh; //first: indices of a center second: its low and high value
			while(c != 5 && c != 7 && long_center_is_reached == false){//As logn as there is a short center
				std::vector<int> center;
				size_t acc;
				in.read(c1,sizeof(size_t));
				size_t *temp = reinterpret_cast<size_t*> (c1);
				acc = *temp;
				std::map<size_t , std::vector<int> >::iterator it = acc_of_long_center.find(acc);
				if(it == acc_of_long_center.end()){
					acc_of_long_center.insert(std::make_pair(acc, std::vector<int>()));
					it = acc_of_long_center.find(acc);
				}
				std::cout << "acc "<< acc <<std::endl;
				it->second.push_back(index);
				center.push_back(index);
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t i = 0; i < binary_high_value.size();i++){
					high_value += binary_high_value.at(i)*(1<<i);
				}
				std::vector<unsigned int> boundry(2,0);
				boundry.at(1)=high_value;
				std::cout << "high "<< high_value << std::endl;
				boundry.at(0)=low;
				std::cout << "low "<< low << std::endl;
				low = boundry.at(1);
				longCenterHigh.push_back(std::make_pair(center, boundry));
				index ++;
				c = in.get();
			}//Now there is no more a short center
			if(c == 5){
				long_center_is_reached =true;
				std::cout << " c is 5 " << index <<std::endl;
				c=in.get();
			}
			while(c != 7 && long_center_is_reached == true){
				//Add index
				char h1;
				h1 = in.get();
			//	size_t *temp = reinterpret_cast<size_t*> (h1);
			//	size_t center_length = *temp;
				std::cout << "size 0f the center "<< int(h1)<<std::endl;
				vector<int>center;
				for(size_t i =0; i < h1 ; i++){
					char c1[sizeof(int)];
					in.read(c1,sizeof(int));
					int *temp = reinterpret_cast<int*> (c1);	
					center.push_back(*temp);
				}
				for(size_t i =0; i < center.size();i++){
					std::cout << " index "<< center.at(i)<<std::endl;
				}
				std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				unsigned int high_value = 0;		
				for(size_t i = 0; i < binary_high_value.size();i++){
					high_value += binary_high_value.at(i)*(1<<i);
				}
				std::vector<unsigned int> boundry(2,0);
				boundry.at(1)=high_value;
				std::cout << "high "<< high_value << std::endl;
				boundry.at(0)=low;
				std::cout << "low "<< low << std::endl;
				low = boundry.at(1);
				longCenterHigh.push_back(std::make_pair(center, boundry));
				c = in.get();
				std::cout<<"c " << int(c) <<std::endl;
			}
			long_centers_high_value.at(i) = longCenterHigh;
			i = i + 1;
			std::cout <<  "i " << i <<std::endl;
			c = in.get();
			std::cout << "c " << int(c) <<std::endl;
		}
		std::cout << "c " << int(c) <<std::endl;
		for(size_t i =0; i < long_centers_high_value.size();i++){
			std::vector<bool> binary_high_value(0);
			size_t bound = bit/8;
			for(size_t j = 0 ; j < bound ; j++){ 
				h=in.get();
				for(size_t k = 0; k < 8 ; k++){
					binary_high_value.push_back(h%2);
					h = h/2;	
				}
			}
			unsigned int high_value = 0;	
			std::cout << binary_high_value.size() <<std::endl;
			for(size_t i = 0; i < binary_high_value.size();i++){
				high_value += binary_high_value.at(i)*(1<<i);
			}
			std::cout << high_value <<std::endl;
			long_center_partitions_high_value.push_back(high_value);
		}*/
	}
	template<typename T>
	void dynamic_decoder<T>::decoding_long_centers(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec){

		size_t total_length = 0;
		set_long_centers_from_stream(in); 
	//	set_long_center_partition_high_from_stream(in); 
		set_center_flags_from_the_stream(in);
		model.calculate_center_high_values(center_flags.at(0));
		std::cout << "all the parameters are set" <<std::endl;
		size_t acc = 0;
		dec.set_stream(in);	
		std::string first_pattern = first_sequence_pattern();
		unsigned int target;
		size_t lengthOfSeq = 0;
		unsigned int total = TOTAL_FOR_ENCODING;
		size_t numCenter = 0;
		size_t counter =0;
		for(size_t i = 0; i < long_centers_high_value.size(); i++){
			for(size_t j = 0; j < long_centers_high_value.at(i).size();j++){
				numCenter = numCenter +1 ;
			}
		}
		std::cout<< "partition size: "<<long_centers_high_value.size()<<std::endl;
		std::cout << " numbe of centers: "<< numCenter <<std::endl;
		while(counter < long_centers_high_value.size()){
			std::cout << "partition " << counter << " " << std::endl;
			size_t m = 0;			
			std::vector<int> centerOfAPartition;
			for(size_t j =0 ; j < long_centers_high_value.at(counter).size();j++){
				if(long_centers_high_value.at(counter).at(j).first.size()==1){
				//	std::cout << " low "<< long_centers_high_value.at(counter).at(j).second.at(0) << " high "<< long_centers_high_value.at(counter).at(j).second.at(1)<<std::endl;
					centerOfAPartition.push_back(long_centers_high_value.at(counter).at(j).first.at(0));
				}

			}
		//	for(std::map<std::vector<std::string>, std::vector<unsigned int> >::iterator it =long_centers_high_value.at(counter).begin(); it != long_centers_high_value.at(counter).end();it ++){
		//		if(it->first.size()==1){
		//			centerOfAPartition.push_back(it->first.at(0));
		//			std::cout<< "center " << it->first.at(0)<< std::endl;
		//		}		
		//	}
			std::cout<< "size of short centers of a partition: "<<centerOfAPartition.size()<<std::endl;
			std::map<int, std::string> intermediate;
			target = dec.get_target(total);
			std::cout<< "target "<<target << std::endl;
			size_t center_counter = 0;
			while(center_counter < centerOfAPartition.size()){
				std::cout << " m "<< m <<std::endl;
				int c_id = centerOfAPartition.at(m);
				std::cout << "c_id " << c_id << std::endl;
				for(std::map<size_t, std::vector<int> >::iterator it = acc_of_long_center.begin();it != acc_of_long_center.end(); it ++){
					for(size_t k = 0; k < it->second.size();k++){		
						if(it ->second.at(k) == c_id){
							acc = it->first;
							std::cout << "acc " << acc <<std::endl;
							break;
						}
					}
				}
		//		std::cout<< "acc in decoding is "<< acc<<std::endl;
				std::string center_seq;
				std::vector<unsigned int>high(6,0);
				std::vector<unsigned int>low(6,0);
				lengthOfSeq = 1;
				std::string pattern = first_pattern;
				std::cout<<"first pattern: " << first_pattern<<std::endl;
				std::map<std::string, std::vector<uint32_t> >::const_iterator it=model.get_center_model_highs(acc).find(first_pattern);
				assert(it != model.get_center_model_highs(acc).end());
				for(size_t j=0; j<5; j++){
					high.at(j) = it->second.at(j);
					if(j!= 0){
						low.at(j)=high.at(j-1);
					}else low.at(j)= 0;
				}
				low.at(5)= high.at(4);
				std::cout<< "high.at(4)" << high.at(4)<<std::endl;
				high.at(5)= center_flags.at(0)+center_flags.at(1);
				std::cout << "high.at(5) " << high.at(5) <<std::endl;
				size_t base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				center_seq +=dnastring::index_to_base(base);
				std::cout << "base "<< base << "l " << low.at(base) << " h " << high.at(base) << " target "<< target <<std::endl;
				dec.decode(low.at(base),high.at(base));
				wrappers.decode(low.at(base),high.at(base),total);
				std::cout<<"first base: "<<base<<std::endl;
			//	std::cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << std::endl;
				target = dec.get_target(total);
				while(target <  center_flags.at(0)){//end of a center	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					std::string current_pattern;
					std::stringstream s;
					if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL > 1){
						for(size_t M=1; M< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; M++){// ye doone function to model vase current pattern ham benevis!
							s<<pattern.at(M);
						}
						s<<p;
						s>>current_pattern;
					}
					if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL == 1){
					//	s<<pattern.at(0);
						s<<p;
						s>>current_pattern;
					}
				//	std::cout<<"current pattern: "<< current_pattern<<std::endl;
					std::map<std::string, std::vector<uint32_t> >::const_iterator it1=model.get_center_model_highs(acc).find(current_pattern);

					assert(it1 !=model.get_center_model_highs(acc).end());
					for(size_t j=0; j<5; j++){
						high.at(j) = it1->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}	
			//		high.at(4)=  model.get_powerOfTwo().at(bit);
			//		low.at(5)= high.at(4);
			//		high.at(5)= model.get_powerOfTwo().at(bit) + 5;
			//		low.at(6)= high.at(5);
			//		high.at(6) = model.get_powerOfTwo().at(bit) + 10;
					base = 10;
					for(size_t n = 0; n < 5 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					center_seq +=dnastring::index_to_base(base);
					dec.decode(low.at(base),high.at(base));
					wrappers.decode(low.at(base),high.at(base),total);
				//	std::cout << "l " << low.at(base) << " h " << high.at(base) <<std::endl;
				//	std::cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << std::endl;
				//	std::cout << "base: "<< base << std::endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				}
				std::cout << "length of seq " << lengthOfSeq << " center is " << center_seq << std::endl;
				base = 10;
				for(size_t n = 0; n < 6 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				assert(base == 5);
				std::cout<<"flag base : "<< base << std::endl;
				dec.decode(low.at(base),high.at(base));
				wrappers.decode(low.at(base),high.at(base),total);
				int center_id = centerOfAPartition.at(m);
				m = m +1;
				std::cout<< "id of cent in dec: " << center_id <<std::endl;
				intermediate.insert(std::make_pair(center_id,center_seq));
				target = dec.get_target(total);	
				std::cout<< "target: "<<target << std::endl;
				center_counter++;
			}
			//Decoding the flag which represents end of a partition
			if(center_counter ==  centerOfAPartition.size()){
				std::cout<< "end of partition target: "<<target << std::endl;					
				decoded_long_center_in_partition.push_back(intermediate);
			}
			counter = counter + 1 ;
		}
	///	for(size_t i =0; i < decoded_center_in_partition.size(); i++){
	//		for(std::map<std::string, string>::iterator it = decoded_center_in_partition.at(i).begin();it != decoded_center_in_partition.at(i).end();it++){
	//			std::string center = it ->first;
	//			std::cout<<"center size in decoded map is: "<<std::endl;				
	//			std::cout<< center.size() << std::endl;
	//		}
	//	}
	}


	template<typename T>
	void dynamic_decoder<T>::al_decode_with_long_center(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out){
		model.set_patterns(in);
		std::cout<< "sequence and alignemnts patterns are retrieved! "<<std::endl;
//		set_long_centers_from_stream(in);
		std::cout << "long centers are set! "<<std::endl;
		decoding_long_centers(in,dec);
		std::cout<< "centers are decoded! "<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		size_t sequence_counter = 0;
		out.write(reinterpret_cast<char*>(&sequence_counter),sizeof(uint32_t));
		unsigned int total = TOTAL_FOR_ENCODING;
//		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		std::string first_al_pattern = first_alignment_pattern();
		std::string first_pattern = first_sequence_pattern();
		std::string al_pattern;
		unsigned int target;
		size_t acc = 0;
		std::cout << acc_of_long_center.size() <<std::endl;
		while(acc<model.get_acc()){
			std::cout<< "accession "<< acc <<std::endl;
			std::vector<uint32_t> al_begin(2);
			std::vector<uint32_t> seq_acc_end(2);
			model.get_seq_flag(acc,al_begin,seq_acc_end);
			std::cout << "al begin low " << al_begin.at(0) << " al begin high " << al_begin.at(1) <<std::endl;
			std::cout << " seq end low " << seq_acc_end.at(0) << " seq end high " << seq_acc_end.at(1) << std::endl;
			target = dec.get_target(total); 
			while(target < seq_acc_end.at(0)){//end of an accession
				std::string pattern;
				size_t test_out = 0;
				if(target<al_begin.at(0)){//If there is no alignment on the first position
					std::cout<< "Sequence doesn't start with an alignment!" <<std::endl;
					std::vector<unsigned int>high(5,0);
					std::vector<unsigned int>low(5,0);
					pattern = first_pattern;
					std::cout<< "pattern size: "<<pattern.size()<<std::endl;
					std::map<std::string, std::vector<unsigned int> >::const_iterator it=model.get_sequence_model_highs(acc).find(first_pattern);
					for(size_t j=0; j<5; j++){
						high.at(j) = it->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}
					assert(high.at(4) == al_begin.at(0));
					base = 12;
					for(size_t n = 0; n < 5 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					assert(base < 5);
					char p = dnastring::index_to_base(base);
					out<<p;	
					test_out++;
					dec.decode(low.at(base),high.at(base));	
					wrappers.decode(low.at(base),high.at(base),total);
					int base_int = base;
					wrappers.decodeContext(base_int);
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}
				else {
					std::cout<<"Sequence starts with an alignment!"<<std::endl;
				}
				size_t counter = 0;
				while(target< seq_acc_end.at(0)){//end of each seq
					while(target< al_begin.at(0)){//If there is no alignment on that position
						char p = dnastring::index_to_base(base);
						std::stringstream s;
						std::string current_pattern;
						if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL > 1){
							for(size_t M=1; M< MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; M++){
								s<<pattern.at(M);
							}
							s<<p;
							current_pattern = s.str();
						}
						if(MAX_SEQUENCE_MARKOV_CHAIN_LEVEL ==1){	
							s<<p;
							s>>current_pattern;
						}
					//	std::cout<< current_pattern<<std::endl;
						std::vector<unsigned int>high(5,0);
						std::vector<unsigned int>low(5,0);
						std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_sequence_model_highs(acc).find(current_pattern);
						assert(it1 !=model.get_sequence_model_highs(acc).end());
						for(size_t j=0; j<5; j++){
							high.at(j) = it1->second.at(j);
							if(j!= 0){
								low.at(j)=high.at(j-1);
							}else low.at(j)= 0;
						}	
						assert(high.at(4)== al_begin.at(0));
						base = 12;
						for(size_t n = 0; n < 5 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						assert(base < 5);
						char p1 = dnastring::index_to_base(base);
						out<<p1;
						test_out++;
						counter = counter + 1;
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						int base_int = base;
						wrappers.decodeContext(base_int);
					//	std::cout << "base " << base << " low " << low.at(base) << " high " << high.at(base) << std::endl;
						pattern = current_pattern;
						target = dec.get_target(total);	
					//	std::cout << "tar1 "<<target <<std::endl;
						if(target>= al_begin.at(0)){
							std::cout<< "starting of an alignemnt's flag!"<<std::endl;
							std::cout<< test_out<<std::endl;
						}
					}
					std::cout<< "counter "<< counter << std::endl;
					counter = 0;
					std::cout << "target "<<target <<std::endl;
					if(al_begin.at(0)<=target && target < seq_acc_end.at(0)){//Beginning of an al but not end of the seq
						if( al_begin.at(0)<=target && target < al_begin.at(1)){
							dec.decode( al_begin.at(0), al_begin.at(1));
							wrappers.decode( al_begin.at(0), al_begin.at(1),total);
						}else{
							std::cout<< "there is a problem in bal flag"<<std::endl;
						}
						target = dec.get_target(total);
						std::vector<int> center; 
						size_t number_of_par = 0;
						std::cout<< "par target: "<< target <<std::endl;
						unsigned int low_par = 0;
						unsigned int high_par = 0;
						for(size_t j = 0; j < long_center_partitions_high_value.size(); j++){
							low_par = 0;
							high_par = long_center_partitions_high_value.at(j);
							if(j != 0){
								low_par = long_center_partitions_high_value.at(j-1);
							}
							if(target >= low_par && target < high_par){
								number_of_par = j;
								std::cout<< "number of par: "<< number_of_par << std::endl;
								std::cout<< "low par: "<< low_par << " high par"<< high_par <<std::endl;
								break;
							}
						}
						dec.decode(low_par,high_par);
						wrappers.decode(low_par,high_par,total);
						target = dec.get_target(total);
						std::cout << "tar2: "<< target << std::endl;
						for(size_t j =0; j < long_centers_high_value.at(number_of_par).size();j++){
							if(long_centers_high_value.at(number_of_par).at(j).second.at(0)<=target && long_centers_high_value.at(number_of_par).at(j).second.at(1)>target){
								std::cout << long_centers_high_value.at(number_of_par).at(j).second.at(0) << " " << long_centers_high_value.at(number_of_par).at(j).second.at(1) <<std::endl;
								center = long_centers_high_value.at(number_of_par).at(j).first;
								dec.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1));
								wrappers.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1),total);
								break;
							}
						} 
						//Decoding modification!	
						size_t cent_acc;
						std::string decodedCenter;
						std::cout<< "center size "<< center.size() << std::endl;
						for(size_t n =0; n < center.size();n++){
							std::cout << center.at(n)<< std::endl;
						}
						target = dec.get_target(total);
						for(std::map<size_t, std::vector<int> >::iterator it_acc = acc_of_long_center.begin(); it_acc != acc_of_long_center.end(); it_acc++){
							for(size_t g= 0; g< it_acc->second.size();g++){
								if(it_acc->second.at(g)== center.at(0)){
									std::cout << "center at 0 "<< center.at(0)<<std::endl;
									cent_acc = it_acc->first;
									std::cout<< "acc of cent: "<<cent_acc<<std::endl;
									break;
								}else continue;
							}
						}
						std::cout << "target before " << target <<std::endl;
						bool reverse_center = false;
						bool reverse_member = false;
						bool reverse_both = false;
						bool forward_both = false;
						decoding_directions(in, dec, reverse_center, reverse_member, reverse_both, forward_both, target);
						std::cout << "target after " << target <<std::endl;

						for(size_t m =0; m < center.size();m++){
							std::cout << "on center " << m <<std::endl;
							if(center.at(m) > 0){
								if(reverse_center == true){
									for(size_t j =0; j <decoded_long_center_in_partition.size() ;j++){
										std::map<int, std::string>::iterator seq = decoded_long_center_in_partition.at(j).find(center.at(center.size()-1-m));
										if(seq != decoded_long_center_in_partition.at(j).end()){
											std::cout << "decoded center size here is " << seq->second.size() << std::endl;
											decodedCenter += seq->second;											
											std::cout<< decodedCenter.size() << std::endl;
											break;
										}
									}


								}else{
									for(size_t j =0; j <decoded_long_center_in_partition.size() ;j++){
										std::map<int, std::string>::iterator seq = decoded_long_center_in_partition.at(j).find(center.at(m));
										if(seq != decoded_long_center_in_partition.at(j).end()){
											std::cout << "decoded center size " << seq->second.size() << std::endl;
											decodedCenter += seq->second;
											std::cout<< decodedCenter.size() << std::endl;
											break;
										}
									}
								}
							}else{//add reverse complement
								std::cout << "rev_comp is added !" <<std::endl;
								for(size_t j =0; j <decoded_long_center_in_partition.size() ;j++){
									std::map<int, std::string>::iterator seq = decoded_long_center_in_partition.at(j).find(-1*center.at(m));
									if(seq != decoded_long_center_in_partition.at(j).end()){
										std::cout << "decoded center size " << seq->second.size() << std::endl;
										for(size_t n = seq->second.size() ; n > 0 ; n--){
											std::cout << seq->second.at(n-1) << " " << dnastring::complement(seq->second.at(n-1)) <<std::endl;
											decodedCenter += dnastring::complement(seq->second.at(n-1));
										}
										break;
									}
								}
							}
						}
						std::cout<< decodedCenter.size()<< " " <<std::endl;
						std::cout << decodedCenter << std::endl;
						al_pattern = first_al_pattern;//decoding an alignment starts here!
						size_t pos = 0;
						size_t modify = 0;
						std::string member;
						std::vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
						std::vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
						size_t last_base = dnastring::base_to_index(decodedCenter.at(decodedCenter.size()-1));
						std::vector<uint32_t> al_end(2);
						model.get_end_al_flag(cent_acc, acc ,al_end);
						std::cout << "al end: " << al_end.at(0) << " "<<al_end.at(1) <<std::endl;			

					//	std::cout << "target before " << target <<std::endl;
					//	bool reverse_center = false;
					//	bool reverse_member = false;
					//	bool reverse_both = false;
					//	bool forward_both = false;
					//	decoding_directions(in, dec, reverse_center, reverse_member, reverse_both, forward_both, target);
						std::cout << "target after " << target <<std::endl;
						while(target < al_end.at(0) && pos < decodedCenter.size()){//while we are on an alignment
							std::vector<unsigned int>al_high(NUM_MODIFICATIONS,0);
							std::vector<unsigned int>al_low(NUM_MODIFICATIONS,0);
							size_t last_base;
							if(reverse_center == true){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							if(reverse_both == true){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							if(forward_both == true){//both references are forward
								std::cout << "both are forward" <<std::endl;
								last_base = dnastring::base_to_index(decodedCenter.at(pos));
							}
							if(reverse_member == true){
								last_base = dnastring::base_to_index(decodedCenter.at(pos));
								std::cout << "member is reversed! " << std::endl;
							}

							std::string al_context = al_pattern;
							al_context += last_base;
							std::cout<< "decoding_context: ";
							for(size_t k =0; k < al_context.size(); k++){
								std::cout<< int(al_context.at(k))<< " ";
								int con =  int(al_context.at(k));
								wrappers.decodeContext(con);
							}
							std::cout<< " " <<std::endl;
							
						std::map<std::string, std::vector<unsigned int> >::const_iterator mod = model.get_alignment_model_highs(cent_acc,acc).find(al_context);//modifiactions are considered from center to the other reference
							assert(mod != model.get_alignment_model_highs(cent_acc,acc).end());
							for(size_t k = 0 ; k < mod->second.size();k ++){
								al_high.at(k) = mod->second.at(k);
								if(k > 0){
									al_low.at(k) = al_high.at(k-1);
								}else{
									al_low.at(k) = 0;
								}				
							}
							assert(al_high.at(mod->second.size()-1)= al_end.at(0));
							size_t modification = NUM_KEEP+NUM_DELETE+20;
							for(size_t n = 0; n < NUM_KEEP+NUM_DELETE+10; n++){
								if(al_low.at(n)<=target && al_high.at(n)> target){
									modification = n;
									break;
								}
							}
							std::cout << "modification "<<modification << " low "<< al_low.at(modification)<< " high " <<al_high.at(modification)<<std::endl;
							dec.decode(al_low.at(modification),al_high.at(modification));
							wrappers.decode(al_low.at(modification),al_high.at(modification),total);
							modify = modification;
							std::cout<< " modify: " << modify << std::endl;
							stringstream s;
							if(MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL > 1){
								for(size_t M=1; M< al_context.size()-1; M++){
									std::cout << int(al_context.at(M)) <<std::endl;
									s << al_context.at(M);
								}
								al_pattern=s.str();
								al_pattern += modify;
							}else{
								al_pattern = modify;
							}
							std::cout<<"al_pattern: ";
							for(size_t al = 0; al < al_pattern.size();al++){
								std::cout<< int(al_pattern.at(al))<< " ";
							}
							std::cout<<" " << std::endl;
							if(reverse_center==true || reverse_both ==true){
								std::string rev_decodedCenter;
								for(size_t rev = 0; rev< decodedCenter.size(); rev++){
									char chr= dnastring::complement(decodedCenter.at(decodedCenter.size()-1-rev));	
									rev_decodedCenter +=chr;
								}
								std::cout << "associatedMem_both reverse "<<std::endl;
								member += associatedMember(rev_decodedCenter,modification,pos);							
							}else{
								std::cout << "associatedMember "<<std::endl;
								member += associatedMember(decodedCenter,modification,pos);
							//	std::cout<< "member "<< member << " modification "<<modification<<std::endl;
							}
							pos += model.modification_length(modification);
							assert(pos <= decodedCenter.length());
							std::cout<< "position on center: " << pos << std::endl;
							std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
							target = dec.get_target(total);	
							std::cout << "tar is "<< target << std::endl;
						}
						std::string temp;
						if(member.length()>0){
							if(reverse_member == true || reverse_both == true){
								std::string rev_member;
								for(size_t rev = 0; rev< member.size(); rev++){
									char chr= dnastring::complement(member.at(member.size()-1-rev));	
									rev_member +=chr;
								}
								out<<rev_member;
								test_out +=rev_member.size();
								std::cout<< "add to the 'out' 1"<<endl;								
								for(size_t mem =0; mem<rev_member.length();mem++){
									std::cout<<rev_member.at(mem);
								}
								std::cout << " " << std::endl;
								for(size_t sl =MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl>0; sl--){
									temp += rev_member.at(rev_member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(rev_member.at(rev_member.length()-1));
							}else{
								std::cout<< "add to the 'out' 2 " <<  " Target " << target <<endl;
								for(size_t mem =0; mem<member.length();mem++){
									std::cout<< member.at(mem);
								}
								out<<member;
								test_out += member.size();
								std::cout << " " << std::endl;
								for(size_t sl =MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl>0; sl--){
									temp += member.at(member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(member.at(member.length()-1));
							}
						}
						else{
							std::cout<< "add to the 'out' when center on th ref"<<endl;
							out << decodedCenter;
							test_out +=decodedCenter.size();
							for(size_t sl=MAX_SEQUENCE_MARKOV_CHAIN_LEVEL; sl > 0; sl--){
								char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
								temp +=center_base;
							}
							pattern = temp; 
							base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
						}	
						//end of an al flag:
						if(al_end.at(0) <= target && target < al_end.at(1)){
							dec.decode(al_end.at(0), al_end.at(1));
							wrappers.decode(al_end.at(0), al_end.at(1),total);
							std::cout << " End of an al! "<<std::endl;
							std::cout << "t " << target << " " <<al_end.at(0)<< " "<< al_end.at(1)<<std::endl;
						}
						else{
							std::cout<<"there is something wrong with end of al flag!"<<std::endl;
							std::cout << "decoded seq length "<< test_out << std::endl;							
							exit(1);
						}
	
						std::cout<< "back to seq target: "<< target << std::endl;
					}
					target = dec.get_target(total);	
					if(target<al_begin.at(0)){
						std::cout<< "back to sequence from an alignment target " << target <<std::endl;
					}
				}
				assert(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1));
				if(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1)){
					dec.decode(seq_acc_end.at(0),seq_acc_end.at(1));
					wrappers.decode(seq_acc_end.at(0),seq_acc_end.at(1),total);
					std::cout << " End of a sequence! "<<std::endl;
					std::cout << "decoded seq length "<< test_out << std::endl;
		//		//	out << std::endl;
				//	save << sequence_counter+ 1;
					sequence_counter ++;
					std::cout << "seq counter "<< sequence_counter <<std::endl;
					out.write(reinterpret_cast<char*>(&sequence_counter),sizeof(uint32_t));

				}else{
					std::cout<<"there is something wrong with end of seq flag!"<<std::endl;
					exit(1);
				}
				target = dec.get_target(total);		
				std::cout << "target: "<<target <<std::endl;
			}	
			assert(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1));
			if(seq_acc_end.at(0) <= target && target < seq_acc_end.at(1)){
				dec.decode(seq_acc_end.at(0), seq_acc_end.at(1));
				wrappers.decode(seq_acc_end.at(0), seq_acc_end.at(1),total);
				std::cout << " End of an accession! "<<std::endl;
			}else{
				std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
				exit(1);
			}
			acc++;
			std::cout << "start of acc " << acc <<std::endl;
		}//end of looping over number of accessions 

	}

#endif
