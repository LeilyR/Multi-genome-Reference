#include "encoder.hpp"

#ifndef ENCODER_CPP
#define ENCODER_CPP

	

//	template<typename T>
//	encoder<T>::encoder( all_data & d, const T & a_model,wrapper & wrap): data(d),model(a_model),wrappers(wrap),upper_bound(d.numSequences()),AlignmentsFromClustering(data.numSequences()){
//
//	}
//	template<typename T>
//	encoder<T>::~encoder(){}
	template<typename T>
	void encoder<T>::arithmetic_encoding_seq(std::ofstream & outs){
		size_t bit =12;
//		std::ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
//		if(outs.is_open()){
			dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
			enc -> set_stream(outs);
			for(size_t i = 0; i< data.numAcc(); i++){
			//	std::cout<< "accession: " << i << std::endl;
			//	std::cout<< "from encoding: "<< dnastring::base_to_index(data.getSequence(data.getAcc(i).at(0)).at(0))<<std::endl;
				for(size_t k = 0; k <data.getAcc(i).size(); k++){
					const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
					//std::cout<< "sequence id: "<< k << std::endl;
			//		std::cout<< "sequence length: "<< sequence.length() << std::endl;
					//std::cout << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << std::endl;
					for(size_t m = 0; m < sequence.length();m++){
						unsigned int l = 0;
					//	std::cout << m << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << std::endl;
					//	std::cout << m << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << std::endl;
						size_t base = dnastring::base_to_index(sequence.at(m));
						unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),m).at(base);
						if(base !=0){
							l = model.get_high_at_position(data.getAcc(i).at(k),m).at(base-1);
						}else l = 0;
						unsigned int t = model.get_powerOfTwo().at(bit) + 10;
					/*	if(i == 0){
							std::cout<< "base "<< base<<std::endl;
							for(size_t H=0;H<5;H++){
								std::cout<< "high at "<< H << " is "<< model.get_high_at_position(data.getAcc(i).at(k),m).at(H)<<std::endl;
							}
						}*/
					/*	if(i == 1){
							std::cout<<"context in encoding function "<<std::endl;
							std::string pattern = model.get_context(m, data.getAcc(i).at(k));
							for(size_t j = 0; j < pattern.length(); j++){
								std::cout<<pattern.at(j)<< std::endl;
							}
						}*/
					//	std::cout<<"base: "<<base<<std::endl;
					//	std::map<std::string, std::vector<unsigned int > >::const_iterator it1 = model.get_high(i).find(pattern);
					//	for(size_t j = 0; j < 5 ; j++){
					//		std::cout << " h " << it1 ->second.at(j) << std::endl;
					//	}
					//	std::cout<< "accession: "<< i<<" low: "<< l << " high: " << h << " total: "<< t <<" base: "<< base<<std::endl;
						if(base == 4){
							h=  model.get_powerOfTwo().at(bit);
						}
						enc -> encode(l,h,t);	
					/*	if(i == 1){
							std::cout<< "accession: "<< i<< " position "<< m <<" low: "<< l << " high: " << h << " total: "<< t <<" base: "<< base<<std::endl;
						}*/

					}

				//	size_t base = dnastring::base_to_index(sequence.at(sequence.length()-1));
					unsigned int l1	= model.get_powerOfTwo().at(bit);
					unsigned int h1 = model.get_powerOfTwo().at(bit) + 5;					
					unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
					enc -> encode(l1,h1,t1);  
				//	std::cout<< "accession: "<< i<<" low: "<< l1 << " high: " << h1 << " total: "<< t1 <<" base: "<< "5" <<std::endl;
				}
			//	const dnastring & sequence = data.getSequence(data.getAcc(i).at(data.getAcc(i).size()-1));
			//	size_t base = dnastring::base_to_index(sequence.at(sequence.length()-1));
				unsigned int l1 =  model.get_powerOfTwo().at(bit) + 5;
				unsigned int h1 = model.get_powerOfTwo().at(bit) + 10;
				unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
				enc -> encode(l1,h1,t1); 
			//	std::cout<< "accession: "<< i<<" low: "<< l1 << " high: " << h1 << " total: "<< t1 <<" base: "<<  "6" <<std::endl;
			}
			delete enc;
	//	}
	//	outs.close();
	}
//	template<typename T>
//	std::string encoder<T>::associatedMember(std::string & center, size_t & modificationPattern, size_t & position)
	template<typename T>
	void encoder<T>::calculating_clusters_high(std::map<std::string, unsigned int> & weight){
		size_t max_bit = 8;
		unsigned int h = 0;
		for(std::map<std::string, unsigned int>::iterator it = weight.begin(); it != weight.end(); it++){
			std::string center = it->first;
			std::map<std::string , std::vector<unsigned int> >::iterator it1 =cluster_high.find(center);
			if(it1 == cluster_high.end()){
				cluster_high.insert(std::make_pair(center,std::vector<unsigned int>(2,0)));
				it1 = cluster_high.find(center);
			}
			it1 ->second.at(0) = h;
			it1 -> second.at(1) = h + it->second;
			h = it1 -> second.at(1);
		}
	//	assert(h < model.get_powerOfTwo().at(31));
	//	for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high.begin(); it !=cluster_high.end(); it++){
	//		std::cout<< "low and high value for  cluster center "<< it->first <<" is "<< it->second.at(0) << " and "<<it ->second.at(1) <<std::endl;			
	//	}
		for(std::map<std::string, std::vector<unsigned int> >::iterator it=cluster_high.begin(); it != cluster_high.end();it++){
		//	std::cout<<"center0: " << it->first << std::endl;
		//	std::cout<< "low: "<< it->second.at(0)<<" high: " << it->second.at(1)<<std::endl;
		}
	}
//	template<typename T>
//	void encoder<T>::partition_centers(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::map<std::string, unsigned int> & weight)
//	template<typename T>
//	void encoder<T>::calculate_high_in_partition(std::map<std::string, unsigned int> & weight, std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters)
//	template<typename T>
//	void encoder<T>::partition_high()
//	template<typename T>
//	void encoder<T>::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters)
	template<typename T>
	const std::multimap<size_t, pw_alignment*> & encoder<T>::get_alignment(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, size_t seq_id){
		setOfAlignments(alignmentsOfClusters);
		return AlignmentsFromClustering.at(seq_id);
	}
	template<typename T>
	void encoder<T>::set_pattern_from_stream(std::ifstream & in){
	/*	char c;
		c=in.get();
		size_t sequence_id=0;
		while(c != 7){
			std::map<size_t,std::vector<std::string> >::iterator it =first_pattern_after_al.find(sequence_id);
			if(it == first_pattern_after_al.end()){
				first_pattern_after_al.insert(std::make_pair(sequence_id,std::vector<std::string>()));
			}
			while(c!=8){
				std::string pattern;
				std::stringstream s;
				while(c!=0){
					s<<c;
					c=in.get();
				}
				s>>pattern;
				std::map<size_t,std::vector<std::string> >::iterator it1 =first_pattern_after_al.find(sequence_id);
				it1->second.push_back(s.str());
				c=in.get();
			}
			sequence_id = sequence_id +1;
			c=in.get();
		}
		for(std::map<size_t,std::vector<std::string> >::iterator it1 = first_pattern_after_al.begin(); it1 != first_pattern_after_al.end(); it1 ++){
			for(size_t i =0; i< it1->second.size();i++){
				std::cout<< it1->second.at(i)<<std::endl;
			}
		}*/
	}
//	template<typename T>
//	void encoder<T>::add_partition_high_to_stream(std::ofstream & outs)
//	template<typename T>	
//	void encoder<T>::set_partition_high_from_stream(std::ifstream & in)
//	template<typename T>
//	void encoder<T>:: al_encoding(std::map<std::string, unsigned int> & weight, std::map<std::string, std::string > & cluster_members, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &outs,dlib::entropy_encoder_kernel_1 & enc )
//	template<typename T>
//	void encoder<T>::al_decoding(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec )
/*	void encoder::test_al_decoding(std::ifstream & in){
		dlib::entropy_decoder_kernel_1  dec;
		arithmetic_decoding_centers(in);
		std::cout<<"decoding the center has been done!"<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		std::string first_al_pattern = model.get_firstAlignmentPattern();
		std::string al_pattern;
		dec.std::set_stream(in);
		unsigned int target;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		std::string center = '2:13223';
		size_t part = 0;
		size_t flag = model.get_powerOfTwo().at(bit);	//8192
		size_t i = 1; 
		target = dec.get_target(total);
		std::cout << "target : "<< target <<std::endl;
		while(target < flag + 5){
			while(target < flag){
			size_t pos = 0;
			std::string center;
			size_t number_of_par = 0;
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
			dec.decode(low_par,high_par);
			target = dec.get_target(total);
			for(std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(number_of_par).begin(); it2!= cluster_high_partition.at(number_of_par).end();it2++){
				if(it2->second.at(0) <= target && it2->second.at(1) > target){
					center = it2->first;
					std::cout <<"center: "<< center<< std::endl;
					std::cout<< "center target: "<< target << "total" << total<<std::endl;
					dec.decode(it2->second.at(0),it2->second.at(1));
					std::cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<std::endl;
					break;
				}	
			}		
			std::map<std::string, string>::iterator seq = decoded_center_in_partition.at(0).find(center);
			std::string decodedCenter =seq->second;
			std::cout<< "decoded center: "<< decodedCenter <<std::endl;
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int cent_ref = atoi(center_parts.at(0).c_str());
			unsigned int cent_left = atoi(center_parts.at(1).c_str());
			std::cout<< "cent left: "<< cent_left << " cent_ref: "<<cent_ref << std::endl;
			size_t cent_acc;
			for(std::map<size_t, std::vector<std::string> >::iterator it_acc = acc_of_center.begin(); it_acc != acc_of_center.end(); it_acc++){
				for(size_t g= 0; g< it_acc->second.size();g++){
					if(it_acc->second.at(g)== center){
						cent_acc = it_acc->first;
						std::cout<< "acc of cent: "<<cent_acc<<std::endl;
						break;
					}else continue;
				}
			}
			al_pattern = first_al_pattern;
			target = dec.get_target(total);
			std::cout<<" target b_al: "<< target << std::endl;
			while(target < flag){//decoding modifications:
				std::vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
				std::vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
				size_t last_base = dnastring::base_to_index(decodedCenter.at(pos));
				std::string al_context = al_pattern;
				al_context += last_base;
				std::cout<< "size of context: "<< al_context.size() << std::endl;
				std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are always from center to the other sample
					assert(mod != model.get_highValue(cent_acc,i).end());
					for(size_t k = 0 ; k < mod->second.size();k ++){
						al_high.at(k) = mod->second.at(k);
						if(k > 0){
							al_low.at(k) = al_high.at(k-1);
						}else{
							al_low.at(k) = 0;
						}				
					}
					al_high.at(mod->second.size()-1)= model.get_powerOfTwo().at(bit);
					size_t modification = NUM_KEEP+NUM_DELETE+20;
					for(size_t n = 0; n < NUM_KEEP+NUM_DELETE+10; n++){
						if(al_low.at(n)<=target && al_high.at(n)> target){
							modification = n;
							break;
						}
					}
					std::cout<< "modification: " << modification << std::endl;
					dec.decode(al_low.at(modification),al_high.at(modification));
					al_pattern = modification;
					std::cout<<"al_mod_pattern: "<< al_pattern << std::endl;
					pos += model.modification_length(modification);
					std::cout<< "position on center: " << pos << std::endl;
					std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
					assert(pos < decodedCenter.size());
					target = dec.get_target(total);						
			}
				target= dec.get_target(total);
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit) +5;
				if(target >= l1 && target < h1){
					dec.decode(l1,h1);
					target = dec.get_target(total);
				}
			}
		}
		unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
		unsigned int h1 = model.get_powerOfTwo().at(bit) +10;
		dec.decode(l1,h1);
	}*/
//	template<typename T>	
//	void encoder<T>::add_center_to_stream(std::ofstream & outs)
//	void encoder::add_acc_to_stream(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters,std::ofstream & outs){
	template<typename T>
	void encoder<T>::add_acc_to_stream(std::ofstream & outs){	
	//	std::ofstream outs("encode",std::ofstream::binary);
			for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high.begin(); it != cluster_high.end(); it++){
			//	assert(it->second.size() != 0);
				std::string center = it ->first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int cent_ref = atoi(center_parts.at(0).c_str());
				unsigned int cent_left = atoi(center_parts.at(1).c_str());
				size_t acc_id = data.accNumber(cent_ref);
				outs << (char) 0;
				outs << acc_id;
				outs << (char) 0;
				outs << cent_ref;
				outs<< (char) 0;
				outs<< cent_left;
				outs << (char)0;
		//		std::map<std::string, std::vector<unsigned int> >::iterator it1 = cluster_high.find(center);
				std::vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
			//	std::cout<< "center1: "<< it->first<< "low : "<< it->second.at(0) << "high: "<< it->second.at(1)<<std::endl;	
				for(size_t i = 0 ; i < 32; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
			/*	for(size_t n = 0 ; n < bit_to_byte.size();n++){
					std::cout<< " "<< bit_to_byte.at(n);
				}
					std::cout<< "" << std::endl;*/
			//	std::cout<<"bit to byte: "<< bit_to_byte.size() <<std::endl;
				size_t count =0;
				for(size_t n = 0; n < bit_to_byte.size(); n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
						count = count +1;
					}
			//		std::cout<< "a: "<< int(a) << std::endl;					
					n= n+7;
					outs<< a;
				}
			//	std::cout << "count: "<< count << std::endl;
			}
			outs << (char) 8 ;
	}
	template<typename T>
	void encoder<T>::set_acc_from_stream(std::ifstream & in){
		size_t bit = 32;
		char c ; 
		unsigned char h;
		c = in.get();
		unsigned int low = 0;
		while(c != 8){
			size_t id;
			in >> id;
			std::map<size_t , std::vector<std::string> >::iterator it = acc_of_center.find(id);
			if(it == acc_of_center.end()){
				acc_of_center.insert(std::make_pair(id, std::vector<std::string>()));
			}
			c = in.get();
			unsigned int cent_ref;
			in >> cent_ref;
			c = in.get();
			unsigned int cent_left;
			in >> cent_left;
			std::stringstream center_id;
			std::string center;
			center_id << cent_ref << ":" << cent_left;
			std::map<size_t , std::vector<std::string> >::iterator it1 = acc_of_center.find(id);
			it1 ->second.push_back(center_id.str());
			c = in.get();
			std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high.find(center_id.str());
			if(it2 == cluster_high.end()){
				cluster_high.insert(std::make_pair(center,std::vector<unsigned int>(2,0)));
				it2 = cluster_high.find(center);	
			}
			std::vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
				//	std::cout<< "H: "<<int(h)<<std::endl;
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
			/*	for(size_t n = 0 ; n < binary_high_value.size();n++){
					std::cout<< " "<< binary_high_value.at(n);
				}
					std::cout<< "" << std::endl;*/
				unsigned int high_value = 0;		
			//	std::cout<<"binary high value size: "<<  binary_high_value.size()<<std::endl;
				for(size_t i = 0; i < binary_high_value.size();i++){
						high_value += binary_high_value.at(i)*model.get_powerOfTwo().at(i);
				}
		/*		for(size_t i = 0; i < binary_high_value.size();i++){
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(j-i);
					}
					i=i+bit-1;
				}*/
			//	std::cout << "high value: "<< high_value << std::endl;
				it2 -> second.at(1)=high_value;
				it2 ->second.at(0)=low;
				low= it2->second.at(1);
				c=in.get();
		}
	/*	for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high.begin(); it3 != cluster_high.end();it3++){
			std::cout<< "center_red_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
		}*/
	}
//	template<typename T>
//	void encoder<T>::set_center_from_stream(std::ifstream & in)
//	template<typename T>
//	void encoder<T>::arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream & outs, dlib::entropy_encoder_kernel_1 & enc)
	template<typename T>
	void encoder<T>::arithmetic_encoding_centId(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream & outs){
	/*	add_center_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		// TODO wrong!
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		for(size_t j =0 ; j < partition.size(); j++){
			for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(j).begin(); it != cluster_high_partition.at(j).end(); it++){
				std::string center = it->first;
				unsigned int l = it ->second.at(0);
				unsigned int h = it->second.at(1);
				unsigned int total =  model.get_powerOfTwo().at(19);
			//	std::cout<< "cent_id_enc"<<center<< " low: "<< it->second.at(0) << " high: "<< it->second.at(1)<< std::endl;		
				enc -> encode(l,h,total);
			}
		}
		delete enc;*/
	}
	template<typename T>
	void encoder<T>::arithmetic_decoding_centId(std::ifstream& in, dlib::entropy_decoder_kernel_1& dec){
	/*	model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);

		// TODO wrong!
		dec.set_stream(in);	
		unsigned int target;
		size_t numCenter = 0;
		size_t counter =0;
		std::string center;
		unsigned int t = model.get_powerOfTwo().at(19);
		for(size_t i = 0; i < cluster_high_partition.size(); i++){
			for(size_t j = 0; j < cluster_high_partition.at(i).size();j++){
				numCenter = numCenter +1 ;
			}
		}
	//	std::cout << " numbe of centers: "<< numCenter <<std::endl;
	//	while(counter < numCenter){
			for(size_t i = 0; i < cluster_high_partition.size(); i++){
				for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end();it++){
					target = dec.get_target(t);
					if(it->second.at(0) <= target && it->second.at(1) > target){
						center = it ->first;
						dec.decode(it->second.at(0), it->second.at(1));
					//	centerId.push_back(center);
					}	
				}
			}
	//	}*/
	}
	template<typename T>
	const std::map<std::string, std::vector<unsigned int> >& encoder<T>:: get_center_high() const{
		return cluster_high;
	}
//	template<typename T>
//	void encoder<T>::arithmetic_decoding_centers(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec)
	template<typename T>
	void encoder<T>::encoding_seq_test(){
	/*	unsigned int low = 0;
		unsigned int high = 0;
		unsigned int total = 0;
		std::vector<double> num(5,0);
		std::string seq = "";
	//	std::string seq1 = "ATAATAAA";
	//	for(size_t n= 0; n< 1000; n++){
	//		seq+=seq1;	
	//	}
		const dnastring seq = data.getSequence(15);
		std::cout<< "seq length: "<< seq.length()<<std::endl;
		for(size_t i = 0; i < seq.length();i++){
			size_t base = dnastring::base_to_index(seq.at(i));
			num.at(base)++;	
		}
		for(size_t j = 0; j < 5; j++){
			total += num.at(j);	
		}
		for(size_t k=0; k<5; k++){
			lower_bound.at(k)=low;
			high += num.at(k);
			low = high;
			upper_bound.at(k)=high;
		}
		std::ofstream outs("test", std::ofstream::binary);
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1() ;
		enc->std::set_stream(outs);
		for(size_t i = 0; i < seq.length(); i++){
			size_t base = dnastring::base_to_index(seq.at(i));	
			unsigned int l = lower_bound.at(base);
			unsigned int h = upper_bound.at(base);
			unsigned int t = total;
			enc->encode(l, h, t);
		}
		for(size_t i=0; i<5; i++) {
			std::cout<< "low " << lower_bound.at(i) << " | high " << upper_bound.at(i) << "  " << total <<std::endl;
			std::cout<<"ratio: "<< (upper_bound.at(i)-lower_bound.at(i))/total <<std::endl;

		}
		std::cout<< "sequence: "<<std::endl;
		for(size_t k =0; k < 50; k++){
			std::cout<< seq.at(k)<<std::endl;
		}	
		delete enc;
		decoding_seq(total);*/		
	}
	template<typename T>
	void encoder<T>::arithmetic_decoding_seq_test(int t){

		dlib::entropy_decoder_kernel_1  dec;
		std::ifstream in("test", std::ifstream::binary);
	/*	dec.std::set_stream(in);

		for(size_t i=0; i<5; i++) {
			size_t target = dec.get_target(t);

			size_t base = 18;
			for(size_t j=0; j<5; j++) {
				if(lower_bound.at(j)<= target && upper_bound.at(j) >= target) {
					base = j;
					break;
				}
			}
			dec.decode(lower_bound.at(base), upper_bound.at(base));
			std::cout << " at " << i << " see " << base << std::endl;

		}*/

	}
	template<typename T>
	void encoder<T>::arithmetic_decoding_seq(){
		std::ifstream in("encode",std::ifstream::binary);
		size_t bit =12;
		model.set_patterns(in);//returns high values we saved in the stream is called "in"
		model.set_alignment_pattern(in);
		set_acc_from_stream(in);//retrieves accession id & high of each cluster center.
		dlib::entropy_decoder_kernel_1  dec;
		// TODO wrong!
		dec.set_stream(in);	
		std::string first_pattern = model.get_firstPattern();
		unsigned int target;
		size_t i = 0;
		size_t lengthOfSeq = 0;
		size_t flag = model.get_powerOfTwo().at(bit);
		unsigned int total = model.get_powerOfTwo().at(bit) + 10;
//		for(size_t i = 0; i < data.numOfAcc(); i++)	
		while(i<data.numOfAcc()){
		//	target = 0;
			target = dec.get_target(total);
			while(target < flag+5){
				std::vector<unsigned int>high(7,0);
				std::vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				std::string pattern = first_pattern;
				std::map<std::string, std::vector<unsigned int> >::const_iterator it=model.get_high(i).find(first_pattern);
			//	std::cout<< "first pattern: " << first_pattern << std::endl;
				for(size_t j=0; j<5; j++){
					high.at(j) = it->second.at(j);
					if(j!= 0){
						low.at(j)=high.at(j-1);
					}else low.at(j)= 0;
				}
				high.at(4)=  model.get_powerOfTwo().at(bit);
				low.at(5)= high.at(4);
				//low.at(5)= it->second.at(4);
				high.at(5)= model.get_powerOfTwo().at(bit) + 5;
				low.at(6)= high.at(5);
				high.at(6) = model.get_powerOfTwo().at(bit) + 10;
			//	target = dec.get_target(total);
			//	std::cout<< "target1: "<< target <<std::endl;
				size_t base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				dec.decode(low.at(base),high.at(base));
				std::cout<<"first base: "<<base<<std::endl;
			//	std::cout<< "base1: " << base << " " << low.at(base) << " " << high.at(base) << std::endl;
		/*		for(size_t H =0 ; H < 7 ; H++){
					std::cout << "low at " << H  << " is " << low.at(H) <<std::endl;
				}*/
			//	while(1)
				target = dec.get_target(total);
				while(target < flag){	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					std::string current_pattern;
					std::stringstream s;
					for(size_t M=1; M< Sequence_level; M++){// ye doone function to model vase current pattern ham benevis!
						s<<pattern.at(Sequence_level-M);
					}
					s<<p;
					s>>current_pattern;
			//		std::cout<<"current pattern: "<< current_pattern<<std::endl;
					std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
					assert(it1 != model.get_high(i).end());
					for(size_t j=0; j<5; j++){
						high.at(j) = it1->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}	
				/*	std::cout<<"high values of current pattern in decoding: "<<std::endl;
					for(size_t k = 0 ; k < 5; k++){
						std::cout<< it1 -> second.at(k)<<std::endl;
					}*/
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
				//	low.at(5)= it1->second.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + 5;
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + 10;
				//	total = model.get_powerOfTwo().at(bit) + 10;
				//	target = dec.get_target(total);
				//	std::cout<< "target: "<< target <<std::endl;
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						std::cout << "low at " << H  << " is " << low.at(H) <<std::endl;
					}*/
					dec.decode(low.at(base),high.at(base));
				//	std::cout<< "base: " << base <<" " << low.at(base) << " " << high.at(base)  << std::endl;
				//	std::cout << "base: "<< base << std::endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				/*	if(base >= 5){
						break;
					}*/
				}

				std::cout << "length of seq" << lengthOfSeq << std::endl;
			//	std::cout<<"here1! "<<std::endl;
			//	std::cout << "target 2 : "<< target <<std::endl;	
				base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						std::cout << "low at " << H  << " is " << low.at(H) <<std::endl;
					}*/
				dec.decode(low.at(base),high.at(base));	
				target = dec.get_target(total);	
			//	std::cout<< "target3: "<< target <<std::endl;
				if(target >= flag){
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						std::cout << "low at " << H  << " is " << low.at(H) <<std::endl;
					}*/
					dec.decode(low.at(base),high.at(base));
			//		std::cout<< "base3: " << base <<" " << low.at(base) << " " << high.at(base)  << std::endl;
				} else continue;
			}	
			std::cout<<"here!"<<std::endl;
			i = i+1;
		}
	}
	template<typename T>
	void encoder<T>::write_to_stream( std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster ,std::ofstream & outs){		
		add_center_to_stream(outs);
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
	}
	template<typename T>
	void encoder<T>::read_from_stream(std::ifstream & in){
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
	}
	template<typename T>
	void encoder<T>::size_calculator( size_t & file_size,std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::map<std::string, std::string> & cluster_members ){
		for(size_t i = 0; i< data.numAcc(); i++){
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
			//+creating cost of center( + modification cost of its correspondence associated member)
						pw_alignment *p = it->second;
						double c1,c2,m1,m2;
						model.cost_function(*p, c1, c2, m1, m2);
						size_t left_1,left_2,right_1,right_2; 
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);		
						std::stringstream mem;
						mem << sequenceId << ":" << n;
						std::map<std::string, std::string>::iterator cl = cluster_members.find(mem.str());
						assert(cl != cluster_members.end());
						std::string center = cl->second;
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(0).c_str());
						unsigned int cent_left = atoi(center_parts.at(1).c_str());
						if(p->getreference1()== sequenceId && cent_ref == sequenceId && cent_left == n){
							file_size += c1;
							n = right_1;
						}
						else if(p->getreference2()==sequenceId && cent_ref == sequenceId && cent_left == n){
							file_size += c2;
							n = right_2;
						}
						else{
							if(p->getreference1()==sequenceId){
								file_size += c2 + m2;
								n = right_1;
							}else{
								file_size +=c1+m1;
								n = right_2;
							}
						}
				
					}else{//If there is no alignment in that position
			//+creating cost of sequnce.at(n)
						size_t base = dnastring::base_to_index(sequence.at(n));
						size_t creating_cost = model.get_create_cost(i).at(base);
						file_size +=creating_cost;
					}
				}
			}
		}

	}


#endif
