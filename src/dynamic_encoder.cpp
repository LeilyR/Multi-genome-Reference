#include "dynamic_encoder.hpp"

#ifndef DYNAMIC_ENCODER
#define DYNAMIC_ENCODER


	template<typename T>
	dynamic_encoder<T>::dynamic_encoder(all_data & d, T & a_model, wrapper & wrap):data(d),model(a_model),wrappers(wrap),AlignmentsFromClustering(data.numSequences()){}

	template<typename T>
	dynamic_encoder<T>::~dynamic_encoder(){}

	template<typename T>
	void dynamic_encoder<T>::partition_centers(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::map<std::string, unsigned int> & weight){
		size_t sum_of_weight = 0;
		for(std::map<std::string, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){
			sum_of_weight = sum_of_weight + it->second ;
		}
		std::cout << "total weight "<< sum_of_weight <<std::endl;
		size_t numberOfPartitions = (sum_of_weight/TOTAL_FOR_ENCODING) + 1;
			std::cout<< "no. of partition:" << numberOfPartitions << std::endl;
		//	std::cout<< "size of al of cluster" << alignmentsOfClusters.size() << std::endl;
		std::vector<std::string> centers;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it= alignmentsOfClusters.begin(); it!= alignmentsOfClusters.end();it++){
			std::string c_id = it->first;
			centers.push_back(c_id);
		}
		size_t number_of_center = 0;
		for(size_t j = 0 ; j < numberOfPartitions ; j ++){	
			partition.insert(std::make_pair(j,std::vector<std::string>( )));
			size_t center_weight = 0;
			for(size_t i = number_of_center; i < centers.size();i++){
				std::map<std::string, unsigned int >::iterator it1 = weight.find(centers.at(i));
				center_weight +=  it1->second;
				if(center_weight <= TOTAL_FOR_ENCODING){
					std::map<size_t , std::vector<std::string> >::iterator it = partition.find(j);
					it->second.push_back(centers.at(i));
				}else{
					std::cout << "center weight " <<center_weight <<std::endl;
					number_of_center = i;
					break;
				}

			}
		/*	for(size_t i =j*alignmentsOfClusters.size()/numberOfPartitions; i < (j+1)*alignmentsOfClusters.size()/numberOfPartitions; i++){	
				std::map<size_t , std::vector<std::string> >::iterator it = partition.find(j);
				it->second.push_back(center.at(i));
			} */
		}
	}
	template<typename T>
	void dynamic_encoder<T>::calculate_high_in_partition(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::map<std::string, unsigned int> & weight){
		partition_centers(alignmentsOfClusters, weight);
		std::cout<<"partition size: "<< partition.size() <<std::endl;
		for(size_t i = 0; i < partition.size(); i ++){
			unsigned int low = 0;
			std::map<std::string, std::vector<unsigned int> > high;
			std::map<size_t , std::vector<std::string> >::iterator it = partition.find(i);
			assert(it != partition.end());
			for(size_t j = 0; j <it->second.size(); j++){
				std::string center = it->second.at(j);
				std::map<std::string, unsigned int>::iterator it2 = weight.find(center);
				assert(it2 != weight.end());
				std::map<std::string , std::vector<unsigned int> >::iterator it1 =high.find(center);
				if(it1 == high.end()){
					high.insert(std::make_pair(center,std::vector<unsigned int>(2,0)));
					it1 = high.find(center);
				}
				it1 ->second.at(0) = low;
				it1 -> second.at(1) = low + it2->second;
				low = it1 -> second.at(1);
				if(j == it->second.size()-1){
				//	std::cout << "h "<< h << "it2 second "<<it2->second << std::endl;
					assert(it1->second.at(1) <= TOTAL_FOR_ENCODING);
					it1->second.at(1) = TOTAL_FOR_ENCODING;
				}
			}
			cluster_high_partition.push_back(high);	
		}
	/*	for(size_t i =0; i < partition.size(); i++){
			for(std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).begin(); it2 !=cluster_high_partition.at(i).end(); it2++){
				std::cout<<"center "<< "in partition " << i << " is "<< it2->first << "high: "<< it2->second.at(1)<<std::endl;
			}
		}*/
		partition_high();
	}
	template<typename T>
	void dynamic_encoder<T>::partition_high(){
		unsigned int high = 0;
		for(size_t i = 0; i < partition.size(); i++){
			high = high + TOTAL_FOR_ENCODING/partition.size();
			partitionHigh.push_back(high);
		}
		partitionHigh.at(partitionHigh.size()-1) = TOTAL_FOR_ENCODING;
	/*	for(size_t i =0; i < partition.size();i++){
			std::cout << partitionHigh.at(i)<<std::endl;
		}*/
	}
	template<typename T>
	void dynamic_encoder<T>::add_center_to_the_stream(std::ofstream & outs){
		size_t bit =32;
		size_t sizeOfPartition = partition.size();
		outs.write(reinterpret_cast<char*>(&sizeOfPartition), sizeof(size_t));
		std::cout << "partition size in encoding " << partition.size() << std::endl;
		for(size_t i = 0 ; i < partition.size(); i ++){
		//	size_t index =0;
			size_t id = 0;
			for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				std::cout<< "id " << id << " center_write: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
				id++;
			}
			for ( std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end(); it ++ ){
				std::string center = it -> first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int cent_ref = atoi(center_parts.at(1).c_str());
				unsigned int cent_left = atoi(center_parts.at(2).c_str());
				size_t acc_id = data.accNumber(cent_ref);
				outs.put((char)0);
				outs.write(reinterpret_cast<char*>(&acc_id),sizeof(size_t));
				std::vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
				for(size_t i = 0 ; i < bit; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
			//	std::cout << "bit to byte "<<std::endl;
				for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						unsigned int position = 1<<(m-n);
						a+= position* bit_to_byte.at(m);
				//		std::cout << bit_to_byte.at(m)<<std::endl;
					}
					n= n+7;
				//	std::cout << "a is " << std::endl;
				//	std::cout << int(a) <<std::endl;
					outs.put(a);
				}
			//	index +=1;
			}
			outs.put((char) 7);
		}
		outs.put((char)8);
	}
	template<typename T>
	void dynamic_encoder<T>::add_partition_high_to_the_stream(std::ofstream & outs){
		for(size_t i = 0; i < partitionHigh.size(); i++){
			outs.put((char)0);
			std::vector<bool> bit_to_byte(0);
			unsigned int high = partitionHigh.at(i);
			std::cout<< "high_write: "<< high << std::endl;
			for(size_t j = 0 ; j < 32; j++){
				bit_to_byte.push_back(high%2);
				high = high/2;
			}
			for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
				unsigned char a = 0;//we can save it as char and read it as char but use it later on as size_t 'cus it is never higher than 2^8
				for(size_t m = n; m <n+8; m++){
					unsigned int position = 1<<(m-n);
					a+= position* bit_to_byte.at(m);
				}
				n= n+7;
				outs.put(a);
			}
			outs.put((char)7);
		}
		outs.put((char)8);
	}
	template<typename T>
	void dynamic_encoder<T>::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				std::string center = it->first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());	
				for(size_t j = 0; j < it->second.size();j++){
					pw_alignment * p = & it->second.at(j);
					size_t left1; 
					size_t left2;
					size_t right1;
					size_t right2;
					p->get_lr1(left1,right1);
					p->get_lr2(left2,right2);
					if(i != center_ref){
						if(p->getreference1()== i && p->getreference2()== center_ref&&left2==center_left){
							std::multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left1,p));
							}
						}else if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
							std::multimap<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
							if(it2 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left2,p));
							}
						}
					}
					if(i == center_ref){
						if(p->getreference1()== i && p->getreference2()== i){
							std::multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(left1);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left1,p));
							}
							std::multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left2);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left2,p));
							}
						}else{
							std::multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(center_left);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(center_left,p));
							}
						}
					}
				}
			}
		}
	}
	template<typename T>
	void dynamic_encoder<T>::add_center_flags_to_the_stream(std::vector<uint32_t> & bits , unsigned char & model_index, std::ofstream &outs){
		outs.put(model_index);
		for(size_t i =0; i < bits.size();i++){	
			size_t bit = bits.at(i);
			outs.write(reinterpret_cast<char*>(&bit),sizeof(uint32_t));
		}
	}

	template<typename T>
	void dynamic_encoder<T>::calculate_end_of_center_flags(std::map<std::string, std::vector<pw_alignment> > & alignmentsOfCluster,std::vector<uint32_t> & widths, std::ofstream & outs){
		vector<size_t> counts(2,0);
		size_t length = 0;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfCluster.begin(); it != alignmentsOfCluster.end(); it++){
			std::string center = it->first;
			std::vector<std::string> split;
			strsep(center, ":" , split);
			size_t cent_ref = atoi(split.at(1).c_str());
			size_t cent_left = atoi(split.at(2).c_str());
			pw_alignment * p = & it->second.at(0);
			size_t left1;
			size_t right1;
			size_t left2;
			size_t right2;
			p->get_lr1(left1,right1);
			p->get_lr2(left2,right2);
			if(cent_ref == p->getreference1()&& cent_left == left1){
				length += right1 - left1 + 1;
			}else{
				length += right2 - left2 + 1;
			}
		}
		counts.at(0) =  length + alignmentsOfCluster.size();
		counts.at(1) = alignmentsOfCluster.size();
		unsigned char model_index;
		std::vector<uint32_t> bits;
		uint32_t target_total = TOTAL_FOR_ENCODING;
		model.calculate_center_flags(counts,target_total, model_index , bits, widths);
		std::cout << "model index "<< int (model_index) << std::endl;
		std::cout << "bits " << bits.size()<<std::endl;
		for(size_t i = 0; i < bits.size(); i++){
			std::cout << bits.at(i) << std::endl;
		}
		//ADD mode_index & bits to the outs stream.
		add_center_flags_to_the_stream(bits,model_index,outs);
		std::cout << "widths "<<std::endl;
		for(size_t i =0;i < widths.size();i++){
			std::cout << widths.at(i) << std::endl;
		}
		std::cout << "target total "<< target_total << std::endl;
		assert(target_total < TOTAL_FOR_ENCODING);
		model.calculate_center_high_values(target_total);//recalculate the high value of bases on a center based on the valueof target_total
	}
	template<typename T>
	void dynamic_encoder<T>::arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &outs,dlib::entropy_encoder_kernel_1 & enc ){
		std::vector<uint32_t> widths(2,0);
		calculate_end_of_center_flags(alignmentOfCluster, widths,outs);
		size_t counter = 0;
		size_t index = 0;
		unsigned int total = TOTAL_FOR_ENCODING; //It includes all the flags. 
		enc.set_stream(outs);
		for(size_t i = 0 ; i < cluster_high_partition.size(); i++){
			for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin();it != cluster_high_partition.at(i).end(); it++){
			//	std::cout << " low " << it->second.at(0) << " high " << it->second.at(1) <<std::endl;
				std::string center = it ->first;
				std::cout << "center " << center << std::endl;
				std::vector<std::string> split;
				strsep(center, ":" , split);
				size_t cent_right = 0;
				size_t cent_ref = atoi(split.at(1).c_str());
				size_t cent_left = atoi(split.at(2).c_str());
				std::map<std::string,std::vector<pw_alignment> >::iterator it1 = alignmentOfCluster.find(center);
				assert(it1 != alignmentOfCluster.end());
				pw_alignment * p = & it1->second.at(0);
				const dnastring & seq = data.getSequence(cent_ref);
				size_t left1;
				size_t right1;
				size_t left2;
				size_t right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
			//	p->print();
				double gain = model.get_the_gain(*p, center);
				std::cout << "gain " << gain <<std::endl;
				std::cout<< "id of cent in enc: " << center <<std::endl;
				if(cent_ref == p->getreference1()&& cent_left == left1){
					if(p->getbegin1() < p->getend1()){//Encode the sequence itself
						std::cout<< "id of cent in enc1: " << center <<std::endl;
						cent_right = p->getend1();
					}else{
						cent_right = p->getbegin1();
					}				
					for(size_t k = cent_left; k < cent_right+1;k++){
						size_t base = dnastring::base_to_index(seq.at(k));				
						unsigned int l = 0;
						unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
						assert(h <= widths.at(0));
						if(base !=0){
							l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
						}else l = 0;
						if(base == 4){
							std::cout << "largest high is "<< h <<std::endl;
						//	h=  model.get_powerOfTwo().at(bit); 
						}
						counter++;
					//	std::cout << "base " <<base<<std::endl;
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
					}
					std::cout << "length " << cent_right+1 - cent_left <<std::endl;
					std::cout << "index " << index << std::endl;
					index++;
				}
				if(cent_ref == p->getreference2()&& cent_left == left2){
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
					}else{
						cent_right = p->getbegin2();
					}
					std::cout<< "id of cent in enc3: " << center  << " its left " << left << " its right " << right <<std::endl;
					for(size_t k = cent_left; k < cent_right+1;k++){
						size_t base = dnastring::base_to_index(seq.at(k));				
						unsigned int l = 0;
						unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
						assert(h <= widths.at(0));
						if(base !=0){
							l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
						}else l = 0;
						if(base == 4){
							std::cout << "largest high is "<< h <<std::endl;
						//	h=  model.get_powerOfTwo().at(bit);
						}
						counter++;
					//	std::cout << "base " <<base<<std::endl;
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
					}
					std::cout<<"cent length: "<< cent_right - cent_left+1 << std::endl;
					std::cout << "index " << index << std::endl;
					index++;
				}
				unsigned int l1	= widths.at(0);
				unsigned int h1 = widths.at(0)+widths.at(1);					
				enc.encode(l1,h1,total); 
				wrappers.encode(l1,h1,total);
				std::cout << "end of a center "<<std::endl;
			}
		//	unsigned int l1	=widths.at(0)+widths.at(1);
		//	unsigned int h1 =widths.at(0)+widths.at(1)+widths.at(2);
		//	enc.encode(l1,h1,total); 
		//	wrappers.encode(l1,h1,total);
			std::cout<< "end of partition: " << i <<std::endl;
			std::cout << "counter "<<counter<<std::endl;
		}
	}
	template<typename T>
	void dynamic_encoder<T>::encoding_directions(unsigned int & cent_ref, unsigned int & cent_left, const pw_alignment & p, std::ofstream & outs, dlib::entropy_encoder_kernel_1 & enc){
		//Reverse flags: //One fourth of TOTAL_FOR_ENCODING for each of them
		std::cout << "direction flag ";
		unsigned int total = TOTAL_FOR_ENCODING;
		size_t left_1; 
		size_t left_2;
		size_t right_1;
		size_t right_2;
		p.get_lr1(left_1,right_1);
		p.get_lr2(left_2,right_2);

		if((cent_ref == p.getreference1() && cent_left == p.getend1()&& left_2==p.getbegin2())||(cent_ref == p.getreference2()&& cent_left ==p.getend2()&&left_1 == p.getbegin1())){
			std::cout << "rev_cent "<<std::endl;
			//Center is reverse
			unsigned int l1 =  0;
			unsigned int h1 = TOTAL_FOR_ENCODING/4;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getbegin1()&& left_2==p.getend2())||(cent_ref == p.getreference2()&& cent_left ==p.getbegin2()&&left_1 == p.getend1())){
			std::cout << "rev_other "<<std::endl;
			//Other ref is reverse
			unsigned int l1 = TOTAL_FOR_ENCODING/4;
			unsigned int h1 = TOTAL_FOR_ENCODING/2;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getend1()&& left_2==p.getend2())||(cent_ref == p.getreference2() && cent_left ==p.getend2()&&left_1 == p.getend1())){
			std::cout << "rev_both "<<std::endl;
			//Both are reverse
			unsigned int l1 = TOTAL_FOR_ENCODING/2;
			unsigned int h1 = 3*TOTAL_FOR_ENCODING/4;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getbegin1()&& left_2==p.getbegin2())||(cent_ref == p.getreference2()&& cent_left ==p.getbegin2()&&left_1 == p.getbegin1())){
			std::cout << "forward_both "<<std::endl;
			//Both are forwards
			unsigned int l1 = 3*TOTAL_FOR_ENCODING/4;
			unsigned int h1 = TOTAL_FOR_ENCODING;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
	}
	template<typename T>
	void dynamic_encoder<T>::al_encode(std::map<std::string, unsigned int> & weight, std::map<std::string, std::string > & cluster_members, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &outs,dlib::entropy_encoder_kernel_1 & enc ){
		calculate_high_in_partition(alignmentOfCluster,weight);
		add_center_to_the_stream(outs);
		add_partition_high_to_the_stream(outs);
		setOfAlignments(alignmentOfCluster);
		arithmetic_encoding_centers(alignmentOfCluster,outs,enc);
		unsigned int total = TOTAL_FOR_ENCODING; //It includes all the flags. 
		std::vector<uint32_t> al_begin;
		std::vector<uint32_t> al_end;
		std::vector<uint32_t> seq_acc_end;
		std::cout << "total members "<<cluster_members.size()<<std::endl;
		for(size_t i = 0; i< data.numAcc(); i++){
			model.get_seq_flag(i , al_begin , seq_acc_end);
			std::cout << "al begin flag low "<< al_begin.at(0) << " al begin flag high " << al_begin.at(1) <<std::endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
				size_t length_between_two_centers = 0;				
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
						std::cout << "length between two centers " << length_between_two_centers <<std::endl;
						length_between_two_centers = 0;
						std::cout<< "al position: "<< n << std::endl;
						size_t left_1; 
						size_t left_2;
						size_t right_1;
						size_t right_2;
						const pw_alignment *p = it->second;
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);
						size_t acc1 = data.accNumber(p->getreference1());
						size_t acc2 = data.accNumber(p->getreference2());
						//A fixed flag before encoding an al:
						unsigned int l1 = al_begin.at(0);
						unsigned int h1 = al_begin.at(1);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::stringstream mem;
						mem << sequenceId << ":" << n;
						std::cout << mem.str() << std::endl;
						std::map<std::string, std::string>::iterator cl = cluster_members.find(mem.str());
						assert(cl != cluster_members.end());
						std::string center = cl->second;
				//		if(center == "0:11:15630"){
				//			p->print();
				//		}
						std::cout << "center " << center << std::endl;
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(1).c_str());
						unsigned int cent_left = atoi(center_parts.at(2).c_str());
						unsigned int center_l;
						unsigned int center_h;
						for(size_t j = 0; j < partition.size(); j++){
							std::map<std::string, std::vector<unsigned int> >::iterator it1 = cluster_high_partition.at(j).find(center);		
							if(it1 != cluster_high_partition.at(j).end()){
								unsigned int par_low = 0;
								if(j != 0){
									par_low = partitionHigh.at(j-1);
								}
								unsigned int par_high = partitionHigh.at(j);
								enc.encode(par_low,par_high,total);
								wrappers.encode(par_low,par_high,total);
								center_l = it1 ->second.at(0);
								center_h = it1->second.at(1);
								break;
							}
						}
						enc.encode(center_l, center_h, total);
						wrappers.encode(center_l,center_h,total);
					//	if(acc1 == data.accNumber(cent_ref) && left_1 ==cent_left ){
					//		model.get_end_al_flag(acc1,acc2, al_end);
					//	}else{
					//		model.get_end_al_flag(acc2,acc1, al_end);							
					//	}
					//	encoding_directions(cent_ref, cent_left, *p, outs, enc);
						// When center is on the sequence:					
						if(cent_ref == sequenceId && cent_left == n){
							if(acc1 == data.accNumber(cent_ref) && left_1 ==cent_left ){
								model.get_end_al_flag(acc1,acc1, al_end);
							}else{
								model.get_end_al_flag(acc2,acc2, al_end);							
							}
							encoding_directions(cent_ref, cent_left, *p, outs, enc);
						}else{//Otherwise: 
							if(acc1 == data.accNumber(cent_ref) && left_1 ==cent_left ){
								model.get_end_al_flag(acc1,acc2, al_end);
							}else{
								model.get_end_al_flag(acc2,acc1, al_end);							
							}
							encoding_directions(cent_ref, cent_left, *p, outs, enc);
							std::vector< std::vector< unsigned int> > low_high;
							model.arith_encode_al(*p,cent_ref, cent_left, low_high);
							assert(low_high.size() < p->alignment_length());
							for(size_t k =0; k < low_high.size(); k++){
								l1 = low_high.at(k).at(0);
								h1 = low_high.at(k).at(1);
								assert(h1<=al_end.at(0));
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							//	if(center == "0:11:15630"){
							//		std::cout << "k " << k << " " << l1 << " " << h1 <<std::endl;
							//		assert(h1<= 16206);
							//	}
							}
						}//end of an al 
						l1 = al_end.at(0);
						h1 = al_end.at(1);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout << l1 << " " << h1 << std::endl;
						std::cout<< "end of an al in enc. "<<std::endl;
					//	p->print();
						if(p->getreference1()== sequenceId && n == left_1){
							n = right_1;
							std::cout<< "n_1 "<< n << std::endl;
						}else{
							n = right_2;
							std::cout<< "n_2 "<< n << std::endl;
						}
					}else{//if there is no alignment in that position 
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
					//	std::cout<< "base "<<base << "position "<< n << std::endl;
						unsigned int h = model.get_seq_high_at_position(sequenceId,n).at(base);
						if(base !=0){
							l = model.get_seq_high_at_position(sequenceId,n).at(base-1);
						}
						if (base == 4){
							assert(h == al_begin.at(0));
						//	h = model.get_powerOfTwo().at(bit);
						}
						std::cout << "base " << base << " low " << l << " high " << h << std::endl;
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
						int base_int = base;
						wrappers.context(base_int);
						assert(h <= al_begin.at(0));
						length_between_two_centers ++;
					}
				}//end of each seq 
				unsigned int l1	= seq_acc_end.at(0);
				unsigned int h1 = seq_acc_end.at(1);
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout << l1 << " " << h1 <<std::endl;
				std::cout<<"end of a seq in enc. "<<std::endl;
			}//end of all sequences of an accsseion 
			unsigned int l1 = seq_acc_end.at(0);
			unsigned int h1 = seq_acc_end.at(1);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout << l1 << " " << h1 <<std::endl;
			std::cout<<"end of an acc"<<std::endl;
		}
		std::cout<< "encoding is finished!"<<std::endl;
	}
	template<typename T>
	void dynamic_encoder<T>::partitioning_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string ,std::vector<pw_alignment> > & alignmentsOfClusters ,std::map<std::vector<std::string>, unsigned int> & weight){//it includes both long and original centers 
		std::vector<std::vector<std::string> >all_centers;
		size_t sumOfWeights = 0;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){//First original centers are added to all centers
			std::vector<std::string> temp_center;
			temp_center.push_back(it->first);
			all_centers.push_back(temp_center);
		}
	//	Afterward long centers are added
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it =new_centers.begin(); it != new_centers.end();it++){
			if(it->first.size() != 1){
				all_centers.push_back(it->first);
			}
		}
		std::cout << "all_centers size is : "<< all_centers.size()<<std::endl;
		for(size_t i = 0 ; i < all_centers.size();i++){
			std::cout << "i " << i<<std::endl;
			for(size_t j =0 ; j < all_centers.at(i).size();j++){
				std::cout<< "j "<< j <<std::endl;
				std::cout<< all_centers.at(i).at(j)<< " ";
			}
			std::cout << "" <<std::endl;
		}
		for(std::map<std::vector<std::string>, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){// All the weights are summed up.
			sumOfWeights += it->second ;
		}
		std::cout << "sum of weight "<< sumOfWeights <<std::endl;
		size_t numberOfPartitions = (sumOfWeights/TOTAL_FOR_ENCODING) + 1;
		std::cout << numberOfPartitions <<std::endl;
		size_t number_of_center = 0;
		for(size_t j = 0 ; j < numberOfPartitions ; j ++){
			long_center_partition.insert(std::make_pair(j,std::vector< std::vector<std::string> >( )));
			size_t center_weight = 0;
			for(size_t i = number_of_center; i < all_centers.size(); i++){
				if(center_weight <= TOTAL_FOR_ENCODING){
					std::map<size_t , std::vector< std::vector<std::string> > >::iterator it = long_center_partition.find(j);
					assert(it != long_center_partition.end());
					it->second.push_back(all_centers.at(i));
					std::map<std::vector<std::string>, unsigned int >::iterator it1 = weight.find(all_centers.at(i));
					center_weight +=  it1->second;
				}else{
					number_of_center = i;
					break;
				}
			}
		}
	}
	template<typename T>
	void dynamic_encoder<T>::calculate_long_center_high(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers,  std::map<std::string ,std::vector<pw_alignment> > & alignmentsOfClusters ,std::map<std::vector<std::string>, unsigned int> & weight){
		std::cout << "calculate long centers high value"<<std::endl;
		partitioning_long_centers(new_centers,alignmentsOfClusters, weight);//Includes both type of centers
		unsigned int partition_high = 0;
		for(size_t i =0; i < long_center_partition.size(); i++){//go through number of partitions
			partition_high = partition_high + (TOTAL_FOR_ENCODING/long_center_partition.size());//calculates the high value of the partion.
			unsigned int h = 0;
			std::vector< std::pair<std::vector<std::string>, std::vector<unsigned int> > >high;
			std::map<size_t , std::vector<std::vector<std::string> > >::iterator it = long_center_partition.find(i);
			assert(it != long_center_partition.end());
			for(size_t j = 0; j <it->second.size(); j++){
				std::vector<std::string> center = it->second.at(j);
				std::map<std::vector<std::string>, unsigned int>::iterator it2 = weight.find(center);
				assert(it2 != weight.end());
				for(size_t k =0; k < center.size();k++){
					std::cout << center.at(k)<< " ";
				}
				std::cout << " "<<std::endl;
				unsigned int Low = h;
				unsigned int High = h + it2->second;
				h = High;
				if(j == it->second.size()-1){
					High = TOTAL_FOR_ENCODING;
				}
				std::vector<unsigned int>temp;
				temp.push_back(Low);
				temp.push_back(High);
				high.push_back(make_pair(center,temp));
				std::cout << Low << " " << High << std::endl;
			}
			long_center_high.push_back(high);//Note that short centers were saved first and then long ones go afterwards.
			HighOfPartition.push_back(partition_high);
		}
		HighOfPartition.at(HighOfPartition.size()-1)= TOTAL_FOR_ENCODING;
		std::cout << "size of the long center high "<< std::endl;
		for(size_t i =0; i < long_center_partition.size(); i++){
			std::cout << long_center_high.at(i).size() <<std::endl;
		}
		
	}
	template<typename T>
	void dynamic_encoder<T>::add_long_center_to_the_stream(std::ofstream & outs){
		std::map<std::string, size_t> index;
		int id =1;
		bool last_short_center = false;
		size_t sizeOfPartition = long_center_partition.size();
		outs.write(reinterpret_cast<char*>(&sizeOfPartition), sizeof(size_t));
		std::cout<< "number of paritions "<< sizeOfPartition << std::endl;
		std::cout << "add all the centers to the stream"<<std::endl;
		for(size_t i = 0 ; i < long_center_partition.size(); i ++){//short centers go first and then come the long ones
			std::cout << "partition "<< i << std::endl;
			std::vector<std::pair<std::vector<std::string>, std::vector<unsigned int> > >center_high = long_center_high.at(i);
			for(size_t j =0; j < center_high.size(); j++){
				std::vector<std::string> center = center_high.at(j).first;
				if(center.size() == 1 ){
					std::cout << "short center " << center.at(0) <<std::endl;
					std::vector<std::string> center_parts;
					strsep(center.at(0), ":" , center_parts);
					unsigned int cent_dir = atoi(center_parts.at(0).c_str());
					unsigned int cent_ref = atoi(center_parts.at(1).c_str());
					unsigned int cent_left = atoi(center_parts.at(2).c_str());
					size_t acc_id = data.accNumber(cent_ref);
					outs.put((char)0);
					outs.write(reinterpret_cast<char*>(&acc_id),sizeof(size_t));
					if(cent_dir == 0){
						index.insert(make_pair(center.at(0),id));
						std::cout << "index is " << id << std::endl;
					}
					else{
						index.insert(make_pair(center.at(0),id));

						std::cout << "index is " << id << std::endl;

					}
					id +=1;
					last_short_center = true;
				}else{
					if(last_short_center == true){
						outs.put(char(5));
						last_short_center = false;					
					}
					outs.put((char)0);
					outs.put((char) center.size());
					for(size_t k = 0 ; k < center.size(); k++){
						std::cout<< center.at(k)<<" ";
						std::vector<std::string> center_parts;
						strsep(center.at(k), ":" , center_parts);
						unsigned int cent_dir = atoi(center_parts.at(0).c_str());
						unsigned int cent_ref = atoi(center_parts.at(1).c_str());
						unsigned int cent_left = atoi(center_parts.at(2).c_str());
						std::map<std::string, size_t>::iterator it = index.find(center.at(k));//write the index of long centers.
						if(it == index.end()){
							unsigned int rev_cent_dir;
							if(cent_dir == 0){
								rev_cent_dir = 1;
							}else{
								rev_cent_dir = 0;
							}
							stringstream rev_cent;
							rev_cent << rev_cent_dir<<":"<< cent_ref<<":"<< cent_left;
							it = index.find(rev_cent.str());
							assert(it != index.end());
						}
						if(cent_dir == 0){
							outs.write(reinterpret_cast<char*>(&it->second),sizeof(int));
							std::cout << "index "<< it->second << " ";
						}else{
							int rev_index = -1*it->second;
							outs.write(reinterpret_cast<char*>(&rev_index),sizeof(int));
							std::cout << "rev_index "<< rev_index<< " " ;
						}
					}
					std::cout << " "<<std::endl;
				}
				std::cout << " " << std::endl;				
				std::vector<bool> bit_to_byte(0);
				int high = center_high.at(j).second.at(1);
				std::cout<< "high " << high << std::endl;
				for(size_t i = 0 ; i < 32; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
				size_t counter = 0;
				for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= (1<<(m-n))* bit_to_byte.at(m);
					}
					n= n+7;
					outs.put(char(a));
					counter = counter + 1;
				}
				std::cout << "counter " << counter << std::endl;
			}
			outs.put(char(7));
		}
		std::cout << "index size is "<< index.size() << std::endl;
		outs.put(char(8));
		std::cout<< "partition size " <<HighOfPartition.size() << std::endl;
		for(size_t i =0; i < HighOfPartition.size(); i++){
		//	outs<<(char)0;
			std::vector<bool> bit_to_byte(0);
			unsigned int high = HighOfPartition.at(i);
			std::cout<<"its high is " << high << std::endl;
			for(size_t j = 0 ; j < 32; j++){
				bit_to_byte.push_back(high%2);
				high = high/2;
			}
			size_t counter = 0;
			for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
				unsigned char a = 0;
				for(size_t m = n; m <n+8; m++){
					a+= (1<<(m-n))* bit_to_byte.at(m);
				}
				n= n+7;
				outs.put(char(a));
				counter = counter + 1;
			}
			std::cout << "counter1 " << counter << std::endl;
		}
	}

	template<typename T>
	void dynamic_encoder<T>::encoding_long_center(std::map<std::string,std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream & outs, dlib::entropy_encoder_kernel_1 & enc){
		setOfAlignments(alignmentOfCluster);
		std::vector<uint32_t> widths(2,0);
		calculate_end_of_center_flags(alignmentOfCluster, widths,outs);
		size_t counter = 0;
		size_t index = 0;
		unsigned int total = TOTAL_FOR_ENCODING;
		enc.set_stream(outs);
		for(size_t i = 0 ; i < long_center_high.size(); i++){
			std::vector<std::pair<std::vector<std::string>, std::vector<unsigned int> > >center_high = long_center_high.at(i);
			for(size_t j =0; j < center_high.size();j++){
				if(center_high.at(j).first.size()==1){//If there is a short center
					std::cout << "low "<< center_high.at(j).second.at(0) << " high "<< center_high.at(j).second.at(1)<<std::endl;
					std::string center = center_high.at(j).first.at(0);
					std::cout << "center at cent_enc is "<< center << std::endl;
					if(center == "1:4:113"){ std::cout << center << std::endl;}
					std::vector<std::string> split;
					strsep(center, ":" , split);
					size_t cent_right = 0;
					size_t cent_ref = atoi(split.at(1).c_str());
					size_t cent_left = atoi(split.at(2).c_str());
					std::map<std::string,std::vector<pw_alignment> >::iterator it1 = alignmentOfCluster.find(center);
					assert(it1 != alignmentOfCluster.end());
					pw_alignment * p = & it1->second.at(0);
					const dnastring & seq = data.getSequence(cent_ref);
					size_t left1;
					size_t right1;
					size_t left2;
					size_t right2;
					p->get_lr1(left1,right1);
					p->get_lr2(left2,right2);
					std::cout<< "id of cent in enc: " << center <<std::endl;
					if(cent_ref == p->getreference1()&& cent_left == left1){
						if(p->getbegin1() < p->getend1()){//Encode the sequence itself
							std::cout<< "id of cent in enc1: " << center <<std::endl;
							cent_right = p->getend1();
							for(size_t k = cent_left; k < cent_right+1;k++){
								size_t base = dnastring::base_to_index(seq.at(k));				
								unsigned int l = 0;
								unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
								if(base !=0){
									l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
								}else l = 0;
								if(base == 4){
							//		h=  TOTAL_FOR_ENCODING;
								}
							//	std::cout<< "l " << l << " h "<< h <<std::endl;
								enc.encode(l,h,total);
								wrappers.encode(l,h,total);
							}
						}else{ //Encode the reverse complement of the sequence
							cent_right = p->getbegin1();
							std::cout<< "id of cent in enc2: " << center <<std::endl;
							for(size_t k = cent_left; k < cent_right+1;k++){
								size_t base = dnastring::base_to_index(seq.at(k));
								size_t com_base = dnastring::base_to_index(base);				
								unsigned int l = 0;
								unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
								if(base !=0){
									l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
								}else l = 0;
								if(base == 4){
							//		h=  TOTAL_FOR_ENCODING;
								}
								enc.encode(l,h,total);
								wrappers.encode(l,h,total);
							//	std::cout<< "l " << l << " h "<< h <<std::endl;

							}
						}
					}
					else if(cent_ref == p->getreference2()&& cent_left == left2){
						if(p->getbegin2()<p->getend2()){//Encode the sequence itself
							cent_right = p->getend2();
							std::cout<< "id of cent in enc3: " << center <<std::endl;
							for(size_t k = cent_left; k < cent_right+1;k++){
								size_t base = dnastring::base_to_index(seq.at(k));				
								unsigned int l = 0;
								unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
								if(base !=0){
									l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
								}else l = 0;
								if(base == 4){
							//		h= TOTAL_FOR_ENCODING;
								}
								enc.encode(l,h,total);
								wrappers.encode(l,h,total);
							//	std::cout<< "l " << l << " h "<< h <<std::endl;

							}
						}else{ //Encode the reverse complement of the sequence
							cent_right = p->getbegin2();
							std::cout<< "id of cent in enc4: " << center <<std::endl;
							for(size_t k = cent_left; k < cent_right+1;k++){
								size_t base = dnastring::base_to_index(seq.at(k));
								size_t com_base = dnastring::base_to_index(base);				
								unsigned int l = 0;
								unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
								if(base !=0){
									l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
								}else l = 0;
								if(base == 4){
							//		h=  model.get_powerOfTwo().at(bit);
								}
								enc.encode(l,h,total);
							//	std::cout<< "l " << l << " h "<< h <<std::endl;
								wrappers.encode(l,h,total);


							}
						}
					}
					unsigned int l1	= widths.at(0);
					unsigned int h1 = widths.at(0)+widths.at(1);					
					enc.encode(l1,h1,total); 
					wrappers.encode(l1,h1,total);
					std::cout << "end of a center "<<std::endl;
					std::cout << "length of the center is " << cent_right+1-cent_left << std::endl;

				}
			}
			std::cout<< "end of a partition"<<std::endl;
		}
	}
	template<typename T>
	void dynamic_encoder<T>::make_fully_reverse_center(std::vector<std::string> & center , std::vector<std::string> & reverse){
		for(size_t i =center.size() ; i >0 ;i--){
			std::vector<std::string> center_parts;
			strsep(center.at(i-1), ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			if(center_dir == 0){
				std::stringstream rev;
				rev << 1 << ":"<< center_ref<<":"<<center_left;
				reverse.push_back(rev.str());
			}else{
				std::stringstream rev;
				rev << 0 << ":" << center_ref << ":" << center_left;
				reverse.push_back(rev.str());
			}
		}

	}
	template<typename T>
	void dynamic_encoder<T>::al_encode_with_long_center(std::vector<std::map<size_t , std::string> >& centerOnSequence, std::map<std::vector<std::string> , unsigned int> & long_center_weight, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::vector<std::map<size_t, std::vector<std::string> > >& centersPositionOnASeq ,std::map<std::string, std::string > & cluster_members,std::ofstream & outs ,dlib::entropy_encoder_kernel_1 & enc, std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
		calculate_long_center_high(new_centers,alignmentOfCluster ,long_center_weight);//Calculates low and high value of all the centers.
		add_long_center_to_the_stream(outs);//Adds all the centers, their accessions and high values to the stream.
		encoding_long_center(alignmentOfCluster,outs,enc);//dont replace it with the one i used for short centers because the order matters!
		unsigned int total = TOTAL_FOR_ENCODING;
		std::vector<uint32_t> al_begin;
		std::vector<uint32_t> al_end;
		std::vector<uint32_t> seq_acc_end;
		std::cout<<"number of acc:" <<data.numAcc()<<std::endl;
		for(size_t i = 0; i< data.numAcc(); i++){
			model.get_seq_flag(i , al_begin , seq_acc_end);
			std::cout << "al begin flag low "<< al_begin.at(0) << " al begin flag high " << al_begin.at(1) <<std::endl;
			std::cout<< " accession "<< i <<std::endl;
			std::cout << "number of seq in this acc : "<< data.getAcc(i).size() <<std::endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
				std::cout << "seq id "<<sequenceId << "length of seq "<< sequence.length() <<std::endl;
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){
						//This flag represents the beginning of an alignment
						//A fixed flag before encoding an al:
						unsigned int l1 = al_begin.at(0);
						unsigned int h1 = al_begin.at(1);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout<< "al position: "<< n << std::endl;
					//	pw_alignment p1 = *it->second;
						std::vector<std::string> current_center;
						std::map<size_t , std::vector<std::string> >::iterator it1=centersPositionOnASeq.at(sequenceId).find(n);
						if(it1 != centersPositionOnASeq.at(sequenceId).end()){//If there is a long center in that position
							std::cout << "LONG CENTER "<<std::endl;
							current_center = it1->second;
						}else{
							std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(sequenceId).find(n);
							current_center.push_back(it1->second);
						} 
						size_t partition = 0; 
						bool CenterIsFound = false;
						for(size_t j = 0; j < long_center_partition.size(); j++){
							std::map<size_t , std::vector<std::vector<std::string> > >::iterator par = long_center_partition.find(j);
							for(size_t k =0; k < par->second.size();k++){
								if(par->second.at(k) == current_center){//If it reaches the end, you should check for the reverse
									partition = j;
									CenterIsFound = true;
									break;
								}
							}
						}
						if(CenterIsFound == false){
							for(size_t j = 0; j < long_center_partition.size(); j++){
								std::map<size_t , std::vector<std::vector<std::string> > >::iterator par = long_center_partition.find(j);
								std::cout << "check for the reverse! " << " center size is " << current_center.size() <<std::endl;
								std::vector<std::string> reverse_center;
								make_fully_reverse_center(current_center,reverse_center);
								for(size_t k =0; k < par->second.size();k++){
									if(par->second.at(k) == reverse_center){
										partition = j;
										current_center = reverse_center;
										CenterIsFound = true;
										break;
									}
								}
							}
						}
						unsigned int low = 0;
						if(partition != 0){
							low = HighOfPartition.at(partition -1);
						}
						unsigned high = HighOfPartition.at(partition);
						enc.encode(low,high,total);//long center partition number is encoded!
						wrappers.encode(low,high,total);
						std::cout << "partition "<< partition << " " << low<< " " << high<<std::endl;
						std::vector< std::pair<std::vector<std::string> , vector<unsigned int> > > cent_high = long_center_high.at(partition);
						for(size_t j = 0; j < cent_high.size();j++){
							if(cent_high.at(j).first == current_center){
								low = cent_high.at(j).second.at(0);
								high = cent_high.at(j).second.at(1);
								break;
							}
						}	
						std::cout << "l "<< low << "h "<< high << std::endl;
						enc.encode(low,high,total);// center id is encoded!
						wrappers.encode(low,high,total);
						unsigned int cent_ref;
						unsigned int cent_left;
					//	size_t gap = 0;
						size_t length = n;//length of alignments and gaps are added to it.
						std::cout << "current center size "<< current_center.size() << std::endl;
						for(size_t j =0; j < current_center.size();j++){
							std::string center = current_center.at(j);
							std::cout<< "center: "<< center << std::endl;
						}
						std::cout << " " <<std::endl;
						std::string center = current_center.at(0);
						std::cout<< "center: "<< center << std::endl;
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						cent_ref = atoi(center_parts.at(1).c_str());
						cent_left = atoi(center_parts.at(2).c_str());
						if(current_center.size() > 1){//long al is encoded. The low/high values of the accession of the first part is used.
							std::cout << "long center occurs in encoding! "<<std::endl;
							std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it2 = new_centers.find(current_center);
							assert(it2 != new_centers.end());
							assert(it2->second.size()>=1);
							for(size_t k =0; k < it2->second.size(); k++){
								size_t left_1,right_1;
								const pw_alignment & p = it2->second.at(k);
								p.get_lr1(left_1,right_1);
								std::cout << " " << std::endl;
								if(p.getreference1()== sequenceId && left_1 ==n){	
									std::cout << " seq id "<< sequenceId << "l1 " << n <<std::endl;
									//The non center reference of the alignment is always forward
									//If center is reverse (50% it is not reverse, 50% it is reverse!!)
									if(p.getbegin2() > p.getend2()){
										std::cout << "reverse long center " << std::endl; 
										//Center is reverse
										unsigned int l1 =  0;
										unsigned int h1 = TOTAL_FOR_ENCODING/4;
										enc.encode(l1,h1,total);
										wrappers.encode(l1,h1,total);
									}else{
										//Both are forwards
										std::cout << "Forward long center " << std::endl; 
										unsigned int l1 = 3*TOTAL_FOR_ENCODING/4;
										unsigned int h1 = TOTAL_FOR_ENCODING;
										enc.encode(l1,h1,total);
										wrappers.encode(l1,h1,total);
									}
									size_t accession = data.accNumber(cent_ref);
									model.get_end_al_flag(accession,i, al_end);// It should be found between the acc of current seq and the acc of the first piece of the center
									std::cout << "al end: " << al_end.at(0) << " "<<al_end.at(1) <<std::endl;
									std::cout<< "flags are set!"<<std::endl;
									std::vector< std::vector< unsigned int> > low_high;
									model.arith_encode_long_al(p, i, accession ,cent_ref, cent_left, low_high);
									assert(low_high.size() < p.alignment_length());
									for(size_t k =0; k < low_high.size(); k++){
										unsigned int l1 = low_high.at(k).at(0);
										unsigned int h1 = low_high.at(k).at(1);
									//	std::cout << "l "<< l1 << " h "<< h1 <<std::endl;
										assert(h1<=al_end.at(0));
										enc.encode(l1,h1,total);
										wrappers.encode(l1,h1,total);
									}
									length = right_1;
									n = right_1;
									break;
								}	
							}		
						}else{//short al
							const pw_alignment & p= *it->second;
							size_t left_1; 
							size_t left_2;
							size_t right_1;
							size_t right_2;
							p.get_lr1(left_1,right_1);
							p.get_lr2(left_2,right_2);
							size_t acc1 = data.accNumber(p.getreference1());
							size_t acc2 = data.accNumber(p.getreference2());
							if(cent_ref == sequenceId && cent_left == n){//If center is on the the seq
								if(acc1 == data.accNumber(cent_ref) && left_1 ==cent_left ){
									model.get_end_al_flag(acc1,acc1, al_end);
								}else{
									model.get_end_al_flag(acc2,acc2, al_end);							
								}
								encoding_directions(cent_ref, cent_left, p, outs, enc);
							}else{//Otherwise
								if(acc1 == data.accNumber(cent_ref) && left_1 ==cent_left ){
									model.get_end_al_flag(acc1,acc2, al_end);
								}else{
									model.get_end_al_flag(acc2,acc1, al_end);							
								}
								encoding_directions(cent_ref, cent_left, p, outs, enc);
								std::vector< std::vector< unsigned int> > low_high;
								model.arith_encode_al(p,cent_ref, cent_left, low_high);
								assert(low_high.size() < p.alignment_length());
								for(size_t k =0; k < low_high.size(); k++){
									l1 = low_high.at(k).at(0);
									h1 = low_high.at(k).at(1);
									assert(h1<=al_end.at(0));
									enc.encode(l1,h1,total);
									wrappers.encode(l1,h1,total);
								}
							}
							std::cout << "al end: " << al_end.at(0) << " "<<al_end.at(1) <<std::endl;
							size_t l1,l2,r1,r2;
							p.get_lr1(l1,r1);
							p.get_lr2(l2,r2);
							if(p.getreference1()== sequenceId && n == l1){
								length = r1;
								n = length;
							}else{
								length = r2;
								n = length;
							}
						}
						std::cout<< "end position of an al " << n <<std::endl;
						//Flag shows end of an alignment
						l1 = al_end.at(0);
						h1 = al_end.at(1);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout << l1 << " " << h1 << std::endl;
						std::cout<< "end of an al in enc. "<<std::endl;
					}else{ //If there is no alignment in that position
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
						unsigned int h = model.get_seq_high_at_position(sequenceId,n).at(base);
						if(base !=0){
							l = model.get_seq_high_at_position(sequenceId,n).at(base-1);
						}
						if (base == 4){
							assert(h == al_begin.at(0));
						}
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
						int base_int = base;
					//	std::cout << "base " << base << " low " << l << " high " << h << std::endl;

						wrappers.context(base_int);
						assert(h <= al_begin.at(0));
					}
				}//End of a sequence
				unsigned int l1	= seq_acc_end.at(0);
				unsigned int h1 = seq_acc_end.at(1);
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout << l1 << " " << h1 <<std::endl;
				std::cout<<"end of a seq in enc. "<<std::endl;
			}//end of all sequences of an accsseion
			unsigned int l1	= seq_acc_end.at(0);
			unsigned int h1 = seq_acc_end.at(1);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout << l1 << " " << h1 <<std::endl;
			std::cout<<"end of an acc"<<std::endl;
		}
	}

#endif
