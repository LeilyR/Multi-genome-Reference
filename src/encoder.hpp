#ifndef ENCODER_HPP
#define ENCODER_HPP

#include "data.hpp"
#include "dynamic_mc.hpp"
#include <iostream>
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

#include <map>
#include <vector>
#include <cassert>


template<typename T>
class encoder{
	public:
	encoder(all_data & d, T & a_model,wrapper & wrap): data(d),model(a_model),wrappers(wrap),upper_bound(d.numSequences()),AlignmentsFromClustering(data.numSequences()){

	}

//	encoder( all_data& , const T &, wrapper &);
	~encoder(){};
	void arithmetic_encoding_seq(std::ofstream&);
	void arithmetic_decoding_seq_test(int t);
	void calculating_clusters_high(std::map<std::string, unsigned int> & weight);
	void encoding_seq_test();
	void arithmetic_decoding_seq();
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//all the alignment that has a certain sequence as at least one of their references.
	//	std::ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				std::string center = it->first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_ref = atoi(center_parts.at(0).c_str());
				unsigned int center_left = atoi(center_parts.at(1).c_str());	
				for(size_t j = 0; j < it->second.size();j++){
					pw_alignment * p = & it->second.at(j);
					size_t left1; 
					size_t left2;
					size_t right1;
					size_t right2;
					p->get_lr1(left1,right1);
					p->get_lr2(left2,right2);
				//	size_t ref1 = p->getreference1();
				//	size_t ref2 = p->getreference2();
					if(i != center_ref){
						if(p->getreference1()== i && p->getreference2()== center_ref&&left2==center_left){
							std::multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left1,p));
							//	std::string pattern = model.get_context(right1,i);
							//	outs << pattern;
							//	outs << char(0);
							}else continue;
						}else if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
							std::multimap<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
							if(it2 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(std::make_pair(left2,p));
							//	std::string pattern = model.get_context(right2,i);
							//	outs<< pattern;
							//	outs<< char(0);
							}else continue;
						}else continue;
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
						//	size_t center_end = center_left + p->alignment_length()-1;
						//	std::string pattern = model.get_context(center_end,i);
						//	outs << pattern;
						//	outs<< char(0);
							}else continue;
						}
					}
				}
			}
		//	outs<<char(8);
		}
	//	for(std::multimap<size_t,pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(0).begin();it2 !=AlignmentsFromClustering.at(0).end();it2++){
	//		pw_alignment *p = it2 -> second;
	//		p->print();
		//	std::cout<<"left on the reference: "<< it2->first<<std::endl;
	//	}
	//	outs << char(7);
	}
	void add_acc_to_stream(std::ofstream &);
	void set_acc_from_stream(std::ifstream &);
	void set_pattern_from_stream(std::ifstream & );
	const std::map<std::string, std::vector<unsigned int> > & get_center_high() const;
	void partition_centers(std::map<std::string, vector<pw_alignment> > & alignmentsOfClusters, std::map<std::string, unsigned int>& weight){ //dividing centeres into different partitions (For facing the problem of using not larger than 2^13)
		size_t bit  =13;
		unsigned int length = model.get_powerOfTwo().at(bit);		
		size_t sum_of_weight = 0;
		for(std::map<std::string, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){
			sum_of_weight = sum_of_weight + it->second ;
		}
		size_t numberOfPartitions = (sum_of_weight/length) + 1;
			std::cout<< "no. of partition:" << numberOfPartitions << std::endl;
		//	std::cout<< "size of al of cluster" << alignmentsOfClusters.size() << std::endl;
		std::vector<std::string> center;
//		for(size_t k = 0; k < data.numAcc(); k++){
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it= alignmentsOfClusters.begin(); it!= alignmentsOfClusters.end();it++){
				std::string c_id = it->first;
				std::vector<std::string> split;
				strsep(c_id, ":" , split);
				size_t c_ref = atoi(split.at(0).c_str());
				size_t c_left = atoi(split.at(1).c_str());
				size_t acc = data.accNumber(c_ref);
			//	if(acc == k){
					center.push_back(c_id);
			//	}
			}
	//	}
		for(size_t j = 0 ; j < numberOfPartitions ; j ++){	
			partition.insert(std::make_pair(j,std::vector<std::string>( )));//check if doesnt make a center with a high value bigger than total in encoding
			for(size_t i =j*alignmentsOfClusters.size()/numberOfPartitions; i < (j+1)*alignmentsOfClusters.size()/numberOfPartitions; i++){	
				std::map<size_t , std::vector<std::string> >::iterator it = partition.find(j);
				it->second.push_back(center.at(i));
			} 
		}
	/*	for(std::map<size_t , std::vector<std::string> >::iterator it = partition.begin(); it != partition.end(); it++){
			std::cout<< "partition "<< it->first<<std::endl;
			for(size_t i = 0; i <it->second.size(); i++){
				std::cout<< it->second.at(i)<<std::endl;
			}
		}*/

	}
	void calculate_high_in_partition(std::map<std::string, unsigned int> & weight, std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){
		partition_centers(alignmentsOfClusters, weight);
		std::cout<<"partition size: "<< partition.size() <<std::endl;
		for(size_t i = 0; i < partition.size(); i ++){
			unsigned int h = 0;
			size_t bit = 13;
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
				it1 ->second.at(0) = h;
				it1 -> second.at(1) = h + it2->second;
				h = it1 -> second.at(1);
				if(j == it->second.size()-1){
					it1->second.at(1) = model.get_powerOfTwo().at(bit);
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
	void partition_high(){
		size_t bit = 13;
		unsigned int weight = model.get_powerOfTwo().at(bit);		
		unsigned int high = 0;
		for(size_t i = 0; i < partition.size(); i++){
			high = high + weight/partition.size();
			partitionHigh.push_back(high);
		}
	}
	void add_center_to_stream(std::ofstream & outs){
		size_t bit =32;
		for(size_t i = 0 ; i < partition.size(); i ++){
		/*	for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				std::cout<< "center_write: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
			}*/
			for ( std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end(); it ++ ){
				std::string center = it -> first;
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
				std::vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
				for(size_t i = 0 ; i < bit; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
			//	for(size_t n = 0 ; n < bit_to_byte.size();n++){
			//		std::cout<< " "<< bit_to_byte.at(n);
			//	}
			//		std::cout<< "" << std::endl;
			//	std::cout<<"bit to byte: "<< bit_to_byte.size() <<std::endl;
				for(size_t n = 0; n < bit_to_byte.size(); n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
					}
					n= n+7;
					outs<< a;
				}
			}
			outs << (char) 7;
		}
		outs << (char) 8;
	}
	void set_center_from_stream (std::ifstream & in){
		size_t bit = 32;
		char c ; 
		unsigned char h;
		c = in.get();
		size_t i = 0;
		while(c != 8){
			unsigned int low = 0;
			while(c != 7){
				std::cout << "here!"<<std::endl;
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
				std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).find(center_id.str());
				if(it2 == cluster_high_partition.at(i).end()){
					cluster_high_partition.at(i).insert(std::make_pair(center,std::vector<unsigned int>(2,0)));
					it2 = cluster_high_partition.at(i).find(center);	
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
			c = in.get();
		/*	for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				std::cout<< "center_red_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
			}*/
			i = i + 1;
			std::cout<< "i: "<< i << std::endl;
		}
	//	std::cout<< "END!"<<std::endl;
	/*	for(std::map<size_t, std::vector<std::string> >::iterator it = acc_of_center.begin();it!=acc_of_center.end();it++){
			std::cout<<"acc std::set in the map: "<<it->first<<std::endl;
		}*/
	}
	void arithmetic_encoding_centId(std::map<std::string, std::vector<pw_alignment> > &, std::ofstream &);
	void arithmetic_decoding_centId(std::ifstream&, dlib::entropy_decoder_kernel_1&);
	void arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster,std::ofstream & outs,dlib::entropy_encoder_kernel_1 & enc){// after partitioning the clusters (New one!)
	//	arithmetic_encoding_centId(alignmentOfCluster,outs);
		add_center_to_stream(outs);//center id, its acc and high
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		size_t bit = 13;
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
	//	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc.set_stream(outs);
	//	std::cout<< "size of  cluster_high_partition: "<<  cluster_high_partition.size()<<std::endl;
		for(size_t i = 0 ; i < cluster_high_partition.size(); i++){
		//	for(size_t j = 0 ; j < data.numAcc(); j++)
			for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin();it != cluster_high_partition.at(i).end(); it++){
				std::string center = it ->first;
				std::vector<std::string> split;
				strsep(center, ":" , split);
				size_t cent_right = 0;
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				//if(data.accNumber(cent_ref) == j)
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
				//	std::cout<< "encoded center1: "<< it->first << " its al length: "<< p->alignment_length()<<std::endl;
					if(p->getbegin1() < p->getend1()){//Encode the sequence itself
						std::cout<< "id of cent in enc1: " << center <<std::endl;
						cent_right = p->getend1();
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	std::cout<<"base _enc: "<<base <<std::endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
						//	std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
					//	std::cout<<"cent length1: "<< cent_right - cent_left << std::endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin1();
					/*	std::cout<< "id of cent in enc2: " << center <<std::endl;
						for(size_t k = cent_right; k >= cent_left;k--){
							size_t base = dnastring::base_to_index(seq.at(k));
							char co_base = dnastring::complement(seq.at(k));
							size_t com_base = dnastring::base_to_index(co_base);	
						//	std::cout<<"base _enc: "<<base <<std::endl;
						//	std::cout<< "com_base_enc: "<< com_base <<std::endl;
							unsigned int l = 0;
							unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
							if(com_base !=0){
								l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
							}else l = 0;
							if(com_base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
						//	std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
					//	std::cout<<"cent length2: "<< cent_right - cent_left << std::endl;*/
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
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
						}

					}
				}
				if(cent_ref == p->getreference2()&& cent_left == left2){
				//	std::cout<< "encoded center2: "<< it->first << " its al length: "<< p->alignment_length()<<std::endl;				
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
						std::cout<< "id of cent in enc3: " << center <<std::endl;
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	std::cout<<"base _enc: "<<base <<std::endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
					//		std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
					//	std::cout<<"cent length3: "<< cent_right - cent_left << std::endl;
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
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
						}
					}
				}
		//		std::cout<<"end of a center!"<<std::endl;
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit)+5;					
				unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
				enc.encode(l1,h1,t1); 
			}
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit)+10;					
			unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
			enc.encode(l1,h1,t1); 
		//	std::cout<< "partition: " << i <<std::endl;
		}
	//	delete enc; 
	//	std::cout<<"end of cent encoding!"<<std::endl;
	}
	void arithmetic_decoding_centers(std::ifstream & in,dlib::entropy_decoder_kernel_1 & dec){// After partitioning the cluster! (new one!)
		std::ofstream save1;
		save1.open("decode.txt");
	//	arithmetic_decoding_centId(in,dec);
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
		size_t bit =13;
		size_t acc = 0;
	//	dlib::entropy_decoder_kernel_1 dec;
		dec.set_stream(in);	
		std::string first_pattern = model.get_firstPattern();
		unsigned int target;
		size_t lengthOfSeq = 0;
		size_t flag = model.get_powerOfTwo().at(bit);
		unsigned int total = model.get_powerOfTwo().at(bit)+10;
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
			size_t m = 0;			
			std::vector<std::string> centerOfAPartition;
			for(std::map<std::string, std::vector<unsigned int> >::iterator it =cluster_high_partition.at(counter).begin(); it != cluster_high_partition.at(counter).end();it ++){
				centerOfAPartition.push_back(it->first);		
			}
			std::cout<< "size of centerOfAPartition: "<<centerOfAPartition.size()<<std::endl;
			std::map<std::string, std::string> intermediate;
		//	while(acc < data.numOfAcc())
			target = dec.get_target(total);
			while(target<flag+5){//end of all sequences of a partition
		//		std::cout << " m "<< m <<std::endl;
				std::string c_id = centerOfAPartition.at(m);
				for(std::map<size_t, std::vector<std::string> >::iterator it = acc_of_center.begin();it != acc_of_center.end(); it ++){	
					for(size_t k = 0; k < it->second.size();k++){		
						if(it ->second.at(k) == c_id){
							acc = it->first;
							break;
						}
					}
				}
		//		std::cout<< "acc in decoding is "<< acc<<std::endl;
				std::string center_seq;
				std::vector<unsigned int>high(7,0);
				std::vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				std::string pattern = first_pattern;
				std::cout<<"first pattern: " << first_pattern<<std::endl;
				std::map<std::string, std::vector<unsigned int> >::const_iterator it=model.get_high(acc).find(first_pattern);
				for(size_t j=0; j<5; j++){
					high.at(j) = it->second.at(j);
					if(j!= 0){
						low.at(j)=high.at(j-1);
					}else low.at(j)= 0;
				}
				high.at(4)=  model.get_powerOfTwo().at(bit);
				low.at(5)= high.at(4);
				high.at(5)= model.get_powerOfTwo().at(bit) + 5;
				low.at(6)= high.at(5);
				high.at(6) = model.get_powerOfTwo().at(bit) + 10;
				size_t base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				center_seq +=dnastring::index_to_base(base);
				dec.decode(low.at(base),high.at(base));
			//	std::cout<<"first base: "<<base<<std::endl;
			//	std::cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << std::endl;
				target = dec.get_target(total);
				while(target < flag){	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					std::string current_pattern;
					std::stringstream s;
					if(Sequence_level > 1){
						for(size_t M=1; M< Sequence_level; M++){// ye doone function to model vase current pattern ham benevis!
							s<<pattern.at(Sequence_level-M);
						}
						s<<p;
						s>>current_pattern;
					}
					if(Sequence_level == 1){
					//	s<<pattern.at(0);
						s<<p;
						s>>current_pattern;
					}
					std::cout<<"current pattern: "<< current_pattern<<std::endl;
					std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_high(acc).find(current_pattern);
					assert(it1 != model.get_high(acc).end());
					for(size_t j=0; j<5; j++){
						high.at(j) = it1->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}	
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + 5;
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + 10;
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					center_seq +=dnastring::index_to_base(base);
					dec.decode(low.at(base),high.at(base));
				//	save1<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " t: "<< total << std::endl;
				//	std::cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << std::endl;
				//	std::cout << "base: "<< base << std::endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				}
			//	std::cout << "length of seq" << lengthOfSeq << std::endl;
			//	save1 << "length "<< lengthOfSeq <<std::endl;
				base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
			//	std::cout<<"flag: "<< base << std::endl;
				dec.decode(low.at(base),high.at(base));
				std::string center_id = centerOfAPartition.at(m);
				m = m +1;
			//	std::cout<< "id of cent in dec: " << center_id <<std::endl;
				intermediate.insert(std::make_pair(center_id,center_seq));
				target = dec.get_target(total);	
				std::cout<< "target: "<<target << std::endl;
				//Decoding the flag which represents end of a partition
				if(target >= flag){
					std::cout<< "end of partition target: "<<target << std::endl;					
					decoded_center_in_partition.push_back(intermediate);
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));
			//		std::cout<< "base3: " << base <<" " << low.at(base) << " " << high.at(base)  << std::endl;
				} else continue;
			}
		//	std::cout << " counter "<< counter <<std::endl;
			counter = counter + 1 ;
		}
	/*	for(size_t i =0; i < decoded_center_in_partition.size(); i++){
			for(std::map<std::string, string>::iterator it = decoded_center_in_partition.at(i).begin();it != decoded_center_in_partition.at(i).end();it++){
				std::string center = it ->first;
				std::cout<<"center size in decoded map is: "<<std::endl;				
				std::cout<< center.size() << std::endl;
			}
		}*/
		save1.close();
	}
	void add_partition_high_to_stream(std::ofstream & outs){
		size_t bit = 32;
		for(size_t i = 0; i < partitionHigh.size(); i++){
			outs<<(char)0;
			std::vector<bool> bit_to_byte(0);
			unsigned int high = partitionHigh.at(i);
		//	std::cout<< "high_write: "<< high << std::endl;
			for(size_t j = 0 ; j < 32; j++){
				bit_to_byte.push_back(high%2);
				high = high/2;
			}
			for(size_t n = 0; n < bit_to_byte.size(); n++){
				unsigned char a = 0;
				for(size_t m = n; m <n+8; m++){
					a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
				}
				n= n+7;
				outs<< a;
			//	std::cout<< int(a)<<std::endl;
			}
			outs << (char) 7;
		}
		outs << (char) 8;
	}
	void set_partition_high_from_stream(std::ifstream& in){
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
					high_value += binary_high_value.at(i)*model.get_powerOfTwo().at(i);
				}
				partitionHigh.push_back(high_value);
			//	std::cout<< "high_read: "<< high_value<<std::endl;
				c=in.get();
			}
			c=in.get();
		}
	}
	void al_encoding(std::map<std::string, unsigned int> & weight, std::map<std::string, std::string > & cluster_members, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &outs,dlib::entropy_encoder_kernel_1 & enc ){
		calculate_high_in_partition(weight,alignmentOfCluster);
		encoding_functor functor(data,&model,wrappers,enc);
		arithmetic_encoding_centers(alignmentOfCluster,outs,enc);
		size_t bit = 13;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
//		size_t i = 1;
		std::cout<<"number of acc:" <<data.numAcc()<<std::endl;
		for(size_t i = 0; i< data.numAcc(); i++){
			std::cout<< " i "<< i <<std::endl;
			std::cout << "number of seq in an acc : "<< data.getAcc(i).size() <<std::endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
		/*		if(sequenceId == 20){
					std::cout<< "sequence 20 at 203 and 204:"<<std::endl;
					for(size_t j = 203; j<=204; j++){
						std::cout<< sequence.at(j);
					}
				}
				std::cout<< " " <<std::endl;*/
				std::cout<< "sequence length is: "<< sequence.length()<< " sequence id: "<< sequenceId <<std::endl;
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
						std::cout<< "al position: "<< n << std::endl;
						size_t left_1; 
						size_t left_2;
						size_t right_1;
						size_t right_2;
						pw_alignment *p = it->second;
						std::cout<< "al length: "<< p->alignment_length()<<std::endl;
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);		
						std::cout<<" l1 "<<left_1 << " l2 " << left_2 << " r1 "<< right_1 << " r2 "<<right_2 <<std::endl;
						std::cout << "begin1 "<< p->getbegin1() << " begin2 "<< p->getbegin2() << " end1 "<< p->getend1() << " end2 " << p->getend2() << std::endl;
						//A fixed flag before encoding a center:
						std::cout<<"enc  bal flag"<<std::endl;
						unsigned int l1 = model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit)+5;
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::stringstream mem;
						mem << sequenceId << ":" << n;
						std::map<std::string, std::string>::iterator cl = cluster_members.find(mem.str());
						assert(cl != cluster_members.end());
						std::string center = cl->second;
						std::cout<< "center: "<< center << std::endl;
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(0).c_str());
						unsigned int cent_left = atoi(center_parts.at(1).c_str());
						unsigned int center_l;
						unsigned int center_h;
						size_t part;
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
								std::cout << "par low: "<< par_low << " par high: "<< par_high << std::endl;									
								center_l = it1 ->second.at(0);
								center_h = it1->second.at(1);
								part = j;
								break;
							}
						}
						std::cout<< "center low "<< center_l << " center high " << center_h<<"center left : "<< cent_left <<std::endl;
						enc.encode(center_l, center_h, total);
						wrappers.encode(center_l,center_h,total);
						if(cent_ref == sequenceId && cent_left == n){
							std::cout<< " center is on the ref "<< std::endl;	
							if(cent_left == p->getend1()||cent_left ==p->getend2()){
								std::cout<<"reverse center!"<<std::endl;
							//	p->print();
							//	unsigned int l1 =  model.get_powerOfTwo().at(bit) + 15;
							//	unsigned int h1 = model.get_powerOfTwo().at(bit) + 20;
							//	enc.encode(l1,h1,total);
							//	wrappers.encode(l1,h1,total);
							}else{
							}					
						}else{
							if((cent_ref == p->getreference1() && cent_left == p->getend1()&& left_2==p->getbegin2())||(cent_ref == p->getreference2()&& cent_left ==p->getend2()&&left_1 == p->getbegin1())){
								std::cout<<"reverse center!"<<std::endl;
								//center is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 0;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 5;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);

							}
							if((cent_ref == p->getreference1()&& cent_left == p->getbegin1()&& left_2==p->getend2())||(cent_ref == p->getreference2()&& cent_left ==p->getbegin2()&&left_1 == p->getend1())){
							//	p->print();
								std::cout<<"other ref is reverse!"<<std::endl;
								//other ref is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 10;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 15;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}
							if((cent_ref == p->getreference1()&& cent_left == p->getend1()&& left_2==p->getend2())||(cent_ref == p->getreference2() && cent_left ==p->getend2()&&left_1 == p->getend1())){
								p->print();
								std::cout<<"reverse_both!"<<std::endl;
								//both are reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 15;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 20;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}

							model.get_encoded_member(*p,cent_ref,cent_left,functor,outs);
							std::cout<< " center is not on the ref "<< std::endl;
						}//end of each al
						l1 = model.get_powerOfTwo().at(bit)+5;
						h1 = model.get_powerOfTwo().at(bit) +10;
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout<< "end of an al in enc. "<<std::endl;
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
						unsigned int h = model.get_high_at_position(sequenceId,n).at(base);
						if(base !=0){
							l = model.get_high_at_position(sequenceId,n).at(base-1);
						}
						if (base == 4){
							h = model.get_powerOfTwo().at(bit);
						}
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
						int base_int = base;
						wrappers.context(base_int);
					}
				}//end of each seq
				unsigned int l1	= model.get_powerOfTwo().at(bit)+10;
				unsigned int h1 = model.get_powerOfTwo().at(bit) +15;
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout<<"end of a seq in enc. "<<std::endl;
			}//end of all sequences of an accsseion
			unsigned int l1 = model.get_powerOfTwo().at(bit)+15;
			unsigned int h1 = model.get_powerOfTwo().at(bit) +20;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout<<"end of an acc"<<std::endl;
		}
		std::cout<< "encoding is finished!"<<std::endl;
	}
//	void al_encoding(std::map<std::string, unsigned int> & , std::map<std::string, std::string > & , std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &, dlib::entropy_encoder_kernel_1 &);
	void al_decoding(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec){
		arithmetic_decoding_centers(in,dec);
		std::cout<<"decoding the center has been done!"<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		std::string first_al_pattern = model.get_firstAlignmentPattern();
		std::string first_pattern = model.get_firstPattern();
		std::string al_pattern;
		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);
		size_t i = 0;
		while(i<data.numOfAcc()){
			std::cout<< "accession "<< i <<std::endl;
			target = dec.get_target(total);
			while(target < flag + 15){//end of an accession
				std::string decodedCenter;
				std::string pattern;
				std::vector<unsigned int>high(9,0);
				std::vector<unsigned int>low(9,0);
				if(target<flag){//if there is no alignment on the first position
					std::cout<< "Sequence doesn't start with an alignment!" <<std::endl;
					pattern = first_pattern;
					std::cout<< "pattern size: "<<pattern.size()<<std::endl;
					std::map<std::string, std::vector<unsigned int> >::const_iterator it=model.get_high(i).find(first_pattern);
					for(size_t j=0; j<5; j++){
						high.at(j) = it->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + 5;
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + 10;
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + 15;
					low.at(8)= high.at(7);
					high.at(8) = model.get_powerOfTwo().at(bit) + 20;
					base = 12;
					for(size_t n = 0; n < 9 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));	
					wrappers.decode(low.at(base),high.at(base),total);
					int base_int = base;
					wrappers.decodeContext(base_int);
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}
				else {
					std::cout<<"Sequence starts with an alignment!"<<std::endl;
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + 5;
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + 10;
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + 15;
					low.at(8)= high.at(7);
					high.at(8) = model.get_powerOfTwo().at(bit) + 20;
				}
				while(target<flag+10){//end of each seq
					size_t reverse_center = 0;
					size_t reverse_member = 0;
					size_t reverse_both = 0;
					while(target<flag){//If there is no alignment on that position (see if i can move it to the else of alignment if)
						char p = dnastring::index_to_base(base);	
						std::stringstream s;
						std::string current_pattern;
						if(Sequence_level > 1){// TODO: Make a function in model class which returns the current pattern!
							for(size_t M=1; M< Sequence_level; M++){
								s<<pattern.at(Sequence_level-M);
							}
							s<<p;
							s>>current_pattern;
						}
						if(Sequence_level ==1){	
						//	s<<pattern.at(0);
							s<<p;
							s>>current_pattern;
						}
					//	std::cout<< "current pattern: " << current_pattern<<std::endl;
						std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
						assert(it1 != model.get_high(i).end());
						for(size_t j=0; j<5; j++){
							high.at(j) = it1->second.at(j);
							if(j!= 0){
								low.at(j)=high.at(j-1);
							}else low.at(j)= 0;
						}	
						high.at(4)=  model.get_powerOfTwo().at(bit);
						low.at(5)= high.at(4);
						high.at(5)= model.get_powerOfTwo().at(bit) + 5;
						low.at(6)= high.at(5);
						high.at(6) = model.get_powerOfTwo().at(bit) + 10;
						low.at(7)= high.at(6);
						high.at(7) = model.get_powerOfTwo().at(bit) + 15;
						low.at(8)= high.at(7);
						high.at(8) = model.get_powerOfTwo().at(bit) + 20;
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						int base_int = base;
						wrappers.decodeContext(base_int);
					//	std::cout<<"base: "<<base<<std::endl;
						pattern = current_pattern;
						target = dec.get_target(total);	
					//	std::cout << "target2: "<<target<<std::endl;
						if(target> flag){
							std::cout<< "starting of an alignemnt's flag!"<<std::endl;
						}
					}
					if(flag<=target && target < flag+10){
						if(flag<=target && target < flag+5){
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
						}else{
							std::cout<< "there is a problem in bal flag"<<std::endl;
						}
						std::string center;
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
						dec.decode(low_par,high_par);
						wrappers.decode(low_par,high_par,total);
						target = dec.get_target(total);
						for(std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(number_of_par).begin(); it2!= cluster_high_partition.at(number_of_par).end();it2++){
							if(it2->second.at(0) <= target && it2->second.at(1) > target){
								center = it2->first;
								std::cout <<"center: "<< center<< std::endl;
								std::cout<< "center target: "<< target << "total" << total<<std::endl;
								dec.decode(it2->second.at(0),it2->second.at(1));
								wrappers.decode(it2->second.at(0),it2->second.at(1),total);
								std::cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<std::endl;
								break;
							}	
						}		
						std::map<std::string, std::string>::iterator seq = decoded_center_in_partition.at(number_of_par).find(center);
						assert(seq != decoded_center_in_partition.at(number_of_par).end());
						decodedCenter =seq->second;
						std::cout<< "decoded center: "<< decodedCenter <<std::endl;
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
						size_t pos = 0;
						//when rev_center
						target = dec.get_target(total);	
						if(target>=flag+0&&target<flag+5){
							std::cout<<"rev_center"<<std::endl;
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
							reverse_center ++;
							target = dec.get_target(total);
						}
						//when rev_member
						if(target>=flag+10&&target<flag+15){
							std::cout<<"rev_member"<<std::endl;
							dec.decode(low.at(7),high.at(7));
							wrappers.decode(low.at(7),high.at(7),total);
							reverse_member ++;
							target = dec.get_target(total);
						}
						//both rev
						if(target>=flag+15&&target<flag+20){
							std::cout<<"rev_both"<<std::endl;
							dec.decode(low.at(8),high.at(8));
							wrappers.decode(low.at(8),high.at(8),total);
							reverse_both ++;
							target = dec.get_target(total);
						}
						size_t modify = 0;
						std::string member;
						while(target < flag){// the way a context is created doesn't work if the alignment_level is greater than 1!!
							std::vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
							std::vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
							size_t last_base;
							if(reverse_center > 0){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							//	std::string rev_decodedCenter;
							//	for(size_t rev = 0; rev< decodedCenter.size(); rev++){
							//		char chr= dnastring::complement(decodedCenter.at(decodedCenter.size()-1-rev));	
							//		rev_decodedCenter +=chr;
							//	}
							//	last_base = dnastring::base_to_index(rev_decodedCenter.at(pos));
						//		last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
						//	}
							//else if(reverse_member > 0){
							//	last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							//	if(pos==0){
							//		std::cout<< "last base: " << last_base << std::endl;
							//	}
						//	
							else if(reverse_both > 0){
							//	std::string rev_decodedCenter;
							//	for(size_t rev = 0; rev< decodedCenter.size(); rev++){
							//		char chr= dnastring::complement(decodedCenter.at(decodedCenter.size()-1-rev));	
							//		rev_decodedCenter +=chr;
							//	}
							//	last_base = dnastring::base_to_index(rev_decodedCenter.at(pos));
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							else{//both references are forward
								last_base = dnastring::base_to_index(decodedCenter.at(pos));
							}
							std::string al_context = al_pattern;
							al_context += last_base;
						//	std::cout<< "size of context: "<< al_context.size() << std::endl;
							std::cout<< "decoding_context: ";
							for(size_t k =0; k < al_context.size(); k++){
								std::cout<< int(al_context.at(k));
								int con =  int(al_context.at(k));
								wrappers.decodeContext(con);
							}
							std::cout<< " " <<std::endl;
							std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are considered from center to the other reference
							assert(mod != model.get_highValue(cent_acc,i).end());
						//	std::cout << "cent_acc " << cent_acc << " other_acc " << i <<std::endl;
							for(size_t k = 0 ; k < mod->second.size();k ++){
								al_high.at(k) = mod->second.at(k);
								if(k > 0){
									al_low.at(k) = al_high.at(k-1);
								}else{
									al_low.at(k) = 0;
								}				
							}
						//	std::cout << "mod->second size: " << mod->second.size() << " target: "<< target<<std::endl;
							al_high.at(mod->second.size()-1)= model.get_powerOfTwo().at(bit);
							size_t modification = NUM_KEEP+NUM_DELETE+20;
							for(size_t n = 0; n < NUM_KEEP+NUM_DELETE+10; n++){
								if(al_low.at(n)<=target && al_high.at(n)> target){
									modification = n;
									break;
								}
							}
							dec.decode(al_low.at(modification),al_high.at(modification));
							wrappers.decode(al_low.at(modification),al_high.at(modification),total);
							modify = modification;
							std::cout<< " modify: " << modify << std::endl;
						//	string cur_pattern;
						//	stringstream s;
						//	if(Alignment_level > 1){
						//		for(size_t M=0; M< Alignment_level; M++){
						//			s<<al_pattern.at(Alignment_level-1-M);
						//		}
						//	}
						//	s<<modification;
						//	s>>cur_pattern;
						//	al_pattern=cur_pattern;
							al_pattern = modification;
						//	std::cout<<"al_pattern: ";
						//	for(size_t al = 0; al < al_pattern.size();al++){
						//		std::cout<< int(al_pattern.at(al));
						//	}
						//	std::cout<<" " << std::endl;
							if(reverse_center>0 || reverse_both > 0){
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

							}
							assert(pos < decodedCenter.size());
							pos += model.modification_length(modification);
							std::cout<< "position on center: " << pos << std::endl;
							std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
							target = dec.get_target(total);	
						}
						if(flag+5 <= target && target < flag+10){
							dec.decode(low.at(6), high.at(6));
							wrappers.decode(low.at(6),high.at(6),total);
							std::cout << " End of an al! "<<std::endl;
						}
						else{
						std::cout<<"there is something wrong!"<<std::endl;
						}
						std::string temp;
						std::cout<< "length of the member: "<<member.length()<<std::endl;
						if(member.length()> 0){
							if(reverse_member > 0 || reverse_both > 0){
								std::string rev_member;
								for(size_t rev = 0; rev< member.size(); rev++){
									char chr= dnastring::complement(member.at(member.size()-1-rev));	
									rev_member +=chr;
								}
								for(size_t mem =0; mem<rev_member.length();mem++){
									std::cout<<rev_member.at(mem);
								}
								std::cout << " " << std::endl;
								for(size_t sl =Sequence_level; sl>0; sl--){
									temp += rev_member.at(rev_member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(rev_member.at(rev_member.length()-1));
							}else{
								for(size_t mem =0; mem<member.length();mem++){
									std::cout<< member.at(mem);
								}
								std::cout << " " << std::endl;
								for(size_t sl =Sequence_level; sl>0; sl--){
									temp += member.at(member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(member.at(member.length()-1));
							}
						}else{
							for(size_t sl=Sequence_level; sl > 0; sl--){
								//reverse complement!
							//	std::cout<<"here!"<<std::endl;
							//	char co_base = dnastring::complement(decodedCenter.at(sl));
							//	temp += co_base;
								char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
								temp +=center_base;
							}
							pattern = temp; 
							base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
						}
//						if(reverse_center ==0 && member.length()==0){
//							for(size_t sl =Sequence_level; sl>0; sl--){
//								temp += decodedCenter.at(decodedCenter.length()-1-sl);
//							}
//							pattern = temp; 
//							base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
//						}
						target = dec.get_target(total);	
						std::cout<< "problematic target: "<< target << std::endl;
						reverse_center = 0;
						reverse_member = 0;
						reverse_both = 0;
					}
				//	target = dec.get_target(total);	
					if(target<flag){
						std::cout<< "back to sequence from an alignment"<<std::endl;
					}
				}
				if(flag+10 <= target && target < flag+15){
				//	unsigned int l = flag+10;
				//	unsigned int h = flag + 15;
					dec.decode(low.at(7),high.at(7));
					wrappers.decode(low.at(7),high.at(7),total);
					std::cout << " End of a sequence! "<<std::endl;
				}else std::cout<<"there is something wrong here!"<<std::endl;
				target = dec.get_target(total);		
			}	
			if(flag+15 <= target && target < flag+20){
				unsigned int l = flag+15;
				unsigned int h = flag + 20;
				dec.decode(l, h);
				wrappers.decode(l,h,total);
				std::cout << " End of an accession! "<<std::endl;
			}else std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
			i = i + 1;
		}
	}
	void write_to_stream(std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &);
	void read_from_stream(std::ifstream&);
	std::string associatedMember(std::string & center , size_t& modificationPattern, size_t& position){
		std::string member = "";
		if(modificationPattern<5){
			char base = dnastring::index_to_base(modificationPattern);
			member += base;	
			std::cout<< "modification"<<std::endl;
		//	std::cout<< "member is "<<member<<std::endl;
			return member;

		}
		if(modificationPattern>=5 && modificationPattern<5+NUM_DELETE){
			std::cout<<"deletion happened!"<<std::endl;
			return member;
		}
		if(modificationPattern>=5+NUM_DELETE && modificationPattern<5+NUM_KEEP+NUM_DELETE){
			size_t length = model.modification_length(modificationPattern);
			std::cout << "center length "<< center.size() << " pattern length: "<< length <<" pattern index : "<< modificationPattern << std::endl;
			for(size_t i = 0 ; i < length;i++){
				assert(position+ i <center.length());
				member += center.at(position+i);
			//	std::cout<< "member is "<<member<<std::endl;
			}
			return member;
		}
		if(modificationPattern>=5+NUM_KEEP+NUM_DELETE){
			char base = dnastring::index_to_base(20-modificationPattern);
			member += base;
			std::cout<< "insertion"<<std::endl;
		//	std::cout<< "member is "<<member<<std::endl;

			return member;

		}
		assert(0);
		return "";
	}
	const std::multimap<size_t, pw_alignment*> & get_alignment(std::map<std::string,std::vector<pw_alignment> > & , size_t);
	void size_calculator( size_t &, std::map<std::string, std::vector<pw_alignment> > &, std::map<std::string, std::string>&);
	private:
	all_data & data;
	T & model;
	wrapper& wrappers;
	std::map<std::string,std::vector<double> > all_the_patterns;
	std::map<std::string, unsigned int> total;
	std::map< size_t, std::vector<unsigned int> > lower_bound;
	std::vector<std::map< size_t, vector<unsigned int> > >upper_bound;
	std::map<std::string, std::vector<unsigned int> >cluster_high;
	std::vector<std::multimap<size_t , pw_alignment *> > AlignmentsFromClustering; //vector---> goes over all the references, size_t --> left of an alignment on the sequence
	std::map<size_t, std::vector<std::string> > acc_of_center;
	std::vector< std::map<std::string, vector<unsigned int> > >cluster_high_partition;//vector ---> partitions , map<> --> center of each partition with their high and low value.
	std::map<std::string,std::string> decoded_centers; // string --->dna sequence, string ---> ref:left(ino estefade nemikonam)
//	std::vector<std::string> centerId;//keep the centers in the order that has been decoded
	std::map<size_t , std::vector<std::string> > partition; // size_t ---> partition number, vector---> centers of a partition
	std::vector<unsigned int>partitionHigh; //high value of each partition.
	std::vector<std::map <std::string, std::string> >decoded_center_in_partition;//string ---> ref:left, string --->dna sequence



};



#endif

