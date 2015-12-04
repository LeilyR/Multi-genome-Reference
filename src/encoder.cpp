#include "encoder.hpp"

#ifndef ENCODER_CPP
#define ENCODER_CPP

	

	template<typename T>
	encoder<T>::encoder( all_data & d, T & a_model, wrapper & wrap): data(d),model(a_model),wrappers(wrap),upper_bound(d.numSequences()),AlignmentsFromClustering(data.numSequences()){

	}
	template<typename T>
	encoder<T>::~encoder(){}
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
	template<typename T>
	void encoder<T>::partition_centers(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::map<std::string, unsigned int> & weight){ //dividing centeres into different partitions
		size_t bit  =13;
		unsigned int length = model.get_powerOfTwo().at(bit);		
		size_t sum_of_weight = 0;
		for(std::map<std::string, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){
			sum_of_weight = sum_of_weight + it->second ;
		}
		std::cout << "total weight "<< sum_of_weight <<std::endl;
		size_t numberOfPartitions = (sum_of_weight/length) + 1;
			std::cout<< "no. of partition:" << numberOfPartitions << std::endl;
		//	std::cout<< "size of al of cluster" << alignmentsOfClusters.size() << std::endl;
		std::vector<std::string> centers;
//		for(size_t k = 0; k < data.numAcc(); k++){
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it= alignmentsOfClusters.begin(); it!= alignmentsOfClusters.end();it++){
			std::string c_id = it->first;
		//	std::vector<std::string> split;
		//	strsep(c_id, ":" , split);
		//	size_t c_ref = atoi(split.at(0).c_str());
		//	size_t c_left = atoi(split.at(1).c_str());
		//	size_t acc = data.accNumber(c_ref);
		//	if(acc == k){
			centers.push_back(c_id);
		//	}
		}
	//	}
		size_t number_of_center = 0;
		for(size_t j = 0 ; j < numberOfPartitions ; j ++){	
			partition.insert(std::make_pair(j,std::vector<std::string>( )));
			size_t center_weight = 0;
			for(size_t i = number_of_center; i < centers.size();i++){
				std::map<std::string, unsigned int >::iterator it1 = weight.find(centers.at(i));
				center_weight +=  it1->second;
				if(center_weight <= length){
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
	void encoder<T>::partitioning_long_centers( std::vector<std::vector<std::string> > & long_centers, std::map<std::string ,std::vector<pw_alignment> > & alignmentsOfClusters ,std::map<std::vector<std::string>, unsigned int> & weight){//it includes both long and original centers 
		size_t bit = 13;
		std::vector<std::vector<std::string> >all_centers;
		unsigned int max_length = model.get_powerOfTwo().at(bit);
		size_t sumOfWeights = 0;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){//original centers are added to all centers
			std::vector<std::string> temp_center;
			temp_center.push_back(it->first);
			all_centers.push_back(temp_center);
		}
		for(size_t i =0; i < long_centers.size(); i++){
			all_centers.push_back(long_centers.at(i));
		}
		for(std::map<std::vector<std::string>, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){// All the weights are summed up.
			sumOfWeights += it->second ;
		}
		std::cout << "sum of weight "<< sumOfWeights <<std::endl;
		size_t numberOfPartitions = (sumOfWeights/max_length) + 1;
		std::cout << numberOfPartitions <<std::endl;
		size_t number_of_center = 0;
		for(size_t j = 0 ; j < numberOfPartitions ; j ++){
			long_center_partition.insert(std::make_pair(j,std::vector< std::vector<std::string> >( )));
			size_t center_weight = 0;
			for(size_t i = number_of_center; i < all_centers.size(); i++){
				if(center_weight <= max_length){
					std::map<size_t , std::vector< std::vector<std::string> > >::iterator it = long_center_partition.find(j);
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
	void encoder<T>::calculate_high_in_partition(std::map<std::string, unsigned int> & weight, std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){
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
				//	std::cout << "h "<< h << "it2 second "<<it2->second << std::endl;
					assert(h+it2->second <= model.get_powerOfTwo().at(bit));
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
	template<typename T>
	void encoder<T>::partition_high(){
		size_t bit = 13;
		unsigned int weight = model.get_powerOfTwo().at(bit);		
		unsigned int high = 0;
		for(size_t i = 0; i < partition.size(); i++){
			high = high + weight/partition.size();
			partitionHigh.push_back(high);
		}
		partitionHigh.at(partitionHigh.size()-1) = weight;
	/*	for(size_t i =0; i < partition.size();i++){
			std::cout << partitionHigh.at(i)<<std::endl;
		}*/
	}
	template<typename T>
	void encoder<T>::calculate_long_center_high(std::vector<std::vector<std::string> > & long_centers,  std::map<std::string ,std::vector<pw_alignment> > & alignmentsOfClusters ,std::map<std::vector<std::string>, unsigned int> & weight){
		std::cout << "calculate long centers "<<std::endl;
		partitioning_long_centers(long_centers,alignmentsOfClusters, weight);//Includes both type of centers
		size_t bit = 13;
		unsigned int partition_high = 0;
		for(size_t i =0; i < long_center_partition.size(); i++){//go through number of partitions
			partition_high = partition_high + (model.get_powerOfTwo().at(bit)/long_center_partition.size());//calculates the high value of the partion.
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
					High = model.get_powerOfTwo().at(bit);
				}
				std::vector<unsigned int>temp;
				temp.push_back(Low);
				temp.push_back(High);
				high.push_back(make_pair(center,temp));
				std::cout << Low << " " << High << std::endl;
			}
			long_center_high.push_back(high);
			HighOfPartition.push_back(partition_high);
		}
	}
	template<typename T>
	void encoder<T>::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//all the alignment that has a certain sequence as at least one of their references.
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
				//	size_t ref1 = p->getreference1();
				//	size_t ref2 = p->getreference2();
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
	//	for(std::multimap<size_t,pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(0).begin();it2 !=AlignmentsFromClustering.at(0).end();it2++){
	//		pw_alignment *p = it2 -> second;
	//		p->print();
		//	std::cout<<"left on the reference: "<< it2->first<<std::endl;
	//	}
	//	outs << char(7);
	}
/*	template<typename T>
	const std::multimap<size_t, pw_alignment*> & encoder<T>::get_alignment(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, size_t seq_id){
		setOfAlignments(alignmentsOfClusters);
		return AlignmentsFromClustering.at(seq_id);
	}*/
	template<typename T>
	void encoder<T>::add_partition_high_to_stream(std::ofstream & outs){
		size_t bit = 32;
		for(size_t i = 0; i < partitionHigh.size(); i++){
			outs<<(char)0;
			std::vector<bool> bit_to_byte(0);
			unsigned int high = partitionHigh.at(i);
			std::cout<< "high_write: "<< high << std::endl;
			for(size_t j = 0 ; j < 32; j++){
				bit_to_byte.push_back(high%2);
				high = high/2;
			}
			for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
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
	template<typename T>	
	void encoder<T>::add_center_to_stream(std::ofstream & outs){
		size_t bit =32;
		outs<< (char)0;
		outs<<partition.size();
		std::cout << "partition size in encoding " << partition.size() << std::endl;
		for(size_t i = 0 ; i < partition.size(); i ++){
			for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				std::cout<< "center_write: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
			}
			for ( std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end(); it ++ ){
				std::string center = it -> first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int cent_ref = atoi(center_parts.at(1).c_str());
				unsigned int cent_left = atoi(center_parts.at(2).c_str());
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
				for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
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
template<typename T>
	void encoder<T>::add_long_center_to_the_stream(std::ofstream& outs){
	size_t bit =32;
		outs<< (char)0;
		outs<<long_center_partition.size();
		std::cout << "add all the centers to the stream"<<std::endl;
		for(size_t i = 0 ; i < long_center_partition.size(); i ++){
			std::cout << "partition "<< i << std::endl;
			std::vector<std::pair<std::vector<std::string>, std::vector<unsigned int> > >center_high = long_center_high.at(i);
			for(size_t j =0; j < center_high.size(); j++){
//			for(std::map<std::vector<std::string>, std::vector<unsigned int> >::iterator it = long_center_high.at(i).begin(); it != long_center_high.at(i).end();it++)
				std::vector<std::string> center = center_high.at(j).first;
				for(size_t k = 0 ; k < center.size(); k++){
					std::cout<< center.at(k)<<" ";
					std::vector<std::string> center_parts;
					strsep(center.at(k), ":" , center_parts);
					unsigned int cent_ref = atoi(center_parts.at(1).c_str());
					unsigned int cent_left = atoi(center_parts.at(2).c_str());
					if(center.size()==1){
						size_t acc_id = data.accNumber(cent_ref);
						outs << (char) 5;
						outs << acc_id;
					}
					outs << (char) 0;
					outs << cent_ref;
					outs<< (char) 0;
					outs<< cent_left;
				}
				outs << (char) 6;
				std::cout << " " << std::endl;				
				std::vector<bool> bit_to_byte(0);
				int high = center_high.at(j).second.at(1);
				std::cout<< "high " << high << std::endl;
				for(size_t i = 0 ; i < bit; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
				size_t counter = 0;
				for(size_t n = 0; n <= bit_to_byte.size()-8; n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
					}
					n= n+7;
					outs<< a;
					counter = counter + 1;
				}
			//	std::cout << "counter " << counter << std::endl;
			}
			outs << (char) 7;
		}
		outs << (char) 8;
		std::cout<< "partition size " <<HighOfPartition.size() << std::endl;
		for(size_t i =0; i < HighOfPartition.size(); i++){
			outs<<(char)0;
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
					a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
				}
				n= n+7;
				outs<< a;
				counter = counter + 1;
			}
			std::cout << "counter1 " << counter << std::endl;
			outs << (char) 7;
		}
		outs << (char) 8;
	}
template<typename T>
	void encoder<T>::arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream & outs, dlib::entropy_encoder_kernel_1 & enc){// after partitioning the clusters (New one!)
	//	add_center_to_stream(outs);//center id, its acc and high
	//	add_partition_high_to_stream(outs);
	//	setOfAlignments(alignmentOfCluster);
		size_t bit = 13;
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
	//	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc.set_stream(outs);
	//	std::cout<< "size of  cluster_high_partition: "<<  cluster_high_partition.size()<<std::endl;
		for(size_t i = 0 ; i < cluster_high_partition.size(); i++){
		//	for(size_t j = 0 ; j < data.numAcc(); j++)
			for(std::map<std::string, std::vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin();it != cluster_high_partition.at(i).end(); it++){
				std::string center = it ->first;
				std::cout << "center " << center << std::endl;
				std::vector<std::string> split;
				strsep(center, ":" , split);
				size_t cent_right = 0;
				size_t cent_ref = atoi(split.at(1).c_str());
				size_t cent_left = atoi(split.at(2).c_str());
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
					}else{
						cent_right = p->getbegin1();
					}				
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
							wrappers.encode(l,h,t);
						//	std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
						std::cout<<"cent length1: "<< cent_right - cent_left << std::endl;
				/*	}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin1();
					//	std::cout<< "id of cent in enc2: " << center <<std::endl;
					//	for(size_t k = cent_right; k >= cent_left;k--){
					//		size_t base = dnastring::base_to_index(seq.at(k));
					//		char co_base = dnastring::complement(seq.at(k));
					//		size_t com_base = dnastring::base_to_index(co_base);	
						//	std::cout<<"base _enc: "<<base <<std::endl;
						//	std::cout<< "com_base_enc: "<< com_base <<std::endl;
					//		unsigned int l = 0;
					//		unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
					//		if(com_base !=0){
					//			l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
					//		}else l = 0;
					//		if(com_base == 4){
					//			h=  model.get_powerOfTwo().at(bit);
					//		}
					//		enc.encode(l,h,t);
						//	std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
					//	std::cout<<"cent length2: "<< cent_right - cent_left << std::endl;
						std::cout<< "id of cent in enc2: " << center <<std::endl;
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));
						//	size_t com_base = dnastring::base_to_index(base);				
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
							wrappers.encode(l,h,t);
						}

					}*/
				}
				if(cent_ref == p->getreference2()&& cent_left == left2){
				//	std::cout<< "encoded center2: "<< it->first << " its al length: "<< p->alignment_length()<<std::endl;				
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
					}else{
						cent_right = p->getbegin2();
					}
						std::cout<< "id of cent in enc3: " << center  << " its left " << left << " its right " << right <<std::endl;
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
							wrappers.encode(l,h,t);
					//		std::cout << "l: " << l << " h: "<< h << " base: " << base << std::endl;
						}
						std::cout<<"cent length3: "<< cent_right - cent_left << std::endl;
				/*	}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin2();
						std::cout<< "id of cent in enc4: " << center <<std::endl;
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc.encode(l,h,t);
							wrappers.encode(l,h,t);
						}
					}*/
				}
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit)+5;					
				enc.encode(l1,h1,t); 
				wrappers.encode(l1,h1,t);
				std::cout << "end of a center "<<std::endl;
			}
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit)+10;					
			enc.encode(l1,h1,t); 
			wrappers.encode(l1,h1,t);
			std::cout<< "partition: " << i <<std::endl;
		}
	}
	template<typename T>
	void encoder<T>::encoding_centers_with_long_center(std::map<std::string,std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream & outs, dlib::entropy_encoder_kernel_1 & enc){
		setOfAlignments(alignmentOfCluster);
		std::ofstream save("encode" , std::ofstream::binary);
		size_t bit = 13;
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
		enc.set_stream(outs);
		for(size_t i = 0 ; i < long_center_high.size(); i++){
			std::vector<std::pair<std::vector<std::string>, std::vector<unsigned int> > >center_high = long_center_high.at(i);
			for(size_t j =0; j < center_high.size();j++){
				if(center_high.at(j).first.size()==1){
					std::string center = center_high.at(j).first.at(0);
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
									h=  model.get_powerOfTwo().at(bit);
								}
								save<< "l" << l << "h"<< h ;								
							//	std::cout<< "l " << l << " h "<< h <<std::endl;
								enc.encode(l,h,t);
							//	wrappers.encode(l,h,t);
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
									h=  model.get_powerOfTwo().at(bit);
								}
								enc.encode(l,h,t);
							//	wrappers.encode(l,h,t);
								save<< "l" << l << "h"<< h ;								
							//	std::cout<< "l " << l << " h "<< h <<std::endl;

							}
						}
					}
					if(cent_ref == p->getreference2()&& cent_left == left2){
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
									h=  model.get_powerOfTwo().at(bit);
								}
								enc.encode(l,h,t);
							//	wrappers.encode(l,h,t);
								save<< "l" << l << "h"<< h ;								
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
									h=  model.get_powerOfTwo().at(bit);
								}
								enc.encode(l,h,t);
								save<< "l" << l << "h"<< h ;
							//	std::cout<< "l " << l << " h "<< h <<std::endl;
							//	wrappers.encode(l,h,t);


							}
						}
					}
					std::cout<<"end of a center!"<<std::endl;
					unsigned int l1	= model.get_powerOfTwo().at(bit);
					unsigned int h1 = model.get_powerOfTwo().at(bit)+5;					
					unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
					enc.encode(l1,h1,t1);
				//	wrappers.encode(l1,h1,t);
 
				}
			}
			std::cout<< "end of a partition"<<std::endl;
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit)+10;					
			unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
			enc.encode(l1,h1,t1); 
		//	wrappers.encode(l1,h1,t);
		}
		save.close();
	}
	template<typename T>
	void encoder<T>::encode_flags_of_reverse_parts(pw_alignment & p, unsigned int & cent_ref, unsigned int & cent_left, dlib::entropy_encoder_kernel_1 & enc){
		size_t bit = 13;
	//	unsigned int total = model.get_powerOfTwo().at(bit)+20;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		size_t left_1; 
		size_t left_2;
		size_t right_1;
		size_t right_2;
		p.get_lr1(left_1,right_1);
		p.get_lr2(left_2,right_2);		
		std::cout<< "beginning of reverse function"<<std::endl;
		if((cent_ref == p.getreference1() && cent_left == p.getend1()&& left_2==p.getbegin2())||(cent_ref == p.getreference2()&& cent_left ==p.getend2()&&left_1 == p.getbegin1())){
			std::cout<<"reverse center!"<<std::endl;
			//center is reverse
			unsigned int l1 =  model.get_powerOfTwo().at(bit);
			unsigned int h1 = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getbegin1()&& left_2==p.getend2())||(cent_ref == p.getreference2()&& cent_left ==p.getbegin2()&&left_1 == p.getend1())){
			std::cout<<"other ref is reverse!"<<std::endl;
			//other ref is reverse
			unsigned int l1 =  model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
			unsigned int h1 = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getend1()&& left_2==p.getend2())||(cent_ref == p.getreference2() && cent_left ==p.getend2()&&left_1 == p.getend1())){
			std::cout<<"reverse_both!"<<std::endl;
			//both are reverse
			unsigned int l1 =  model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
			unsigned int h1 = 2*model.get_powerOfTwo().at(bit);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		std::cout<< "end of reverse function"<<std::endl;
	}
	template<typename T>
	void encoder<T>::encode_optimized_flags_of_reverse_parts(pw_alignment & p, unsigned int & cent_ref ,unsigned int & cent_left, dlib::entropy_encoder_kernel_1 & enc){
		size_t bit = 13;
//		unsigned int total = model.get_powerOfTwo().at(bit)+high_of_flags.at(7);
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		size_t left_1; 
		size_t left_2;
		size_t right_1;
		size_t right_2;
		p.get_lr1(left_1,right_1);
		p.get_lr2(left_2,right_2);		
		std::cout<< "beginning of reverse function"<<std::endl;
		if((cent_ref == p.getreference1() && cent_left == p.getend1()&& left_2==p.getbegin2())||(cent_ref == p.getreference2()&& cent_left ==p.getend2()&&left_1 == p.getbegin1())){
			std::cout<<"reverse center!"<<std::endl;
			//center is reverse
			unsigned int l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(0);
			unsigned int h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(1);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getbegin1()&& left_2==p.getend2())||(cent_ref == p.getreference2()&& cent_left ==p.getbegin2()&&left_1 == p.getend1())){
			std::cout<<"other ref is reverse!"<<std::endl;
			//other ref is reverse
			unsigned int l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(1);
			unsigned int h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(2);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		if((cent_ref == p.getreference1()&& cent_left == p.getend1()&& left_2==p.getend2())||(cent_ref == p.getreference2() && cent_left ==p.getend2()&&left_1 == p.getend1())){
			std::cout<<"reverse_both!"<<std::endl;
			//both are reverse
			unsigned int l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(2);
			unsigned int h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(3);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
		}
		std::cout<< "end of reverse function with optimized flags"<<std::endl;
	}
	template<typename T>
	void encoder<T>:: al_encoding(std::map<std::string, unsigned int> & weight, std::map<std::string, std::string > & cluster_members, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &outs,dlib::entropy_encoder_kernel_1 & enc ){
		calculate_high_in_partition(weight,alignmentOfCluster);
		encoding_functor functor(data,&model,wrappers,enc);
		add_center_to_stream(outs);
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		arithmetic_encoding_centers(alignmentOfCluster,outs,enc);
		size_t bit = 13;
//		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
//		size_t i = 1;
//		std::cout<<"number of acc:" <<data.numAcc()<<std::endl;
		size_t al_number =0;
		for(size_t i = 0; i< data.numAcc(); i++){
//			std::cout<< " i "<< i <<std::endl;
//			std::cout << "number of seq in an acc : "<< data.getAcc(i).size() <<std::endl;
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
//				std::cout<< "sequence length is: "<< sequence.length()<< " sequence id: "<< sequenceId <<std::endl;
				if(sequenceId == 2){
					for (std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).begin() ; it != AlignmentsFromClustering.at(sequenceId).end(); it++){
						std::cout << "it ->first "<< it->first <<std::endl;
						pw_alignment *p = it->second;
						p->print();
						
					}
				}
				al_number +=AlignmentsFromClustering.at(sequenceId).size();
				std::cout << "no. of als on a seq "<< al_number <<std::endl;
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
						std::cout<< "al position: "<< n << std::endl;
						size_t left_1; 
						size_t left_2;
						size_t right_1;
						size_t right_2;
						pw_alignment *p = it->second;
//						std::cout<< "al length: "<< p->alignment_length()<<std::endl;
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);		
//						std::cout<<" l1 "<<left_1 << " l2 " << left_2 << " r1 "<< right_1 << " r2 "<<right_2 <<std::endl;
//						std::cout << "begin1 "<< p->getbegin1() << " begin2 "<< p->getbegin2() << " end1 "<< p->getend1() << " end2 " << p->getend2() << std::endl;
						//A fixed flag before encoding a center:
//						std::cout<<"enc  bal flag"<<std::endl;
						unsigned int l1 = model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/4);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::stringstream mem;
						mem << sequenceId << ":" << n;
						std::cout << mem.str() << std::endl;
						std::map<std::string, std::string>::iterator cl = cluster_members.find(mem.str());
						assert(cl != cluster_members.end());
						std::string center = cl->second;
//						std::cout<< "center from cluster member: "<< center << std::endl;
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(1).c_str());
						unsigned int cent_left = atoi(center_parts.at(2).c_str());
						unsigned int center_l;
						unsigned int center_h;
					//	size_t part;
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
//								std::cout << "par low: "<< par_low << " par high: "<< par_high << std::endl;									
								center_l = it1 ->second.at(0);
								center_h = it1->second.at(1);
							//	part = j;
								break;
							}
						}
//						std::cout<< "center low "<< center_l << " center high " << center_h<<"center left : "<< cent_left <<std::endl;
						enc.encode(center_l, center_h, total);
						wrappers.encode(center_l,center_h,total);
						if(cent_ref == sequenceId && cent_left == n){
//							std::cout<< " center is on the ref "<< std::endl;	
							if(cent_left == p->getend1()||cent_left ==p->getend2()){
//								std::cout<<"reverse center!"<<std::endl;
							//	p->print();
							//	unsigned int l1 =  model.get_powerOfTwo().at(bit) + 15;
							//	unsigned int h1 = model.get_powerOfTwo().at(bit) + 20;
							//	enc.encode(l1,h1,total);
							//	wrappers.encode(l1,h1,total);
							}else{
							}					
						}else{
							if((cent_ref == p->getreference1() && cent_left == p->getend1()&& left_2==p->getbegin2())||(cent_ref == p->getreference2()&& cent_left ==p->getend2()&&left_1 == p->getbegin1())){
//								std::cout<<"reverse center!"<<std::endl;
								//center is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 0;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);

							}
							if((cent_ref == p->getreference1()&& cent_left == p->getbegin1()&& left_2==p->getend2())||(cent_ref == p->getreference2()&& cent_left ==p->getbegin2()&&left_1 == p->getend1())){
							//	p->print();
//								std::cout<<"other ref is reverse!"<<std::endl;
								//other ref is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
								unsigned int h1 = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}
							if((cent_ref == p->getreference1()&& cent_left == p->getend1()&& left_2==p->getend2())||(cent_ref == p->getreference2() && cent_left ==p->getend2()&&left_1 == p->getend1())){
								p->print();
//								std::cout<<"reverse_both!"<<std::endl;
								//both are reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
								unsigned int h1 = 2*model.get_powerOfTwo().at(bit);
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}

							model.get_encoded_member(*p,cent_ref,cent_left,functor,outs);
//							std::cout<< " center is not on the ref "<< std::endl;
						}//end of each al
						l1 = model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/4);
						h1 = model.get_powerOfTwo().at(bit) +(model.get_powerOfTwo().at(bit)/2);
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
				unsigned int l1	= model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/2);
				unsigned int h1 = model.get_powerOfTwo().at(bit) +(3*model.get_powerOfTwo().at(bit)/4);
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout<<"end of a seq in enc. "<<std::endl;
			}//end of all sequences of an accsseion
			unsigned int l1 = model.get_powerOfTwo().at(bit)+(3*model.get_powerOfTwo().at(bit)/4);
			unsigned int h1 = 2*model.get_powerOfTwo().at(bit);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout<<"end of an acc"<<std::endl;
		}
		std::cout<< "encoding is finished!"<<std::endl;
	}
	template<typename T>
	void encoder<T>::al_encode_with_long_center(std::vector<std::map<size_t , std::string> >& centerOnSequence, std::map<std::vector<std::string> , unsigned int> & long_center_weight, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::vector<std::map<size_t, std::vector<std::string> > >& centersPositionOnASeq ,std::map<std::string, std::string > & cluster_members, std::vector<std::vector<std::string> > & long_centers,std::ofstream & outs ,dlib::entropy_encoder_kernel_1 & enc){
		calculate_long_center_high(long_centers,alignmentOfCluster ,long_center_weight);//Calculates low and high value of all the centers.
		encoding_functor functor(data,&model,wrappers,enc);
		add_long_center_to_the_stream(outs);//Adds all the centers, their accessions and high values to the stream.
		encoding_centers_with_long_center(alignmentOfCluster,outs,enc);
		size_t bit = 13;
//		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		std::cout<<"number of acc:" <<data.numAcc()<<std::endl;
		for(size_t i = 0; i< data.numAcc(); i++){
			std::cout<< " accession "<< i <<std::endl;
			std::cout << "number of seq in this acc : "<< data.getAcc(i).size() <<std::endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
			/*	if(sequenceId == 0){
					std::cout<< "all the centers on seq 0"<<std::endl;
					for(std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(0).begin(); it1 != centerOnSequence.at(0).end(); it1++){
						std::cout<< it1->second << " , " ;
						std::cout << it1->first << " , ";

					}
					std::cout<< "" <<std::endl;
				}*/
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){
						//This flag represents the beginning of an alignment
						unsigned int l1	= model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout<< "al position: "<< n << std::endl;
					//	pw_alignment p1 = *it->second;
						std::vector<std::string> current_center;
						std::map<size_t , std::vector<std::string> >::iterator it1=centersPositionOnASeq.at(sequenceId).find(n);
						if(it1 != centersPositionOnASeq.at(sequenceId).end()){//If there is a long center in that position
							current_center = it1->second;
						}else{
							std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(sequenceId).find(n);
							current_center.push_back(it1->second);
						} 
						size_t partition = 0; 
						for(size_t j = 0; j < long_center_partition.size(); j++){
							std::map<size_t , std::vector<std::vector<std::string> > >::iterator par = long_center_partition.find(j);
							for(size_t k =0; k < par->second.size();k++){
								if(par->second.at(k) == current_center){
									partition = j;
								break;
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
						enc.encode(low,high,total);// center itself is encoded!
						wrappers.encode(low,high,total);
						unsigned int cent_ref;
						unsigned int cent_left;
						size_t gap = 0;
						size_t length = n;//length of alignments and gaps are added to it.
						std::cout << "current center size "<< current_center.size() << std::endl;
						for(size_t j =0; j < current_center.size();j++){//find alignment of each center in the long center
							std::string final_context;
							pw_alignment p;
							size_t l_1,r_1,l_2,r_2;
							//break the center here! check its ref and left
							std::string center = current_center.at(j);
							std::cout<< "center: "<< center << std::endl;
							std::vector<std::string> center_parts;
							strsep(center, ":" , center_parts);
							cent_ref = atoi(center_parts.at(1).c_str());
							cent_left = atoi(center_parts.at(2).c_str());
							std::cout<<"length "<< length << std::endl;
							bool center_on_ref = false;
							if(cent_ref == sequenceId && cent_left == length){
								std::cout<<"center is on ref"<<std::endl;
								std::map<std::string, std::vector<pw_alignment> >::iterator find_al = alignmentOfCluster.find(current_center.at(j));
								pw_alignment p1 = find_al->second.at(0);
								p1.get_lr1(l_1,r_1);
								p1.get_lr2(l_2,r_2);	
								if(cent_ref == p1.getreference1() && cent_left == l_1){
									length = r_1+1;
								}else{
									length = r_2+1;
								}
								center_on_ref = true;
							}else{
								std::multimap<size_t, pw_alignment* >::iterator al =AlignmentsFromClustering.at(sequenceId).find(length);
								p = *al->second;
								p.get_lr1(l_1,r_1);
								p.get_lr2(l_2,r_2);	
								if(sequenceId == p.getreference1() && length == l_1){
									length = r_1+1;
								}else{
									length = r_2+1;
								}
								//encoding modification!
								encode_flags_of_reverse_parts(p,cent_ref,cent_left,enc);
								model.get_encoded_member_long_center(p,cent_ref,cent_left,final_context,functor,outs);//final_context is needed for the gap encoding
							}
							if(j != current_center.size()-1){//Find the gap between two successive centers!
								std::cout<< "finding gap between two centers "<<std::endl;
								pw_alignment p1;
								size_t left1,right1,left2,right2;
								std::map<std::string, std::vector<pw_alignment> >::iterator find_al = alignmentOfCluster.find(current_center.at(j+1));
								for(size_t m =0; m < find_al->second.size();m++){
									p1 = find_al->second.at(m);
									p1.get_lr1(left1,right1);
									p1.get_lr2(left2,right2);	
									if((p.getreference1()== sequenceId && left1 >= length && left1 < length + 5)|| (p.getreference2()== sequenceId && left2 >= length && left2 < length +5)){
										for(std::multimap<size_t, pw_alignment* >::iterator pos =AlignmentsFromClustering.at(sequenceId).begin(); pos !=AlignmentsFromClustering.at(sequenceId).end();pos++){
											if(pos->second == &p1){
												gap = pos->first - length;
												break;
											}
										}
										break;
									}
								}
								std::cout << "final context "<< final_context << std::endl;
								std::cout << "gap "<< gap<<std::endl;
								for(size_t m = 0; m < gap ; m++){//TODO it is hard to see the errors in this part if there is any. because there is no gap in this data set
									unsigned int l1;
									unsigned int h1;
									char seq_base = sequence.at(length+m);
									char s1,s2;
									if(center_on_ref == false){
										p.alignment_col(p.alignment_length()-1, s1, s2);
										if(cent_ref == p.getreference1() && cent_left == l_1){
											model.get_insertion_high_between_centers(sequenceId , seq_base , s1 , cent_ref, final_context ,h1, l1);
											enc.encode(l1,h1,total);
											wrappers.encode(l1,h1,total);
										}else{
											model.get_insertion_high_between_centers(sequenceId , seq_base , s2 , cent_ref, final_context,h1, l1);
											enc.encode(l1,h1,total);
											wrappers.encode(l1,h1,total);
										}
									}else{
										model.get_insertion_high_between_centers(sequenceId , seq_base , seq_base , cent_ref, final_context ,h1, l1);
										enc.encode(l1,h1,total);
										wrappers.encode(l1,h1,total);
									}
								}
								l1 = model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/4);
								h1 = model.get_powerOfTwo().at(bit) +(model.get_powerOfTwo().at(bit)/2);
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
								length = length + gap;
								gap = 0;
							}	
						}
						n = length-1;
						std::cout<< "end position of long center " << n <<std::endl;
						//Flag shows end of an alignment
						l1 = model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/4);
						h1 = model.get_powerOfTwo().at(bit) +(model.get_powerOfTwo().at(bit)/2);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
					}else{ //If there is no alignment in that position
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
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
				}//End of a sequence
				unsigned int l1	= model.get_powerOfTwo().at(bit)+(model.get_powerOfTwo().at(bit)/2);
				unsigned int h1 = model.get_powerOfTwo().at(bit) +(3*model.get_powerOfTwo().at(bit)/4);
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout<<"end of a seq in enc. "<<std::endl;

			}//end of all sequences of an accsseion
			unsigned int l1 = model.get_powerOfTwo().at(bit)+(3*model.get_powerOfTwo().at(bit)/4);
			unsigned int h1 = 2*model.get_powerOfTwo().at(bit);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout<<"end of an acc"<<std::endl;
		}

	}
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
	void encoder<T>::write_to_stream( std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster ,std::ofstream & outs){		
		add_center_to_stream(outs);
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
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
	template<typename T>
	void encoder<T>::count_flags(std::map<std::vector<std::string>, std::vector<pw_alignment> >& new_centers, std::map<std::string, std::vector<pw_alignment> >& al_in_cluster,std::vector<size_t> &flagsCounts ){
		size_t al_begin = 0;
		size_t end_of_each_al_part = 0;
		size_t end_of_al = 0;
		size_t end_of_seq = data.numSequences();
		size_t end_of_acc = data.numAcc();
		size_t reverse_center = 0; 
		size_t reverse_other_ref = 0;
		size_t reverse_both = 0;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end(); it ++){
			if(it->first.size() == 1){
				for(size_t i = 0 ; i < it->second.size(); i++){
					al_begin ++;
					end_of_al ++;
				}

			}else{	
				for(size_t i =0; i < it->first.size()-1; i++){
					end_of_each_al_part++;
				}
				for(size_t i = 0; i < it->second.size();i++){
					al_begin ++;
					end_of_al ++;
				}				
			}
		}
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = al_in_cluster.begin(); it != al_in_cluster.end(); it ++){
			std::string center = it->first;
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int cent_ref = atoi(center_parts.at(0).c_str());
			unsigned int cent_left = atoi(center_parts.at(1).c_str());
			for(size_t i = 0 ; i < it->second.size(); i++){
				al_begin ++;
				pw_alignment p = it->second.at(i);
				size_t left_1; 
				size_t left_2;
				size_t right_1;
				size_t right_2;
				p.get_lr1(left_1,right_1);
				p.get_lr2(left_2,right_2);
				if((cent_ref == p.getreference1() && cent_left == p.getend1()&& left_2==p.getbegin2())||(cent_ref == p.getreference2()&& cent_left ==p.getend2()&&left_1 == p.getbegin1())){
					//center is reverse
					reverse_center ++;
				}
				if((cent_ref == p.getreference1()&& cent_left == p.getbegin1()&& left_2==p.getend2())||(cent_ref == p.getreference2()&& cent_left ==p.getbegin2()&&left_1 == p.getend1())){
					//other ref is reverse
					reverse_other_ref ++;
				}
				if((cent_ref == p.getreference1()&& cent_left == p.getend1()&& left_2==p.getend2())||(cent_ref == p.getreference2() && cent_left ==p.getend2()&&left_1 == p.getend1())){
					//both are reverse
					reverse_both ++;
				}
			}
		}
		flagsCounts.push_back(al_begin);
		flagsCounts.push_back(reverse_center);
		flagsCounts.push_back(reverse_other_ref);
		flagsCounts.push_back(reverse_both);
		flagsCounts.push_back(end_of_each_al_part);
		flagsCounts.push_back(end_of_al);
		flagsCounts.push_back(end_of_seq);
		flagsCounts.push_back(end_of_acc); 
		std::cout << " flags counts " << std::endl;
		for(size_t i =0; i < flagsCounts.size(); i++){
			std::cout << flagsCounts.at(i)<<std::endl;
		}
	}
	template<typename T>
	void encoder<T>::weight_flags(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, std::vector<pw_alignment> > & al_in_cluster, std::ofstream & outs){//TODO One way could be calculate the weight for two diffente data set and fix them based on what we got from those two sets. The other way could be saving them on the stream. Should try and see which one returns better result!
		std::vector<size_t> counts;
		size_t bit = 13;
		count_flags(new_centers, al_in_cluster,counts);
		size_t sum = 0;
		for(size_t i =0; i < counts.size();i++){
			sum += counts.at(i);
		}
		for(size_t i =0; i < counts.size(); i++){
			counts.at(i)= counts.at(i)*(model.get_powerOfTwo().at(bit))/sum;
			if(counts.at(i)== 0){
				counts.at(i)=1;
			}
		}
		std::cout << " weight " << std::endl;
		for(size_t i =0; i < counts.size(); i++){
			std::cout << counts.at(i)<<std::endl;
		}
		unsigned int h = 0;
		for(size_t i =0; i< counts.size();i++){
			h +=counts.at(i);
			std::cout << "h: "<< h<<std::endl;
			high_of_flags.push_back(h);
		}
		high_of_flags.at(high_of_flags.size()-1) = model.get_powerOfTwo().at(bit); //It may reduce the efficiency of picking the optimal high value
		add_flag_to_stream(outs);
	}
	template<typename T>
	void encoder<T>::add_flag_to_stream(std::ofstream & outs){
		size_t bit = 16;
		for(size_t i = 0; i < high_of_flags.size(); i++){
			outs << char(0);
			unsigned int high =high_of_flags.at(i);
			std::cout << "high "<<high <<std::endl;
			vector<bool> bit_to_byte;
			for(size_t j = 0 ; j < bit; j++){
				bit_to_byte.push_back(high%2);
				high = high/2;
			}
			std::cout << "  bit_to_byte size " << bit_to_byte.size() <<std::endl;
			for(size_t j =0; j < bit_to_byte.size();j++){
				std::cout << bit_to_byte.at(j) << " ";
			}
			std::cout << " " <<std::endl;
			for(size_t j = 0 ; j<= bit_to_byte.size()-8; j++){
				std::cout << "j "<< j <<std::endl;
				unsigned char a = 0;
				for(size_t k =j; k < j + 8; k++){
					a += model.get_powerOfTwo().at(k-j)*bit_to_byte.at(k);
				}
				outs << a;
				j = j + 7;
				std::cout << "a "<<int(a)<<std::endl;
			}
			outs << char(7);
		}
		outs << char(8);
	}
	template<typename T>
	void encoder<T>::al_encode_long_center_optimized_flags(std::vector<std::map<size_t , std::string> >& centerOnSequence, std::map<std::vector<std::string> , unsigned int> & long_center_weight, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::vector<std::map<size_t, std::vector<std::string> > >& centersPositionOnASeq ,std::map<std::string, std::string > & cluster_members, std::vector<std::vector<std::string> > & long_centers,std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers,std::ofstream & outs ,dlib::entropy_encoder_kernel_1 & enc){
		calculate_long_center_high(long_centers,alignmentOfCluster ,long_center_weight);//Calculates low and high value of all the centers.
		encoding_functor functor(data,&model,wrappers,enc);
		add_long_center_to_the_stream(outs);//Adds all the centers, their accessions and high values to the stream.
		weight_flags(new_centers,alignmentOfCluster,outs);
		encoding_centers_with_long_center(alignmentOfCluster,outs,enc);
		size_t bit = 13;
	//	unsigned int total = model.get_powerOfTwo().at(bit)+high_of_flags.at(7);
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		std::cout << "total "<< total <<std::endl;		
		std::cout<<"number of acc:" <<data.numAcc()<<std::endl;
		size_t i =0;
		size_t k = 0;
		for(size_t i = 0; i< data.numAcc(); i++){
			std::cout<< " accession "<< i <<std::endl;
			std::cout << "number of seq in this acc : "<< data.getAcc(i).size() <<std::endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
				if(sequenceId == 0){
					std::cout<< "all the centers on seq 0"<<std::endl;
					for(std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(0).begin(); it1 != centerOnSequence.at(0).end(); it1++){
						std::cout<< it1->second << " , " ;
						std::cout << it1->first << " , ";

					}
					std::cout<< "" <<std::endl;
				}
				for(size_t n= 0; n < sequence.length(); n++){
					std::multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){
						//This flag represents the beginning of an alignment
						unsigned int l1	= model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit) + high_of_flags.at(0);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						std::cout<< "al position: "<< n << std::endl;
					//	pw_alignment p1 = *it->second;
						std::vector<std::string> current_center;
						std::map<size_t , std::vector<std::string> >::iterator it1=centersPositionOnASeq.at(sequenceId).find(n);
						if(it1 != centersPositionOnASeq.at(sequenceId).end()){//If there is a long center in that position
							current_center = it1->second;
						}else{
							std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(sequenceId).find(n);
							current_center.push_back(it1->second);
						} 
						size_t partition = 0; 
						for(size_t j = 0; j < long_center_partition.size(); j++){
							std::map<size_t , std::vector<std::vector<std::string> > >::iterator par = long_center_partition.find(j);
							for(size_t k =0; k < par->second.size();k++){
								if(par->second.at(k) == current_center){
									partition = j;
									break;
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
						enc.encode(low,high,total);// center itself is encoded!
						wrappers.encode(low,high,total);
						unsigned int cent_ref;
						unsigned int cent_left;
						size_t gap = 0;
						size_t length = n;//length of alignments and gaps are added to it.
						std::cout << "current center size "<< current_center.size() << std::endl;
						for(size_t j =0; j < current_center.size();j++){//find alignment of each center in the long center
							std::string final_context;
							pw_alignment p;
							size_t l_1,r_1,l_2,r_2;
							//break the center here! check its ref and left
							std::string center = current_center.at(j);
							std::cout<< "center: "<< center << std::endl;
							std::vector<std::string> center_parts;
							strsep(center, ":" , center_parts);
							cent_ref = atoi(center_parts.at(1).c_str());
							cent_left = atoi(center_parts.at(2).c_str());
							std::cout<<"length "<< length << std::endl;
							bool center_on_ref = false;
							if(cent_ref == sequenceId && cent_left == length){
								std::cout<<"center is on ref"<<std::endl;
								std::map<std::string, std::vector<pw_alignment> >::iterator find_al = alignmentOfCluster.find(current_center.at(j));
								pw_alignment p1 = find_al->second.at(0);
								p1.get_lr1(l_1,r_1);
								p1.get_lr2(l_2,r_2);	
								if(cent_ref == p1.getreference1() && cent_left == l_1){
									length = r_1+1;
								}else{
									length = r_2+1;
								}
								center_on_ref = true;
							}else{
								std::multimap<size_t, pw_alignment* >::iterator al =AlignmentsFromClustering.at(sequenceId).find(length);
								p = *al->second;
								p.get_lr1(l_1,r_1);
								p.get_lr2(l_2,r_2);	
								if(sequenceId == p.getreference1() && length == l_1){
									length = r_1+1;
								}else{
									length = r_2+1;
								}
								//encoding modification!
								encode_optimized_flags_of_reverse_parts(p,cent_ref,cent_left,enc);
								model.get_encoded_member_long_center(p,cent_ref,cent_left,final_context,functor,outs);//final_context is needed for the gap encoding
							}
							if(j != current_center.size()-1){//Find the gap between two successive centers!
								std::cout<< "finding gap between two centers "<<std::endl;
								pw_alignment p1;
								size_t left1,right1,left2,right2;
								std::map<std::string, std::vector<pw_alignment> >::iterator find_al = alignmentOfCluster.find(current_center.at(j+1));
								for(size_t m =0; m < find_al->second.size();m++){
									p1 = find_al->second.at(m);
									p1.get_lr1(left1,right1);
									p1.get_lr2(left2,right2);	
									if((p.getreference1()== sequenceId && left1 >= length && left1 < length + 5)|| (p.getreference2()== sequenceId && left2 >= length && left2 < length +5)){
										for(std::multimap<size_t, pw_alignment* >::iterator pos =AlignmentsFromClustering.at(sequenceId).begin(); pos !=AlignmentsFromClustering.at(sequenceId).end();pos++){
											if(pos->second == &p1){
												gap = pos->first - length;
												break;
											}
										}
										break;
									}
								}
								std::cout << "final context "<< final_context << std::endl;
								std::cout << "gap "<< gap<<std::endl;
								for(size_t m = 0; m < gap ; m++){//TODO it is hard to see the errors in this part if there is any. because there is no gap in this data set
									unsigned int l1;
									unsigned int h1;
									char seq_base = sequence.at(length+m);
									char s1,s2;
									if(center_on_ref == false){
										p.alignment_col(p.alignment_length()-1, s1, s2);
										if(cent_ref == p.getreference1() && cent_left == l_1){
											model.get_insertion_high_between_centers(sequenceId , seq_base , s1 , cent_ref, final_context ,h1, l1);
											enc.encode(l1,h1,total);
											wrappers.encode(l1,h1,total);
										}else{
											model.get_insertion_high_between_centers(sequenceId , seq_base , s2 , cent_ref, final_context,h1, l1);
											enc.encode(l1,h1,total);
											wrappers.encode(l1,h1,total);
										}
									}else{
										model.get_insertion_high_between_centers(sequenceId , seq_base , seq_base , cent_ref, final_context ,h1, l1);
										enc.encode(l1,h1,total);
										wrappers.encode(l1,h1,total);
									}
								}
								l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(3);
								h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(4);
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
								length = length + gap;
								gap = 0;
							}	
						}
						n = length-1;
						std::cout<< "end position of long center " << n <<std::endl;
						//Flag shows end of an alignment
						l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(4);
						h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(5);
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
					}else{ //If there is no alignment in that position
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
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
				}//End of a sequence
				unsigned int l1	= model.get_powerOfTwo().at(bit)+high_of_flags.at(5);
				unsigned int h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(6);
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				std::cout<<"end of a seq in enc. "<<std::endl;

			}//end of all sequences of an accsseion
			unsigned int l1 = model.get_powerOfTwo().at(bit)+high_of_flags.at(6);
			unsigned int h1 = model.get_powerOfTwo().at(bit) +high_of_flags.at(7);
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			std::cout<<"end of an acc"<<std::endl;
		}

	}
	template<typename T>
	decoder<T>::decoder(all_data & d, T & a_model, decoding_wrapper & wrap):data(d), model(a_model),wrappers(wrap){}
	template<typename T>
	decoder<T>::~decoder(){}
template<typename T>
	std::string decoder<T>::associatedMember(std::string & center, size_t & modificationPattern, size_t & position){
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
	void decoder<T>::set_partition_high_from_stream(std::ifstream & in){
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
				std::cout<< "high_read: "<< high_value<<std::endl;
				c=in.get();
			}
			c=in.get();
		}
	}
	template<typename T>
	void decoder<T>::set_center_from_stream(std::ifstream & in){
		size_t bit = 32;
		size_t partition_size = 0;
		char c ; 
		c = in.get();
		in >> partition_size;
		std::vector<std::pair<std::string, vector<unsigned int> > >value;
		for(size_t i = 0; i < partition_size; i++){
			cluster_high_partition.push_back(value);
		}
		unsigned char h;
		c = in.get();
		size_t i = 0;
		while(c != 8){
			unsigned int low = 0;
			while(c != 7){
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
				std::vector<unsigned int> temp(2,0) ;
				std::pair<std::string, std::vector<unsigned int> >Centerhigh;
				Centerhigh.second = temp;
				Centerhigh.first = center_id.str();				
//				std::pair<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).find(center_id.str());
//				if(it2 == cluster_high_partition.at(i).end()){
				std::cout << "center: " << center_id.str() <<std::endl;
//				cluster_high_partition.at(i).push_back(std::make_pair(center_id.str(),std::vector<unsigned int>(2,0)));
//				it2 = cluster_high_partition.at(i).find(center_id.str());	
//				}
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
					for(size_t j = 0; j < binary_high_value.size();j++){
							high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(j);
					}
		/*		for(size_t i = 0; i < binary_high_value.size();i++){
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(j-i);
					}
					i=i+bit-1;
				}*/
			//	std::cout << "high value: "<< high_value << std::endl;
			//		it2 -> second.at(1)=high_value;
			//		it2 ->second.at(0)=low;
			//		low= it2->second.at(1);
					Centerhigh.second.at(0) = low;
					Centerhigh.second.at(1) = high_value;
					low = Centerhigh.second.at(1);
					c=in.get();
					cluster_high_partition.at(i).push_back(Centerhigh);
			}
			c = in.get();
		/*	for(std::map<std::string, std::vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				std::cout<< "center_read_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<std::endl;	
			}*/
			i = i + 1;
			std::cout<< "i: "<< i << std::endl;
		}
	//	std::cout<< "END!"<<std::endl;
	/*	for(std::map<size_t, std::vector<std::string> >::iterator it = acc_of_center.begin();it!=acc_of_center.end();it++){
			std::cout<<"acc std::set in the map: "<<it->first<<std::endl;
		}*/
	/*	std::cout << "cluster high parttion "<< cluster_high_partition.size() << std::endl;
		for(size_t i =0; i <  cluster_high_partition.size(); i++){
			for(std::map<std::string, vector<unsigned int> >::iterator it =  cluster_high_partition.at(i).begin(); it !=  cluster_high_partition.at(i).end(); it++){
				std::cout<< "center" <<it->first <<std::endl;
				std::cout << it->second.at(0) << " " << it->second.at(1) <<std::endl;
			}
		}*/
	}
	template<typename T>
	void decoder<T>::set_acc_from_stream(std::ifstream & in){
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
template<typename T>
	void decoder<T>::set_long_centers_from_stream(std::ifstream & in){
		size_t bit = 32;
		size_t partition_size = 0;
		char c ; 
		c = in.get();
		in >> partition_size;
		std::cout << "par size: " << partition_size << std::endl;
		std::vector<std::pair<std::vector<std::string>, std::vector<unsigned int> > >value;
		for(size_t i = 0; i < partition_size; i++){
			long_centers_high_value.push_back(value);
		}
		std::cout << "long center partition size "<<long_centers_high_value.size()<<std::endl;
		unsigned char h;
		c = in.get();
		size_t i = 0;
		while(c != 8){
			size_t low = 0;
			std::vector<std::pair<std::vector<std::string> , std::vector<unsigned int> > >longCenterHigh;
			while(c != 7){//This loop goes over the partitions
				std::vector<std::string> center;
				std::string cent;
				bool short_center = false;
				size_t accession;
				while(c != 6){//each member of a partition
					if(c == 5){
						unsigned int cent_acc;
						in>> cent_acc;
						std::map<size_t , std::vector<std::string> >::iterator it = acc_of_center.find(cent_acc);
						if(it == acc_of_center.end()){
							acc_of_center.insert(std::make_pair(cent_acc, std::vector<std::string>()));
						}
						accession = cent_acc;
						short_center = true;
						std::cout << "short center!"<<std::endl;
						c= in.get();
					}
					unsigned int cent_ref;
					in >> cent_ref;
					c = in.get();
					unsigned int cent_left;
					in >> cent_left;
					std::stringstream center_id;
					center_id << cent_ref << ":" << cent_left;
					center.push_back(center_id.str());
					cent = center_id.str();
					std::cout<< center_id.str() << " " ;
					c = in.get();
				}
				std::cout << " " << std::endl;
				if(short_center == true){
					std::map<size_t , std::vector<std::string> >::iterator it = acc_of_center.find(accession);
					it->second.push_back(cent);
				}
				std::cout << "a long center: " ;
				for(size_t i =0; i < center.size();i++){
					std::cout << center.at(i) << " ";
				}
				std::cout << " " <<std::endl;
//				std::map<std::pair<std::vector<std::string>, std::vector<unsigned int> >::iterator it2 = long_centers_high_value.at(i).find(center);
//				if(it2 == long_centers_high_value.at(i).end()){
//					long_centers_high_value.at(i).insert(std::make_pair(center,std::vector<unsigned int>(2,0)));
//					it2 = long_centers_high_value.at(i).find(center);	
//				}
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
				std::vector<unsigned int> temp(2,0);
				temp.at(1)=high_value;
//				it2 -> second.at(1)=high_value;
				std::cout << "high "<< high_value << std::endl;
//				it2 ->second.at(0)=low;
				temp.at(0)=low;
				std::cout << "low "<< low << std::endl;
				low = temp.at(1);
//				low= it2->second.at(1);
				longCenterHigh.push_back(std::make_pair(center, temp));
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
		c = in.get();
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
				long_center_partitions_high_value.push_back(high_value);
				c=in.get();
			}
			c=in.get();
		}
	}
	template<typename T>
	void decoder<T>::set_flag_from_stream(std::ifstream & in){
		size_t bit = 16;
		char c;
		unsigned char h;
		c = in.get();
		std::cout << "c1 " << c << std::endl;
		std::cout << flag_value.size()<<std::endl;
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
				flag_value.push_back(high_value);
				c=in.get();
				std::cout << "c " << c << std::endl;
			}
			c=in.get();
		}
		for(size_t i =0 ; i < flag_value.size(); i++){
			std::cout << "high of flag " << i << " is " << flag_value.at(i)<<std::endl;
		}
	}
	template<typename T>
	void decoder<T>::decode_flag_of_reverse_parts(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out, bool & reverse_center, bool &reverse_member, bool & reverse_both, unsigned int& target){
		size_t bit = 13;
//		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
//		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);
		unsigned int high=0;
		unsigned int low=0;
		//when rev_center
		if(target>=flag&&target<flag+(flag/4)){
			low = model.get_powerOfTwo().at(bit);
			high =  model.get_powerOfTwo().at(bit) + (flag/4);
			std::cout<<"rev_center"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_center = true;
			target = dec.get_target(total);
			std::cout<< "target in reverse "<< target <<std::endl;
		}
		//when rev_member
		if(target>=flag+(flag/2)&&target<flag+(3*flag/4)){
			low =  model.get_powerOfTwo().at(bit) + (flag/2);
			high =  model.get_powerOfTwo().at(bit) + (3*flag/4);
			std::cout<<"rev_member"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_member = true;
			target = dec.get_target(total);
		}
		//both rev
		if(target>=flag+(3*flag/4)&&target<2*flag){
			low =  model.get_powerOfTwo().at(bit) + (3*flag/4);
			high =  model.get_powerOfTwo().at(bit) + 2*flag;
			std::cout<<"rev_both"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_both = true;
			target = dec.get_target(total);
		}
	}
	template<typename T>
	void decoder<T>::decode_optimized_flag_of_reverse_parts(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out, bool & reverse_center, bool &reverse_member, bool & reverse_both, unsigned int& target){
		size_t bit = 13;
	//	unsigned int total = flag_value.at(7);
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		std::cout << "total in rev function "<< total <<std::endl;
		unsigned int high=0;
		unsigned int low=0;
		//when rev_center
		if(target>= flag_value.at(0)&& target< flag_value.at(1)){
			low = flag_value.at(0);
			high = flag_value.at(1);
			std::cout<<"rev_center"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_center = true;
			target = dec.get_target(total);
			std::cout<< "target in reverse "<< target <<std::endl;
		}
		//when rev_member
		else if(target>= flag_value.at(1) && target< flag_value.at(2)){
			low = flag_value.at(1);
			high = flag_value.at(2);
			std::cout<<"rev_member"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_member = true;
			target = dec.get_target(total);
		}
		//both rev
		else if(target>= flag_value.at(2) && target< flag_value.at(3)){
			low =  flag_value.at(2);
			high = flag_value.at(3);
			std::cout<<"rev_both"<<std::endl;
			dec.decode(low,high);
			wrappers.decode(low,high,total);
			reverse_both = true;
			target = dec.get_target(total);
		}
	}
	template<typename T>
	void decoder<T>::al_decoding(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out){
		model.set_patterns(in);
		std::cout<< "sequence patterna are retrieved! "<<std::endl;
		model.set_alignment_pattern(in);
		std::cout<< "alignment patterns are retrieved! "<<std::endl;
		arithmetic_decoding_centers(in,dec);
		std::cout<<"decoding the center has been done!"<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		size_t sequence_counter = 0;
		out << 0;
	//	unsigned int total = model.get_powerOfTwo().at(bit)+20;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		std::string first_al_pattern = model.get_firstAlignmentPattern();
		std::string first_pattern = model.get_firstPattern();
		std::string al_pattern;
		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);
		size_t i = 0;
		while(i<data.numOfAcc()){
			std::cout<< "accession "<< i <<std::endl;
			target = dec.get_target(total);
			while(target < flag + (3*flag/4)){//end of an accession
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
					high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
					low.at(8)= high.at(7);
					high.at(8) = 2*model.get_powerOfTwo().at(bit);
					base = 12;
					for(size_t n = 0; n < 9 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					char p = dnastring::index_to_base(base);
					out<<p;	
					dec.decode(low.at(base),high.at(base));	
					wrappers.decode(low.at(base),high.at(base),total);
					int base_int = base;
					wrappers.decodeContext(base_int);
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}
				else {
					
					std::cout<<"Sequence starts with an alignment!" << target <<std::endl;
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
					low.at(8)= high.at(7);
					high.at(8) = 2*model.get_powerOfTwo().at(bit);
				}
				while(target<low.at(7)){//end of each seq
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
						high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
						low.at(6)= high.at(5);
						high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
						low.at(7)= high.at(6);
						high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
						low.at(8)= high.at(7);
						high.at(8) = 2*model.get_powerOfTwo().at(bit);
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						char p1 = dnastring::index_to_base(base);
						out<<p1;
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
						//	out<<std::endl;
						//	out<< "start of an al" <<std::endl;
						}
					}
					std::cout << "target "<<target <<std::endl;
					if(low.at(5)<=target && target < low.at(7)){//Beginning of an al but not end of the seq
						if(low.at(5)<=target && target < high.at(5)){
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
						dec.decode(low_par,high_par);//Dont need to be added to the 'out'
						wrappers.decode(low_par,high_par,total);
						target = dec.get_target(total);
						for(size_t j =0; j < cluster_high_partition.at(number_of_par).size();j++){
							if(cluster_high_partition.at(number_of_par).at(j).second.at(0) <= target && cluster_high_partition.at(number_of_par).at(j).second.at(1)> target){
								center = cluster_high_partition.at(number_of_par).at(j).first;
								dec.decode(cluster_high_partition.at(number_of_par).at(j).second.at(0),cluster_high_partition.at(number_of_par).at(j).second.at(1));
								wrappers.decode(cluster_high_partition.at(number_of_par).at(j).second.at(0),cluster_high_partition.at(number_of_par).at(j).second.at(1),total);
								break;
							}
						}
						std::cout <<"center: "<< center<< std::endl;
					//	for(std::map<std::string, std::vector<unsigned int> >::iterator it2 = cluster_high_partition.at(number_of_par).begin(); it2!= cluster_high_partition.at(number_of_par).end();it2++){
					//		if(it2->second.at(0) <= target && it2->second.at(1) > target){
					//			center = it2->first;
					//			std::cout <<"center: "<< center<< std::endl;
					//			std::cout<< "center target: "<< target << "total" << total<<std::endl;
					//			dec.decode(it2->second.at(0),it2->second.at(1));
					//			wrappers.decode(it2->second.at(0),it2->second.at(1),total);
					//			std::cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<std::endl;
					//			break;
					//		}	
					//	}		
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
						if(target>=low.at(5)&&target<high.at(5)){
							std::cout<<"rev_center"<<std::endl;
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
							reverse_center ++;
							target = dec.get_target(total);
						}
						//when rev_member
						if(target>=low.at(7)&&target<high.at(7)){
							std::cout<<"rev_member"<<std::endl;
							dec.decode(low.at(7),high.at(7));
							wrappers.decode(low.at(7),high.at(7),total);
							reverse_member ++;
							target = dec.get_target(total);
						}
						//both rev
						if(target>=low.at(8)&&target<high.at(8)){
							std::cout<<"rev_both"<<std::endl;
							dec.decode(low.at(8),high.at(8));
							wrappers.decode(low.at(8),high.at(8),total);
							reverse_both ++;
							target = dec.get_target(total);
						}
						size_t modify = 0;
						std::string member;
						while(target < flag){// the way a context is created doesn't work if the alignment_level is greater than 1!!(TODO Can't remember why I wrote this comment. Check it!)
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
								std::cout<< "member "<< member << " modification "<<modification<<std::endl;

							}
							assert(pos < decodedCenter.size());
							pos += model.modification_length(modification);
							std::cout<< "position on center: " << pos << std::endl;
							std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
							target = dec.get_target(total);	
						}
						if(low.at(6) <= target && target < high.at(6)){
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
								out<<rev_member;
								std::cout<< "add to the 'out' 1"<<endl;								
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
								std::cout<< "add to the 'out' 2"<<endl;
								for(size_t mem =0; mem<member.length();mem++){
									std::cout<< member.at(mem);
								}
								out<<member;
								std::cout << " " << std::endl;
								for(size_t sl =Sequence_level; sl>0; sl--){
									temp += member.at(member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(member.at(member.length()-1));
							}
						}else{
							std::cout<< "add to the 'out' 3"<<endl;
							out << decodedCenter;
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
				if(low.at(7) <= target && target < high.at(7)){
				//	unsigned int l = flag+10;
				//	unsigned int h = flag + 15;
					dec.decode(low.at(7),high.at(7));
					wrappers.decode(low.at(7),high.at(7),total);
					std::cout << " End of a sequence! "<<std::endl;
				//	out << std::endl;
					out << sequence_counter + 1;
					sequence_counter = sequence_counter + 1;
				}else std::cout<<"there is something wrong here!"<<std::endl;
				target = dec.get_target(total);		
			}	
			if(flag+(3*flag/4) <= target && target < 2*flag){
				unsigned int l = flag+(3*flag/4);
				unsigned int h = 2*flag;
				dec.decode(l, h);
				wrappers.decode(l,h,total);
				std::cout << " End of an accession! "<<std::endl;
			}else std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
			i = i + 1;
		}
	}
	template<typename T>
	void decoder<T>::al_decode_with_long_center(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out){
		model.set_patterns(in);
		std::cout<< "sequence patterna are retrieved! "<<std::endl;
		model.set_alignment_pattern(in);
		std::cout<< "alignment patterns are retrieved! "<<std::endl;
		set_long_centers_from_stream(in);
		std::cout << "long centers are set! "<<std::endl;
		decoding_centers_with_long_center(in,dec);
		std::cout<< "centers are decoded! "<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		size_t sequence_counter = 0;
		out << 0;
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
//		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		std::string first_al_pattern = model.get_firstAlignmentPattern();
		std::string first_pattern = model.get_firstPattern();
		std::string al_pattern;
		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);
		size_t i = 0;
		while(i<data.numOfAcc()){
			std::cout<< "accession "<< i <<std::endl;
			target = dec.get_target(total);
			while(target < flag + (3*flag/4)){//end of an accession
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
					high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
					low.at(8)= high.at(7);
					high.at(8) = 2*model.get_powerOfTwo().at(bit);
					base = 12;
					for(size_t n = 0; n < 9 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					char p = dnastring::index_to_base(base);
					out<<p;	
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
					high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
					low.at(7)= high.at(6);
					high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
					low.at(8)= high.at(7);
					high.at(8) = 2*model.get_powerOfTwo().at(bit);
				}
				size_t counter = 0;
				bool reverse_center = false;
				bool reverse_member = false;
				bool reverse_both = false;
				while(target<low.at(7)){//end of each seq
					std::cout << "here "<<std::endl;
					while(target<flag){//If there is no alignment on that position
						char p = dnastring::index_to_base(base);
						std::stringstream s;
						std::string current_pattern;
						if(Sequence_level > 1){// TODO: Make a function in model class which returns the current pattern!
							for(size_t M=1; M< Sequence_level; M++){
								s<<pattern.at(Sequence_level-M);
							//	std::cout << pattern.at(Sequence_level-M) << " ";
							}
						//	std::cout<< "" <<std::endl;
							s<<p;
							s>>current_pattern;
					//		std::cout << current_pattern<< std::endl;

						}
						if(Sequence_level ==1){	
							s<<p;
							s>>current_pattern;
						}
				//		std::cout << "i "<< i <<std::endl;
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
						high.at(5)= model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/4);
						low.at(6)= high.at(5);
						high.at(6) = model.get_powerOfTwo().at(bit) + (model.get_powerOfTwo().at(bit)/2);
						low.at(7)= high.at(6);
						high.at(7) = model.get_powerOfTwo().at(bit) + (3*model.get_powerOfTwo().at(bit)/4);
						low.at(8)= high.at(7);
						high.at(8) = 2*model.get_powerOfTwo().at(bit);
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						char p1 = dnastring::index_to_base(base);
						out<<p1;
						counter = counter + 1;
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						int base_int = base;
						wrappers.decodeContext(base_int);
						pattern = current_pattern;
						target = dec.get_target(total);	
						if(target> flag){
							std::cout<< "starting of an alignemnt's flag!"<<std::endl;
						}
					}
					std::cout<< "target at beginning of an al "<<target <<std::endl;
					std::cout<< "counter "<< counter << std::endl;
					counter = 0;
					if(low.at(5)<=target && target < low.at(7)){//beginning of an alignment but not end of a sequence
						if(low.at(5)<=target && target < high.at(5)){//decoding the flag of beginning of an alignment
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
						}else{
							std::cout<< "counter "<< counter << std::endl;				
							std::cout<< "there is something wrong in b_al flag point!"<<std::endl;
						}
						target = dec.get_target(total);
						std::vector<std::string> center; 
						size_t number_of_par = 0;
						target=dec.get_target(total);
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
						for(size_t j =0; j < long_centers_high_value.at(number_of_par).size();j++){
							if(long_centers_high_value.at(number_of_par).at(j).second.at(0)<=target && long_centers_high_value.at(number_of_par).at(j).second.at(1)>target){
								center = long_centers_high_value.at(number_of_par).at(j).first;
								dec.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1));
								wrappers.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1),total);
								break;
							}
						}
				/*		for(std::map<std::vector<std::string>,std::vector<unsigned int> >::iterator it2 =long_centers_high_value.at(number_of_par).begin(); it2!= long_centers_high_value.at(number_of_par).end();it2++){
							std::cout << "center " << it2->first.at(0) << " l "<< it2->second.at(0)<< " h "<< it2->second.at(1) << "target " << target <<std::endl;
							if(it2->second.at(0) <= target && it2->second.at(1) > target){
								center = it2->first;
								for(size_t i =0; i < center.size();i++){
									std::cout <<"center: "<< center.at(i)<< std::endl;
								}
								std::cout<< "center target: "<< target << "total" << total<<std::endl;
								dec.decode(it2->second.at(0),it2->second.at(1));
								wrappers.decode(it2->second.at(0),it2->second.at(1),total);
								std::cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<std::endl;
								break;
							}	
						}*///Decoding modification!	
						size_t cent_acc;
						std::cout<< "center size "<< center.size() << std::endl;
						for(size_t m =0; m < center.size();m++){
							std::cout << "m "<< m << std::endl;
							for(std::map<size_t, std::vector<std::string> >::iterator it_acc = acc_of_center.begin(); it_acc != acc_of_center.end(); it_acc++){
								for(size_t g= 0; g< it_acc->second.size();g++){
									if(it_acc->second.at(g)== center.at(m)){
										cent_acc = it_acc->first;
										std::cout<< "acc of cent: "<<cent_acc<<std::endl;
										break;
									}else continue;
								}
							}
							for(size_t j =0; j <decoded_center_in_partition.size() ;j++){
								std::map<std::string, std::string>::iterator seq = decoded_center_in_partition.at(j).find(center.at(m));
								if(seq != decoded_center_in_partition.at(j).end()){
									decodedCenter = seq->second;
									break;
								}
							}
							std::cout<< decodedCenter.size()<< " "<< decodedCenter <<std::endl;
							al_pattern = first_al_pattern;//decoding an alignment starts here!
							size_t pos = 0;
							size_t modify = 0;
							std::string member;
							std::vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
							std::vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
							size_t last_base = dnastring::base_to_index(decodedCenter.at(decodedCenter.size()-1));
							target = dec.get_target(total);
							std::cout << "reverse target "<< target <<std::endl;
							decode_flag_of_reverse_parts(in, dec, out, reverse_center, reverse_member, reverse_both,target);
							std::cout << "returned target "<< target <<std::endl;
							while(pos < decodedCenter.length() && target<flag){//we need target< flag for the case of single center that happens on the center
							//	while(target<flag)
									if(reverse_center == true){
										last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
									}
									else if(reverse_both == true){
										last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
									}
									else{//Either both references are forward or member is reverse!??
										last_base = dnastring::base_to_index(decodedCenter.at(pos));
										std::cout<<"both are forward"<<std::endl;
									}
									std::string al_context = al_pattern;
									al_context += last_base;
									std::cout<< "decoding_context: ";
									for(size_t k =0; k < al_context.size(); k++){
										std::cout<< int(al_context.at(k));
										int con =  int(al_context.at(k));
										wrappers.decodeContext(con);
									}
									std::cout<< " " <<std::endl;
									std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are considered from center to the other reference
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
									dec.decode(al_low.at(modification),al_high.at(modification));
									wrappers.decode(al_low.at(modification),al_high.at(modification),total);
									modify = modification;
									std::cout<< " modify: " << modify << std::endl;
									al_pattern = modification;
									if(reverse_center==true || reverse_both ==true){
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
										std::cout<< "member "<< member << " modification "<<modification<<std::endl;
									}
									pos += model.modification_length(modification);
									assert(pos <= decodedCenter.length());
									std::cout<< "position on center: " << pos << std::endl;
									std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
									target = dec.get_target(total);	
								}
								std::cout << " end of al " << target <<std::endl;
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
										std::cout<< "add to the 'out' 1"<<endl;								
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
										std::cout<< "add to the 'out' 2"<<endl;
										for(size_t mem =0; mem<member.length();mem++){
											std::cout<< member.at(mem);
										}
										out<<member;
										std::cout << " " << std::endl;
										for(size_t sl =Sequence_level; sl>0; sl--){
											temp += member.at(member.length()-1-sl);
										}
										pattern = temp; 
										base = dnastring::base_to_index(member.at(member.length()-1));
									}
								}else{
									std::cout<< "add to the 'out' when center on th ref"<<endl;
									out << decodedCenter;
									for(size_t sl=Sequence_level; sl > 0; sl--){
										char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
										temp +=center_base;
									}
									pattern = temp; 
									base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
								}		
								//End of decoding the alignment part
								//Decode few insertions!/
								if(m != center.size()-1){
									al_pattern = first_al_pattern;
									while(target < flag){
										al_pattern += base;
										std::cout<< "decoding_context: ";
										for(size_t k =0; k < al_pattern.size(); k++){
											std::cout<< int(al_pattern.at(k));
											int con =  int(al_pattern.at(k));
											wrappers.decodeContext(con);
										}
										std::cout<< " " <<std::endl;
										std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_pattern);//modifiactions are considered from center to the other reference
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
										//insert_character(modification, base);
										dec.decode(al_low.at(modification),al_high.at(modification));
										wrappers.decode(al_low.at(modification),al_high.at(modification),total);
										modify = modification;
										std::cout<< " modify: " << modify << std::endl;
										out<< dnastring::index_to_base(modification);//Add insertions to the sequence in out!
										al_pattern = modification;
										target = dec.get_target(total);	
									}
									std::cout<< "tar "<<target <<std::endl;
									assert(target >= low.at(6) && target < high.at(6));
									dec.decode(low.at(6), high.at(6));
									wrappers.decode(low.at(6),high.at(6),total);
									std::cout << " End of a piece of an al! "<<std::endl;
									target = dec.get_target(total);
								}
								reverse_center = false;
								reverse_member = false;
								reverse_both = false;
							}//end of looping over the length of "long center"
							//end of an al flag:
						//	target = dec.get_target(total);	
							if(low.at(6) <= target && target < high.at(6)){
								dec.decode(low.at(6), high.at(6));
								wrappers.decode(low.at(6),high.at(6),total);
								std::cout << " End of an al! "<<std::endl;
							}
							target = dec.get_target(total);	
							std::cout<< "back to seq target: "<< target << std::endl;
					}
					if(target<flag){
						std::cout<< "back to sequence from an alignment target " << target <<std::endl;
					}
				}
				if(low.at(7) <= target && target < high.at(7)){
					dec.decode(low.at(7),high.at(7));
					wrappers.decode(low.at(7),high.at(7),total);
					std::cout << " End of a sequence! "<<std::endl;
					out << sequence_counter + 1;
					sequence_counter = sequence_counter + 1;
				}else{
					std::cout<<"there is something wrong here!"<<std::endl;
				}
				target = dec.get_target(total);		
			}	
			if(flag+(3*flag/4) <= target && target < 2*flag){
				unsigned int l = flag+(3*flag/4);
				unsigned int h = 2*flag;
				dec.decode(l, h);
				wrappers.decode(l,h,total);
				std::cout << " End of an accession! "<<std::endl;
			}else std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
			i = i + 1;
		}//end of looping over number of accessions

	}
	template<typename T>
	void decoder<T>::al_decode_long_center_optimized_flag(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec, std::ofstream & out){
		model.set_patterns(in);
		std::cout<< "sequence patterna are retrieved! "<<std::endl;
		model.set_alignment_pattern(in);
		std::cout<< "alignment patterns are retrieved! "<<std::endl;
		set_long_centers_from_stream(in);
		std::cout << "long centers are set! "<<std::endl;
		set_flag_from_stream(in);
		std::cout<< "flags are set! "<<std::endl;
		decoding_centers_with_long_center(in,dec);
		std::cout<< "centers are decoded! "<<std::endl;
		size_t bit = 13;
		size_t base = 0;
		size_t sequence_counter = 0;
		out << 0;
		std::string first_al_pattern = model.get_firstAlignmentPattern();
		std::string first_pattern = model.get_firstPattern();
		std::string al_pattern;
		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);//low value of the first flag
		for(size_t i =0; i < flag_value.size();i++){
			flag_value.at(i)+=flag;
		}
		for(size_t i = 0; i < flag_value.size();i++){
			std::cout << flag_value.at(i) << std::endl;
		}
	//	unsigned int total = flag_value.at(7);
		unsigned int total = 2*model.get_powerOfTwo().at(bit);
		std::cout << "total "<< total <<std::endl;
		size_t i = 0;
		while(i<data.numOfAcc()){
			std::cout<< "accession "<< i <<std::endl;
			target = dec.get_target(total);
			while(target < flag_value.at(6)){//end of an accession
				std::string decodedCenter;
				std::string pattern;
				std::vector<unsigned int>high(5,0);
				std::vector<unsigned int>low(5,0);
				size_t counter = 0;
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
					base = 12;
					for(size_t n = 0; n < 5 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					char p = dnastring::index_to_base(base);
					out<<p;	
					dec.decode(low.at(base),high.at(base));	
					wrappers.decode(low.at(base),high.at(base),total);
					int base_int = base;
					wrappers.decodeContext(base_int);
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}
				bool reverse_center = false;
				bool reverse_member = false;
				bool reverse_both = false;
				while(target<flag_value.at(5)){//end of each seq
					std::cout << "here "<<std::endl;
					while(target<flag){//If there is no alignment on that position
						char p = dnastring::index_to_base(base);
						std::stringstream s;
						std::string current_pattern;
						if(Sequence_level > 1){// TODO: Make a function in model class which returns the current pattern!
							for(size_t M=1; M< Sequence_level; M++){
								s<<pattern.at(Sequence_level-M);
							//	std::cout << pattern.at(Sequence_level-M) << " ";
							}
						//	std::cout<< "" <<std::endl;
							s<<p;
							s>>current_pattern;
						//	std::cout << current_pattern<< std::endl;

						}
						if(Sequence_level ==1){	
							s<<p;
							s>>current_pattern;
						}
				//		std::cout << "i "<< i <<std::endl;
						std::map<std::string, std::vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
						assert(it1 != model.get_high(i).end());
						for(size_t j=0; j<5; j++){
							high.at(j) = it1->second.at(j);
							if(j!= 0){
								low.at(j)=high.at(j-1);
							}else low.at(j)= 0;
						}	
						high.at(4)=  model.get_powerOfTwo().at(bit);
						base = 12;
						for(size_t n = 0; n < 5 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						char p1 = dnastring::index_to_base(base);
						out<<p1;
						counter +=1;
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						int base_int = base;
						wrappers.decodeContext(base_int);
						pattern = current_pattern;
						target = dec.get_target(total);	
					//	std::cout<< "counter1 "<< counter<< " " << base_int << " " << high.at(base)<< std::endl;
						if(target> flag){
							std::cout<< "starting of an alignemnt's flag!"<<std::endl;
						}
					}
					std::cout<< "target at beginning of an al "<<target <<std::endl;
					std::cout<< "counter "<< counter << std::endl;
					counter = 0;
					if(flag<=target && target < flag_value.at(5)){//beginning of an alignment but not end of a sequence
						if(flag<=target && target < flag_value.at(0)){//decoding the flag of beginning of an alignment
							dec.decode(flag,flag_value.at(0));
							wrappers.decode(flag,flag_value.at(0),total);
						}else{
							std::cout<< "counter "<< counter << std::endl;				
							std::cout<< "there is something wrong in b_al flag point!"<<std::endl;
						}
						target = dec.get_target(total);
						std::vector<std::string> center; 
						size_t number_of_par = 0;
						target=dec.get_target(total);
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
						for(size_t j =0; j < long_centers_high_value.at(number_of_par).size();j++){
							if(long_centers_high_value.at(number_of_par).at(j).second.at(0) <=target && long_centers_high_value.at(number_of_par).at(j).second.at(1)>target){
								center = long_centers_high_value.at(number_of_par).at(j).first;
								dec.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1));
								wrappers.decode(long_centers_high_value.at(number_of_par).at(j).second.at(0),long_centers_high_value.at(number_of_par).at(j).second.at(1),total);
								break;
							}
						}
					/*	for(std::map<std::vector<std::string>,std::vector<unsigned int> >::iterator it2 =long_centers_high_value.at(number_of_par).begin(); it2!= long_centers_high_value.at(number_of_par).end();it2++){
							std::cout << "center " << it2->first.at(0) << " l "<< it2->second.at(0)<< " h "<< it2->second.at(1) << "target " << target <<std::endl;
							if(it2->second.at(0) <= target && it2->second.at(1) > target){
								center = it2->first;
								for(size_t i =0; i < center.size();i++){
									std::cout <<"center: "<< center.at(i)<< std::endl;
								}
								std::cout<< "center target: "<< target << "total" << total<<std::endl;
								dec.decode(it2->second.at(0),it2->second.at(1));
								wrappers.decode(it2->second.at(0),it2->second.at(1),total);
								std::cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<std::endl;
								break;
							}	
						}*///Decoding modification!	
						size_t cent_acc;
						std::cout<< "center size "<< center.size() << std::endl;
						for(size_t m =0; m < center.size();m++){
							std::cout << "m "<< m << std::endl;
							for(std::map<size_t, std::vector<std::string> >::iterator it_acc = acc_of_center.begin(); it_acc != acc_of_center.end(); it_acc++){
								for(size_t g= 0; g< it_acc->second.size();g++){
									if(it_acc->second.at(g)== center.at(m)){
										cent_acc = it_acc->first;
										std::cout<< "acc of cent: "<<cent_acc<<std::endl;
										break;
									}else continue;
								}
							}
							for(size_t j =0; j <decoded_center_in_partition.size() ;j++){
								std::map<std::string, std::string>::iterator seq = decoded_center_in_partition.at(j).find(center.at(m));
								if(seq != decoded_center_in_partition.at(j).end()){
									decodedCenter = seq->second;
									break;
								}
							}
							std::cout<< decodedCenter.size()<< " "<< decodedCenter <<std::endl;
							al_pattern = first_al_pattern;//decoding an alignment starts here!
							size_t pos = 0;
							size_t modify = 0;
							std::string member;
							std::vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
							std::vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
							size_t last_base = dnastring::base_to_index(decodedCenter.at(decodedCenter.size()-1));//????
							target = dec.get_target(total);
							std::cout << "reverse target "<< target <<std::endl;
							decode_optimized_flag_of_reverse_parts(in, dec, out, reverse_center, reverse_member, reverse_both,target);
							std::cout << "returned target "<< target <<std::endl;
							while(pos < decodedCenter.length() && target<flag){//we need target< flag for the case of single center that happens on the center
								if(reverse_center == true){
									last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
								}
								else if(reverse_both == true){
									last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
								}
								else{//Either both references are forward or member is reverse!??
									last_base = dnastring::base_to_index(decodedCenter.at(pos));
									std::cout<<"both are forward"<<std::endl;
								}
								std::string al_context = al_pattern;
								al_context += last_base;
								std::cout<< "decoding_context: ";
								for(size_t k =0; k < al_context.size(); k++){
									std::cout<< int(al_context.at(k));
									int con =  int(al_context.at(k));
									wrappers.decodeContext(con);
								}
								std::cout<< " " <<std::endl;
								std::cout << "cent acc "<< cent_acc << "i "<< i <<std::endl;
								std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are considered from center to the other reference
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
										std::cout << "n "<< n << std::endl;
										break;
									}
								}
								dec.decode(al_low.at(modification),al_high.at(modification));
								wrappers.decode(al_low.at(modification),al_high.at(modification),total);
								modify = modification;
								std::cout<< " modify: " << modify << std::endl;
								al_pattern = modification;
								if(reverse_center==true || reverse_both ==true){
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
									std::cout<< "member "<< member << " modification "<<modification<<std::endl;
								}
								pos += model.modification_length(modification);
								assert(pos <= decodedCenter.length());
								std::cout<< "position on center: " << pos << std::endl;
								std::cout<< "center size: "<< decodedCenter.size()<<std::endl;
								target = dec.get_target(total);	
								std::cout << "ta "<< target << std::endl;
							}
							std::cout << " end of a piece of an al " << target <<std::endl;
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
									std::cout<< "add to the 'out' 1"<<endl;								
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
									std::cout<< "add to the 'out' 2"<<endl;
									for(size_t mem =0; mem<member.length();mem++){
										std::cout<< member.at(mem);
									}
									out<<member;
									std::cout << " " << std::endl;
									for(size_t sl =Sequence_level; sl>0; sl--){
										temp += member.at(member.length()-1-sl);
									}
									pattern = temp; 
									base = dnastring::base_to_index(member.at(member.length()-1));
								}
							}else{
								std::cout<< "add to the 'out' when center on th ref"<<endl;
								out << decodedCenter;
								for(size_t sl=Sequence_level; sl > 0; sl--){
									char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
									temp +=center_base;
								}
								pattern = temp; 
								base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
							}		
							//End of decoding the alignment part
							//Decode few insertions!/
							if(m != center.size()-1){
								al_pattern = first_al_pattern;
								while(target < flag){
									al_pattern += base;
									std::cout<< "decoding_context: ";
									for(size_t k =0; k < al_pattern.size(); k++){
										std::cout<< int(al_pattern.at(k));
										int con =  int(al_pattern.at(k));
										wrappers.decodeContext(con);
									}
									std::cout<< " " <<std::endl;
									std::map<std::string, std::vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_pattern);//modifiactions are considered from center to the other reference
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
									//insert_character(modification, base);
									dec.decode(al_low.at(modification),al_high.at(modification));
									wrappers.decode(al_low.at(modification),al_high.at(modification),total);
									modify = modification;
									std::cout<< " modify: " << modify << std::endl;
									out<< dnastring::index_to_base(modification);//Add insertions to the sequence in out!
									al_pattern = modification;
									target = dec.get_target(total);	
								}
								std::cout<< "tar "<<target <<std::endl;
								assert(target >= flag_value.at(3) && target < flag_value.at(4));
								dec.decode(flag_value.at(3), flag_value.at(4));
								wrappers.decode(flag_value.at(3),flag_value.at(4),total);
								std::cout << " End of a gap! "<<std::endl;
								target = dec.get_target(total);
							}
							reverse_center = false;
							reverse_member = false;
							reverse_both = false;
						}//end of looping over the length of "long center"
						//end of an al flag:
						//	target = dec.get_target(total);	
						std::cout << "target "<< target<< " " << flag_value.at(4)<< " " << flag_value.at(5)<<std::endl;
						if(flag_value.at(4) <= target && target < flag_value.at(5)){
							dec.decode(flag_value.at(4), flag_value.at(5));
							wrappers.decode(flag_value.at(4),flag_value.at(5),total);
							std::cout << " End of an al! "<<std::endl;
						}
						target = dec.get_target(total);	
						std::cout<< "back to seq target: "<< target << std::endl;
					}
					if(target<flag){
						std::cout<< "back to sequence from an alignment target " << target <<std::endl;
					}
				}
				if(flag_value.at(5)<= target && target < flag_value.at(6)){
					unsigned int l = flag_value.at(5);
					unsigned int h = flag_value.at(6);
					dec.decode(l,h);
					wrappers.decode(l,h,total);
					std::cout << " End of a sequence! "<<std::endl;
					out << sequence_counter + 1;
					sequence_counter = sequence_counter + 1;
				}else{
					std::cout<<"there is something wrong here!"<<std::endl;
				}
				target = dec.get_target(total);		
			}	
			if(flag_value.at(6) <= target && target < flag_value.at(7)){
				unsigned int l = flag_value.at(6);
				unsigned int h = flag_value.at(7);
				dec.decode(l, h);
				wrappers.decode(l,h,total);
				std::cout << " End of an accession! "<<std::endl;
			}else std::cout<<"there is something wrong at the end of an accession!"<<std::endl;
			i = i + 1;
		}
	}
	template<typename T>	
	void decoder<T>::insert_character(size_t & modification, size_t & base){
		if(modification == NUM_KEEP+NUM_DELETE+5){
			base = 0;
		}
		if(modification == NUM_KEEP+NUM_DELETE+6){
			base = 1;
		}
		if(modification == NUM_KEEP+NUM_DELETE+7){
			base = 2;
		}
		if(modification == NUM_KEEP+NUM_DELETE+8){
			base = 3;
		}
	}
	template<typename T>
	void decoder<T>::arithmetic_decoding_centers(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec){// After partitioning the clusters!
		std::ofstream save1;
		save1.open("decode.txt");
	//	model.set_patterns(in);
	//	std::cout<< "sequence patterna are retrieved! "<<std::endl;
	//	model.set_alignment_pattern(in);
	//	std::cout<< "alignment patterns are retrieved! "<<std::endl;
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
		std::cout << "all the parameters are set" <<std::endl;
		size_t bit =13;
		size_t acc = 0;
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
			for(size_t i =0; i < cluster_high_partition.at(counter).size();i++){
				centerOfAPartition.push_back(cluster_high_partition.at(counter).at(i).first);
			}
	//		for(std::pair<std::string, std::vector<unsigned int> >::iterator it =cluster_high_partition.at(counter).begin(); it != cluster_high_partition.at(counter).end();it ++){
	//			centerOfAPartition.push_back(it->first);		
	//		}
			std::cout<< "size of centerOfAPartition: "<<centerOfAPartition.size()<<std::endl;
			std::map<std::string, std::string> intermediate;
		//	while(acc < data.numOfAcc())
			target = dec.get_target(total);
			while(target<flag+5){//end of all sequences of a partition
			//	std::cout << " m "<< m <<std::endl;
				std::string c_id = centerOfAPartition.at(m);
				std::cout << "center that is supposed to be decoded " << c_id <<std::endl;
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
			//	std::cout<<"first pattern: " << first_pattern<<std::endl;
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
				wrappers.decode(low.at(base),high.at(base),total);
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
				//	std::cout<<"current pattern: "<< current_pattern<<std::endl;
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
					wrappers.decode(low.at(base),high.at(base),total);
				//	save1<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " t: "<< total << std::endl;
				//	std::cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << std::endl;
				//	std::cout << "base: "<< base << std::endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				}
				std::cout << "length of center " << lengthOfSeq << std::endl;
			//	save1 << "length "<< lengthOfSeq <<std::endl;
				base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				std::cout<<"flag base : "<< base << std::endl;
				dec.decode(low.at(base),high.at(base));
				std::string center_id = centerOfAPartition.at(m);
				m = m +1;
			//	std::cout<< "id of cent in dec: " << center_id <<std::endl;
				intermediate.insert(std::make_pair(center_id,center_seq));
				target = dec.get_target(total);	
			//	std::cout<< "target: "<<target << std::endl;
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
					wrappers.decode(low.at(base),high.at(base),total);
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
	template<typename T>
	void decoder<T>::decoding_centers_with_long_center(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec){
		size_t bit =13;
		size_t acc = 0;
		std::ofstream save1("decode", std::ofstream::binary);
		dec.set_stream(in);	
		std::string first_pattern = model.get_firstPattern();
		unsigned int target;
		size_t lengthOfSeq = 0;
		size_t flag = model.get_powerOfTwo().at(bit);
		std::cout << "flag "<<flag <<std::endl;
		unsigned int total = model.get_powerOfTwo().at(bit)+10;
		size_t numCenter = 0;
		size_t counter =0;
		for(size_t i = 0; i < long_centers_high_value.size(); i++){
			for(size_t j = 0; j < long_centers_high_value.at(i).size();j++){
				numCenter = numCenter +1 ;
			}
		}
		std::cout<< "partition size: "<<cluster_high_partition.size()<<std::endl;
		std::cout << " numbe of centers: "<< numCenter <<std::endl;
		while(counter < long_centers_high_value.size()){
			size_t m = 0;			
			std::vector<std::string> centerOfAPartition;
			for(size_t j =0 ; j < long_centers_high_value.at(counter).size();j++){
				if(long_centers_high_value.at(counter).at(j).first.size()==1){
					centerOfAPartition.push_back(long_centers_high_value.at(counter).at(j).first.at(0));
				}

			}
		/*	for(std::map<std::vector<std::string>, std::vector<unsigned int> >::iterator it =long_centers_high_value.at(counter).begin(); it != long_centers_high_value.at(counter).end();it ++){
				if(it->first.size()==1){
					centerOfAPartition.push_back(it->first.at(0));
					std::cout<< "center " << it->first.at(0)<< std::endl;
				}		
			}*/
			std::cout<< "size of centerOfAPartition: "<<centerOfAPartition.size()<<std::endl;
			std::map<std::string, std::string> intermediate;
		//	while(acc < data.numOfAcc())
			target = dec.get_target(total);
			while(target<flag+5){//end of all sequences of a partition 
				std::cout << " m "<< m <<std::endl;
				std::string c_id = centerOfAPartition.at(m);
				std::cout << "c_id " << c_id << std::endl;
				for(std::map<size_t, std::vector<std::string> >::iterator it = acc_of_center.begin();it != acc_of_center.end(); it ++){	
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
				std::vector<unsigned int>high(7,0);
				std::vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				std::string pattern = first_pattern;
			//	std::cout<<"first pattern: " << first_pattern<<std::endl;
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
			//	std::cout << "l " << low.at(base) << " h " << high.at(base) <<std::endl;
				save1<<"l"<< low.at(base)<<"h"<< high.at(base);
				dec.decode(low.at(base),high.at(base));
			//	wrappers.decode(low.at(base),high.at(base),total);
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
				//	std::cout<<"current pattern: "<< current_pattern<<std::endl;
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
				//	wrappers.decode(low.at(base),high.at(base),total);
				//	std::cout << "l " << low.at(base) << " h " << high.at(base) <<std::endl;
					save1<<"l"<< low.at(base)<<"h"<< high.at(base);
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
				std::cout<< "target1: "<<target << std::endl;
				std::cout<<"end of a centerflag: "<< base << std::endl;
				dec.decode(low.at(base),high.at(base));
			//	wrappers.decode(low.at(base),high.at(base),total);
				std::string center_id = centerOfAPartition.at(m);
				m = m +1;
			//	std::cout<< "id of cent in dec: " << center_id <<std::endl;
				intermediate.insert(std::make_pair(center_id,center_seq));
				target = dec.get_target(total);	
				std::cout<< "target: "<<target << std::endl;
				//Decoding the flag which represents end of a partition
				if(target >= flag+5){
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
			//		wrappers.decode(low.at(base),high.at(base),total);
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
	template<typename T>
	void decoder<T>::arithmetic_decoding_seq_test(int t){

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
	void decoder<T>::arithmetic_decoding_seq(){
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
	void decoder<T>::read_from_stream(std::ifstream & in){
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
	}
	

#endif
