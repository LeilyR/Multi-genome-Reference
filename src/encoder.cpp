#include "encoder.hpp"

#ifndef ENCODER_CPP
#define ENCODER_CPP

	


	encoder::encoder( all_data & d, mc_model & mc,wrapper & wrap): data(d),model(mc),wrappers(wrap),upper_bound(d.numSequences()),AlignmentsFromClustering(data.numSequences()){

	}
	encoder::~encoder(){}
	void encoder::arithmetic_encoding_seq(ofstream & outs){
		size_t bit =12;
//		ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
//		if(outs.is_open()){
			dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
			enc -> set_stream(outs);
			for(size_t i = 0; i< data.numAcc(); i++){
			//	cout<< "accession: " << i << endl;
			//	cout<< "from encoding: "<< dnastring::base_to_index(data.getSequence(data.getAcc(i).at(0)).at(0))<<endl;
				for(size_t k = 0; k <data.getAcc(i).size(); k++){
					const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
					//cout<< "sequence id: "<< k << endl;
			//		cout<< "sequence length: "<< sequence.length() << endl;
					//cout << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << endl;
					for(size_t m = 0; m < sequence.length();m++){
						unsigned int l = 0;
					//	cout << m << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << endl;
					//	cout << m << " accession "<< i << " seq " << k << " " << data.getAcc(i).at(k) << " has length " << sequence.length() <<  " " <<data.getSequence(data.getAcc(i).at(k)).length() << endl;
						size_t base = dnastring::base_to_index(sequence.at(m));
						unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),m).at(base);
						if(base !=0){
							l = model.get_high_at_position(data.getAcc(i).at(k),m).at(base-1);
						}else l = 0;
						unsigned int t = model.get_powerOfTwo().at(bit) + 10;
					/*	if(i == 0){
							cout<< "base "<< base<<endl;
							for(size_t H=0;H<5;H++){
								cout<< "high at "<< H << " is "<< model.get_high_at_position(data.getAcc(i).at(k),m).at(H)<<endl;
							}
						}*/
					/*	if(i == 1){
							cout<<"context in encoding function "<<endl;
							string pattern = model.get_context(m, data.getAcc(i).at(k));
							for(size_t j = 0; j < pattern.length(); j++){
								cout<<pattern.at(j)<< endl;
							}
						}*/
					//	cout<<"base: "<<base<<endl;
					//	map<string, vector<unsigned int > >::const_iterator it1 = model.get_high(i).find(pattern);
					//	for(size_t j = 0; j < 5 ; j++){
					//		cout << " h " << it1 ->second.at(j) << endl;
					//	}
					//	cout<< "accession: "<< i<<" low: "<< l << " high: " << h << " total: "<< t <<" base: "<< base<<endl;
						if(base == 4){
							h=  model.get_powerOfTwo().at(bit);
						}
						enc -> encode(l,h,t);	
					/*	if(i == 1){
							cout<< "accession: "<< i<< " position "<< m <<" low: "<< l << " high: " << h << " total: "<< t <<" base: "<< base<<endl;
						}*/

					}

				//	size_t base = dnastring::base_to_index(sequence.at(sequence.length()-1));
					unsigned int l1	= model.get_powerOfTwo().at(bit);
					unsigned int h1 = model.get_powerOfTwo().at(bit) + 5;					
					unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
					enc -> encode(l1,h1,t1);  
				//	cout<< "accession: "<< i<<" low: "<< l1 << " high: " << h1 << " total: "<< t1 <<" base: "<< "5" <<endl;
				}
			//	const dnastring & sequence = data.getSequence(data.getAcc(i).at(data.getAcc(i).size()-1));
			//	size_t base = dnastring::base_to_index(sequence.at(sequence.length()-1));
				unsigned int l1 =  model.get_powerOfTwo().at(bit) + 5;
				unsigned int h1 = model.get_powerOfTwo().at(bit) + 10;
				unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
				enc -> encode(l1,h1,t1); 
			//	cout<< "accession: "<< i<<" low: "<< l1 << " high: " << h1 << " total: "<< t1 <<" base: "<<  "6" <<endl;
			}
			delete enc;
	//	}
	//	outs.close();
	}
	void encoder::calculating_clusters_high(map<string, unsigned int> & weight){
		size_t max_bit = 8;
		unsigned int h = 0;
		for(map<string, unsigned int>::iterator it = weight.begin(); it != weight.end(); it++){
			string center = it->first;
			map<string , vector<unsigned int> >::iterator it1 =cluster_high.find(center);
			if(it1 == cluster_high.end()){
				cluster_high.insert(make_pair(center,vector<unsigned int>(2,0)));
				it1 = cluster_high.find(center);
			}
			it1 ->second.at(0) = h;
			it1 -> second.at(1) = h + it->second;
			h = it1 -> second.at(1);
		}
	//	assert(h < model.get_powerOfTwo().at(31));
	//	for(map<string, vector<unsigned int> >::iterator it = cluster_high.begin(); it !=cluster_high.end(); it++){
	//		cout<< "low and high value for  cluster center "<< it->first <<" is "<< it->second.at(0) << " and "<<it ->second.at(1) <<endl;			
	//	}
		for(map<string, vector<unsigned int> >::iterator it=cluster_high.begin(); it != cluster_high.end();it++){
		//	cout<<"center0: " << it->first << endl;
		//	cout<< "low: "<< it->second.at(0)<<" high: " << it->second.at(1)<<endl;
		}
	}
	void encoder::partition_centers(map<string,vector<pw_alignment> > & alignmentsOfClusters, map<string, unsigned int> & weight){ //dividing centeres into different partitions (For facing the problem of using not larger than 2^13)
		size_t bit  =13;
		unsigned int length = model.get_powerOfTwo().at(bit);		
		size_t sum_of_weight = 0;
		for(map<string, unsigned int >::iterator it = weight.begin(); it != weight.end(); it ++){
			sum_of_weight = sum_of_weight + it->second ;
		}
		size_t numberOfPartitions = (sum_of_weight/length) + 1;
			cout<< "no. of partition:" << numberOfPartitions << endl;
		//	cout<< "size of al of cluster" << alignmentsOfClusters.size() << endl;
		vector<string> center;
//		for(size_t k = 0; k < data.numAcc(); k++){
			for(map<string, vector<pw_alignment> >::iterator it= alignmentsOfClusters.begin(); it!= alignmentsOfClusters.end();it++){
				string c_id = it->first;
				vector<string> split;
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
			partition.insert(make_pair(j,vector<string>( )));//check if doesnt make a center with a high value bigger than total in encoding
			for(size_t i =j*alignmentsOfClusters.size()/numberOfPartitions; i < (j+1)*alignmentsOfClusters.size()/numberOfPartitions; i++){	
				map<size_t , vector<string> >::iterator it = partition.find(j);
				it->second.push_back(center.at(i));
			} 
		}
	/*	for(map<size_t , vector<string> >::iterator it = partition.begin(); it != partition.end(); it++){
			cout<< "partition "<< it->first<<endl;
			for(size_t i = 0; i <it->second.size(); i++){
				cout<< it->second.at(i)<<endl;
			}
		}*/

	}
	void encoder::calculate_high_in_partition(map<string, unsigned int> & weight, map<string,vector<pw_alignment> > & alignmentsOfClusters){
		partition_centers(alignmentsOfClusters, weight);
		cout<<"partition size: "<< partition.size() <<endl;
		for(size_t i = 0; i < partition.size(); i ++){
			unsigned int h = 0;
			size_t bit = 13;
			map<string, vector<unsigned int> > high;
			map<size_t , vector<string> >::iterator it = partition.find(i);
			assert(it != partition.end());
			for(size_t j = 0; j <it->second.size(); j++){
				string center = it->second.at(j);
				map<string, unsigned int>::iterator it2 = weight.find(center);
				assert(it2 != weight.end());
				map<string , vector<unsigned int> >::iterator it1 =high.find(center);
				if(it1 == high.end()){
					high.insert(make_pair(center,vector<unsigned int>(2,0)));
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
		for(size_t i =0; i < partition.size(); i++){
			for(map<string, vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).begin(); it2 !=cluster_high_partition.at(i).end(); it2++){
				cout<<"center "<< "in partition " << i << " is "<< it2->first << "high: "<< it2->second.at(1)<<endl;
			}
		}
		partition_high();
	}
	void encoder::partition_high(){
		size_t bit = 13;
		unsigned int weight = model.get_powerOfTwo().at(bit);		
		unsigned int high = 0;
		for(size_t i = 0; i < partition.size(); i++){
			high = high + weight/partition.size();
			partitionHigh.push_back(high);
		}
	}
	void encoder::setOfAlignments(map<string,vector<pw_alignment> > & alignmentsOfClusters){//all the alignment that has a certain sequence as at least one of their references.
	//	ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (map<string,vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				string center = it->first;
				vector<string> center_parts;
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
							multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							//	string pattern = model.get_context(right1,i);
							//	outs << pattern;
							//	outs << char(0);
							}else continue;
						}else if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
							multimap<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
							if(it2 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							//	string pattern = model.get_context(right2,i);
							//	outs<< pattern;
							//	outs<< char(0);
							}else continue;
						}else continue;
					}
					if(i == center_ref){
						if(p->getreference1()== i && p->getreference2()== i){
							multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(left1);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							}
							multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left2);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							}
						}else{
							multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(center_left);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(center_left,p));
						//	size_t center_end = center_left + p->alignment_length()-1;
						//	string pattern = model.get_context(center_end,i);
						//	outs << pattern;
						//	outs<< char(0);
							}else continue;
						}
					}
				}
			}
		//	outs<<char(8);
		}
	//	for(multimap<size_t,pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(0).begin();it2 !=AlignmentsFromClustering.at(0).end();it2++){
	//		pw_alignment *p = it2 -> second;
	//		p->print();
		//	cout<<"left on the reference: "<< it2->first<<endl;
	//	}
	//	outs << char(7);
	}
	void encoder::set_pattern_from_stream(ifstream & in){
	/*	char c;
		c=in.get();
		size_t sequence_id=0;
		while(c != 7){
			map<size_t,vector<string> >::iterator it =first_pattern_after_al.find(sequence_id);
			if(it == first_pattern_after_al.end()){
				first_pattern_after_al.insert(make_pair(sequence_id,vector<string>()));
			}
			while(c!=8){
				string pattern;
				stringstream s;
				while(c!=0){
					s<<c;
					c=in.get();
				}
				s>>pattern;
				map<size_t,vector<string> >::iterator it1 =first_pattern_after_al.find(sequence_id);
				it1->second.push_back(s.str());
				c=in.get();
			}
			sequence_id = sequence_id +1;
			c=in.get();
		}
		for(map<size_t,vector<string> >::iterator it1 = first_pattern_after_al.begin(); it1 != first_pattern_after_al.end(); it1 ++){
			for(size_t i =0; i< it1->second.size();i++){
				cout<< it1->second.at(i)<<endl;
			}
		}*/
	}
	void encoder::add_partition_high_to_stream(ofstream & outs){
		size_t bit = 32;
		for(size_t i = 0; i < partitionHigh.size(); i++){
			outs<<(char)0;
			vector<bool> bit_to_byte(0);
			unsigned int high = partitionHigh.at(i);
			cout<< "high_write: "<< high << endl;
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
			//	cout<< int(a)<<endl;
			}
			outs << (char) 7;
		}
		outs << (char) 8;
	}
	
	void encoder::set_partition_high_from_stream(ifstream & in){
		size_t bit = 32;
		unsigned char h;
		char c;
		c = in.get();
		unsigned int low = 0;
		while(c != 8){
			while(c!=7){
				vector<bool> binary_high_value(0);
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
				cout<< "high_read: "<< high_value<<endl;
				c=in.get();
			}
			c=in.get();
		}
	}
	void encoder::arithmetic_encoding_alignment(map<string, unsigned int> & weight, map<string, string > & cluster_members, map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs){
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		ofstream save;
		save.open("values.txt");
		calculate_high_in_partition(weight,alignmentOfCluster);
		encoding_functor functor(data,&model,wrappers);//enc o bezar inja va az jahay ezafe pakesh kon!!
		arithmetic_encoding_centers(alignmentOfCluster,outs);
		string first_pattern = model.get_firstPattern();
		size_t bit = 13;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
	//	for(size_t i = 0; i< data.numAcc(); i++){
			size_t i = 1;
			size_t k = 0;
		//	for(size_t k = 0; k <data.getAcc(i).size(); k++){//retruns all sequences of an acc
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
				unsigned test_counter = 0;
				cout<< "sequence length is: "<< sequence.length()<<endl;
			/*	for(size_t n= 29005; n < 29093; n++){
					cout<< sequence.at(n);
				}
				cout << " " << endl;*/
				for(size_t n= 0; n < sequence.length(); n++){
					multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
						cout<< "al position: "<< n << endl;
						size_t left_1; 
						size_t left_2;
						size_t right_1;
						size_t right_2;
						pw_alignment *p = it->second;
					//	p->print();
						cout<< "al length: "<< p->alignment_length()<<endl;
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);
						//A fixed flag before encoding a center:
						cout<<"enc  bal flag"<<endl;
						unsigned int l1 = model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit)+5;
						enc-> encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						cout << "low: "<< l1 << " high: "<< h1 << endl;					
						save << "low: "<< l1 << " high: "<< h1 << endl;
						if(left_1 == n && p->getreference1()==sequenceId){//curret sequence is its first reference
							cout<< "we are at ref1"<<endl;
							stringstream member;
							member << sequenceId << ":" << n;
							map<string, string>::iterator cl = cluster_members.find(member.str());
							string center = cl->second;
						//	for(map<string,vector<string> >::iterator cl = cluster_members.begin(); cl !=cluster_members.end();cl++){
						//		for(size_t j = 0; j < cl->second.size(); j++){
						//			string member = cl ->second.at(j);
						//			vector<string> member_parts;
						//			strsep(member, ":" , member_parts);
						//			unsigned int ref = atoi(member_parts.at(0).c_str());
						//			unsigned int left = atoi(member_parts.at(1).c_str());
						//			if(ref == sequenceId && left == n){
						//				string center = cl->first;
										cout<< "center1: "<< center << endl;
										vector<string> center_parts;
										strsep(center, ":" , center_parts);
										unsigned int cent_ref = atoi(center_parts.at(0).c_str());
										unsigned int cent_left = atoi(center_parts.at(1).c_str());
										//Encode the center, just as a flag:
										unsigned int center_l;
										unsigned int center_h;
										for(size_t j = 0; j < partition.size(); j++){
											map<string, vector<unsigned int> >::iterator it1 = cluster_high_partition.at(j).find(center);		
											if(it1 != cluster_high_partition.at(j).end()){
												unsigned int par_low = 0;
												if(j != 0){
													par_low = partitionHigh.at(j-1);
												}
												unsigned int par_high = partitionHigh.at(j);
												enc->encode(par_low,par_high,total);
												save << "low: "<< par_low << " high: "<< par_high << endl;
												cout << "par low: "<< par_low << " par high: "<< par_high << endl;									
												wrappers.encode(par_low,par_high,total);
												center_l = it1 ->second.at(0);
												center_h = it1->second.at(1);
												break;
											}
										}
										cout<< "center low: "<< center_l << " center high: "<< center_h << endl;
										enc -> encode(center_l, center_h, total);
										save << "low: "<< center_l << " high: "<< center_h << endl;
										wrappers.encode(center_l,center_h,total);
										if(cent_ref == sequenceId && cent_left == n){ //when center itself is on the sequence
											cout<< "center is on the sequence"<<endl;
										} else { //Otherwise 
											model.get_encoded_member(*p,cent_ref,functor,outs);
											cout<<"modencode ref1"<<endl;
										}
										//encode another flag shows alignment's end.
										unsigned int l2 = model.get_powerOfTwo().at(bit)+5;
										unsigned int h2 = model.get_powerOfTwo().at(bit)+10;
										enc-> encode(l2,h2,total);
										save << "low: "<< l2 << " high: "<< h2 << endl;
										wrappers.encode(l2,h2,total);
									//	break;
								//	}else continue;
							//	}
						//	}
							n = right_1;
							cout << " n: " << n << endl;							
						}
						if (left_2 == n && p->getreference2()==sequenceId){//Just the reference has been changed!
							cout<< "we are at ref2"<<endl;
							stringstream member;
							member << sequenceId << ":" << n;
							map<string, string>::iterator cl = cluster_members.find(member.str());
							string center = cl->second;
					//		for(map<string,vector<string> >::iterator cl = cluster_members.begin(); cl !=cluster_members.end();cl++){
					//			for(size_t j = 0; j < cl->second.size();j++){
					//				string member = cl ->second.at(j);
					//				vector<string> member_parts;
					//				strsep(member, ":" , member_parts);
					//				unsigned int ref = atoi(member_parts.at(0).c_str());
					//				unsigned int left = atoi(member_parts.at(1).c_str());
					//				if(ref == sequenceId && left == n){
					//					string center = cl->first;
										cout<< "center2: "<< center << endl;
										vector<string> center_parts;
										strsep(center, ":" , center_parts);
										unsigned int cent_ref = atoi(center_parts.at(0).c_str());
										unsigned int cent_left = atoi(center_parts.at(1).c_str());
										//Encode the center just as a flag:
										unsigned int center_l;
										unsigned int center_h;
										for(size_t j = 0; j < partition.size(); j++){
											map<string, vector<unsigned int> >::iterator it1 = cluster_high_partition.at(j).find(center);
										//	assert(it1 != cluster_high_partition.at(j).end()); 
											if(it1 != cluster_high_partition.at(j).end()){
												unsigned int par_low = 0;
												if(j != 0){
													par_low = partitionHigh.at(j-1);
												}
												unsigned int par_high = partitionHigh.at(j);
												enc->encode(par_low,par_high,total);
												save << "low: "<< par_low << " high: "<< par_high << endl;
												cout << "low: "<< par_low << " high: "<< par_high << endl;
												wrappers.encode(par_low,par_high,total);
												center_l = it1 ->second.at(0);
												center_h = it1->second.at(1);
												break;
											}
										}
										cout<< "center low: "<< center_l << " center high: "<< center_h << endl;
										enc -> encode(center_l, center_h, total);
										save << "low: "<< center_l << " high: "<< center_h << endl;
										wrappers.encode(center_l,center_h,total);
										if(cent_ref == sequenceId && cent_left == n){ //when cluster center itself is on the sequence
											cout<<"center is on ref"<<endl;
										} else { //Otherwise
											model.get_encoded_member(*p,cent_ref,functor,outs);
											cout<<"modencode ref2"<<endl;
										}
										//encode another flag that shows alignment's end.
										cout<< "end of al flag in encode! "<<endl;
										unsigned int l2 = model.get_powerOfTwo().at(bit)+5;
										unsigned int h2 = model.get_powerOfTwo().at(bit)+10;
										enc-> encode(l2,h2,total);
										save << "low: "<< l2 << " high: "<< h2 << endl;
										wrappers.encode(l2,h2,total);
									//	break;
								//	}else continue;
							//	}
						//	}
							n = right_2;	
							cout << " n: " << n << endl;
						}
					}else{//If there is no alignment on that position
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
						cout<< "seq base: "<<base<< " at position: "<< n << endl;
						unsigned int h = model.get_high_at_position(sequenceId,n).at(base);
						if(base !=0){
							l = model.get_high_at_position(sequenceId,n).at(base-1);
						}
						if (base == 4){
							h = model.get_powerOfTwo().at(bit);
						}
						enc -> encode(l,h,total);
						save << "low: "<< l << " high: "<< h << endl;
						wrappers.encode(l,h,total);
						test_counter = test_counter + 1;
						cout<<"base "<< base <<" l "<< l << " h "<< h<<" position "<< n <<endl;		
					//	save<<"base "<< base <<" l "<< l << " h "<< h<<" position "<< n <<endl;
					}
					cout<< "test counter: "<< test_counter <<endl;//starts at 1
				}
				unsigned int l1	= model.get_powerOfTwo().at(bit)+ 10;
				unsigned int h1 = model.get_powerOfTwo().at(bit) + 15;					
				enc -> encode(l1,h1,total); 
				save << "low: "<< l1 << " high: "<< h1 << endl;
				wrappers.encode(l1,h1,total);
				cout<<"end of a enc seq" << endl;
		//	}
			 l1 =  model.get_powerOfTwo().at(bit) + 15;
			 h1 = model.get_powerOfTwo().at(bit) + 20;
			enc -> encode(l1,h1,total);
			save << "low: "<< l1 << " high: "<< h1 << endl;
			wrappers.encode(l1,h1,total); 
	//	}
		delete enc;
		save.close();// Doesn't have alignements' modifications
	}
	void encoder::arithmetic_decoding_alignment(ifstream& in){
		dlib::entropy_decoder_kernel_1  dec;
		arithmetic_decoding_centers(in);
		cout<<"decoding the center has been done!"<<endl;
		size_t bit = 13;
		size_t base = 0;
		string first_pattern = model.get_firstPattern();
		string first_al_pattern = model.get_firstAlignmentPattern();
		dec.set_stream(in);
		unsigned int target;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		size_t flag = model.get_powerOfTwo().at(bit);	//8192
		size_t i = 1; //Accession number
		size_t seq_id = 0; //Sequence number
		string pattern;
		string al_pattern;
	//	while(i<data.numOfAcc()){
			target = dec.get_target(total);
			cout << "target : "<< target <<endl;
			size_t lengthOfSeq = 0;
			while(target < flag + 15){//End of all the sequences of an accession
				string decodedCenter;
				vector<unsigned int>high(9,0);
				vector<unsigned int>low(9,0);
				if(target<flag){
					cout<< "milestone 1" <<endl;
					lengthOfSeq = 1;
					pattern = first_pattern;
					map<string, vector<unsigned int> >::const_iterator it=model.get_high(i).find(first_pattern);
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
					size_t position = 0;
					cout<<"first base: "<< base << endl;// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}else{// when a sequence starts with an alignment
					cout<< "milestone 2"<<endl;
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
				while(target < flag + 10){//End of a sequence!
					cout<< "target1: "<<target <<endl;
					while(target<flag){//If there is no alignment on that position
						lengthOfSeq = lengthOfSeq + 1;
						char p = dnastring::index_to_base(base);	
						stringstream s;
						string current_pattern;
						for(size_t M=1; M< Sequence_level; M++){// Improvement: Make a function in model class which returns the current pattern!
							s<<pattern.at(Sequence_level-M);
						}
						s<<p;
						s>>current_pattern;
						cout<< "current pattern: " << current_pattern<<endl;
						map<string, vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
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
					/*	for(size_t k =0; k <5; k++){
								cout<< high.at(k)<< " ; ";
							}
							cout<< " "<< endl;*/
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						cout<<"base: "<<base<<endl;
						pattern = current_pattern;
						target = dec.get_target(total);	
						cout << "target2: "<<target<<endl;
					}
					cout<< "length of sequence: " << lengthOfSeq << endl;
					if(target < flag+ 10){
						//Decoding the beginning flag of an alignment:
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						cout<<"bflag: "<<base << endl;//should be 5
						//Decoding the center of a cluster in its own partition:
						string center;
						size_t number_of_par = 0;
						target = dec.get_target(total);	
						cout<< "par target: "<< target <<endl;
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
								cout<< "number of par: "<< number_of_par << endl;
								cout<< "low par: "<< low_par << " high par"<< high_par <<endl;
								break;
							}
						}
						dec.decode(low_par,high_par);
						wrappers.decode(low_par,high_par,total);
						target = dec.get_target(total);
						cout<<"tar: "<< target << endl;
						for(map<string, vector<unsigned int> >::iterator it2 = cluster_high_partition.at(number_of_par).begin(); it2!= cluster_high_partition.at(number_of_par).end();it2++){
							if(it2->second.at(0) <= target && it2->second.at(1) > target){
								center = it2->first;
								cout <<"center: "<< center<< endl;
								cout<< "center target: "<< target << "total" << total<<endl;
								dec.decode(it2->second.at(0),it2->second.at(1));
								wrappers.decode(it2->second.at(0),it2->second.at(1),total);
								cout << " l " << it2->second.at(0) << " h " << it2->second.at(1) <<endl;
								break;
							}	
						}		
						map<string, string>::iterator seq = decoded_center_in_partition.at(number_of_par).find(center);
						decodedCenter =seq->second;
						cout<< "decoded center: "<< decodedCenter <<endl;
						vector<string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(0).c_str());
						unsigned int cent_left = atoi(center_parts.at(1).c_str());
						cout<< "cent left: "<< cent_left << " cent_ref: "<<cent_ref << endl;
						size_t cent_acc;
						for(map<size_t, vector<string> >::iterator it_acc = acc_of_center.begin(); it_acc != acc_of_center.end(); it_acc++){
							for(size_t g= 0; g< it_acc->second.size();g++){
								if(it_acc->second.at(g)== center){
									cent_acc = it_acc->first;
									cout<< "acc of cent: "<<cent_acc<<endl;
									break;
								}else continue;
							}
						}
						al_pattern = first_al_pattern;
						target = dec.get_target(total);
						cout<<" target b_al: "<< target << endl;
						size_t pos = 0;
						while(target < flag+5){//decoding modifications:
							vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
							vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
							cout<< "here!"<<endl;
							size_t last_base = dnastring::base_to_index(decodedCenter.at(pos));
						//	string last_char =  model.print_modification_character(last_base);
							string al_context = al_pattern;
							al_context += last_base;
							cout<< "al context: " << al_context << endl;
							cout<< "size of context: "<< al_context.size() << endl;
							for(size_t m =0; m < al_context.size(); m++){
								cout<< "al_context: "<< int(al_context.at(m)) <<endl;
							}
							map<string, vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are always from center to the other sample
							assert(mod != model.get_highValue(cent_acc,i).end());
							for(size_t k = 0 ; k < mod->second.size();k ++){
								al_high.at(k) = mod->second.at(k);
								if(k > 0){
									al_low.at(k) = al_high.at(k-1);
								}else{
									al_low.at(k) = 0;
								}				
							}	
							size_t modification = NUM_KEEP+NUM_DELETE+20;
							for(size_t n = 0; n < NUM_KEEP+NUM_DELETE+10; n++){
								if(al_low.at(n)<=target && al_high.at(n)> target){
									modification = n;
									break;
								}
							}
							cout<< "modification: " << modification << endl;
							dec.decode(al_low.at(modification),al_high.at(modification));
							wrappers.decode(al_low.at(modification),al_high.at(modification),total);
						//	al_pattern = model.print_modification_character(modification);// return string
							al_pattern = modification;
							cout<<"al_mod_pattern: "<< al_pattern << endl;
							pos += model.modification_length(modification);
							cout<< "position on center: " << pos << endl;
							cout<< "center size: "<< decodedCenter.size()<<endl;
							assert(pos < decodedCenter.size());
							target = dec.get_target(total);		
						}//Decoding the end of alignment's flag!
						base = 12;
						for(size_t n = 0; n < 9 ; n ++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;//should be 6
								cout<< "end of an al flag in dec: " << base << endl;
								break;
							}
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
						string temp;
						for(size_t sl =Sequence_level; sl>0; sl--){
							temp += decodedCenter.at(decodedCenter.length()-1-sl);
						}
						pattern = temp;
						cout<< "first pat of seq after an al: "<<pattern<<endl;
						base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
						cout<<"last base on a center: "<< base <<endl;
						target = dec.get_target(total);
						cout << "target_back to sequence: "<< target << endl;
					}
				}//Decoding the flag shows the end of sequece!
				base = 12;
				for(size_t n = 0; n < 9 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;//should be 7
						cout<< "end of a seq flag in dec: " << base << endl;						
						break;
					}
				}
				dec.decode(low.at(base),high.at(base));
				wrappers.decode(low.at(base),high.at(base),total);
				seq_id = seq_id + 1;
				target =dec.get_target(total);
				if(target >= flag + 5){
				//Decoding the flag shows end of the accession!
					base = 12;
					for(size_t n = 0; n < 9 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;//should be 8
							cout<< "end of an acc flag in dec: " << base << endl;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));
					wrappers.decode(low.at(base),high.at(base),total);
				//	target =dec.get_target(total);
				}		
			}
			
		//	i = i+1;
	//	} 

	}
	void encoder::add_center_to_stream(ofstream & outs){
		size_t bit =32;
		for(size_t i = 0 ; i < partition.size(); i ++){
		/*	for(map<string, vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				cout<< "center_write: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<endl;	
			}*/
			for ( map<string, vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end(); it ++ ){
				string center = it -> first;
				vector<string> center_parts;
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
				vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
				for(size_t i = 0 ; i < bit; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
			//	for(size_t n = 0 ; n < bit_to_byte.size();n++){
			//		cout<< " "<< bit_to_byte.at(n);
			//	}
			//		cout<< "" << endl;
			//	cout<<"bit to byte: "<< bit_to_byte.size() <<endl;
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
//	void encoder::add_acc_to_stream(map<string,vector<pw_alignment> > & alignmentsOfClusters,ofstream & outs){
	void encoder::add_acc_to_stream(ofstream & outs){	
	//	ofstream outs("encode",std::ofstream::binary);
			for(map<string, vector<unsigned int> >::iterator it = cluster_high.begin(); it != cluster_high.end(); it++){
			//	assert(it->second.size() != 0);
				string center = it ->first;
				vector<string> center_parts;
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
		//		map<string, vector<unsigned int> >::iterator it1 = cluster_high.find(center);
				vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
			//	cout<< "center1: "<< it->first<< "low : "<< it->second.at(0) << "high: "<< it->second.at(1)<<endl;	
				for(size_t i = 0 ; i < 32; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
			/*	for(size_t n = 0 ; n < bit_to_byte.size();n++){
					cout<< " "<< bit_to_byte.at(n);
				}
					cout<< "" << endl;*/
			//	cout<<"bit to byte: "<< bit_to_byte.size() <<endl;
				size_t count =0;
				for(size_t n = 0; n < bit_to_byte.size(); n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= model.get_powerOfTwo().at(m-n)* bit_to_byte.at(m);
						count = count +1;
					}
			//		cout<< "a: "<< int(a) << endl;					
					n= n+7;
					outs<< a;
				}
			//	cout << "count: "<< count << endl;
			}
			outs << (char) 8 ;
	}
	void encoder::set_acc_from_stream(ifstream & in){
		size_t bit = 32;
		char c ; 
		unsigned char h;
		c = in.get();
		unsigned int low = 0;
		while(c != 8){
			size_t id;
			in >> id;
			map<size_t , vector<string> >::iterator it = acc_of_center.find(id);
			if(it == acc_of_center.end()){
				acc_of_center.insert(make_pair(id, vector<string>()));
			}
			c = in.get();
			unsigned int cent_ref;
			in >> cent_ref;
			c = in.get();
			unsigned int cent_left;
			in >> cent_left;
			stringstream center_id;
			string center;
			center_id << cent_ref << ":" << cent_left;
			map<size_t , vector<string> >::iterator it1 = acc_of_center.find(id);
			it1 ->second.push_back(center_id.str());
			c = in.get();
			map<string, vector<unsigned int> >::iterator it2 = cluster_high.find(center_id.str());
			if(it2 == cluster_high.end()){
				cluster_high.insert(make_pair(center,vector<unsigned int>(2,0)));
				it2 = cluster_high.find(center);	
			}
			vector<bool> binary_high_value(0);
				size_t bound = bit/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
				//	cout<< "H: "<<int(h)<<endl;
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
			/*	for(size_t n = 0 ; n < binary_high_value.size();n++){
					cout<< " "<< binary_high_value.at(n);
				}
					cout<< "" << endl;*/
				unsigned int high_value = 0;		
			//	cout<<"binary high value size: "<<  binary_high_value.size()<<endl;
				for(size_t i = 0; i < binary_high_value.size();i++){
						high_value += binary_high_value.at(i)*model.get_powerOfTwo().at(i);
				}
		/*		for(size_t i = 0; i < binary_high_value.size();i++){
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(j-i);
					}
					i=i+bit-1;
				}*/
			//	cout << "high value: "<< high_value << endl;
				it2 -> second.at(1)=high_value;
				it2 ->second.at(0)=low;
				low= it2->second.at(1);
				c=in.get();
		}
	/*	for(map<string, vector<unsigned int> >::iterator it3 = cluster_high.begin(); it3 != cluster_high.end();it3++){
			cout<< "center_red_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<endl;	
		}*/
	}
	void encoder::set_center_from_stream(ifstream & in){
		size_t bit = 32;
		char c ; 
		unsigned char h;
		c = in.get();
		size_t i = 0;
		while(c != 8){
			unsigned int low = 0;
			while(c != 7){
			//	cout << "here!"<<endl;
				size_t id;
				in >> id;
				map<size_t , vector<string> >::iterator it = acc_of_center.find(id);
				if(it == acc_of_center.end()){
					acc_of_center.insert(make_pair(id, vector<string>()));
				}
				c = in.get();
				unsigned int cent_ref;
				in >> cent_ref;
				c = in.get();
				unsigned int cent_left;
				in >> cent_left;
				stringstream center_id;
				string center;
				center_id << cent_ref << ":" << cent_left;
				map<size_t , vector<string> >::iterator it1 = acc_of_center.find(id);
				it1 ->second.push_back(center_id.str());
				c = in.get();
				map<string, vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).find(center_id.str());
				if(it2 == cluster_high_partition.at(i).end()){
					cluster_high_partition.at(i).insert(make_pair(center,vector<unsigned int>(2,0)));
					it2 = cluster_high_partition.at(i).find(center);	
				}
				vector<bool> binary_high_value(0);
					size_t bound = bit/8;
					for(size_t j = 0 ; j < bound ; j++){ 
						h=in.get();
				//	cout<< "H: "<<int(h)<<endl;
						for(size_t k = 0; k < 8 ; k++){
							binary_high_value.push_back(h%2);
							h = h/2;	
						}
					}
			/*	for(size_t n = 0 ; n < binary_high_value.size();n++){
					cout<< " "<< binary_high_value.at(n);
				}
					cout<< "" << endl;*/
					unsigned int high_value = 0;		
			//	cout<<"binary high value size: "<<  binary_high_value.size()<<endl;
					for(size_t i = 0; i < binary_high_value.size();i++){
							high_value += binary_high_value.at(i)*model.get_powerOfTwo().at(i);
					}
		/*		for(size_t i = 0; i < binary_high_value.size();i++){
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(j-i);
					}
					i=i+bit-1;
				}*/
			//	cout << "high value: "<< high_value << endl;
					it2 -> second.at(1)=high_value;
					it2 ->second.at(0)=low;
					low= it2->second.at(1);
					c=in.get();
			}
			c = in.get();
			for(map<string, vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				cout<< "center_red_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<endl;	
			}
			i = i + 1;
			cout<< "i: "<< i << endl;
		}
	//	cout<< "END!"<<endl;
	/*	for(map<size_t, vector<string> >::iterator it = acc_of_center.begin();it!=acc_of_center.end();it++){
			cout<<"acc set in the map: "<<it->first<<endl;
		}*/
	}
	void encoder::arithmetic_encoding_centers(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs){// after partitioning the clusters (New one!)
	//	arithmetic_encoding_centId(alignmentOfCluster,outs);
		add_center_to_stream(outs);//center id, its acc and high
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		size_t bit = 13;
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc->set_stream(outs);
	//	cout<< "size of  cluster_high_partition: "<<  cluster_high_partition.size()<<endl;
		for(size_t i = 0 ; i < cluster_high_partition.size(); i++){
		//	for(size_t j = 0 ; j < data.numAcc(); j++)
			for(map<string, vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin();it != cluster_high_partition.at(i).end(); it++){
				string center = it ->first;
				vector<string> split;
				strsep(center, ":" , split);
				size_t cent_right = 0;
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				//if(data.accNumber(cent_ref) == j)
				map<string,vector<pw_alignment> >::iterator it1 = alignmentOfCluster.find(center);
				assert(it1 != alignmentOfCluster.end());
				pw_alignment * p = & it1->second.at(0);
				const dnastring & seq = data.getSequence(cent_ref);
				size_t left1;
				size_t right1;
				size_t left2;
				size_t right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
			//	cout<< "id of cent in enc: " << center <<endl;
				if(cent_ref == p->getreference1()&& cent_left == left1){
					cout<< "encoded center1: "<< it->first << " its al length: "<< p->alignment_length()<<endl;
					if(p->getbegin1() < p->getend1()){//Encode the sequence itself
						cent_right = p->getend1();
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	cout<<"base _enc: "<<base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc-> encode(l,h,t);
						//	cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length1: "<< cent_right - cent_left << endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin1();
						for(size_t k = cent_right; k >= cent_left;k--){
							size_t base = dnastring::base_to_index(seq.at(k));
							char co_base = dnastring::complement(seq.at(k));
							size_t com_base = dnastring::base_to_index(co_base);	
						//	cout<<"base _enc: "<<base <<endl;
							cout<< "com_base_enc: "<< com_base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
							if(com_base !=0){
								l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
							}else l = 0;
							if(com_base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc->encode(l,h,t);
						//	cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length2: "<< cent_right - cent_left << endl;
					}
				}
				if(cent_ref == p->getreference2()&& cent_left == left2){
					cout<< "encoded center2: "<< it->first << " its al length: "<< p->alignment_length()<<endl;				
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
						for(size_t k = cent_left; k < cent_right+1;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	cout<<"base _enc: "<<base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc->encode(l,h,t);
					//		cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length3: "<< cent_right - cent_left << endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin2();
						for(size_t k = cent_right; k >= cent_left;k--){
							size_t base = dnastring::base_to_index(seq.at(k));
							cout<<"original base: "<< base <<endl;
							char co_base = dnastring::complement(seq.at(k));
							size_t com_base = dnastring::base_to_index(co_base);				
						//	cout<<"com_base _enc: "<<com_base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
							if(com_base !=0){
								l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
							}else l = 0;
							if(com_base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc->encode(l,h,t);
						//	cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length4: "<< cent_right - cent_left << endl;
					}
				}
		//		cout<<"end of a center!"<<endl;
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit)+5;					
				unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
				enc->encode(l1,h1,t1); 
			}
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit)+10;					
			unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
			enc->encode(l1,h1,t1); 
		//	cout<< "partition: " << i <<endl;
		}
		delete enc; 
	//	cout<<"end of cent encoding!"<<endl;
	}
	void encoder::arithmetic_enc_centers(map<string, vector<pw_alignment> > & alignmentOfCluster,ofstream & outs){//first sort them by accession! (Old one!)
		arithmetic_enc_centId(alignmentOfCluster, outs);
		size_t bit =12;
		vector<vector<string> >ordered_center(data.numAcc(),vector<string>());
		for(size_t i = 0; i <data.numAcc() ; i++){
			for(map<string,vector<pw_alignment> >::iterator it = alignmentOfCluster.begin(); it != alignmentOfCluster.end();it++){
				string center = it->first;
				vector<string> split;
				strsep(center, ":" , split);
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				size_t acc = data.accNumber(cent_ref);
				if(acc == i){
					ordered_center.at(i).push_back(center);
				}
			}	
		}
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		for(size_t j =0 ; j < data.numAcc(); j++){
			for(size_t i = 0; i < ordered_center.at(j).size();i++){
				string center = ordered_center.at(j).at(i);
				map<string , vector<pw_alignment> >::iterator it= alignmentOfCluster.find(center);
				assert(it != alignmentOfCluster.end());
				pw_alignment * p = & it->second.at(0);
				vector<string> split;
				strsep(center, ":" , split);
				size_t cent_right = 0;
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				const dnastring & seq = data.getSequence(cent_ref);
				size_t left1;
				size_t right1;
				size_t left2;
				size_t right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
		//		cout<< "id of cent in enc: " << center <<endl;
				if(cent_ref == p->getreference1()&& cent_left == left1){
					if(p->getbegin1() < p->getend1()){//Encode the sequence itself
						cent_right = p->getend1();
						for(size_t k = cent_left; k < cent_right;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	cout<<"base _enc: "<<base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc -> encode(l,h,t);
						}
					//	cout<<"cent length1: "<< cent_right - cent_left << endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin1();
						for(size_t k = cent_right; k > cent_left;k--){
							size_t base = dnastring::base_to_index(seq.at(k));
							char co_base = dnastring::complement(seq.at(k));
							size_t com_base = dnastring::base_to_index(co_base);	
						//	cout<<"base _enc: "<<base <<endl;
						//	cout<< "com_base_enc: "<< com_base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
							if(com_base !=0){
								l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
							}else l = 0;
							if(com_base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc -> encode(l,h,t);
						}
					//	cout<<"cent length2: "<< cent_right - cent_left << endl;
					}
				}
				if(cent_ref == p->getreference2()&& cent_left == left2){				
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
						for(size_t k = cent_left; k < cent_right;k++){
							size_t base = dnastring::base_to_index(seq.at(k));				
						//	cout<<"base _enc: "<<base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_center_high_at_position(cent_ref,cent_left,k).at(base);
							if(base !=0){
								l = model.get_center_high_at_position(cent_ref,cent_left,k).at(base-1);
							}else l = 0;
							if(base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc -> encode(l,h,t);
						}
					//	cout<<"cent length3: "<< cent_right - cent_left << endl;

					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin2();
						for(size_t k = cent_right; k > cent_left;k--){
							size_t base = dnastring::base_to_index(seq.at(k));
					//		cout<<"original base: "<< base <<endl;
							char co_base = dnastring::complement(seq.at(k));
							size_t com_base = dnastring::base_to_index(co_base);				
					//		cout<<"com_base _enc: "<<com_base <<endl;
							unsigned int l = 0;
							unsigned int h = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base);
							if(com_base !=0){
								l = model.get_reverse_center_high_at_position(cent_ref,cent_right,k).at(com_base-1);
							}else l = 0;
							if(com_base == 4){
								h=  model.get_powerOfTwo().at(bit);
							}
							enc -> encode(l,h,t);
						}
					//	cout<<"cent length4: "<< cent_right - cent_left << endl;
	
					}
				}
			//	cout<<"end of a center!"<<endl;
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit) + 5;					
				unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
				enc -> encode(l1,h1,t1); 
			//	map<string, vector<unsigned int> >::iterator it1 = cluster_high.find(center);
			//	unsigned int l = it1 ->second.at(0);
			//	unsigned int h = it1->second.at(1);
			//	unsigned int total =  model.get_powerOfTwo().at(10);
			//	cout<< "cent_id_enc"<<center<< " low: "<< it1->second.at(0) << " high: "<< it1->second.at(1)<< endl;		
			//	enc -> encode(l,h,total);
			} 
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit) +10;					
			unsigned int t1 = model.get_powerOfTwo().at(bit) + 10;
			enc -> encode(l1,h1,t1); 
		//	cout<< "acc: " << j <<endl;
		}
		delete enc; 
	//	cout<<"end of cent encoding!"<<endl;
	}
	void encoder::arithmetic_encoding_centId(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs){//new one!
		add_center_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		for(size_t j =0 ; j < partition.size(); j++){
			for(map<string, vector<unsigned int> >::iterator it = cluster_high_partition.at(j).begin(); it != cluster_high_partition.at(j).end(); it++){
				string center = it->first;
				unsigned int l = it ->second.at(0);
				unsigned int h = it->second.at(1);
				unsigned int total =  model.get_powerOfTwo().at(19);
			//	cout<< "cent_id_enc"<<center<< " low: "<< it->second.at(0) << " high: "<< it->second.at(1)<< endl;		
				enc -> encode(l,h,total);
			}
		}
		delete enc;
	}
	void encoder::arithmetic_enc_centId(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs){//old one!
		add_acc_to_stream(outs);//Add cluster high and id to the stream
		setOfAlignments(alignmentOfCluster);
		vector<vector<string> >ordered_center(data.numAcc(),vector<string>());
		for(size_t i = 0; i <data.numAcc() ; i++){
			for(map<string,vector<pw_alignment> >::iterator it = alignmentOfCluster.begin(); it != alignmentOfCluster.end();it++){
				string center = it->first;
				vector<string> split;
				strsep(center, ":" , split);
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				size_t acc = data.accNumber(cent_ref);
				if(acc == i){
					ordered_center.at(i).push_back(center);
				}
			}
		//	cout<< "size: "<< ordered_center.at(i).size() << endl;	
		}
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		for(size_t j =0 ; j < data.numAcc(); j++){
			for(size_t i = 0; i < ordered_center.at(j).size();i++){
				string center = ordered_center.at(j).at(i);
				map<string, vector<unsigned int> >::iterator it1 = cluster_high.find(center);
				unsigned int l = it1 ->second.at(0);
				unsigned int h = it1->second.at(1);
				unsigned int total =  model.get_powerOfTwo().at(19);
		//		cout<< "cent_id_enc"<<center<< " low: "<< it1->second.at(0) << " high: "<< it1->second.at(1)<< endl;		
				enc -> encode(l,h,total);
			}
		}
		delete enc;
	}
	void encoder::arithmetic_decoding_centId(ifstream& in, dlib::entropy_decoder_kernel_1& dec){//new one!
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		dec.set_stream(in);	
		unsigned int target;
		size_t numCenter = 0;
		size_t counter =0;
		string center;
		unsigned int t = model.get_powerOfTwo().at(19);
		for(size_t i = 0; i < cluster_high_partition.size(); i++){
			for(size_t j = 0; j < cluster_high_partition.at(i).size();j++){
				numCenter = numCenter +1 ;
			}
		}
	//	cout << " numbe of centers: "<< numCenter <<endl;
	//	while(counter < numCenter){
			for(size_t i = 0; i < cluster_high_partition.size(); i++){
				for(map<string, vector<unsigned int> >::iterator it = cluster_high_partition.at(i).begin(); it != cluster_high_partition.at(i).end();it++){
					target = dec.get_target(t);
					if(it->second.at(0) <= target && it->second.at(1) > target){
						center = it ->first;
						dec.decode(it->second.at(0), it->second.at(1));
					//	centerId.push_back(center);
					}	
				}
			}
	//	}
	}
	void encoder::arithmetic_dec_centId(ifstream& in, dlib::entropy_decoder_kernel_1& dec){//old one!
		model.set_patterns(in);//returns high values we saved in the stream is called "in"
		model.set_alignment_pattern(in);
		set_acc_from_stream(in);//retrieves accession id & high of each cluster center.
	//	dlib::entropy_decoder_kernel_1  dec;
		dec.set_stream(in);	
		unsigned int target;
		unsigned int t = model.get_powerOfTwo().at(19);
		size_t i = 0;
	//	cout << "cluster high size: "<<cluster_high.size()<<endl;
		while(i<169){//ino bayad avaz koni koli beshe! age kar kard ye flag bayad bezaram- na flag lazem nist bekhatere inke already tedadesshoon ro az rooy set_acc mitooni be dast biary
			target = dec.get_target(t);	
			string center_id;
			for(map<string , vector<unsigned int> >::iterator it2 = cluster_high.begin();it2 != cluster_high.end();it2++){
				if(it2->second.at(0) <= target && it2->second.at(1) > target){
				//	cout<<"target" << target <<endl;
					center_id = it2 ->first;
					dec.decode(it2->second.at(0),it2->second.at(1));
				//	centerId.push_back(center_id);
			//		cout<< "id of cent: " << center_id <<endl;
					break;
				}
			}
			i = i+1;
		}
	}
	const map<string, vector<unsigned int> >& encoder:: get_center_high() const{
		return cluster_high;
	}
	void encoder::arithmetic_decoding_centers(ifstream & in){// After partitioning the cluster! (new one!)
		ofstream save1;
		save1.open("decode.txt");
	//	arithmetic_decoding_centId(in,dec);
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
		size_t bit =13;
		size_t acc = 0;
		dlib::entropy_decoder_kernel_1 dec;
		dec.set_stream(in);	
		string first_pattern = model.get_firstPattern();
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
		cout<< "partition size: "<<cluster_high_partition.size()<<endl;
		cout << " numbe of centers: "<< numCenter <<endl;
		while(counter < cluster_high_partition.size()){
			size_t m = 0;			
			vector<string> centerOfAPartition;
			for(map<string, vector<unsigned int> >::iterator it =cluster_high_partition.at(counter).begin(); it != cluster_high_partition.at(counter).end();it ++){
				centerOfAPartition.push_back(it->first);		
			}
		//	cout<< "size of centerOfAPartition: "<<centerOfAPartition.size()<<endl;
			map<string, string> intermediate;
		//	while(acc < data.numOfAcc())
			target = dec.get_target(total);
			while(target<flag+5){//end of all sequences of a partition
				string c_id = centerOfAPartition.at(m);
				for(map<size_t, vector<string> >::iterator it = acc_of_center.begin();it != acc_of_center.end(); it ++){	
					for(size_t k = 0; k < it->second.size();k++){		
						if(it ->second.at(k) == c_id){
							acc = it->first;
							break;
						}
					}
				}
		//		cout<< "acc in decoding is "<< acc<<endl;
				string center_seq;
				vector<unsigned int>high(7,0);
				vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				string pattern = first_pattern;
		//		cout<<"first pattern: " << first_pattern<<endl;
				map<string, vector<unsigned int> >::const_iterator it=model.get_high(acc).find(first_pattern);
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
			//	cout<<"first base: "<<base<<endl;
			//	cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << endl;
				target = dec.get_target(total);
				while(target < flag){	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					string current_pattern;
					stringstream s;
					for(size_t M=1; M< Sequence_level; M++){// ye doone function to model vase current pattern ham benevis!
						s<<pattern.at(Sequence_level-M);
					}
					s<<p;
				s>>current_pattern;
		//		cout<<"current pattern: "<< current_pattern<<endl;
				map<string, vector<unsigned int> >::const_iterator it1=model.get_high(acc).find(current_pattern);
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
			//	save1<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " t: "<< total << endl;
			//	cout<<" l: "<< low.at(base)<<" h: "<< high.at(base)<< " base: "<< base << endl;
			//	cout << "base: "<< base << endl;
				pattern = current_pattern;
				target = dec.get_target(total);
				}
			//	cout << "length of seq" << lengthOfSeq << endl;
			//	save1 << "length "<< lengthOfSeq <<endl;
				base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
			//	cout<<"flag: "<< base << endl;
				dec.decode(low.at(base),high.at(base));
				string center_id = centerOfAPartition.at(m);
				m = m +1;
			//	cout<< "id of cent in dec: " << center_id <<endl;
				intermediate.insert(make_pair(center_id,center_seq));
				target = dec.get_target(total);	
				//cout<< "target: "<<target << endl;
				//Decoding the flag which represents end of a partition
				if(target >= flag){
					decoded_center_in_partition.push_back(intermediate);
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));
			//		cout<< "base3: " << base <<" " << low.at(base) << " " << high.at(base)  << endl;
				} else continue;
			}
		//	cout << " counter "<< counter <<endl;
			counter = counter + 1 ;
		}
	/*	for(size_t i =0; i < decoded_center_in_partition.size(); i++){
			for(map<string, string>::iterator it = decoded_center_in_partition.at(i).begin();it != decoded_center_in_partition.at(i).end();it++){
				string center = it ->first;
				cout<<"center size in decoded map is: "<<endl;				
				cout<< center.size() << endl;
			}
		}*/
		save1.close();
	}
	void encoder::arithmetic_dec_centers(ifstream & in, dlib::entropy_decoder_kernel_1& dec ){//old one!
		arithmetic_dec_centId(in, dec);
		size_t bit =12;
//		dlib::entropy_decoder_kernel_1  dec;
		dec.set_stream(in);	
		string first_pattern = model.get_firstPattern();
		unsigned int target;
		size_t i = 0;
		size_t m = 0;
		size_t lengthOfSeq = 0;
		size_t flag = model.get_powerOfTwo().at(bit);
		unsigned int total = model.get_powerOfTwo().at(bit) + 10;
	//	unsigned int t = model.get_powerOfTwo().at(10);
		while(i<data.numOfAcc()){
			target = dec.get_target(total);
			while(target < flag+5){
				string center_seq;
				vector<unsigned int>high(7,0);
				vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				string pattern = first_pattern;
				map<string, vector<unsigned int> >::const_iterator it=model.get_high(i).find(first_pattern);
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
				dec.decode(low.at(base),high.at(base));
			//	cout<<"first base: "<<base<<endl;
				target = dec.get_target(total);
				while(target < flag){	
					center_seq +=base;
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					string current_pattern;
					stringstream s;
					for(size_t M=1; M< Sequence_level; M++){// ye doone function to model vase current pattern ham benevis!
						s<<pattern.at(Sequence_level-M);
					}
					s<<p;
					s>>current_pattern;
			//		cout<<"current pattern: "<< current_pattern<<endl;
					map<string, vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
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
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));
				//	cout << "base: "<< base << endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				}
			//	cout << "length of seq" << lengthOfSeq << endl;
				base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
		//		cout<<"flag: "<< base << endl;
				dec.decode(low.at(base),high.at(base));
			//	string center_id = centerId.at(m); ino badan doros kon dige centerId o nadarim
				m = m +1;
		//		cout<< "id of cent in dec: " << center_id <<endl;
		//		decoded_centers.insert(make_pair(center_seq,center_id));
				target = dec.get_target(total);	
				//Decoding the flag which represents end of an accession
				if(target >= flag){
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
					dec.decode(low.at(base),high.at(base));
		//			cout<< "base3: " << base <<" " << low.at(base) << " " << high.at(base)  << endl;
				} else continue;
			}	
	//		cout<<"here!"<<endl;
			i = i+1;
		}
	}
	void encoder::encoding_seq_test(){
	/*	unsigned int low = 0;
		unsigned int high = 0;
		unsigned int total = 0;
		vector<double> num(5,0);
		string seq = "";
	//	string seq1 = "ATAATAAA";
	//	for(size_t n= 0; n< 1000; n++){
	//		seq+=seq1;	
	//	}
		const dnastring seq = data.getSequence(15);
		cout<< "seq length: "<< seq.length()<<endl;
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
		ofstream outs("test", std::ofstream::binary);
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1() ;
		enc->set_stream(outs);
		for(size_t i = 0; i < seq.length(); i++){
			size_t base = dnastring::base_to_index(seq.at(i));	
			unsigned int l = lower_bound.at(base);
			unsigned int h = upper_bound.at(base);
			unsigned int t = total;
			enc->encode(l, h, t);
		}
		for(size_t i=0; i<5; i++) {
			cout<< "low " << lower_bound.at(i) << " | high " << upper_bound.at(i) << "  " << total <<endl;
			cout<<"ratio: "<< (upper_bound.at(i)-lower_bound.at(i))/total <<endl;

		}
		cout<< "sequence: "<<endl;
		for(size_t k =0; k < 50; k++){
			cout<< seq.at(k)<<endl;
		}	
		delete enc;
		decoding_seq(total);*/		
	}
	void encoder::arithmetic_decoding_seq_test(int t){

		dlib::entropy_decoder_kernel_1  dec;
		ifstream in("test", std::ifstream::binary);
	/*	dec.set_stream(in);

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
			cout << " at " << i << " see " << base << endl;

		}*/

	}
	void encoder::arithmetic_decoding_seq(){
		ifstream in("encode",std::ifstream::binary);
		size_t bit =12;
		model.set_patterns(in);//returns high values we saved in the stream is called "in"
		model.set_alignment_pattern(in);
		set_acc_from_stream(in);//retrieves accession id & high of each cluster center.
		dlib::entropy_decoder_kernel_1  dec;
		dec.set_stream(in);	
		string first_pattern = model.get_firstPattern();
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
				vector<unsigned int>high(7,0);
				vector<unsigned int>low(7,0);
				lengthOfSeq = 1;
				string pattern = first_pattern;
				map<string, vector<unsigned int> >::const_iterator it=model.get_high(i).find(first_pattern);
			//	cout<< "first pattern: " << first_pattern << endl;
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
			//	cout<< "target1: "<< target <<endl;
				size_t base = 10;
				for(size_t n = 0; n < 7 ; n ++){
					if(low.at(n) <= target && high.at(n) > target){
						base = n;
						break;
					}
				}
				dec.decode(low.at(base),high.at(base));
				cout<<"first base: "<<base<<endl;
			//	cout<< "base1: " << base << " " << low.at(base) << " " << high.at(base) << endl;
		/*		for(size_t H =0 ; H < 7 ; H++){
					cout << "low at " << H  << " is " << low.at(H) <<endl;
				}*/
			//	while(1)
				target = dec.get_target(total);
				while(target < flag){	
					lengthOfSeq = lengthOfSeq + 1;
					char p = dnastring::index_to_base(base);
					string current_pattern;
					stringstream s;
					for(size_t M=1; M< Sequence_level; M++){// ye doone function to model vase current pattern ham benevis!
						s<<pattern.at(Sequence_level-M);
					}
					s<<p;
					s>>current_pattern;
			//		cout<<"current pattern: "<< current_pattern<<endl;
					map<string, vector<unsigned int> >::const_iterator it1=model.get_high(i).find(current_pattern);
					assert(it1 != model.get_high(i).end());
					for(size_t j=0; j<5; j++){
						high.at(j) = it1->second.at(j);
						if(j!= 0){
							low.at(j)=high.at(j-1);
						}else low.at(j)= 0;
					}	
				/*	cout<<"high values of current pattern in decoding: "<<endl;
					for(size_t k = 0 ; k < 5; k++){
						cout<< it1 -> second.at(k)<<endl;
					}*/
					high.at(4)=  model.get_powerOfTwo().at(bit);
					low.at(5)= high.at(4);
				//	low.at(5)= it1->second.at(4);
					high.at(5)= model.get_powerOfTwo().at(bit) + 5;
					low.at(6)= high.at(5);
					high.at(6) = model.get_powerOfTwo().at(bit) + 10;
				//	total = model.get_powerOfTwo().at(bit) + 10;
				//	target = dec.get_target(total);
				//	cout<< "target: "<< target <<endl;
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						cout << "low at " << H  << " is " << low.at(H) <<endl;
					}*/
					dec.decode(low.at(base),high.at(base));
				//	cout<< "base: " << base <<" " << low.at(base) << " " << high.at(base)  << endl;
					cout << "base: "<< base << endl;
					pattern = current_pattern;
					target = dec.get_target(total);
				/*	if(base >= 5){
						break;
					}*/
				}

				cout << "length of seq" << lengthOfSeq << endl;
			//	cout<<"here1! "<<endl;
			//	cout << "target 2 : "<< target <<endl;	
				base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						cout << "low at " << H  << " is " << low.at(H) <<endl;
					}*/
				dec.decode(low.at(base),high.at(base));	
				target = dec.get_target(total);	
			//	cout<< "target3: "<< target <<endl;
				if(target >= flag){
					base = 10;
					for(size_t n = 0; n < 7 ; n ++){
						if(low.at(n) <= target && high.at(n) > target){
							base = n;
							break;
						}
					}
				/*	for(size_t H =0 ; H < 7 ; H++){
						cout << "low at " << H  << " is " << low.at(H) <<endl;
					}*/
					dec.decode(low.at(base),high.at(base));
			//		cout<< "base3: " << base <<" " << low.at(base) << " " << high.at(base)  << endl;
				} else continue;
			}	
			cout<<"here!"<<endl;
			i = i+1;
		}
	}


#endif
