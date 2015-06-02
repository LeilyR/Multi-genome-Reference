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
	string encoder::associatedMember(string & center, size_t & modificationPattern, size_t & position){
		string member;
		if(modificationPattern<5){
			char base = dnastring::index_to_base(modificationPattern);
			member += base;
		//	cout<< "member is "<<member<<endl;
			return member;

		}
		if(modificationPattern>=5 && modificationPattern<5+NUM_DELETE){

		}
		if(modificationPattern>=5+NUM_DELETE && modificationPattern<5+NUM_KEEP+NUM_DELETE){
			size_t length = model.modification_length(modificationPattern);
			cout << "center length "<< center.size() << " pattern length: "<< length <<" pattern index : "<< modificationPattern << endl;
			for(size_t i = 0 ; i < length;i++){
				member += center.at(position+i);
			//	cout<< "member is "<<member<<endl;
			}
			return member;
		}
		if(modificationPattern>=5+NUM_KEEP+NUM_DELETE){
			char base = dnastring::index_to_base(20-modificationPattern);
			member += base;
		//	cout<< "member is "<<member<<endl;

			return member;

		}
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
	/*	for(size_t i =0; i < partition.size(); i++){
			for(map<string, vector<unsigned int> >::iterator it2 = cluster_high_partition.at(i).begin(); it2 !=cluster_high_partition.at(i).end(); it2++){
				cout<<"center "<< "in partition " << i << " is "<< it2->first << "high: "<< it2->second.at(1)<<endl;
			}
		}*/
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
	const multimap<size_t, pw_alignment*> & encoder::get_alignment(map<string,vector<pw_alignment> > & alignmentsOfClusters, size_t seq_id){
		setOfAlignments(alignmentsOfClusters);
		return AlignmentsFromClustering.at(seq_id);
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
		//	cout<< "high_write: "<< high << endl;
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
			//	cout<< "high_read: "<< high_value<<endl;
				c=in.get();
			}
			c=in.get();
		}
	}
	void encoder:: al_encoding(map<string, unsigned int> & weight, map<string, string > & cluster_members, map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream &outs,dlib::entropy_encoder_kernel_1 & enc ){
		calculate_high_in_partition(weight,alignmentOfCluster);
		encoding_functor functor(data,&model,wrappers,enc);
		arithmetic_encoding_centers(alignmentOfCluster,outs,enc);
		size_t bit = 13;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
//		size_t i = 1;
		cout<<"number of acc:" <<data.numAcc()<<endl;
		for(size_t i = 0; i< data.numAcc(); i++){
			cout<< " i "<< i <<endl;
			for(size_t k = 0; k <data.getAcc(i).size(); k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
				size_t sequenceId = data.getAcc(i).at(k);
			/*	if(sequenceId == 20){
					cout<< "sequence 20 at 203 and 204:"<<endl;
					for(size_t j = 203; j<=204; j++){
						cout<< sequence.at(j);
					}
				}
				cout<< " " <<endl;*/
				cout<< "sequence length is: "<< sequence.length()<< " sequence id: "<< k <<endl;
				for(size_t n= 0; n < sequence.length(); n++){
					multimap<size_t, pw_alignment* >::iterator it = AlignmentsFromClustering.at(sequenceId).find(n);
					if(it != AlignmentsFromClustering.at(sequenceId).end()){//If there is an alignment on that position
						cout<< "al position: "<< n << endl;
						size_t left_1; 
						size_t left_2;
						size_t right_1;
						size_t right_2;
						pw_alignment *p = it->second;
						cout<< "al length: "<< p->alignment_length()<<endl;
						p->get_lr1(left_1,right_1);
						p->get_lr2(left_2,right_2);		
						cout<<" l1 "<<left_1 << " l2 " << left_2 << " r1 "<< right_1 << " r2 "<<right_2 <<endl;
						cout << "begin1 "<< p->getbegin1() << " begin2 "<< p->getbegin2() << " end1 "<< p->getend1() << " end2 " << p->getend2() << endl;
						//A fixed flag before encoding a center:
						cout<<"enc  bal flag"<<endl;
						unsigned int l1 = model.get_powerOfTwo().at(bit);
						unsigned int h1 = model.get_powerOfTwo().at(bit)+5;
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						stringstream mem;
						mem << sequenceId << ":" << n;
						map<string, string>::iterator cl = cluster_members.find(mem.str());
						assert(cl != cluster_members.end());
						string center = cl->second;
						cout<< "center: "<< center << endl;
						vector<string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int cent_ref = atoi(center_parts.at(0).c_str());
						unsigned int cent_left = atoi(center_parts.at(1).c_str());
						unsigned int center_l;
						unsigned int center_h;
						size_t part;
						for(size_t j = 0; j < partition.size(); j++){
							map<string, vector<unsigned int> >::iterator it1 = cluster_high_partition.at(j).find(center);		
							if(it1 != cluster_high_partition.at(j).end()){
								unsigned int par_low = 0;
								if(j != 0){
									par_low = partitionHigh.at(j-1);
								}
								unsigned int par_high = partitionHigh.at(j);
								enc.encode(par_low,par_high,total);
								wrappers.encode(par_low,par_high,total);
								cout << "par low: "<< par_low << " par high: "<< par_high << endl;									
								center_l = it1 ->second.at(0);
								center_h = it1->second.at(1);
								part = j;
								break;
							}
						}
						cout<< "center low "<< center_l << " center high " << center_h<<endl;
						enc.encode(center_l, center_h, total);
						wrappers.encode(center_l,center_h,total);
						if(cent_ref == sequenceId && cent_left == n){
							cout<< " center is on the ref "<< endl;	
							if(cent_left == p->getend1()||cent_left ==p->getend2()){
								cout<<"reverse center!"<<endl;
							}else{
							}					
						}else{
							if((cent_left == p->getend1()&& left_2==p->getbegin2())||(cent_left ==p->getend2()&&left_1 == p->getbegin1())){
								cout<<"revers_center !"<<endl;
								//center is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 0;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 5;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);

							}
							if((cent_left == p->getbegin1()&& left_2==p->getend2())||(cent_left ==p->getbegin2()&&left_1 == p->getend1())){
							//	p->print();
								cout<<"reverse_member!"<<endl;
								//other ref is reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 10;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 15;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}
							if((cent_left == p->getend1()&& left_2==p->getend2())||(cent_left ==p->getend2()&&left_1 == p->getend1())){
								p->print();
								cout<<"reverse_both!"<<endl;
								//both are reverse
								unsigned int l1 =  model.get_powerOfTwo().at(bit) + 15;
								unsigned int h1 = model.get_powerOfTwo().at(bit) + 20;
								enc.encode(l1,h1,total);
								wrappers.encode(l1,h1,total);
							}

							model.get_encoded_member(*p,cent_ref,cent_left,functor,outs);
							cout<< " center is not on the ref "<< endl;
						}//end of each al
						l1 = model.get_powerOfTwo().at(bit)+5;
						h1 = model.get_powerOfTwo().at(bit) +10;
						enc.encode(l1,h1,total);
						wrappers.encode(l1,h1,total);
						cout<< "end of an al in enc. "<<endl;
						if(p->getreference1()== sequenceId && n == left_1){
							n = right_1;
							cout<< "n_1 "<< n << endl;
						}else{
							n = right_2;
							cout<< "n_2 "<< n << endl;
						}
					}else{//if there is no alignment in that position
						unsigned int l = 0;
						size_t base = dnastring::base_to_index(sequence.at(n));
					//	cout<< "base "<<base << "position "<< n << endl;
						unsigned int h = model.get_high_at_position(sequenceId,n).at(base);
						if(base !=0){
							l = model.get_high_at_position(sequenceId,n).at(base-1);
						}
						if (base == 4){
							h = model.get_powerOfTwo().at(bit);
						}
						enc.encode(l,h,total);
						wrappers.encode(l,h,total);
					}
				}//end of each seq
				unsigned int l1	= model.get_powerOfTwo().at(bit)+10;
				unsigned int h1 = model.get_powerOfTwo().at(bit) +15;
				enc.encode(l1,h1,total);
				wrappers.encode(l1,h1,total);
				cout<<"end of a seq in enc. "<<endl;
			}//end of all sequences of an accsseion
			unsigned int l1 = model.get_powerOfTwo().at(bit)+15;
			unsigned int h1 = model.get_powerOfTwo().at(bit) +20;
			enc.encode(l1,h1,total);
			wrappers.encode(l1,h1,total);
			cout<<"end of an acc"<<endl;
		}
		cout<< "encoding is finished!"<<endl;
	}
	void encoder::al_decoding(ifstream & in, dlib::entropy_decoder_kernel_1 & dec ){
		arithmetic_decoding_centers(in,dec);
		cout<<"decoding the center has been done!"<<endl;
		size_t bit = 13;
		size_t base = 0;
		unsigned int total = model.get_powerOfTwo().at(bit)+20;
		string first_al_pattern = model.get_firstAlignmentPattern();
		string first_pattern = model.get_firstPattern();
		string al_pattern;
		unsigned int target;
		unsigned int flag =  model.get_powerOfTwo().at(bit);
		size_t i = 0;
		while(i<data.numOfAcc()){
			cout<< "accession "<< i <<endl;
			target = dec.get_target(total);
			while(target < flag + 15){//end of an accession
				string decodedCenter;
				string pattern;
				vector<unsigned int>high(9,0);
				vector<unsigned int>low(9,0);
				if(target<flag){//if there is no alignment on the first position
					cout<< "Sequence doesn't start with an alignment!" <<endl;
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
					// By now we just decoded the first base on a sequnce, this base happens after an artificial pattern(in this case "A A")
					target = dec.get_target(total);	
				}
				else {
					cout<<"Sequence starts with an alignment!"<<endl;
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
						stringstream s;
						string current_pattern;
						for(size_t M=1; M< Sequence_level; M++){// Improvement: Make a function in model class which returns the current pattern!
							s<<pattern.at(Sequence_level-M);
						}
						s<<p;
						s>>current_pattern;
					//	cout<< "current pattern: " << current_pattern<<endl;
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
						base = 12;
						for(size_t n = 0; n < 9 ; n++){
							if(low.at(n) <= target && high.at(n) > target){
								base = n;
								break;
							}
						}
						dec.decode(low.at(base),high.at(base));
						wrappers.decode(low.at(base),high.at(base),total);
					//	cout<<"base: "<<base<<endl;
						pattern = current_pattern;
						target = dec.get_target(total);	
					//	cout << "target2: "<<target<<endl;
						if(target> flag){
							cout<< "target2 is bigger than flag!"<<endl;
						}
					}
					if(flag<=target && target < flag+10){
						if(flag<=target && target < flag+5){
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
						}else{
							cout<< "there is a problem in bal flag"<<endl;
						}
						string center;
						size_t number_of_par = 0;
						target=dec.get_target(total);
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
						assert(seq != decoded_center_in_partition.at(number_of_par).end());
						decodedCenter =seq->second;
						cout<< "decoded center: "<< decodedCenter <<endl;
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
						size_t pos = 0;
						//when rev_center
						target = dec.get_target(total);	
						if(target>=flag+0&&target<flag+5){
							cout<<"rev_center"<<endl;
							dec.decode(low.at(5),high.at(5));
							wrappers.decode(low.at(5),high.at(5),total);
							reverse_center ++;
							target = dec.get_target(total);
						}
						//when rev_member
						if(target>=flag+10&&target<flag+15){
							cout<<"rev_member"<<endl;
							dec.decode(low.at(7),high.at(7));
							wrappers.decode(low.at(7),high.at(7),total);
							reverse_member ++;
							target = dec.get_target(total);
						}
						//both rev
						if(target>=flag+15&&target<flag+20){
							cout<<"rev_both"<<endl;
							dec.decode(low.at(8),high.at(8));
							wrappers.decode(low.at(8),high.at(8),total);
							reverse_both ++;
							target = dec.get_target(total);
						}
						size_t modify = 0;
						string member;
						while(target < flag){// the way a context is created doesn't work if the alignment_level is greater than 1!!
							vector<unsigned int>al_high(NUM_KEEP+NUM_DELETE+10,0);
							vector<unsigned int>al_low(NUM_KEEP+NUM_DELETE+10,0);
							size_t last_base;
							if(reverse_center > 0){
								last_base = dnastring::base_to_index(dnastring::complement(decodedCenter.at(decodedCenter.length()-1-pos)));
							}
							else if(reverse_both > 0){
							//	string rev_decodedCenter;
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
							string al_context = al_pattern;
							al_context += last_base;
							cout<< "decoding_context: ";
							for(size_t k =0; k < al_context.size(); k++){
								cout<< int(al_context.at(k));
								int con =  int(al_context.at(k));
								wrappers.decodeContext(pos,con);
							}
							cout<< " " <<endl;
							map<string, vector<unsigned int> >::const_iterator mod =model.get_highValue(cent_acc,i).find(al_context);//modifiactions are considered from center to the other reference
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
							cout<< " modify: " << modify << endl;
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
							if(reverse_center>0 || reverse_both > 0){
								string rev_decodedCenter;
								for(size_t rev = 0; rev< decodedCenter.size(); rev++){
									char chr= dnastring::complement(decodedCenter.at(decodedCenter.size()-1-rev));	
									rev_decodedCenter +=chr;
								}
								member += associatedMember(rev_decodedCenter,modification,pos);							
							}else{
								member += associatedMember(decodedCenter,modification,pos);
							}
							cout<<"al_pattern: ";
							for(size_t al = 0; al < al_pattern.size();al++){
								cout<< int(al_pattern.at(al));
							}
							cout<<" " << endl;
							assert(pos < decodedCenter.size());
							pos += model.modification_length(modification);
							cout<< "position on center: " << pos << endl;
							cout<< "center size: "<< decodedCenter.size()<<endl;
							target = dec.get_target(total);	
						}
						if(flag+5 <= target && target < flag+10){
							dec.decode(low.at(6), high.at(6));
							wrappers.decode(low.at(6),high.at(6),total);
							cout << " End of an al! "<<endl;
						}
						else{
						cout<<"there is something wrong!"<<endl;
						}
						string temp;
						cout<< "length of the member: "<<member.length()<<endl;
						if(member.length()> 0){
							if(reverse_member > 0 || reverse_both > 0){
								string rev_member;
								for(size_t rev = 0; rev< member.size(); rev++){
									char chr= dnastring::complement(member.at(member.size()-1-rev));	
									rev_member +=chr;
								}
								for(size_t mem =0; mem<rev_member.length();mem++){
									cout<<rev_member.at(mem);
								}
								cout << " " << endl;
								for(size_t sl =Sequence_level; sl>0; sl--){
									temp += rev_member.at(rev_member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(rev_member.at(rev_member.length()-1));
							}else{
								for(size_t mem =0; mem<member.length();mem++){
									cout<< member.at(mem);
								}
								cout << " " << endl;
								for(size_t sl =Sequence_level; sl>0; sl--){
									temp += member.at(member.length()-1-sl);
								}
								pattern = temp; 
								base = dnastring::base_to_index(member.at(member.length()-1));
							}
						}else{
							for(size_t sl=Sequence_level; sl > 0; sl--){
								char center_base = decodedCenter.at(decodedCenter.length()-sl-1);
								temp +=center_base;
							}
							pattern = temp; 
							base = dnastring::base_to_index(decodedCenter.at(decodedCenter.length()-1));
						}
						target = dec.get_target(total);	
						reverse_center = 0;
						reverse_member = 0;
						reverse_both = 0;
					}
					if(target<flag){
						cout<< "back to sequence from an alignment"<<endl;
					}
				}
				if(flag+10 <= target && target < flag+15){
					dec.decode(low.at(7),high.at(7));
					wrappers.decode(low.at(7),high.at(7),total);
					cout << " End of a sequence! "<<endl;
				}else cout<<"there is something wrong here!"<<endl;
				target = dec.get_target(total);		
			}	
			if(flag+15 <= target && target < flag+20){
				unsigned int l1 = model.get_powerOfTwo().at(bit)+15;
				unsigned int h1 = model.get_powerOfTwo().at(bit) +20;
				dec.decode(l1, h1);
				wrappers.decode(l1,h1,total);
				cout << " End of an accession! "<<endl;
			}else cout<<"there is something wrong at the end of an accession!"<<endl;
			i = i + 1;
		}
	}
	void encoder::add_center_to_stream(ofstream & outs){
		size_t bit =32;
		for(size_t i = 0 ; i < partition.size(); i ++){
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
		/*	for(map<string, vector<unsigned int> >::iterator it3 = cluster_high_partition.at(i).begin(); it3 != cluster_high_partition.at(i).end();it3++){
				cout<< "center_red_check: "<< it3->first<< "low : "<< it3->second.at(0) << "high: "<< it3->second.at(1)<<endl;	
			}*/
			i = i + 1;
			cout<< "i: "<< i << endl;
		}
	//	cout<< "END!"<<endl;
	/*	for(map<size_t, vector<string> >::iterator it = acc_of_center.begin();it!=acc_of_center.end();it++){
			cout<<"acc set in the map: "<<it->first<<endl;
		}*/
	}
	void encoder::arithmetic_encoding_centers(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs, dlib::entropy_encoder_kernel_1 & enc){// after partitioning the clusters (New one!)
	//	arithmetic_encoding_centId(alignmentOfCluster,outs);
		add_center_to_stream(outs);//center id, its acc and high
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
		size_t bit = 13;
		unsigned int t = model.get_powerOfTwo().at(bit) + 10;
	//	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc.set_stream(outs);
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
				cout<< "id of cent in enc: " << center <<endl;
				if(cent_ref == p->getreference1()&& cent_left == left1){
				//	cout<< "encoded center1: "<< it->first << " its al length: "<< p->alignment_length()<<endl;
					if(p->getbegin1() < p->getend1()){//Encode the sequence itself
						cout<< "id of cent in enc1: " << center <<endl;
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
							enc.encode(l,h,t);
						//	cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length1: "<< cent_right - cent_left << endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin1();
					/*	cout<< "id of cent in enc2: " << center <<endl;
						for(size_t k = cent_right; k >= cent_left;k--){
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
							enc.encode(l,h,t);
						//	cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length2: "<< cent_right - cent_left << endl;*/
						cout<< "id of cent in enc2: " << center <<endl;
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
				//	cout<< "encoded center2: "<< it->first << " its al length: "<< p->alignment_length()<<endl;				
					if(p->getbegin2()<p->getend2()){//Encode the sequence itself
						cent_right = p->getend2();
						cout<< "id of cent in enc3: " << center <<endl;
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
							enc.encode(l,h,t);
					//		cout << "l: " << l << " h: "<< h << " base: " << base << endl;
						}
					//	cout<<"cent length3: "<< cent_right - cent_left << endl;
					}else{ //Encode the reverse complement of the sequence
						cent_right = p->getbegin2();
						cout<< "id of cent in enc4: " << center <<endl;
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
		//		cout<<"end of a center!"<<endl;
				unsigned int l1	= model.get_powerOfTwo().at(bit);
				unsigned int h1 = model.get_powerOfTwo().at(bit)+5;					
				unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
				enc.encode(l1,h1,t1); 
			}
			unsigned int l1	= model.get_powerOfTwo().at(bit)+5;
			unsigned int h1 = model.get_powerOfTwo().at(bit)+10;					
			unsigned int t1 = model.get_powerOfTwo().at(bit)+10;
			enc.encode(l1,h1,t1); 
		//	cout<< "partition: " << i <<endl;
		}
	//	delete enc; 
	//	cout<<"end of cent encoding!"<<endl;
	}
	void encoder::arithmetic_encoding_centId(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream & outs){
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
	const map<string, vector<unsigned int> >& encoder:: get_center_high() const{
		return cluster_high;
	}
	void encoder::arithmetic_decoding_centers(ifstream & in, dlib::entropy_decoder_kernel_1 & dec){
		ofstream save1;
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
	void encoder::write_to_stream( map<string, vector<pw_alignment> > & alignmentOfCluster ,ofstream & outs){		
		add_center_to_stream(outs);
		add_partition_high_to_stream(outs);
		setOfAlignments(alignmentOfCluster);
	}
	void encoder::read_from_stream(ifstream & in){
		model.set_patterns(in);
		model.set_alignment_pattern(in);
		set_center_from_stream(in);
		set_partition_high_from_stream(in);
	}


#endif
