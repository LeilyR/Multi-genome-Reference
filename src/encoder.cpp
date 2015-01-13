#include "encoder.hpp"

#ifndef ENCODER_CPP
#define ENCODER_CPP

	


	encoder::encoder( all_data & d, mc_model & mc): data(d),model(mc),upper_bound(d.numSequences()){

	}
	encoder::~encoder(){}
	void encoder::arithmetic_encoding_seq(){
		size_t bit =12;
		ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
		if(outs.is_open()){
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
		}
		outs.close();

		
	}
	void encoder::calculating_clusters_high(map<string, unsigned int> & weight){
		size_t max_bit = 3;
		unsigned int h = 0;
		for(map<string, unsigned int>::iterator it = weight.begin(); it != weight.end(); it++){
			string center = it->first;
			map<string , vector<unsigned int> >::iterator it1 =cluster_high.find(center);
			if(it1 == cluster_high.end()){
				cluster_high.insert(make_pair(center,vector<unsigned int>(0)));
				it1 = cluster_high.find(center);
			}
			it1 ->second.at(0) = h;
			it1 -> second.at(1) = h + it->second;
			h = it1 -> second.at(1);
		}
		for(map<string, vector<unsigned int > >::iterator it = cluster_high.begin(); it !=cluster_high.end(); it++){
			if(it->second.at(0) >= model.get_powerOfTwo().at(max_bit)){
				it -> second.at(0) = model.get_powerOfTwo().at(max_bit)-2;
			}
			if(it->second.at(1) >= model.get_powerOfTwo().at(max_bit)){
				it -> second.at(1) = model.get_powerOfTwo().at(max_bit)-1;
			}
		}
		for(map<string, vector<unsigned int> >::iterator it = cluster_high.begin(); it !=cluster_high.end(); it++){
			cout<< "high value for each cluster center: "<< it ->second.at(1) <<endl;			
		}
	}
//High values of each center will be saved in he stream file.
	void encoder::add_high_to_stream(){
		ofstream outs("encode",std::ofstream::binary|std::ofstream::app);
		if(outs.is_open()){
			for(map<string, vector<unsigned int> >::iterator it = cluster_high.begin(); it != cluster_high.end(); it++){
				outs << it->first;
				outs << (char)0;
				vector<bool> bit_to_byte(0);
				int high = it->second.at(1);
				for(size_t i = 0 ; i < 32; i++){
					bit_to_byte.push_back(high%2);
					high = high/2;
				}
				for(size_t n = 0; n < bit_to_byte.size()-8; n++){
					unsigned char a = 0;
					for(size_t m = n; m <n+8; m++){
						a+= model.get_powerOfTwo().at(n-m+7)* bit_to_byte.at(m);
					}
					n= n+7;
					outs<< a;
				}
				outs<<(char) 8;
			}
		}

	}
	void encoder::set_center_high(ifstream& in){
		//use high values were saved in stream and fill in a data structer is called "cluster_high"
		size_t bit = 12;
		char c;
		char h;
		c=in.get();
		while(c != 8){
			string center;
			stringstream s;
			while(c!=0){
				s<<c;
				c=in.get();
			}
			s>>center;
			unsigned int low = 0;
			map<string, vector<unsigned int> >::iterator it = cluster_high.find(center);
			if(it == cluster_high.end()){
				cluster_high.insert(make_pair(center,vector<unsigned int>(0)));
				it = cluster_high.find(center);	
			}
			vector<bool> binary_high_value(0);
				size_t bound = (5*bit)/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(h%2);
						h = h/2;	
					}
				}
				for(size_t i = 0; i < binary_high_value.size()-bit;i++){
					unsigned int high_value = 0;					
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*model.get_powerOfTwo().at(i-j+bit-1);
					}
					i=i+bit;
					it -> second.at(1)=high_value;
				}

		}
	}
	void encoder::arithmetic_encoding_alignment(map<string, string> & cluster_members, map<string, vector<pw_alignment> > & alignmentOfCluster){//We replace the above function with this one for more compression!
		encoding_functor functor(data);
		add_high_to_stream();
		size_t bit = 12;
		ofstream outs("encode",std::ofstream::binary | std::ofstream::app);
		if(outs.is_open()){
			dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
			enc -> set_stream(outs);
			unsigned int t = model.get_powerOfTwo().at(bit)+10;
			for(size_t i = 0; i< data.numAcc(); i++){
				cout<< "accession: " << i << endl;
				for(size_t k = 0; k <data.getAcc(i).size(); k++){//retrun all sequences of the acc
					const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
					cout<< "length of seq: "<<sequence.length()<<endl;
					size_t position = 0;
					//while(position < sequence.length())
					for(size_t j = 0; j < data.numAlignments(); j++){
						const pw_alignment & al = data.getAlignment(j);
						size_t left1; 
						size_t left2;
						size_t right1;
						size_t right2;
						al.get_lr1(left1,right1);
						al.get_lr2(left2,right2);
						size_t center_id;
						pw_alignment alignment;
						if(al.getreference1()== k){
							for(size_t n = position; n < left1 ; n++){
								unsigned int l = 0;
								size_t base = dnastring::base_to_index(sequence.at(n));
								unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
								if(base !=0){
									l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
								}else l = 0;
								if (base == 4){
									h = model.get_powerOfTwo().at(bit);
								}
								enc -> encode(l,h,t);
							}
							//Cluster center is encoded here!(cluster center can be a piece of original one)
							//The first string is an associated one, the second string is the correspondig center.
							for(map<string,string>::iterator it = cluster_members.begin(); it !=cluster_members.end();it++){
								string member = it ->first;
								vector<string> member_parts;
								strsep(member, ":" , member_parts);
								unsigned int ref = atoi(member_parts.at(0).c_str());
								unsigned int left = atoi(member_parts.at(1).c_str());
								if(ref == k && left1 <= left && right1 >= left){
									for(size_t n = left1; n < left ; n++){
										unsigned int l = 0;
										size_t base = dnastring::base_to_index(sequence.at(n));
										unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
										if(base !=0){
											l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
										}else l = 0;
										if (base == 4){
											h = model.get_powerOfTwo().at(bit);
										}
										enc -> encode(l,h,t);
									}
									string center = it ->second;
									map<string, vector<unsigned int> >::iterator it1 = cluster_high.find(center);
									unsigned int l = it1 ->second.at(0);
									unsigned int h = it1->second.at(1);
									unsigned int t = model.get_powerOfTwo().at(32);
									//Just ad a flag:
									enc ->encode(l,h,t);
									map<string, vector<pw_alignment> >::iterator it2 =  alignmentOfCluster.find(center);
									vector<string> center_parts;
									strsep(member, ":" , center_parts);
									unsigned int center_ref = atoi(center_parts.at(0).c_str());
									unsigned int center_left = atoi(center_parts.at(1).c_str());
									for(size_t k = 0 ; k < it2->second.size(); k++){
										pw_alignment & p = it2->second.at(k);
										if(p.getreference1() == center_ref && p.getreference2()== k){
											center_id = center_ref;
											alignment = p;
										}else if(p.getreference2() == center_ref && p.getreference1() == k){
											center_id = center_ref;
											alignment = p;
										}else continue;
									}
									for(map<string, vector<double> > ::const_iterator it3= model.get_cluster_member_context(alignment,center_id,functor).begin(); it3!= model.get_cluster_member_context(alignment,center_id,functor).end(); it3++){
										string context = it3->first;
										size_t accessionOne = data.accNumber(alignment.getreference1());
										size_t accessionTwo = data.accNumber(alignment.getreference2());
										map<string, vector<unsigned int> > ::const_iterator it4 = model.get_highValue(accessionOne,accessionTwo).find(context);
										vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
										vector<unsigned int> high(NUM_DELETE+NUM_KEEP+10,0);
										for(size_t n = 0; n < NUM_DELETE+NUM_KEEP+10; n++){
											if(n ==0){
												low.at(0) = 0;
											}else low.at(n) = high.at(n-1);
												high.at(n) = it4->second.at(n);
										}	
										for(size_t n = 0; n < NUM_DELETE+NUM_KEEP+10; n++){
								//	cout<<"low: " <<low.at(n) << "high: " << high.at(n)<<endl;
											enc -> encode(low.at(n),high.at(n),t);
										}
										position = left+(alignment.alignment_length());
									}
								}
							}
							//az position ta right1 on encode kon!	
							for(size_t n = position; n < right1 ; n++){
								unsigned int l = 0;
								size_t base = dnastring::base_to_index(sequence.at(n));
								unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
								if(base !=0){
									l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
								}else l = 0;
								if (base == 4){
									h = model.get_powerOfTwo().at(bit);
								}
								enc -> encode(l,h,t);
							}
							position = right1;
//reference 2:
						}else if(al.getreference2()== k){
							for(size_t n = position; n < left2 ; n++){
								unsigned int l = 0;
								size_t base = dnastring::base_to_index(sequence.at(n));
								unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
								if(base !=0){
									l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
								}else l = 0;
								if (base == 4){
									h = model.get_powerOfTwo().at(bit);
								}
								enc -> encode(l,h,t);

							}
							for(map<string,string>::iterator it = cluster_members.begin(); it !=cluster_members.end();it++){
								string member = it->first;
								vector<string> member_parts;
								strsep(member, ":" , member_parts);
								unsigned int ref = atoi(member_parts.at(0).c_str());
								unsigned int left = atoi(member_parts.at(1).c_str());
								if(ref == k && left2 <= left && right2 >= left){
									for(size_t n = left2; n < left ; n++){
										unsigned int l = 0;
										size_t base = dnastring::base_to_index(sequence.at(n));
										unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
										if(base !=0){
											l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
										}else l = 0;
										if (base == 4){
											h = model.get_powerOfTwo().at(bit);
										}
										enc -> encode(l,h,t);
									}
									map<string, vector<unsigned int> >::iterator it1 = cluster_high.find(it->second);
									unsigned int l = it1 ->second.at(0);
									unsigned int h = it1->second.at(1);
									unsigned int t = model.get_powerOfTwo().at(32);
									enc ->encode(l,h,t);
									string center = it ->second;
									map<string, vector<pw_alignment> >::iterator it2 =  alignmentOfCluster.find(center);
									vector<string> center_parts;
									strsep(member, ":" , center_parts);
									unsigned int center_ref = atoi(center_parts.at(0).c_str());
									unsigned int center_left = atoi(center_parts.at(1).c_str());
									for(size_t k = 0 ; k < it2->second.size(); k++){
										pw_alignment & p = it2->second.at(k);
										if(p.getreference1() == center_ref && p.getreference2()== k){
											center_id = center_ref;
											alignment = p;
										}else if(p.getreference2() == center_ref && p.getreference1() == k){
											center_id = center_ref;
											alignment = p;
										}else continue;
									}
									for(map<string, vector<double> > ::const_iterator it3= model.get_cluster_member_context(alignment,center_id,functor).begin(); it3!= model.get_cluster_member_context(alignment,center_id,functor).end(); it3++){
										string context = it3->first;
										size_t accessionOne = data.accNumber(alignment.getreference1());
										size_t accessionTwo = data.accNumber(alignment.getreference2());
										map<string, vector<unsigned int> > ::const_iterator it4 = model.get_highValue(accessionOne,accessionTwo).find(context);
										vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
										vector<unsigned int> high(NUM_DELETE+NUM_KEEP+10,0);
										for(size_t n = 0; n < NUM_DELETE+NUM_KEEP+10; n++){
											if(n ==0){
												low.at(0) = 0;
											}else low.at(n) = high.at(n-1);
												high.at(n) = it4->second.at(n);
										}	
										for(size_t n = 0; n < NUM_DELETE+NUM_KEEP+10; n++){
								//	cout<<"low: " <<low.at(n) << "high: " << high.at(n)<<endl;
											enc -> encode(low.at(n),high.at(n),t);
										}
									}
								}
								position = left+(alignment.alignment_length());
							}
							for(size_t n = position; n < right2 ; n++){
								unsigned int l = 0;
								size_t base = dnastring::base_to_index(sequence.at(n));
								unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
								if(base !=0){
									l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
								}else l = 0;
								if (base == 4){
									h = model.get_powerOfTwo().at(bit);
								}
								enc -> encode(l,h,t);
							}
							position = right2;
						}else{
							continue;
						}
					}
					for(size_t n = position; n < sequence.length();n++){
							unsigned int l = 0;
							unsigned int t = model.get_powerOfTwo().at(bit)+10;
							size_t base = dnastring::base_to_index(sequence.at(n));
							unsigned int h = model.get_high_at_position(data.getAcc(i).at(k),n).at(base);
							if(base !=0){
								l = model.get_high_at_position(data.getAcc(i).at(k),n).at(base-1);
							}else l = 0;
							if (base == 4){
								h = model.get_powerOfTwo().at(bit);
							}
							enc -> encode(l,h,t);
					}
				}
			}
			delete enc;
		}
		outs.close();
	}

	void encoder::arithmetic_decoding_alignment(){
		ifstream in("encode",std::ifstream::binary);
		model.set_patterns(in);//retrieve high values of all patterns for all the sequences of each accession.
		set_center_high(in); //retrieve high values of all centers from the stream file.
		dlib::entropy_decoder_kernel_1  dec;
		dec.set_stream(in);
		unsigned int target;	

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
		dlib::entropy_decoder_kernel_1  dec;
		dec.set_stream(in);	
		string first_pattern = model.get_firstPattern();
		unsigned int target;
		size_t i = 0;
		size_t lengthOfSeq = 0;
		size_t flag = model.get_powerOfTwo().at(bit);
		unsigned int total = model.get_powerOfTwo().at(bit) + 10;
//		for(size_t i = 0; i < data.numOfAcc(); i++){	
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
					pattern = current_pattern;
					target = dec.get_target(total);
				/*	if(base >= 5){
						break;
					}*/
				}

			//	cout << "length of seq" << lengthOfSeq << endl;
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
		//	cout<<"here!"<<endl;
			i = i+1;
		}
	}


#endif
