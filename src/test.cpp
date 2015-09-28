#include "test.hpp"

#ifndef TEST_CPP
#define TEST_CPP


	test_encoder::test_encoder(){}
	test_encoder::~test_encoder(){}
	void test_encoder::encode(){
		low.clear();
		high.clear();
		unsigned int total = 8212;//2^13 +20
		ofstream outs("test",std::ofstream::binary);
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		ifstream read;	
		char c;
		read.open("enc1.txt");
		while(read.good()){
			c = read.get();
			unsigned int l;
			read >> l;
		//	cout <<"low: "<< l << endl;
			low.push_back(l);
			c = read.get();
			unsigned int h;
			read >> h;
		//	cout<<"high: "<<h << endl;
			high.push_back(h);
		}
		cout<< "size of low vector: "<< low.size()<< endl;
		for(size_t m =0; m < low.size(); m++){
			enc->encode(low.at(m),high.at(m),total);
		}
	//	cout<< "low value of test: "<<endl;
	//	for(size_t j=0; j < low.size(); j++){
		//	cout<< "  low: "<< low.at(j)<< " high: " << high.at(j)<<endl;
	//	}
		delete enc;
	}	
	void test_encoder::decode(){
		unsigned int total = 8212;
		unsigned int target;
		size_t i = 0;
		dlib::entropy_decoder_kernel_1  dec;
		ifstream in("test",std::ofstream::binary);
		dec.set_stream(in);
		while(i < low.size()){
			target = dec.get_target(total);
		//	cout<< "target: "<< target << endl;
			if(low.at(i) <= target && high.at(i) > target){
		//		cout<< "base in test: " << i <<"low: " << low.at(i) << "high: " << high .at(i)<< endl;
				dec.decode(low.at(i),high.at(i));
			}
			else{
				cout<< "target is not in a right range!"<<endl;
				cout << "target: "<< target << " low: "<< low.at(i) << " high: "<< high.at(i)<<endl;
				exit(1);
			//	break;
			}
			i = i + 1;
		}
	}
	void test_encoder::compare(){//TODO  is wrong and should be fixed.
		ifstream read;	
		char c;
		low_dec.clear();
		high_dec.clear();
		read.open("dec1.txt");
		while(read.good()){
			c = read.get();
			unsigned int l;
			read >> l;
			low_dec.push_back(l);
			c = read.get();
			unsigned int hi;
			read >> hi;
			high_dec.push_back(hi);
		}
		read.close();
		char h;
		ifstream read1;
		low.clear();
		high.clear();
		read1.open("enc1.txt");
		while(read1.good()){
			h = read1.get();
			unsigned int l;
			read1 >> l;
			low.push_back(l);
			h = read1.get();
			unsigned int hi;
			read1 >> hi;
			high.push_back(hi);


		}
		
				
		std::cout<< "low_dec size: "<< low_dec.size() << "low size "<< low.size()<< std::endl;
		std::cout<< "high_dec size: "<< high_dec.size() << "high size "<< high.size()<< std::endl;
		for(size_t i =0 ; i < high_dec.size();i++){
			cout << "high at " << i << " is " << high.at(i) << " " << high_dec.at(i) << endl;
		}
		for(size_t i = 0; i < low_dec.size(); i++){
			if(low.at(i)!=low_dec.at(i)){
				cout<< "low values at " << i << " are different" << low.at(i) << " " << low_dec.at(i) <<endl;
				exit(1);
			//	break;
			}
			if(high.at(i)!=high_dec.at(i)){
				cout<< "high values at " << i << " are different"<<endl;
				cout<< "high at "<< i << " are " << high.at(i) << " and "<< high_dec.at(i)<<endl;
				exit(1);
			//	break;
			}

		}

		std::cout << "end of control!"<<std::endl;
		
	}
	void test_encoder::context_compare(){//TODO is wrong and should be fixed
		char c;
		vector<size_t> enc_context;
		vector<size_t> dec_context;
		ifstream read;
		read.open("al_encode1.txt");
		while(read.good()){
			c = read.get();
			int h;
			read >> h;
			enc_context.push_back(h);
		}
	//	std::cout<< "encode context: "<< std::endl;
	//	for(size_t i =0; i < enc_context.size(); i++){
	//		std::cout << enc_context.at(i)<<std::endl;
	//	}
		read.close();
		ifstream read1;
		read1.open("al_decode1.txt");
		while(read1.good()){
			c = read1.get();
			int h;
			read1 >> h;
			dec_context.push_back(h);
		}
		std::cout<< "enc size: "<< enc_context.size() << "dec size "<< dec_context.size()<< endl;
	//	assert(enc_context.size()==dec_context.size());
		for(size_t i = 0; i < enc_context.size(); i++){
			if(enc_context.at(i)!= dec_context.at(i)){
				cout<< "contexts at " << i << " are different" << enc_context.at(i) << " " << dec_context.at(i) <<endl;
				exit(1);
			}
		}

		std::cout << "end of context_control!"<<std::endl;



	}
#endif


