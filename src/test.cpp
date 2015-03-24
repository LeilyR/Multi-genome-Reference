#include "test.hpp"

#ifndef TEST_CPP
#define TEST_CPP


	test_encoder::test_encoder(){}
	test_encoder::~test_encoder(){}
	void test_encoder::encode(){
		unsigned int total = 8212;//2^13 +20
		ofstream outs("test",std::ofstream::binary);
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc -> set_stream(outs);
		ifstream read;	
		char c;
		read.open("enc.txt");
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
		cout<< "low value of test: "<<endl;
		for(size_t j=0; j < low.size(); j++){
			cout<< "  low: "<< low.at(j)<< " high: " << high.at(j)<<endl;
		}
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
			if(low.at(i) <= target && high.at(i) > target){
				cout<< "base in test: " << i <<"low: " << low.at(i) << "high: " << high .at(i)<< endl;
				dec.decode(low.at(i),high.at(i));
			}
			else{
				cout<< "target is not in a right range!"<<endl;
				cout << "target: "<< target << " low: "<< low.at(i) << " high: "<< high.at(i)<<endl;
				break;
			}
			i = i + 1;
		}
	}
#endif


