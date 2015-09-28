#ifndef TEST_HPP
#define TEST_HPP

	#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
	#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

	#include<fstream>	
	#include<iostream>
	#include<string>
	#include<sstream>
	#include <map>
	#include <vector>

	using namespace std;

	class test_encoder{
		public:	
		test_encoder();
		~test_encoder();
		void encode();
		void decode();
		void compare();
		void context_compare();
		private:
		std::vector<unsigned int>low;
		std::vector<unsigned int>high;
		std::vector<unsigned int>low_dec;
		std::vector<unsigned int>high_dec;

	



};

#endif


