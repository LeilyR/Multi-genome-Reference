/*
 * needlman.cpp
 *
 *  Created on: May 12, 2015
 *      Author: mdubarry
 */
#include "needleman.hpp"
bool verbose = true;

int** runNeedleman(std::string ref, std::string read){
	if(verbose) std::cout << "runNeedleman " << std::endl;

	//Init score
	int match = 1;
	int mismatch = -1;
	//int gap = -1;

	//int sizeMatrix =  (ref.size()+1)*(read.size()+1);
	int** matrix = new int*[ref.size()+1];
	for(unsigned int i=0; i <= ref.size() ; ++i){ // <= I have a matrix size ref.size()+1
		matrix[i] = new int[read.size()];
		for(unsigned int j=0; j <= read.size(); ++j){
			matrix[0][j] = j*-1;
			matrix[i][0] = i*-1;

			//take max around cell
			//up = matrix[i-1][j];
			// left = matrix[i][j-1];
			// diag = matrix[i-1][j-1];
			if(i >0 && j>0){


				int listNeighbor[] = {matrix[i-1][j],matrix[i][j-1],matrix[i-1][j-1]};
				int maxScore = *(std::max_element(listNeighbor,listNeighbor+3));
				//update cell
				if(ref[i-1] == read[j-1])
					matrix[i][j] = match + matrix[i-1][j-1];
				else
					matrix[i][j] = mismatch + maxScore;
			}

		}
	}

	printMatrix(matrix,ref.size(),read.size());
	return matrix;
}

void printMatrix(int** matrix, int sizeRef, int sizeRead){
	for(int i =0 ; i <= sizeRef ; ++i ){

		for(int j=0 ; j <= sizeRead ; ++j){
			std::cout << matrix[i][j] << "\t";
		}

		std::cout << std::endl;
	}
}

void traceback(int** matrix, std::string ref, std::string read){
	int i = ref.size();
	int j = read.size();
	std::string al1;
	std::string al2;
	int Sij;
	while(i>0 || j>0){
		if(ref[i-1] == read[j-1])
			Sij = 1;
		else
			Sij = -1;
		if(i>0 && j>0 && matrix[i][j] == matrix[i-1][j-1]+Sij){ // 1 = match
			al1 = ref[i-1]+ al1;
			al2 = read[j-1] + al2;
			--i;
			--j;
		}
		else if(i>0 && matrix[i][j] == matrix[i-1][j]-1){
			al2 = "-"+ al2;
			al1 = ref[i-1] + al1;
			--i;
		}
		else{
			al1 = "-" + al1;
			al2 = read[j-1] + al2;
			--j;
		}
	}
	std::cout << al1 << "\n" << al2 <<std::endl;
}
