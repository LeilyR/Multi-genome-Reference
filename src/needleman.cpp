/*
 * needlman.cpp
 *
 *  Created on: May 12, 2015
 *      Author: mdubarry
 */
#include "needleman.hpp"

/*
 * There is 3 different alignments :
 * 		- 1/ Global Alignment when we look at the middle of the read
 * 		- 2/ Semi-Global Alignment with start free when we look at the start of the read
 * 		- 3/ Semi-Global Alignment with end free when we look at the end of the read
 */

// Perspective : change to Hirschberg's algorithm -> divide and conquer strategy or The Four-Russian Algo
bool VERBOSE = false;

void initMatrix(int**& matrix){

}

std::pair<std::string,std::string> runNeedleman(std::string read, std::string ref,int typeOfAlignment){
	if(VERBOSE) std::cout << "runNeedleman " << std::endl;
	int match = 1;
	int mismatch = -1;
	int** matrix = new int*[read.size()+1];
	for(size_t i=0; i <= read.size() ; i++){
		matrix[i] = new int[ref.size()+1];
		for(size_t j =0; j<=ref.size(); ++j){
			if(typeOfAlignment == 1 || typeOfAlignment == 3 ){// Global alignment, semi-global alignment with end free
				matrix[0][j] = j*-1;
				matrix[i][0] = i*-1;
				if(i >0 && j>0){
					int listNeighbor[] = {matrix[i-1][j],matrix[i][j-1],matrix[i-1][j-1]};
					int maxScore = *(std::max_element(listNeighbor,listNeighbor+3));
					//update cell
					if(ref[j-1] == read[i-1])
						matrix[i][j] = match + matrix[i-1][j-1];
					else
						matrix[i][j] = mismatch + maxScore;
				}
			}
			else{// Semi-Global ALignment : Start of The read, allow gap at the beginning
				//init matrix
				matrix[i][0] = i*-1;
				matrix[0][j] = 0;
				if(i >0 && j>0){
					int listNeighbor[] = {matrix[i-1][j],matrix[i][j-1],matrix[i-1][j-1]};
					int maxScore = *(std::max_element(listNeighbor,listNeighbor+3));
					//update cell
					if(ref[j-1] == read[i-1])
						matrix[i][j] = match + matrix[i-1][j-1];
					else
						matrix[i][j] = mismatch + maxScore;
				}
			}
		}
	}
	std::string al1 = "";
	std::string al2 = "";
	if(VERBOSE) printMatrix(matrix,ref.size(),read.size());
	traceback(matrix,ref,read,typeOfAlignment,al1,al2);
	return make_pair(al1,al2);
}

void printMatrix(int** matrix, int sizeRef, int sizeRead){
	for(int i =0 ; i <= sizeRead ; ++i ){
		for(int j=0 ; j <= sizeRef ; ++j){
			std::cout << matrix[i][j] << "\t";
		}

		std::cout << std::endl;
	}
}

void traceback(int** matrix, std::string ref, std::string read, int typeOfAlignment,std::string& al1,std::string& al2){
	int i, j;
	int Sij;
	if(typeOfAlignment == 1 || typeOfAlignment == 2){
		i = read.size();
		j = ref.size();
	}
	else{
		int liste[ref.size()];
		for(unsigned int l=0; l<ref.size();++l){
			liste[l]=matrix[read.size()][l];
		}
		i = read.size();
		j = lastLargestIndex(liste,ref.size());
		for(int tmp = ref.size(); tmp > j; --tmp){
			al1 = ref[tmp-1] + al1;
			al2 += "-";
		}
	}
	while(i>0 || j>0){
		if(ref[j-1] == read[i-1])
			Sij = 1;
		else
			Sij = -1;
		if(i>0 && j>0 && matrix[i][j] == matrix[i-1][j-1]+Sij){ // 1 = match
			al1 = ref[j-1]+ al1;
			al2 = read[i-1] + al2;
			--i;
			--j;
		}
		else if(i>0 && matrix[i][j] == matrix[i-1][j]-1){
			al1 = "-" + al1;
			al2 = read[i-1] + al2;
			--i;
		}
		else{
			al1 = ref[j-1] + al1;
			al2 = "-"+ al2;
			--j;
		}
	}
	//std::cout << al1 << "\n" << al2 <<std::endl;
}

int lastLargestIndex( int arr[], int size )
{
   int maxIndex = 0;
   for( int i = 0; i < size; i++ )
   {
      if ( arr[maxIndex]  <= arr[i] )
      {
         maxIndex = i;
      }
   }
   return maxIndex;
}
