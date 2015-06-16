/*
 * needleman.hpp
 *
 *  Created on: May 12, 2015
 *      Author: mdubarry
 */

#ifndef NEEDLEMAN_HPP_
#define NEEDLEMAN_HPP_
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

std::pair<std::string,std::string> runNeedleman(std::string ref, std::string read, int typeOfAlignment);
void printMatrix(int** matrix, int sizeRef, int sizeRead);
void traceback(int** matrix,std::string ref, std::string read, int typeOfAlignment,std::string& al1,std::string& al2);
int lastLargestIndex( int arr[], int size );



#endif /* NEEDLEMAN_HPP_ */
