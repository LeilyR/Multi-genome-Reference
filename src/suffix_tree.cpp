#include "suffix_tree.hpp"


//Finds all the centers the happen on a sequence:
	finding_centers::finding_centers(all_data & d):data(d),AlignmentsFromClustering(data.numSequences()),centersOnSequence(data.numSequences()), centersOfASequence(data.numSequences()),initial_suffixes(data.numSequences()){
	
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		AlignmentsFromClustering.resize(data.numSequences());
		for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
			assert(it->second.size() != 0);
			for(size_t j = 0; j < it->second.size();j++){
				pw_alignment * p = & it->second.at(j);
				size_t left1; 
				size_t left2;
				size_t right1;
				size_t right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
				AlignmentsFromClustering.at(p->getreference1()).insert(std::make_pair(left1,p));
				AlignmentsFromClustering.at(p->getreference2()).insert(std::make_pair(left2,p));
			}
		}
		for(size_t i = 0 ; i < data.numSequences();i++){
			std::cout << "on sequence " << i << ":" << std::endl;
			for(std::map<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).begin(); it1 != AlignmentsFromClustering.at(i).end(); it1++){
				pw_alignment * p = it1->second;
				p->print();
				std::cout << "position is  "<< it1->first << std::endl;
			}
		}
	}
	void finding_centers::findMemberOfClusters(std::map<std::string,vector<pw_alignment> > & alignmentsOfClusters){//fills in a map with centers and their members and removes the direction from members
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
			assert(it->second.size() != 0);
			std::cout << "no of als " << it->second.size()  << " " << it->first << std::endl;
			for(size_t i = 0; i < it->second.size(); i++){
				size_t ref1;
				size_t ref2;
				size_t left1;
				size_t left2;
				size_t right1;
				size_t right2;
				pw_alignment p =it->second.at(i);
				p.print();
				p.get_lr1(left1,right1);
				p.get_lr2(left2,right2);
				ref1 = p.getreference1();
				ref2 = p.getreference2();
				std::stringstream sample1;
				std::stringstream sample2;
				sample1<< ref1 << ":" << left1;
				sample2<< ref2 << ":" << left2;
				std::stringstream member1;
				std::stringstream member2;
				member1 << ref1 << ":"<<left1;
				member2 << ref2 << ":" << left2;
				if(sample1.str()==it->first){
					std::map<std::string,std::string>::iterator mem = memberOfCluster.find(member2.str());
					assert(mem == memberOfCluster.end());
					memberOfCluster.insert(make_pair(member2.str(),it->first));
				}else{
					std::map<std::string,std::string>::iterator mem = memberOfCluster.find(member1.str());
					assert(mem == memberOfCluster.end());

					memberOfCluster.insert(make_pair(member1.str(),it->first));
				}
			}
			std::map<std::string,std::string>::iterator mem = memberOfCluster.find(it->first);
			assert(mem == memberOfCluster.end());

			memberOfCluster.insert(std::make_pair(it->first,it->first));
		}
		std::cout << "memberOfCluster size "<< memberOfCluster.size() << std::endl;
	}
	void finding_centers::center_frequency(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::vector<std::multimap<size_t, std::string> > & centerOnseq){//it basically returns indices of centers on each sequence
		setOfAlignments(alignmentsOfClusters);//set all the alignments on each sequence
		findMemberOfClusters(alignmentsOfClusters);// returns strings that are member of a cluster in memberOfCluster
		//First: fill in the "center_index" vector
		size_t index_counter = 1;
		std::map<int,bool> if_center;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it2=alignmentsOfClusters.begin(); it2 != alignmentsOfClusters.end();it2++){ //Create index for both directions of each center by adding direction to their name.
			std::stringstream str1;
			std::stringstream str2;
			std::string center = it2->first;
			str1<<0<<":"<<center;
			str2<<1<<":"<<center;
			center_id.insert(std::make_pair(str1.str(),index_counter));
			center_index.insert(std::make_pair(index_counter, str1.str()));
			if_center.insert(std::make_pair(index_counter,false));
			center_id.insert(std::make_pair(str2.str(),(-1*index_counter)));
			center_index.insert(std::make_pair((-1*index_counter),str2.str()));
			if_center.insert(std::make_pair(-1*index_counter,false));
			index_counter++;	
		}
		std::cout<< "center index size before " << center_index.size()<<std::endl;
		for(size_t i = 0; i< data.numSequences(); i++){
			std::cout << "sequence: " << i << std::endl;
			//Second: Find all the centers on a sequence
			find_seq_centers(if_center,i,centerOnseq);
			//Third: Find potential candidates for being a long center. 
			//We are going to find those centers that have less than ALLOWED_GAP base pairs  distances. Their indices are saved in a vector.
			find_long_centers(i);
			std::cout << "check for fully reversed ones" <<std::endl;
			//At this stage i need to check for fully reversed centers and merge two clusters in to one.
			std::map<size_t , std::vector<int> > intermediate;
			for(std::map<size_t,std::vector<int> >::iterator it = centersOfASequence.at(i).begin() ; it != centersOfASequence.at(i).end();it++){
				std::map<size_t,std::vector<int> >::iterator rev_pos = intermediate.find(it->first);
				if(rev_pos == intermediate.end()){
					std::map<int,std::string>::iterator centid = center_index.find(it->second.at(0));
					std::string center = centid->second;
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int dir = atoi(center_parts.at(0).c_str());
					unsigned int ref = atoi(center_parts.at(1).c_str());
					unsigned int left = atoi(center_parts.at(2).c_str());
					for(std::map<size_t , std::vector<int> >::iterator it1 = centersOfASequence.at(i).begin() ; it1 != centersOfASequence.at(i).end();it1++){
						centid = center_index.find(it1->second.at(it1->second.size()-1));
						std::string cent = centid->second;
						std::vector<std::string> cent_parts;
						strsep(cent, ":" , cent_parts);
						unsigned int dir1 = atoi(center_parts.at(0).c_str());
						unsigned int ref1 = atoi(center_parts.at(1).c_str());
						unsigned int left1 = atoi(center_parts.at(2).c_str());
						bool Fullyreverse = true;
						if(ref == ref1 && left == left1 && dir != dir1 && it1->second.size() == it->second.size()){
							for(size_t i = 0; i < it1->second.size();i++){//We are sure that the length is higher than 1
								// All the next left and ref should be checked!
								centid = center_index.find(it1->second.at(i));
								std::string center1 = centid->second;
								std::vector<std::string> center1_parts;
								strsep(cent, ":" , center1_parts);
								unsigned int dir_c1 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c1 = atoi(center_parts.at(1).c_str());
								unsigned int left_c1 = atoi(center_parts.at(2).c_str());
								centid = center_index.find(it->second.at(it1->second.size()-1+i));
								std::string center2 = centid->second;
								std::vector<std::string> center2_parts;
								strsep(cent, ":" , center2_parts);
								unsigned int dir_c2 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c2 = atoi(center_parts.at(1).c_str());
								unsigned int left_c2 = atoi(center_parts.at(2).c_str());
								if(ref_c1 == ref_c2 && left_c1 == left_c2 && dir_c1 != dir_c2){
								
								}else{
									Fullyreverse = false;
									break;
								}
							}
						}else{
							Fullyreverse = false;
						}
						if(Fullyreverse == true){//The fully reverse exists and we keep one of them
							std::cout << "fully reversed! "<<std::endl;
							intermediate.insert(std::make_pair(it1->first,it1->second));
							break;
						}
					}
				}

			}
			//Remove the intermediate from centersOfASequence
			for(std::map< size_t , std::vector<int> >::iterator it = intermediate.begin();it != intermediate.end();it++){
				centersOfASequence.at(i).erase(it);
				std::map<size_t , std::vector<std::vector< int> > >::iterator it1 =long_centers_of_a_sequence.find(i);
				for(size_t i = 0; i < it1->second.size();i++){
					if(it1->second.at(i)==it->second){
						it1->second.erase(it1->second.begin()+(i-1));
						break;
					}
				}
			}
			//Print the potential long centers:
			for(std::map< size_t , std::vector<int> >::iterator it = centersOfASequence.at(i).begin();it != centersOfASequence.at(i).end();it++){
				std::cout << "position "<< it->first <<std::endl;
				for(size_t i =0; i < it->second.size();i++){
					std::cout << it->second.at(i)<<  " ";
				}
				std::cout << ""<<std::endl;
			}
		}
		for(std::map<int, bool>::iterator it = if_center.begin(); it != if_center.end(); it++){
			if(it->second == false){
				std::map<int,std::string>::iterator index = center_index.find(it->first);
				assert(index != center_index.end());
				std::string center = index->second;
				std::map<std::string,int>::iterator id=center_id.find(center);
				assert(id!=center_id.end());
				center_index.erase(index);
				center_id.erase(id);
			}
		}
		std::cout<< "center index size after " << center_index.size() << " " <<center_id.size()<<std::endl;

	}
	void finding_centers::find_seq_centers(std::map<int, bool> & if_center, size_t & seq_id , std::vector<std::multimap<size_t, std::string> > & centerOnseq){
		const dnastring & sequence = data.getSequence(seq_id);
		for(size_t n= 0; n < sequence.length(); n++){
			std::map<size_t, pw_alignment*>::iterator it=AlignmentsFromClustering.at(seq_id).find(n);
			if(it != AlignmentsFromClustering.at(seq_id).end()){
				pw_alignment * p = it->second;
				size_t l1,r1,l2,r2;
				p->get_lr1(l1,r1);
				p->get_lr2(l2,r2);
				std::stringstream member;
				member<< seq_id << ":" << n;
				std::cout << "mem" << member.str() <<std::endl;
				std::map<std::string,std::string>::iterator it1 = memberOfCluster.find(member.str());
				assert(it1 != memberOfCluster.end());
				std::string center = it1->second;				
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int cent_ref = atoi(center_parts.at(0).c_str());
				unsigned int cent_left = atoi(center_parts.at(1).c_str());
				if(p->getreference1()==cent_ref && l1==cent_left){
					if(p->getbegin1()<p->getend1()){
						std::stringstream str;
						str<<0<<":"<<center;
						center = str.str();	
					}else{
						std::stringstream str;
						str<<1<<":"<<center;
						center = str.str();	
					}
				}else{
					assert(p->getreference2()==cent_ref && l2==cent_left);
					if(p->getbegin2()<p->getend2()){
						std::stringstream str;
						str<<0<<":"<<center;
						center = str.str();	
					}else{
						std::stringstream str;
						str<<1<<":"<<center;
						center = str.str();	
					}
				}
				std::cout<< "center: "<< center<<std::endl;//Direction was added to center
				centersOnSequence.at(seq_id).insert(std::make_pair(n,center));
				centerOnseq.at(seq_id).insert(std::make_pair(n,center));//This is the global one. It is initialized in main function
				std::map<std::string,int>::iterator id=center_id.find(center);
				assert(id != center_id.end());
				std::map<int,bool>::iterator it = if_center.find(id->second);
				assert(it!= if_center.end());
				it->second = true;
			}
		}
	}
	void finding_centers::find_long_centers(size_t & seq_id){
		size_t last_position;
		size_t first_position;
		std::map<size_t, std::vector<int> > AllConnectedOnes;
		vector<int> vectorOfcenters;
		std::vector<size_t> position_of_each_suffix; //First one is position, second of is center
		bool direction = true;
		size_t dir;
		size_t dir1;
		for(std::map<size_t, std::string>::iterator it = centersOnSequence.at(seq_id).begin(); it != centersOnSequence.at(seq_id).end();it++){
			std::map<size_t, pw_alignment*>::iterator it1=AlignmentsFromClustering.at(seq_id).find(it->first);
			pw_alignment * al = it1->second;
			size_t l1,l2,r1,r2;
			al->get_lr1(l1,r1);
			al->get_lr2(l2,r2);
			if(l1 == it->first && al->getreference1()== seq_id){
				if(al->getbegin1()<al->getend1()){
					dir1 = 1;
				}else{
					dir1 = 0;
				}
			}
			if(l2 == it->first && al->getreference2() == seq_id){
				if(al->getbegin2()< al->getend2()){
					dir1 = 1;
				}else{
					dir1 = 0;
				}
			}
			if(vectorOfcenters.size()==0){
				last_position = it->first;
				first_position = it->first;
				dir = dir1;
			}
			if(dir ==dir1){
				direction = true;
			}else{
				direction = false;
			}
			std::cout << "position " << it->first << " allowed gap "<< ALLOWED_GAP<<std::endl;
			if((it->first - last_position) < ALLOWED_GAP && direction == true){
				if(vectorOfcenters.size()!=0){
					std::cout << "smaller than ALLOWED_GAP! "<<std::endl;
				}
				std::map<std::string,int>::iterator id=center_id.find(it->second);
				assert(id != center_id.end());
				size_t index = id->second;
				vectorOfcenters.push_back(index);
				position_of_each_suffix.push_back(it->first);
				std::cout << "index "<< index << std::endl;		
			}else{
				std::cout<< "bigger than ALLOWED_GAP: " <<std::endl;
				for(size_t i =0; i < vectorOfcenters.size();i++){
					std::cout << vectorOfcenters.at(i)<< " ";
				}
				std::cout << "-------"<<std::endl;
				AllConnectedOnes.insert(std::make_pair(first_position, vectorOfcenters));
				//Creates suffix and their related position , then clear the 'position_of_each_suffix' 
				if(vectorOfcenters.size()>1){
					for(size_t j =0; j < vectorOfcenters.size(); j++){
						std::vector<int> suffix;
						for(size_t k =j; k < vectorOfcenters.size();k++){
							suffix.push_back(vectorOfcenters.at(k));
						}
						initial_suffixes.at(seq_id).insert(std::make_pair(suffix,position_of_each_suffix.at(j)));
					}
				}
				position_of_each_suffix.clear();
				vectorOfcenters.clear();
			}
			if(l1 == it->first && al->getreference1()== seq_id){
				last_position = r1;
				if(al->getbegin1()<al->getend1()){
					dir = 1;
				}else{
					dir = 0;
				}
			}
			if(l2 == it->first && al->getreference2() == seq_id){
				last_position = r2;
				if(al->getbegin2()< al->getend2()){
					dir = 1;
				}else{
					dir = 0;
				}
			}
			std::cout << "last position "<< last_position << std::endl;
		}
		//If the last center had ALLOWED_GAP with its previous one and there is no more left to go through the else
		AllConnectedOnes.insert(make_pair(first_position, vectorOfcenters));
		if(vectorOfcenters.size()>1){
			for(size_t j =0; j < vectorOfcenters.size(); j++){
				std::vector<int> suffix;
				for(size_t k =j; k < vectorOfcenters.size();k++){
					suffix.push_back(vectorOfcenters.at(k));
				}
				initial_suffixes.at(seq_id).insert(std::make_pair(suffix,position_of_each_suffix.at(j)));
			}
		}
		std::vector<std::vector<int> > temp;
		for(std::map<size_t, std::vector<int> >::iterator it = AllConnectedOnes.begin(); it != AllConnectedOnes.end(); it++){
			if(it->second.size() != 1 && it->second.size() != 0){
				std::cout<< "pos "<< it->first <<std::endl;
				centersOfASequence.at(seq_id).insert(std::make_pair(it->first, it->second));//connected centers that can be potential long centers
				temp.push_back(it->second);
				for(size_t i =0; i < it->second.size();i++){
					std::cout << it->second.at(i)<< " ";
				}
				std::cout << "-------"<<std::endl;
			}
		}
		long_centers_of_a_sequence.insert(std::make_pair(seq_id, temp));
		if(seq_id == 28){
			std::cout << "initial suffix"<<std::endl;
			for(std::multimap<std::vector<int>,size_t >::iterator it = initial_suffixes.at(seq_id).begin(); it != initial_suffixes.at(seq_id).end(); it++){
				std::cout<< it->first << std::endl;
			}
		}
	}
	std::multimap< size_t, std::string> finding_centers::get_sequence_centers(size_t& id)const{
		return centersOnSequence.at(id);
	}
	std::string finding_centers::find_center_name(int & centerIndex)const{
		std::map<int,std::string>::const_iterator it= center_index.find(centerIndex);
		assert(it != center_index.end());
		return it->second;
	}
	std::vector<size_t> finding_centers::get_long_center_position(size_t & seq_id , std::vector<std::string> & long_center){
		std::vector<size_t> position;
		std::vector<int> centers;
		std::cout << "current long cent " <<std::endl;

		for(size_t i =0; i < long_center.size();i++){
			std::map<std::string, int>::iterator it = center_id.find(long_center.at(i));
			assert(it != center_id.end());
			std::cout << long_center.at(i)<< " "<< it->second <<std::endl;
			centers.push_back(it->second);
		}
		std::multimap<std::vector<int>,size_t>::iterator suf = initial_suffixes.at(seq_id).find(centers);
		if(suf != initial_suffixes.at(seq_id).end()){
			std::pair<std::multimap<std::vector<int>,size_t>::iterator , std::multimap<std::vector<int>,size_t>::iterator > it = initial_suffixes.at(seq_id).equal_range(centers);
			for(std::multimap<std::vector<int>,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				position.push_back(it1->second);
			}	
			return position;
		}else{// Need to loop over the initial suffixes multimap which is not so efficient TODO think of a better way
			for(std::multimap<std::vector<int>,size_t>::iterator it = initial_suffixes.at(seq_id).begin() ; it != initial_suffixes.at(seq_id).end();it++){
				if(it->first.at(0)!=centers.at(0) || it->first.size()<centers.size()){
					continue;
				}
				if(it->first.at(0)==centers.at(0) && it->first.size()>centers.size()){
					bool isNot = false;
					for(size_t i =0; i < centers.size();i++){
						if(centers.at(i) != it->first.at(i)){
							isNot = true;
							break;
						}						
					}
					if(isNot == false){
						position.push_back(it->second);
					}
				}
			}
		}
		return position;
	}
	std::map<size_t, std::vector<int> >finding_centers::get_center(size_t & seq_id)const{
		return centersOfASequence.at(seq_id);
	}
	const std::vector<std::vector<int> > & finding_centers::get_centers(size_t & seq_id)const{
		std::vector<std::vector< int> > empty;
		std::map<size_t , std::vector<std::vector<int> > >::const_iterator it =long_centers_of_a_sequence.find(seq_id);
	//	assert(it != long_centers_of_a_sequence.end());
		if(it != long_centers_of_a_sequence.end()){
			return it->second;
		}else return empty;
	}
	size_t finding_centers::get_number_of_centers()const{
		return center_index.size();
	}
	const std::map<std::string, int> finding_centers::get_index()const{
		return center_id;
	}




	suffixTree::suffixTree(size_t & num_seq, finding_centers & cent, std::vector<std::vector<std::vector<int> > > & words):centers(cent){
		this->num_seq = num_seq;
		this->words =words;
		last_node_index=0;
		edges.insert(make_pair(last_node_index,root));
		last_node_index ++;
	}

	suffixTree::~suffixTree(){}

	void suffixTree::add_leaves(size_t & parent , std::vector<int> & new_suffix){
	//	double inf = std::numeric_limits<double>::infinity();
		std::cout<< "last node index "<< last_node_index << " new suffix size "<< new_suffix.size() << std::endl;
		edges_relations.insert(std::make_pair(parent, last_node_index));
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index++;
		std::cout << "relation size "<<edges_relations.size()<< std::endl;
	}
	void suffixTree::add_internal_node(size_t & parent_index , std::vector<int> & new_suffix){//equal range for parent
		std::multimap<size_t,size_t>::iterator it2 = edges_relations.find(parent_index);
		assert(it2 != edges_relations.end());
		std::cout << last_node_index << " " << parent_index <<std::endl;
		std::vector<size_t> intermediate;
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			intermediate.push_back(it1->second);
		}
		edges_relations.erase(parent_index);
		it2 = edges_relations.find(parent_index);
		assert(it2 == edges_relations.end());
		std::cout << "child nodes were deleted"<<std::endl;
		edges_relations.insert(std::make_pair(parent_index, last_node_index));
		std::cout<< "intermediate size "<< intermediate.size()<<std::endl;
		for(size_t i =0; i < intermediate.size();i++){
			edges_relations.insert(std::make_pair(last_node_index, intermediate.at(i)));						
			std::cout << "test at "<< i << " is " << intermediate.at(i)<<std::endl;
		}
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index ++;

	}
	void suffixTree::update_node(size_t & node_index, std::vector<int> & new_content){
		std::cout<< "node index "<< node_index << " new content "<< new_content << std::endl;
		std::map<size_t, std::vector<int> >::iterator it = edges.find(node_index);
		assert(it != edges.end());
		it->second = new_content;

	}
	void suffixTree::split(size_t & parent_index, std::vector<int> & parent , std::vector<int> & new_suffix, std::vector<int> & common_part, bool & make_a_tree){
		assert(common_part.size() <= parent.size());
		std::vector<int> parent_non_common;
		std::cout << "common part size "<< common_part.size() << " parent size "<< parent.size()<<std::endl;
		for(size_t i =common_part.size(); i < parent.size();i++){
			parent_non_common.push_back(parent.at(i));
		}
		if(parent_non_common.size() != 0){//common part shorter than parent's length
			std::cout << "if parent_non_common.size() != 0 "<<std::endl;
			std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
			if(it == edges_relations.end()){
				std::cout<< "add leaves "<<std::endl;
				add_leaves(parent_index,parent_non_common); 
			}else{
				std::cout<< "add internal "<<std::endl;
				add_internal_node(parent_index,parent_non_common);
			}
			std::cout<< "parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() << " " << parent_non_common <<std::endl;
			update_node(parent_index,common_part);//update its content
			std::cout<< "1parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() <<std::endl;
			if(new_suffix.size()>= common_part.size()){
				std::vector<int> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
					std::cout << new_suffix.at(i) << " " ;
					current_non_common.push_back(new_suffix.at(i));
				}
				std::cout<< "current non common size is "<< current_non_common <<std::endl;
				if(current_non_common.size() != 0){
					std::cout<< "add leaves1 "<<std::endl;
					add_leaves(parent_index,current_non_common); 
				}
			}
			make_a_tree = false;
		}else{//common part is equal to parent's length - the new suffix is longer than the parent and fully covers it.
			std::cout << "i am here!" <<std::endl;
			assert(parent_non_common.size() == 0);
			if(common_part.size() < new_suffix.size()){ //TRICKY ONE!!!!
				std::vector<int> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
					std::cout<< new_suffix.at(i) << " ";
					current_non_common.push_back(new_suffix.at(i));
				}
				std::cout << " "<<std::endl;
				assert(current_non_common.size()!=0);
				std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
				if(it != edges_relations.end()){//check on the childrens for a common part
					std::cout<< "has child node"<<std::endl;
					size_t new_parent_index = edges.size();
					find_next_parent(parent_index, current_non_common, common_part, new_parent_index, parent);
					std::cout << new_parent_index <<std::endl;
					if(new_parent_index != edges.size()){
						new_suffix = current_non_common;
						parent_index = new_parent_index;
						std::cout << "here "<< std::endl;
					}else{
						add_leaves(parent_index,current_non_common);
						make_a_tree = false;
					}
				}else{//add it as a leave
					assert(it == edges_relations.end());
					add_leaves(parent_index,current_non_common);
					make_a_tree = false;
				}

			}else{
				assert(common_part.size() == new_suffix.size());
				make_a_tree = false;
			}

		}
	}
	void suffixTree::make_suffix(std::vector<int> & word){
		suffixes.clear();
		for(size_t i =0; i < word.size();i++){
			std::cout << word.at(i) << " ";
			std::vector<int> suffix;
			for(size_t k =i; k < word.size();k++){
				suffix.push_back(word.at(k));
			}//For convenience add 1<<31 at the end of all the suffixes
			suffix.push_back(1<<31);
			suffixes.push_back(suffix);
		} 		std::cout << " "<<std::endl;	
	
	}
	void suffixTree::find_first_edges_after_root(size_t & deepest_parent, std::vector<int> & new_suffix){
		std::map<size_t, size_t >::iterator it = first_edges.find(new_suffix.at(0));
		if(it!= first_edges.end()){
			deepest_parent = it->second;
		}
	
	}
	void suffixTree::find_next_parent(size_t & parent_index, std::vector<int> & current_non_common, std::vector<int> & common_part, size_t & new_parent_index, std::vector<int> & parent){
		common_part.clear();
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			std::cout << "child "<< it1->second << std::endl;
			std::map<size_t, std::vector<int> >::iterator edge=edges.find(it1->second);
			assert(edge != edges.end());
			std::cout << edge->second.at(0)<< " "<<current_non_common.at(0)<<std::endl;
			if(edge->second.at(0)==current_non_common.at(0)){
				std::cout<< "got it"<< std::endl;
				find_common_part(edge->second, current_non_common,common_part);
				new_parent_index = edge->first;
				parent = edge->second;
				std::cout << "new parent index "<< new_parent_index << std::endl;
				break;
			}
		}


	}
	void suffixTree::find_children(std::vector<int> & parent ,std::vector<std::vector<int> >& children){


	}
	void suffixTree::find_common_part(std::vector<int> & parent , std::vector<int> & current, std::vector<int> & common_part){
		size_t length = current.size();
		if(length > parent.size()) length = parent.size();
		for(size_t i =0; i < length; i++){
			if(parent.at(i)==current.at(i)){
				common_part.push_back(parent.at(i));
			}else{
				break;
			}
		}
	}
	void suffixTree::make_tree_for_a_seq(size_t & seq){
		std::cout << "at seq "<< seq <<std::endl;
		for(size_t i = 0; i < words.at(seq).size(); i ++){
			std::vector<int> this_word = words.at(seq).at(i);//current string of centers on a seq
			make_suffix(this_word);
			std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
			for(size_t j =0; j< suffixes.size(); j++){
				size_t deepest_parent =edges.size();
				find_first_edges_after_root(deepest_parent, suffixes.at(j));
				if(deepest_parent == edges.size()){
					std::cout << "no first parent! "<<std::endl;
					first_edges.insert(std::make_pair(suffixes.at(j).at(0),last_node_index));
					size_t root_index =0;
					add_leaves(root_index,suffixes.at(j));
				}else{
					std::map<size_t, std::vector<int> >::iterator it = edges.find(deepest_parent);
					assert(it != edges.end());
					std::cout<< "has a parent"<<std::endl;
					std::vector<int> common_part;
					std::vector<int> parent = it->second;
					std::vector<int> current = suffixes.at(j);
					find_common_part(parent, current,common_part);
					for(size_t n =0; n < common_part.size();n++){
						std::cout<< common_part.at(n) << " ";
					}
					std::cout << " "<<std::endl;
					bool make_a_tree = true;
					size_t parent_index = deepest_parent;
					while(make_a_tree == true){
						split(parent_index,parent,current, common_part, make_a_tree);
					}
				}
			}
		}
//		print_tree();

	}
	void suffixTree::make_tree(std::vector<int> & center_with_highest_gain, int & highest_index){	
		std::cout << "words size is "<< words.size() <<std::endl;
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			if(center_with_highest_gain.size() == 0){
				std::cout<< "seq id "<<seq_id <<std::endl;
			//	words.at(seq_id) = centers.get_centers(seq_id);
				std::cout << "its words size is "<< words.at(seq_id).size() <<std::endl;
			}else{
				update_centers(seq_id, center_with_highest_gain, highest_index);
			}
			make_tree_for_a_seq(seq_id);
		}
		count_branches();
	}
	void suffixTree::update_centers(size_t & seq_id, std::vector<int> & edge_with_highest_gain, int & highest_index){
		for(size_t i = 0; i < words.at(seq_id).size(); i++){
			std::vector<int> centers = words.at(seq_id).at(i);
			for(size_t k =0; k < centers.size();k++){
				if(centers.at(k) == edge_with_highest_gain.at(0)&& ((centers.size()-k)>= edge_with_highest_gain.size())){
					size_t first_common_index;
					std::vector<int> common;
					for(size_t i = 0; i < edge_with_highest_gain.size();i++){
						first_common_index = k;
						if(centers.at(i+k)== edge_with_highest_gain.at(i)){
							common.push_back( edge_with_highest_gain.at(i));
						}else break;
					}
					if(common.size() == edge_with_highest_gain.size()){
						std::vector<int> new_center;
						for(size_t m = 0; m < first_common_index ;m++){
							new_center.push_back(centers.at(m));
						}
						new_center.push_back(highest_index);
						for(size_t m = first_common_index + common.size(); m < centers.size();m++){
							new_center.push_back(centers.at(m));
						}
						centers=new_center;
					}
				}
			}
			words.at(seq_id).at(i) = centers;
		}

	}
	void suffixTree::count_branches(){
		for(size_t seq =0; seq < num_seq; seq++){
			for(size_t w = 0; w< words.at(seq).size(); w++){
				std::vector<int> this_word = words.at(seq).at(w);//current string of centers on a seq
				make_suffix(this_word);
				std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
				for(size_t j =0; j< suffixes.size(); j++){
					std::cout << "j "<<j <<std::endl;
					std::vector<size_t> branch;
					std::vector<int> current = suffixes.at(j);
					size_t deepest_parent =edges.size();
					find_first_edges_after_root(deepest_parent,current);
					assert(deepest_parent != edges.size());
					std::map<size_t, std::vector<int> >::iterator it=edges.find(deepest_parent);
					assert(it != edges.end());
					if(it->second.size()==current.size()){
						branch.push_back(deepest_parent);
						std::cout << "all the suffix is on one edge!"<<std::endl;
					}else{//Note that first parent can not be longer than current.
						std::vector<int> updated_current;
						std::cout<<"updated current: "<<std::endl;
						for(size_t k =it->second.size() ; k < current.size(); k++){
							std::cout << current.at(k)<< " ";
							updated_current.push_back(current.at(k));
						}
						std::cout << " "<<std::endl;
						size_t parent_index = deepest_parent;
						branch.push_back(parent_index);
						current = updated_current;
						while(current.size()>0){
							traverse_the_tree(parent_index,current);
							std::cout << "current size "<< current.size() <<std::endl;
							branch.push_back(parent_index);							
						}
					}
					std::vector<size_t> sub_branch;
					std::cout<<"branch size "<< branch.size()<<std::endl;
					for(size_t i = 0 ; i < branch.size(); i++){
						std::cout << "i "<< i << " " << branch.at(i) <<std::endl;
						sub_branch.push_back(branch.at(i));
						std::map<std::vector<size_t>,size_t>::iterator it1 = branches.find(sub_branch);
						if(it1 == branches.end()){
							branches.insert(make_pair(sub_branch,0));
							it1 = branches.find(sub_branch);
						}
						it1->second = it1->second +1;
					}

				}
			}
		}
	}
	void suffixTree::traverse_the_tree(size_t & parent_index, std::vector<int> & new_suffix){
		std::cout<< "parent index is "<< parent_index <<std::endl;
		std::multimap<size_t, size_t>::iterator check = edges_relations.find(parent_index);
		assert(check != edges_relations.end());
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator > it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t, size_t>::iterator it1 = it.first; it1 !=it.second;it1++){
			std::map<size_t,std::vector<int> >::iterator it2 = edges.find(it1->second);
			std::cout<< it1->second <<std::endl;
			assert(it2 != edges.end());
			if(it2->second.at(0)== new_suffix.at(0)){
				parent_index = it1->second;
				std::cout<< "new parent index " << parent_index<<std::endl;
				std::vector<int> updated_current;
				std::cout<< "updated:"<<std::endl;
				for(size_t i = it2->second.size(); i < new_suffix.size(); i++){
					std::cout<< new_suffix.at(i)<< " ";
					updated_current.push_back(new_suffix.at(i));
				}		
				std::cout << " "<<std::endl;
				new_suffix = updated_current;
				break;
			}
		}
		


	}
	const std::map<size_t, std::vector<int> > suffixTree::get_edges()const{
		return edges;
	}
	const std::map<std::vector<size_t>, size_t > suffixTree::get_branches()const{
		return branches;
	}

	const std::vector<std::vector<int> > suffixTree::get_current_centers(size_t & seqid)const{
		return words.at(seqid);
	}

	void suffixTree::print_tree(){
		for(std::map<size_t, std::vector<int> >::iterator it = edges.begin();it !=edges.end(); it++){
			std::cout<< " node is "<<it->first << " its edge content is ";
			for(size_t i = 0 ; i < it->second.size();i++){
				std::cout << it->second.at(i)<< " ";
			}
			std::cout<< " "<<std::endl;
			std::cout << "its children are "<<std::endl;
			std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it2 = edges_relations.equal_range(it->first);
			for(std::multimap<size_t,size_t>::iterator it1 = it2.first ; it1 != it2.second; it1++){
				std::cout << it1->second << " ";
			}
			std::cout << " " <<std::endl;
		}
	}



