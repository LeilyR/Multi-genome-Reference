#include "map_read.hpp"
void ref_graph::read_dot_file(std::string & refgraph, std::string & refacc){
	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
	std::ifstream dotin(refgraph.c_str());
	size_t counter =0;
	if(!dotin){
		std::cerr << "Error : Cannot open " << refgraph.c_str() << std::endl;
		exit(1);
	}
	else{
		std::string line;
		while(getline(dotin,line)){
			if(line[0] == '/'){
				continue;
			}
			else{
				std::vector<std::string> parts;
				strsep(line, " ", parts);
				std::string from = parts.at(0);
				std::string to;
				if(parts.size()>2){
					to = parts.at(2);
					std::cout<<"from " << from << " to " << to << std::endl;
				}
				std::cout << "ref acc "<< refacc << std::endl;
				std::string dir1(from.end()-1,from.end());
				std::string name1(from.begin(),from.end()-1);
				std::string longname1 = refacc + ":"+name1;
				std::cout<<longname1<<std::endl;
				std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(longname1);
				std::cout << "seq id is "<< findseq1->second << std::endl;
				if(findseq1==longname2seqidx.end()) {
					std::cerr << "Error: unknown sequence in dot File: " << longname1 << std::endl;
					exit(1);
				}
				size_t seqid1 = findseq1->second;
				if(parts.size() == 2){
					graph.set_vertex(seqid1, name1);
					counter ++;
					add_adjacencies(dir1,name1);
				}
				else{
					assert(to != "");
					std::string dir2(to.end()-1,to.end());
					std::string name2(to.begin(),to.end()-1);
					std::string longname2 = refacc + ":"+name2;
					add_adjacencies(dir1,dir2,name1,name2);
					add_adjacencies(dir2,name2);
					std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find(longname2);
					if(findseq2==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << longname2 << std::endl;
						exit(1);
					}
					size_t seqid2 = findseq2->second;
					//case 1 : RtoF
					if((dir1 == "-") && (dir2 == "+")){
						graph.addEdge(1,seqid1,name1,seqid2,name2);
						counter ++;
						continue;
					}
					//case 2 : RtoR
					if((dir1 == "-") && (dir2 == "-")){
						graph.addEdge(2,seqid1,name1,seqid2,name2);
						counter ++;
						continue;
					}
					//case 3 : FtoF
					if((dir1 == "+") && (dir2 == "+")){
						graph.addEdge(3,seqid1,name1,seqid2,name2);
						counter ++;
						continue;
					}
					//case 4 : FtoR
					if((dir1 == "+") && (dir2 == "-")){
						graph.addEdge(4,seqid1,name1,seqid2,name2);
						counter ++;
						continue;
					}
				}
			}
		}
	}
	std::cout << "number of aedges is "<< counter << std::endl;
}

void ref_graph::add_adjacencies(std::string & dir1 , std::string & dir2 , std::string & name1 , std::string & name2){
	int node1 = std::stoi(name1);
	int node2 = std::stoi(name2);
	if((dir1 == "+") && (dir2 == "+")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(node1);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(node1, std::set<int>()));
			adj = adjacencies.find(node1);
		}
		adj->second.insert(node2);
	}
	if((dir1 == "-") && (dir2 == "+")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*node1);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));
			adj = adjacencies.find(-1*node1);
		}
		adj->second.insert(node2);
	}
	if((dir1 == "+") && (dir2 == "-")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(node1);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(node1, std::set<int>()));
			adj = adjacencies.find(node1);
		}
		adj->second.insert(-1*node2);
	}
	if((dir1 == "-") && (dir2 == "-")){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*node1);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(-1*node1, std::set<int>()));
			adj = adjacencies.find(-1*node1);
		}
		adj->second.insert(-1*node2);
	}
}

void ref_graph::add_adjacencies(std::string & dir , std::string & name){
	int cent_id = std::stoi(name);
	if(dir == "+"){
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(cent_id);
		if(adj == adjacencies.end()){
			adjacencies.insert(std::make_pair(cent_id, std::set<int>()));
		}
	}else{
		std::map<int, std::set<int> >::iterator adj = adjacencies.find(-1*cent_id);
		if(adj == adjacencies.end()){
				adjacencies.insert(std::make_pair((-1*cent_id), std::set<int>()));
		}
	}
}
const std::map<int , std::set<int> > & ref_graph::get_adjacencies()const{
	return adjacencies;
}
size_t deep_first::seq_length(std::string & seq_name, std::string & accname){
	std::string longname = accname +":"+seq_name;
	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
	std::map<std::string, size_t>::iterator findseq = longname2seqidx.find(longname);
	assert(findseq != longname2seqidx.end());
	return data.get_seq_size(findseq->second);
	


}
std::string deep_first::seqname(int & node){
	std::string temp;
	std::stringstream ss;
	if(node< 0){
		ss << (-1*node);
		temp = ss.str();
	}else{
		ss << node;
		temp = ss.str();
	}
	return temp;
}

void deep_first::deep_first_search(int & startnode, std::string & refacc){
	paths.clear();
	nodes_on_paths.clear();
	path_length.clear();
	std::string seq_name = seqname(startnode);
	size_t seqlength = seq_length(seq_name, refacc);

	int accu_length = seqlength;
	std::map<int, bool> visited;
	for(std::map<int, std::set<int> >::iterator it = adjacencies.begin() ; it != adjacencies.end() ;it++){
      		visited.insert(std::make_pair(it->first,false));
 	}
	//Recurssive dfs: 
	std::vector<int> apath;
	size_t parent_length = 0;

    	look_for_neighbors(startnode, visited, refacc,accu_length,apath,parent_length);
}
void deep_first::look_for_neighbors(int & node, std::map<int,bool> & visited , std::string & refacc, int & accu_length, std::vector<int> & apath, size_t & parent_length){
	std::cout << "this node is "<< node << std::endl;
	std::map<int,bool>::iterator it = visited.find(node);
	assert(it != visited.end());
	it->second = true;
	std::string seq_name = seqname(node);
	size_t seqlength = seq_length(seq_name, refacc);
	std::cout << "its length is "<< seqlength<<std::endl;
	if(accu_length<=MAXGAP+1000){
		apath.push_back(node);	
		nodes_on_paths.insert(node);
	}
	std::cout << "apath size "<< apath.size() <<std::endl;
	//Look for the adjacent vertices: 
	std::map<int, std::set<int> >::iterator it1 = adjacencies.find(node);
	assert(it1 != adjacencies.end());
	std::cout<< "it has "<< it1->second.size() << " adjacents" <<std::endl;
	if(it1->second.size()==0){
		if(accu_length <= MAXGAP+1000){
			paths.push_back(apath);
			apath.pop_back();
			nodes_on_paths.erase(node);
		}
		accu_length -=seqlength;
	}
	std::cout << "accu length is "<< accu_length << std::endl; 
	std::cout<< "all the adjs "<<std::endl;
	for (std::set<int>::iterator adj = it1->second.begin(); adj != it1->second.end(); adj++){
		std::cout << *adj << " ";
	}
	std::cout<<std::endl;
	for (std::set<int>::iterator adj = it1->second.begin(); adj != it1->second.end(); adj++){
		int adjacent = *adj;
		std::string adj_name = seqname(adjacent);
		size_t adjlength = seq_length(adj_name, refacc);
		accu_length += adjlength;
		if(adj == it1->second.begin()){
			parent_length = seqlength;
			std::cout << "parent length "<< parent_length <<std::endl;
		}
		std::cout<< "adj is "<< *adj << std::endl;
		std::set<int>::iterator adj1 =it1->second.end();
		std::map<int,bool>::iterator it2 = visited.find(*adj);
		assert(it2 != visited.end());
		if (it2->second == false && accu_length <= MAXGAP+1000){
			std::cout << "smaller"<<std::endl;
			int current_node = *adj;
          		look_for_neighbors(current_node, visited,refacc,accu_length,apath,parent_length);
		}
		if(it2->second == false && accu_length > MAXGAP+1000){
			std::cout << "here! " <<std::endl;
			int current_node = *adj;
			it2->second = true; 
			if(apath.size()>1){
				paths.push_back(apath);
				std::cout<< "apath size: "<< apath.size()<<std::endl;
				for(size_t i =0; i < apath.size();i++){
					std::cout << apath.at(i)<< " ";
				}
				std::cout<<" "<<std::endl;
			}
			accu_length -=adjlength;
			std::set<int>::iterator adj2 =it1->second.end();
			if(--adj2 == adj){
				std::cout<<"end of adjs"<<std::endl;
				std::cout << "accu before "<< accu_length <<std::endl;
				accu_length -= parent_length;
				std::cout <<"parent length "<<parent_length<<std::endl;
				std::cout<< "accu after "<< accu_length << std::endl;
			}

		}
	}
	apath.pop_back();
	std::cout<<"paths size is "<<paths.size()<<std::endl;
}

const std::vector<std::vector<int> > deep_first::get_paths()const{
	return paths;
}
const std::set<int> deep_first::get_nodes()const{
	return nodes_on_paths;
}
void als_components::find_als_on_paths(){
	for(size_t i =0; i < alignments.size(); i++){
		std::set<pw_alignment,sort_pw_alignment> copy_al = alignments.at(i);
		std::cout << "copy al size "<< copy_al.size() << std::endl;
		for(std::set<pw_alignment,sort_pw_alignment>::iterator it = alignments.at(i).begin() ; it != alignments.at(i).end(); it++){
			std::set<pw_alignment,sort_pw_alignment> temp;
			pw_alignment p = *it;
			p.print();
			size_t l,r;
			p.get_lr2(l,r);
			for(std::set<pw_alignment,sort_pw_alignment>::iterator it1 = it ; it1 != alignments.at(i).end(); it1++){
				pw_alignment al = *it1;
				size_t l2,r2;
				al.get_lr2(l2,r2);
				al.print();
				if((l2<r)||(l2-r <=MAXGAP)){
					temp.insert(al);
				}else{
					break;
				}
			}
			std::cout<<"temp size is "<< temp.size()<<std::endl;
			std::vector<std::vector<int> > paths;
			get_paths(p,paths);
			std::set<int> nodes = dfirst.get_nodes();
			look_for_neighbors_on_paths(nodes,temp);
		}

	}

}

void als_components::get_paths(const pw_alignment & p, std::vector<std::vector<int> > & paths){
	size_t acc = data.accNumber(p.getreference1());
	std::string ref_accession = data.get_acc(acc);
	std::string seqname = data.get_seq_name(p.getreference1());
	int name = std::stoi(seqname);
	if(p.getbegin1()>p.getend1()){
		name = -1*name;
	}
	std::cout << "name is " << name <<std::endl;
	dfirst.deep_first_search(name, ref_accession);
	paths = dfirst.get_paths();
	std::cout<< "path size is "<< paths.size()<<std::endl;
	for(size_t i =0; i < paths.size(); i++){
		for(size_t j =0; j < paths.at(i).size();j++){
			std::cout << paths.at(i).at(j)<< " ";
		}
		std::cout<< " "<<std::endl;
	}
}
void als_components::look_for_neighbors_on_paths(std::set<int>& nodes , std::set<pw_alignment,sort_pw_alignment> & neighbors){
	std::cout<<"look for neighbors on paths "<<std::endl;
//should check if the first reference of the neighbors is one of the ints in the paths
	std::set<pw_alignment,sort_pw_alignment> temp;
	for(std::set<pw_alignment, sort_pw_alignment>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
		pw_alignment p = *it;
		std::string seqname = data.get_seq_name(p.getreference1());
		int name = std::stoi(seqname);
		if(p.getbegin1()>p.getend1()){
			name = -1*name;
		}
		std::cout << "name is " << name <<std::endl;
		std::set<int>::iterator it1 = nodes.find(name);
		if(it1 == nodes.end()){//if the neighbor is not aligned on any of nodes
			temp.insert(p);
		}else{//make an al using needleman wunsch

		}

	}





}
