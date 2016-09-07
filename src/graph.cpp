/*
 * needleman.hpp
 *
 *  Created on: May 12, 2015
 *      Author: Marion Dubarry
 */

#include "graph.hpp"

//bool VERBOSE = false;

// ---------------------------
// 		Vertex
// ---------------------------

// -------------------------------------------Functions

Graph::Vertex::Vertex(size_t id):idVertex(id){
	name = "";
	distance = std::numeric_limits<double>::infinity();
//	startOnRead = -1;
//	endOnRead = -1 ;
//	startOnNode = -1;
//	endOnNode = -1;
}
Graph::Vertex::Vertex(const Graph::Vertex& v){
	idVertex = v.idVertex;
	name = v.name;
	costScore = v.costScore;
	startOnRead = v.startOnRead;
	endOnRead = v.endOnRead;
	startOnNode = v.startOnNode;
	endOnNode = v.endOnNode;
	distance = v.distance;
	FtoF = v.FtoF;
	RtoR = v.RtoR;
	FtoR = v.FtoR;
	RtoF = v.RtoF;
	previousFtoF = v.previousFtoF;
	previousRtoR = v.previousRtoR;
	previousFtoR = v.previousFtoR;
	previousRtoF = v.previousRtoF;
}

//-------------
//Split a string
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
 }
//--------------------------------------------------------- Set and Get

void Graph::Vertex::setRtoF(const size_t& v){
	if((std::find(RtoF.begin(),RtoF.end(),v)) == RtoF.end())
		RtoF.push_back(v);
}

void Graph::Vertex::setRtoR(const size_t& v){
	if((std::find(RtoR.begin(),RtoR.end(),v) == RtoR.end()))
		RtoR.push_back(v);
}

void Graph::Vertex::setFtoF(const size_t& v){
	if((std::find(FtoF.begin(),FtoF.end(),v) == FtoF.end()))
		FtoF.push_back(v);
}

void Graph::Vertex::setFtoR(const size_t& v){
	if((std::find(FtoR.begin(),FtoR.end(),v) == FtoR.end()))
		FtoR.push_back(v);
}
void Graph::Vertex::setPreviousRtoF(const int& v){
	if((std::find(previousRtoF.begin(),previousRtoF.end(),v)) == previousRtoF.end())
		previousRtoF.push_back(v);
}

void Graph::Vertex::setPreviousRtoR(const int& v){
	if((std::find(previousRtoR.begin(),previousRtoR.end(),v) == previousRtoR.end()))
		previousRtoR.push_back(v);
}

void Graph::Vertex::setPreviousFtoF(const int& v){
	if((std::find(previousFtoF.begin(),previousFtoF.end(),v) == previousFtoF.end()))
		previousFtoF.push_back(v);
}

void Graph::Vertex::setPreviousFtoR(const int& v){
	if((std::find(previousFtoR.begin(),previousFtoR.end(),v) == previousFtoR.end()))
		previousFtoR.push_back(v);
}
void Graph::Vertex::setCostScore(const int& id, const double& score){
	costScore.insert(std::make_pair(id,score));
}

void Graph::Vertex::setDistance(const double& d){distance = d;}
void Graph::Vertex::setStartOnRead(const int& start){startOnRead = start;}
void Graph::Vertex::setEndOnRead(const int& end){endOnRead = end;}
void Graph::Vertex::setStartOnNode(const int& start){startOnNode = start;}
void Graph::Vertex::setEndOnNode(const int& end){endOnNode = end;}
void Graph::Vertex::setName(const std::string& n){name =n;}
const size_t& Graph::Vertex::getIdVertex()const{return idVertex;}
const std::vector<int>& Graph::Vertex::getRtoF()const{return RtoF;}
const std::vector<int>& Graph::Vertex::getRtoR()const{return RtoR;}
const std::vector<int>& Graph::Vertex::getFtoF()const{return FtoF;}
const std::vector<int>& Graph::Vertex::getFtoR()const{return FtoR;}
const std::vector<int>& Graph::Vertex::getPreviousRtoF()const{return previousRtoF;}
const std::vector<int>& Graph::Vertex::getPreviousRtoR()const{return previousRtoR;}
const std::vector<int>& Graph::Vertex::getPreviousFtoF()const{return previousFtoF;}
const std::vector<int>& Graph::Vertex::getPreviousFtoR()const{return previousFtoR;}

double Graph::Vertex::getCostScore(const int& id)const{
	//First try to get the score of the edge, if does not exist get the default score ; which mean the node was a seed
	std::map<int, double>::const_iterator it = costScore.find(id);
	if(it == costScore.end()){
		std::cout << "error cannot find id " << id << std::endl;
		std::map<int, double>::const_iterator it2 = costScore.find(idVertex);
		if(it2 == costScore.end())
			exit(1);
		else{
			std::cout << " score " << it->second << " " << it2->second<<std::endl;
			return it2->second;
		}
	}
	else
		return it->second;
}
double Graph::Vertex::getCostScore()const{return (costScore.find(idVertex))->second;}
const double& Graph::Vertex::getDistance()const{return distance;}
const size_t& Graph::Vertex::getStartOnRead()const{	return startOnRead;}
const size_t& Graph::Vertex::getEndOnRead()const{return endOnRead;}
const size_t& Graph::Vertex::getStartOnNode()const{return startOnNode;}
const size_t& Graph::Vertex::getEndOnNode()const{return endOnNode;}
const std::string& Graph::Vertex::getName()const{return name;}

//-------------------------------
//		oriented_vertex
//-------------------------------
Graph::OrientedVertex::OrientedVertex(const Graph::Vertex& v, const bool& b): vertex(v),isForward(b) {
}

Graph::OrientedVertex::OrientedVertex(const Graph::Vertex& v, const pw_alignment& p): vertex(v.getIdVertex()) {
	//not so good copy of object
	vertex = v;
	if(p.getbegin1() < p.getend1())
		isForward = true;
	else
		isForward = false;
}

void Graph::OrientedVertex::getSuccessors(std::vector<Graph::OrientedVertex>& oriented,Graph g) const{
	if(isForward == true){
		for( uint it =0; it < vertex.getFtoF().size(); ++it){
			size_t id = vertex.getFtoF()[it];
			std::map<size_t,Vertex>::const_iterator it2 = g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cout << "error to find vertex in getSuccessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,true);
			oriented.push_back(orientedObject);
		}
		for( uint it =0; it < vertex.getFtoR().size(); ++it){
			size_t id = vertex.getFtoR()[it];
			std::map<size_t,Vertex>::const_iterator it2 = g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cout << "error to find vertex in getSuccessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,false);
			oriented.push_back(orientedObject);
		}
	}
	else{
		for( uint it =0; it < vertex.getRtoR().size(); ++it){
			size_t id = vertex.getRtoR()[it];
			std::map<size_t,Vertex>::const_iterator it2 =g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cout << "error to find vertex in getSuccessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,false);
			oriented.push_back(orientedObject);
				}
		for( uint it =0; it < vertex.getRtoF().size(); ++it){
			size_t id = vertex.getRtoF()[it];
			std::map<size_t,Vertex>::const_iterator it2 = g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cerr << "error to find vertex in getSuccessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,true);
			oriented.push_back(orientedObject);
		}
	}
}

void Graph::OrientedVertex::getPredecessors(std::vector<OrientedVertex>& oriented,Graph g) const{
	if(isForward == true){
		for( uint it =0; it < vertex.getPreviousFtoF().size(); ++it){
			size_t id = vertex.getPreviousFtoF()[it];
			std::map<size_t,Vertex>::const_iterator it2 =g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cerr << "error to find vertex in getPredecessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,true);
			oriented.push_back(orientedObject);
		}
		for( uint it =0; it < vertex.getPreviousRtoF().size(); ++it){
			size_t id = vertex.getPreviousRtoF()[it];
			std::map<size_t,Vertex>::const_iterator it2 =g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cerr << "error to find vertex in getPredecessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,false);
			oriented.push_back(orientedObject);
		}
	}
	else{
		for( uint it =0; it < vertex.getPreviousRtoR().size(); ++it){
			size_t id = vertex.getPreviousRtoR()[it];
			std::map<size_t,Vertex>::const_iterator it2 =g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cerr << "error to find vertex in getPredecessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,false);
			oriented.push_back(orientedObject);
		}
		for( uint it =0; it < vertex.getPreviousFtoR().size(); ++it){
			size_t id =vertex.getPreviousFtoR()[it];
			std::map<size_t,Vertex>::const_iterator it2 = g.getVertices().find(id);
			if(it2 == g.getVertices().end()){
				std::cerr << "error to find vertex in getPredecessors " <<id<<std::endl;
				exit(1);
			}
			OrientedVertex orientedObject(it2->second,true);
			oriented.push_back(orientedObject);
		}
	}
}

Graph::Vertex& Graph::OrientedVertex::getVertex(){return vertex;}
const bool& Graph::OrientedVertex::getIsForward()const{return isForward;}

std::string Graph::OrientedVertex::getDna(all_data & data, size_t & idOV) const{
	std::string sequence;
//	std::string tmpAccName1 = "noAcc:"+ vertex.getName();
//	std::map<std::string, size_t> longname2seqidx = data.getLongname2seqidx();
//	std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(tmpAccName1);
//	if(findseq1==data.getLongname2seqidx().end()){
//		std::cerr << "Error: unknown sequence in getDna(): " << tmpAccName1 << std::endl;
//		exit(1);
//	}
//	size_t idSeq = findseq1->second;
	dnastring tmpSequence = data.getSequence(idOV);
	if(isForward == true){
		//return sequence = convertDnastringToString(tmpSequence); //cannot call member
		for(uint i = 0 ; i < data.get_seq_size(idOV); ++i){
			sequence += tmpSequence.at(i);
		}
		return sequence;
	}
	else{
		for(uint i = 0 ; i < data.get_seq_size(idOV); ++i){
			sequence += dnastring::complement(tmpSequence.at(i));
		}
		return sequence;
	}
}
// ---------------------------
// 		Graph
// ---------------------------

// -------------------------------------------Functions
void Graph::setStartOnRead(const size_t&id, const size_t&start){
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id);
	if(it == vertices.end()){
			std::cerr << "cannot find node in setstart " << std::endl;
			exit(0);
	}
	it->second.setStartOnRead(start);
}
void Graph::setStartOnNode(const size_t&id, const size_t&start){
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id);
	if(it == vertices.end()){
			std::cerr << "cannot find node in setstart " << std::endl;
			exit(0);
	}
	it->second.setStartOnNode(start);
}
void Graph::setEndOnNode(const size_t&id, const size_t&end){
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id);
		assert(it != vertices.end());
		it->second.setEndOnNode(end);
}
void Graph::setEndOnRead(const size_t&id, const size_t&end){
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id);
		assert(it != vertices.end());

		it->second.setEndOnRead(end);
}
const std::map<size_t,Graph::Vertex>& Graph::getVertices()const{return vertices;}

void Graph::set_vertex(const size_t & seqid , const std::string & seqname){
	std::cout << "set vertex! "<<std::endl;
	Vertex v1= Vertex(seqid);
	v1.setName(seqname);
	addVertex(v1);
}


void Graph::setScore(const size_t& name, const double& score,const size_t& idOrigin){

	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(name);
	if(it == vertices.end()){
		std::cerr << "cannot find node in setScore " << std::endl;
		exit(0);
	}
	it->second.setCostScore(idOrigin,score);
}
const Graph::Vertex& Graph::getVertex(const size_t& id){
	std::cout << "id " << id << std::endl;
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id);
	if(it == vertices.end()){
			std::cerr << "cannot find node in getVertex " << std::endl;
			exit(0);
	}
	return it->second;
}

void Graph::read_dot_file(std::string & dotfile, std::map< std::string, size_t> & longname2seqidx, std::string & ref_accession){
	std::ifstream dotin(dotfile.c_str());
	size_t counter =0;
	if(!dotin){
		std::cerr << "Error : Cannot open " << dotfile.c_str() << std::endl;
		exit(1);
	}
	else{
		std::string line;
		while(getline(dotin,line)){
			if(line[0] == '/'){
				continue;
			}
			else{
				size_t pos = line.find(";");
				if(pos != std::string::npos){
					size_t origin,destination;
					std::string arrow, tmpOrigin,tmpDestination;
					std::istringstream iss(line);
					iss >> tmpOrigin >> arrow >>tmpDestination;
					std::cout << tmpOrigin << " " << tmpDestination << std::endl;
					std::cout << "ref acc "<< ref_accession << std::endl;
					std::string tmp1(tmpOrigin.end()-1,tmpOrigin.end());
					std::string name1(tmpOrigin.begin(),tmpOrigin.end()-1);
					std::string tmpAccName1 = ref_accession + ":"+name1;
					std::cout<<tmpAccName1<<std::endl;
					std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(tmpAccName1);
					std::cout << "seq id is "<< findseq1->second << std::endl;
					if(findseq1==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << tmpAccName1 << std::endl;
						exit(1);
					}
					origin = findseq1->second;
					if(tmpDestination.empty()){
						Graph::Vertex v1= Graph::Vertex(origin);
						v1.setName(name1);
						addVertex(v1);
						counter ++;
					}
					else{
						std::string tmp2(tmpDestination.end()-1,tmpDestination.end());
						std::string name2(tmpDestination.begin(),tmpDestination.end()-1);
						std::string tmpAccName2 = ref_accession + ":"+name2;
						std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find(tmpAccName2);
						if(findseq2==longname2seqidx.end()) {
							std::cerr << "Error: unknown sequence in dot File: " << tmpAccName2 << std::endl;
							exit(1);
						}
						destination = findseq2->second;
						//case 1 : RtoF
						if((tmp1 == "-") && (tmp2 == "+")){
							addEdge(1,origin,name1,destination,name2);
							counter ++;
							continue;
						}
						//case 2 : RtoR
						if((tmp1 == "-") && (tmp2 == "-")){
							addEdge(2,origin,name1,destination,name2);
							counter ++;
							continue;
						}
						//case 3 : FtoF
						if((tmp1 == "+") && (tmp2 == "+")){
							addEdge(3,origin,name1,destination,name2);
							counter ++;
							continue;
						}
						//case 4 : FtoR
						if((tmp1 == "+") && (tmp2 == "-")){
							addEdge(4,origin,name1,destination,name2);
							counter ++;
							continue;
						}
					}
				}
			}
		}
	}
	std::cout << "number of aedges is "<< counter << std::endl;
}
void Graph::writeDotFile(std::string file,std::vector<Graph::OrientedVertex> path){
	std::ofstream dotOut(file.c_str());
	dotOut << "digraph adj {\n";
	if(path.size() == 1){
		dotOut<<path[0].getVertex().getName();
		if(path[0].getIsForward() == true)
			dotOut <<"+;\n";
		else
			dotOut << "-;\n";
	}
	else{
		for(uint i = 0 ; i < path.size(); ++i){
			for(uint j = i+1 ; j < path.size() ; ++j){
				dotOut << path[i].getVertex().getName();
				if(path[i].getIsForward() == true)
					dotOut <<"+ -> ";
				else
					dotOut << "- -> ";
				dotOut << path[j].getVertex().getName();
				if(path[j].getIsForward() == true)
					dotOut <<"+;";
				else
					dotOut << "-;";
			}
		}
	}
	dotOut << "}";
	dotOut.close();

}
void Graph::writeDotFile(std::string file, Graph graph){
	std::ofstream dotOut(file.c_str());
	dotOut << "digraph adj {\n";
	for (std::map<size_t,Vertex>::iterator i = graph.vertices.begin(); i != graph.vertices.end(); ++i){
		if(i->second.getFtoF().empty() && i->second.getFtoR().empty() && i->second.getRtoF().empty() && i->second.getRtoR().empty() )
			dotOut << i->second.getName() << "+ ;\n";
		for(uint j = 0 ; j < i->second.getFtoF().size(); ++j){
			dotOut << i->second.getName() << "+ -> " << graph.vertices.find(i->second.getFtoF()[j])->second.getName() << "+ ;\n";
		}
		for(uint j = 0 ; j < i->second.getFtoR().size(); ++j){
			dotOut << i->second.getName() << "+ -> " << graph.vertices.find(i->second.getFtoR()[j])->second.getName() << "- ;\n";
		}
		for(uint j = 0 ; j < i->second.getRtoF().size(); ++j){
			dotOut << i->second.getName() << "- -> " << graph.vertices.find(i->second.getRtoF()[j])->second.getName() << "+ ;\n";
		}
		for(uint j = 0 ; j < i->second.getRtoR().size(); ++j){
			dotOut << i->second.getName() << "- -> " << graph.vertices.find(i->second.getRtoR()[j])->second.getName() << "- ;\n";
		}
	}
	dotOut << "}";
	dotOut.close();
}
void Graph::alignmentOut(std::string & read_name, std::vector<OrientedVertex> path, all_data data,std::ostream & output){//Add the best path for each read to a fasta file
	//remove virtual node from path
	path.erase(path.begin());
	path.erase(path.end());
	std::string seq, read;
	output <<">" <<  read_name << "\n"; 
	read =  convertDnastringToString(read_name) ;
	for (unsigned i = 0; i < read.size(); i += 60) {
	    output << read.substr(i, 60) << "\n";
	}
	for(std::vector<OrientedVertex>::iterator i = path.begin(); i != path.end(); ++i){
			output <<">" << i->getVertex().getName()<< "\n" ;
			output << data.getAlignment(i->getVertex().getIdVertex()).get_al_ref1() << "\n";
	}
}

void Graph::addVertex(Graph::Vertex v){
	std::cout << "adding here! " << vertices.size() <<std::endl;
	vertices.insert(std::make_pair(v.getIdVertex(), v));
}
template<class InputIterator, class T>
  InputIterator findVertex (InputIterator first, InputIterator last, const T& val)
{
  while (first!=last) {
    if (first->getName()==val.getName() && first->getStart() == val.getStart() && first->getEnd() == val.getEnd()) return first;
    ++first;
  }
  return last;
}

void Graph::addEdge(int typeOfEdge, const size_t & id1, const std::string & name1, const size_t & id2, const std::string & name2){
	// create nodes and add edges
	//TODO change this and use Oriented vertex
	std::cout << "id1 " << id1 << " id2 " << id2 <<std::endl;
	std::map<size_t,Graph::Vertex>::iterator it = vertices.find(id1);
	if(it == vertices.end()){ // If vertex does not exist
		Graph::Vertex v1= Graph::Vertex(id1);
		switch(typeOfEdge){
			case 1 :
				v1.setRtoF(id2);
				break;
			case 2 :
				v1.setRtoR(id2);
				break;
			case 3 :
				v1.setFtoF(id2);
				break;
			case 4 :
				v1.setFtoR(id2);
				break;
		}
		v1.setName(name1);
		addVertex(v1);

	}
	else{ // if exist, update
		switch(typeOfEdge){
			case 1 :
				it->second.setRtoF(id2);
				break;
			case 2 :
				it->second.setRtoR(id2);
				break;
			case 3 :
				it->second.setFtoF(id2);
				break;
			case 4 :
				it->second.setFtoR(id2);
				break;
		}
	}
	std::map<size_t,Graph::Vertex>::iterator it2 = vertices.find(id2);
	if(it2 == vertices.end()){// if v2 does not exist, creation
		Graph::Vertex v2 = Graph::Vertex(id2);
		switch(typeOfEdge){
			case 1 :
				v2.setPreviousRtoF(id1);
				break;
			case 2 :
				v2.setPreviousRtoR(id1);
				break;
			case 3 :
				v2.setPreviousFtoF(id1);
				break;
			case 4 :
				v2.setPreviousFtoR(id1);
				break;
		}
		v2.setName(name2);
		addVertex(v2);
	}
	else{ // else update v2
		switch(typeOfEdge){
			case 1 :
				it2->second.setPreviousRtoF(id1);
				break;
			case 2 :
				it2->second.setPreviousRtoR(id1);
				break;
			case 3 :
				it2->second.setPreviousFtoF(id1);
				break;
			case 4 :
				it2->second.setPreviousFtoR(id1);
				break;
		}
	}
}
int Graph::OrientedVertex::typeOfEdge(OrientedVertex ov2){
		if(isForward == true){
			if(ov2.isForward == true)
				return 3;
			else
				return 4;
		}
		else{
			if (ov2.isForward == true)
				return 1;
			else
				return 2;
		}
}
std::ostream& operator<<(std::ostream &os, Graph & g){
	g.printGraph(os);
	return os;
}

void Graph::printGraph(std::ostream &os){
	os << "There is " << vertices.size() << " nodes" << std::endl;
	for (std::map<size_t,Graph::Vertex>::iterator it=vertices.begin(); it!=vertices.end(); ++it)
		(it->second).printVertex(os) ;
}
void Graph::print(){
	std::cout << "Graph contains "<<vertices.size() << " nodes " <<  std::endl;
	for (std::map<size_t,Graph::Vertex>::iterator it=vertices.begin(); it!=vertices.end(); ++it){
		it->second.printVertex();
	}

}
std::ostream& operator<<(std::ostream &os, Graph::Vertex const& v)
{
	v.printVertex(os);
	return os;
}

void Graph::Vertex::printVertex(std::ostream &os)const{
	os <<"\t---------\n";
	os << " vertex : " << idVertex << " pos["<<startOnRead<<"-"<<endOnRead<<"] name "<< name<<"\n";
	os <<"\t";
	if(RtoF.size() != 0)
		os <<"RtoF: ";
	for(unsigned int i=0 ; i < RtoF.size(); ++i){
		 os << RtoF[i] <<" ";
	}
	os <<"\t";
	if(RtoR.size() != 0)
		os <<"RtoR: ";
	for(unsigned int i=0 ; i < RtoR.size(); ++i){
		os << RtoR[i] <<" ";
	}
	os <<"\t";
	if(FtoF.size() != 0)
		os <<"FtoF: ";
	for(unsigned int i=0 ; i < FtoF.size(); ++i){
		os << FtoF[i] <<" ";
	}
	os <<"\t";
	if(FtoR.size() != 0)
		os <<"FtoR: ";
	for(unsigned int i=0 ; i < FtoR.size(); ++i){
		os << FtoR[i] <<" ";
	}
	os << "\n";

	if(previousRtoF.size() != 0)
		os <<"previousRtoF: ";
	for(unsigned int i=0 ; i < previousRtoF.size(); ++i){
		 os << previousRtoF[i] <<" ";
	}
	os <<"\t";
	if(previousRtoR.size() != 0)
		os <<"previousRtoR: ";
	for(unsigned int i=0 ; i < previousRtoR.size(); ++i){
		os << previousRtoR[i] <<" ";
	}
	os <<"\t";
	if(previousFtoF.size() != 0)
		os <<"previousFtoF: ";
	for(unsigned int i=0 ; i < previousFtoF.size(); ++i){
		os << previousFtoF[i] <<" ";
	}
	os <<"\t";
	if(previousFtoR.size() != 0)
		os <<"previousFtoR: ";
	for(unsigned int i=0 ; i < previousFtoR.size(); ++i){
		os << previousFtoR[i] <<" ";
	}
	os << "\n";
	os<<" score ";
	for(std::map<int,double>::const_iterator it = costScore.begin(); it != costScore.end(); ++it){
		os << it->first <<" " << it->second << "\n";
	}
}
void Graph::Vertex::printVertex()const{
	std::cout <<"\t---------\n";
	std::cout << " vertex : " << idVertex << "  name "<< name<<"\n";
	std::cout <<"\t";
	if(RtoF.size() != 0)
		std::cout <<"RtoF: ";
	for(unsigned int i=0 ; i < RtoF.size(); ++i){
		 std::cout << RtoF[i] <<" ";
	}
	std::cout <<"\t";
	if(RtoR.size() != 0)
		std::cout <<"RtoR: ";
	for(unsigned int i=0 ; i < RtoR.size(); ++i){
		std::cout << RtoR[i] <<" ";
	}
	std::cout <<"\t";
	if(FtoF.size() != 0)
		std::cout <<"FtoF: ";
	for(unsigned int i=0 ; i < FtoF.size(); ++i){
		std::cout << FtoF[i] <<" ";
	}
	std::cout <<"\t";
	if(FtoR.size() != 0)
		std::cout <<"FtoR: ";
	for(unsigned int i=0 ; i < FtoR.size(); ++i){
		std::cout << FtoR[i] <<" ";
	}
	std::cout << "\n";

	if(previousRtoF.size() != 0)
		std::cout <<"previousRtoF: ";
	for(unsigned int i=0 ; i < previousRtoF.size(); ++i){
		 std::cout << previousRtoF[i] <<" ";
	}
	std::cout <<"\t";
	if(previousRtoR.size() != 0)
		std::cout <<"previousRtoR: ";
	for(unsigned int i=0 ; i < previousRtoR.size(); ++i){
		std::cout << previousRtoR[i] <<" ";
	}
	std::cout <<"\t";
	if(previousFtoF.size() != 0)
		std::cout <<"previousFtoF: ";
	for(unsigned int i=0 ; i < previousFtoF.size(); ++i){
		std::cout << previousFtoF[i] <<" ";
	}
	std::cout <<"\t";
	if(previousFtoR.size() != 0)
		std::cout <<"previousFtoR: ";
	for(unsigned int i=0 ; i < previousFtoR.size(); ++i){
		std::cout << previousFtoR[i] <<" ";
	}
	std::cout << "\n";
	std::cout <<" score ";
	for(std::map<int,double>::const_iterator it = costScore.begin(); it != costScore.end(); ++it){
		std::cout << it->first <<" " << it->second << "\n";
	}
}

bool operator==( Graph::OrientedVertex ov1, Graph::OrientedVertex ov2){
	return (ov1.getIsForward() == ov2.getIsForward() &&  ov1.getVertex().getIdVertex() == ov2.getVertex().getIdVertex()) ;//FIXME was getName() -> now getIdVertex() what is better ? ov1.getVertex().getName() == ov2.getVertex().getName() &&
}
bool customComp(Graph::OrientedVertex ov1, Graph::OrientedVertex ov2){
	return (ov1.getIsForward() == ov2.getIsForward() && ov1.getVertex().getName() == ov2.getVertex().getName() && ov1.getVertex().getIdVertex() == ov2.getVertex().getIdVertex()) ;
}
void Graph::updateDistance(Graph::OrientedVertex& ov1, Graph::OrientedVertex& ov2, std::map<int,std::vector<Graph::OrientedVertex> >& previous){
	double weight = ov2.getVertex().getCostScore(ov1.getVertex().getIdVertex());
	if(VERBOSE) std::cout << "updateDistance " << ov1.getVertex() <<"v1.getDistance() "<<ov1.getVertex().getDistance()  << " weight "<< weight <<" "<<ov2.getVertex() <<" v2.getDistance()" << ov2.getVertex().getDistance()<< " START " << ov1.getVertex().getStartOnNode()<<std::endl;
	if(ov1.getVertex().getDistance()+ weight < ov2.getVertex().getDistance()){ //TODO <= ?
		ov2.getVertex().setDistance(ov1.getVertex().getDistance()+weight);
		std::vector<Graph::OrientedVertex> tmpVector = previous[ov1.getVertex().getIdVertex()];
		tmpVector.push_back(ov1);
		previous[ov2.getVertex().getIdVertex()] = tmpVector;
	}
}

std::vector<Graph::OrientedVertex> Graph::dijkstra(OrientedVertex orientedStart,OrientedVertex orientedEnd, size_t & read_id ,std::string initPartOfRead,all_data data, use_model model){
	std::cout << "dijkstra gap " << std::endl; //XXX What gap??
	clock_t startTime = clock();
	dnastring seqOfRead= data.getSequence(read_id);
	std::vector<OrientedVertex> notVisited;
	std::vector<OrientedVertex> visited;
	std::map<int,std::vector<OrientedVertex> > previous;

	for( std::map<size_t,Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ){//XXX WHY?!
		it->second.setDistance(std::numeric_limits<double>::infinity());
	}
	//Initialization depending to the start
	orientedStart.getVertex().setDistance(0); //Start of the gap, distance 0
	notVisited.push_back(orientedStart);
	while (!notVisited.empty()){

		std::multimap<double,OrientedVertex>sortedNotVisited;
		for( uint it = 0; it < notVisited.size(); ++it){
			sortedNotVisited.insert(std::pair<double,OrientedVertex>(notVisited[it].getVertex().getDistance(),notVisited[it]));
		}
		OrientedVertex firstNotVisited = sortedNotVisited.begin()->second;
		size_t startOnRead = firstNotVisited.getVertex().getEndOnRead() +1;
		size_t endOnRead = orientedEnd.getVertex().getStartOnRead()-1;
		std::cout << "firstNotVisited " << firstNotVisited.getVertex().getIdVertex() << " " << firstNotVisited.getVertex().getName() << " startOnRead " << startOnRead << " endOnRead " << endOnRead << std::endl;
		std::cout << " on node " << firstNotVisited.getVertex().getStartOnNode()<<std::endl;
		if(firstNotVisited == orientedEnd )
			break;
		else if(startOnRead >= endOnRead){ //find a way to allow some gap
			if(VERBOSE) std::cout << "reach the end of the read, but we don't reach the endNode "<< std::endl;
			visited.push_back(firstNotVisited);
			std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
			std::vector<OrientedVertex> tmp = previous[firstNotVisited.getVertex().getIdVertex()];
			if(VERBOSE) std::cout << " end on the node " <<firstNotVisited.getVertex().getStartOnRead()<< " - " <<firstNotVisited.getVertex().getEndOnRead() << " " << firstNotVisited.getVertex().getName() << std::endl;
			tmp.push_back(firstNotVisited);
			previous[orientedEnd.getVertex().getIdVertex()] = tmp;
			notVisited.erase(itToErase);
			continue;
		}
		else {
			std::vector<OrientedVertex> successors;
			firstNotVisited.getSuccessors(successors, *this); // I need to sort neighbor to find the next one with the smaller cost
			std::string partOfRead = extractPartOfSeq(seqOfRead,startOnRead, endOnRead);
			std::cout<< "part of read size "<< partOfRead.length() << " successors size "<< successors.size() << std::endl;
			for (uint it = 0; it < successors.size(); ++it){
				if(successors[it] == orientedEnd) //we don't look at the next one
					continue;
				else
					this->lookGap(read_id, partOfRead, successors[it],data, startOnRead,model,1,firstNotVisited.getVertex().getIdVertex());
				updateDistance(firstNotVisited,successors[it],previous);
				if(std::find(visited.begin(),visited.end(),successors[it]) == visited.end()){
						notVisited.push_back(successors[it]);
				}
			}
		}
		visited.push_back(firstNotVisited);
		std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
		notVisited.erase(itToErase);
	}
	previous[orientedEnd.getVertex().getIdVertex()].push_back(orientedEnd); //add the last node
	clock_t endTime = clock()- startTime;
	if(VERBOSE) std::cout <<"djikstra gap time " << ((float)endTime)/CLOCKS_PER_SEC<< std::endl;
	if(VERBOSE)std::cout << "path length "<< previous[orientedEnd.getVertex().getIdVertex()].size()<<std::endl;
	return previous[orientedEnd.getVertex().getIdVertex()];
}
std::vector<Graph::OrientedVertex> Graph::dijkstra(OrientedVertex orientedStart,std::string partOfRead,all_data data ,use_model model, std::string option,Graph newGraph, size_t & read_id, size_t & als_size){
	std::cout << "Dijkstra "<<std::endl;
	dnastring seqOfRead= data.getSequence(read_id);
	std::vector<OrientedVertex> notVisited;
	std::vector<OrientedVertex> visited;
	std::map<int,std::vector<OrientedVertex> > previous;
	OrientedVertex startOV = OrientedVertex(newGraph.getVertex(als_size+1),true);
	OrientedVertex endOV = OrientedVertex(newGraph.getVertex(als_size+2),true);
	for( std::map<size_t,Graph::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
	    	it->second.setDistance(std::numeric_limits<double>::infinity());
	}
	//Initialization depending to the start
	orientedStart.getVertex().setDistance(orientedStart.getVertex().getCostScore());
	orientedStart.getVertex().setDistance(0); //Start of the gap, distance 0
	notVisited.push_back(orientedStart);
	while (!notVisited.empty()){

		std::multimap<double,OrientedVertex> orderOrientedVertex;
		std::multimap<double,OrientedVertex>sortedNotVisited;
		for( uint it = 0; it < notVisited.size(); ++it){
			sortedNotVisited.insert(std::pair<double,OrientedVertex>(notVisited[it].getVertex().getDistance(),notVisited[it]));
		}
		OrientedVertex firstNotVisited = sortedNotVisited.begin()->second;
		if(option == "next"){
			std::vector<OrientedVertex> successors;
			firstNotVisited.getSuccessors(successors, *this); // I need to sort neighbor to find the next one with the smaller cost
			size_t startOnRead = firstNotVisited.getVertex().getEndOnRead() +1;
			size_t endOnRead = seqOfRead.length() -1 ; //startOnRead + partOfRead.size()-1;
			if (startOnRead == seqOfRead.length() ){
				visited.push_back(firstNotVisited);
				std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
				std::vector<OrientedVertex> tmp = previous[firstNotVisited.getVertex().getIdVertex()];
				tmp.push_back(firstNotVisited);
				previous[endOV.getVertex().getIdVertex()] = tmp;
				notVisited.erase(itToErase);
				break; //TODO break or continue ?? //XXX ??!!
			}
			partOfRead= extractPartOfSeq(seqOfRead, firstNotVisited.getVertex().getEndOnRead()+1, endOnRead); //XXX why is it replaced by a new value without using the old value??!!
			// find score of node !
			for (uint it = 0; it < successors.size(); ++it){
				this->lookNext(partOfRead, successors[it],data, startOnRead,model,3,firstNotVisited.getVertex().getIdVertex());
				orderOrientedVertex.insert(std::pair<double,OrientedVertex>(successors[it].getVertex().getDistance(),successors[it]));
			}

		}
		else{ //option == previous
			std::vector<OrientedVertex> predecessors;
			firstNotVisited.getPredecessors(predecessors, *this); // I need to sort neighbor to find the next one with the smaller cost
			int endOnRead = firstNotVisited.getVertex().getStartOnRead() -1;
			size_t startOnRead = 0;
			if (endOnRead <= 0){
				if(VERBOSE) std::cout << "reach the start of the read"<< firstNotVisited.getVertex() << std::endl;
				visited.push_back(firstNotVisited);
				std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
				notVisited.erase(itToErase);
				std::vector<OrientedVertex> tmp = previous[firstNotVisited.getVertex().getIdVertex()];
				tmp.push_back(firstNotVisited);
				previous[startOV.getVertex().getIdVertex()] = tmp;
				break; //TODO break or continue ?? //XXX ??!!
			}
			partOfRead= extractPartOfSeq(seqOfRead, startOnRead, endOnRead);
			// find score of node !
			for (uint it = 0; it < predecessors.size(); ++it){
				this->lookPrevious(partOfRead, predecessors[it],data, endOnRead,model,2,firstNotVisited.getVertex().getIdVertex(), read_id);
				orderOrientedVertex.insert(std::pair<double,OrientedVertex>(predecessors[it].getVertex().getDistance(),predecessors[it]));
			}
		}
		for(std::multimap<double,OrientedVertex>::iterator itOfSuccessor =orderOrientedVertex.begin(); itOfSuccessor != orderOrientedVertex.end(); ++itOfSuccessor ){
			if(option == "next")
				updateDistance(firstNotVisited,itOfSuccessor->second,previous);
			else
				updateDistance(firstNotVisited,itOfSuccessor->second,previous);
			// add node to the queue
			if(std::find(visited.begin(),visited.end(),itOfSuccessor->second) == visited.end() ){
				notVisited.push_back(itOfSuccessor->second);
			}
			//Before erase neighbor, update vertices
			vertices.find(itOfSuccessor->second.getVertex().getIdVertex())->second.setDistance(itOfSuccessor->second.getVertex().getDistance());

		}
		visited.push_back(firstNotVisited);
		std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
		if(itToErase == notVisited.end()){
			std::cerr << "error cannot find vertex in not visited" << std::endl;
			exit(1);
		}
		notVisited.erase(itToErase);
	}
	if(option == "next"){
		previous[endOV.getVertex().getIdVertex()].push_back(endOV);
		return previous[endOV.getVertex().getIdVertex()];
	}
	else{ // if previous
		previous[startOV.getVertex().getIdVertex()].push_back(startOV);
		std::vector<OrientedVertex> tmp2 =  previous[startOV.getVertex().getIdVertex()];
		std::reverse(tmp2.begin(),tmp2.end());
		return tmp2;
	}
}

std::vector<Graph::OrientedVertex> Graph::dijkstra(OrientedVertex orientedStart,OrientedVertex orientedEnd){
	std::cout << "Final dijkstra "<< std::endl;
	int loop = 0 ;
	std::vector<OrientedVertex> notVisited;
	std::vector<OrientedVertex> visited;
	std::map<int,std::vector<OrientedVertex> > previous;
  	orientedStart.getVertex().setDistance(orientedStart.getVertex().getCostScore());
	notVisited.push_back(orientedStart);
	orientedStart.getVertex().setDistance(0);//notVisited.begin()->getVertex().getCostScore()
	while (!notVisited.empty()){
		std::multimap<double,OrientedVertex>sortedNotVisited;
		for( uint it = 0; it < notVisited.size(); ++it){ //TODO just replace notVisited by multimap
			sortedNotVisited.insert(std::pair<double,OrientedVertex>(notVisited[it].getVertex().getDistance(),notVisited[it]));
		}
		OrientedVertex firstNotVisited = sortedNotVisited.begin()->second;
		if(firstNotVisited == orientedEnd){
			break;
		}
		else {
			std::vector<OrientedVertex> successors;
			firstNotVisited.getSuccessors(successors, *this);
			std::multimap<double,OrientedVertex> testSuccessor;
			for( uint it = 0; it < successors.size(); ++it){
				testSuccessor.insert(std::pair<double,OrientedVertex>(successors[it].getVertex().getDistance(),successors[it]));
			}
			for(std::multimap<double,OrientedVertex>::iterator itOfSuccessor =testSuccessor.begin(); itOfSuccessor != testSuccessor.end(); ++itOfSuccessor ){
				updateDistance(firstNotVisited,itOfSuccessor->second,previous);
				if(std::find_if(visited.begin(), visited.end(), special_compare(itOfSuccessor->second)) == visited.end() && std::find_if(notVisited.begin(), notVisited.end(), special_compare(itOfSuccessor->second)) == notVisited.end() ){
					notVisited.push_back(itOfSuccessor->second);//TODO if not in notVisited !!
				}
				std::cout << std::endl;
				//Before erase neighbor, update vertices
				vertices.find(itOfSuccessor->second.getVertex().getIdVertex())->second.setDistance(itOfSuccessor->second.getVertex().getDistance());
			}
		}
		visited.push_back(firstNotVisited);
		std::vector<OrientedVertex>::iterator itToErase = std::find(notVisited.begin(),notVisited.end(),firstNotVisited);
		if(itToErase == notVisited.end()){
			std::cerr << "error cannot find remove vertex in notVisited" << std::endl;
			exit(1);
		}
		notVisited.erase(itToErase);
		++loop;
	}
		previous[orientedEnd.getVertex().getIdVertex()].push_back(orientedEnd);
		printFinalPath(previous[orientedEnd.getVertex().getIdVertex()]);
//		for(std::map<int,vector<OrientedVertex> >::iterator it = previous.begin(); it != previous.end();++it){
//			std::cout << "-------------"<< it->first << std::endl;
//			printFinalPath(it->second);
//		}
	return previous[orientedEnd.getVertex().getIdVertex()];
}

double Graph::countCostScore(use_model m,pw_alignment p){
	clock_t timeCost = clock();
	m.cost_function(p);
	timeCost = clock() - timeCost;
	if(VERBOSE)std::cout <<"countCostScore time " <<((float)timeCost)/CLOCKS_PER_SEC<< std::endl;
	return p.get_modify1();
}

void Graph::lookPrevious(std::string partOfRead,Graph::OrientedVertex& previousNode, all_data data, size_t endOnRead, use_model model, int optionNeedleman,int idOrigin, size_t & read_id){
	std::string partOfPreviousNode, newPartOfRead, restPartOfRead ;
	int startToExtract = 0;
	int endToExtract = partOfRead.size()-1;
	size_t lengthToLook = partOfRead.size();
	size_t lengthPreviousNode = data.get_seq_size(previousNode.getVertex().getIdVertex());

	if(lengthToLook > lengthPreviousNode){
		lengthToLook = lengthPreviousNode; //reduce the length
		startToExtract = endToExtract - (lengthToLook -1) ;
		partOfRead = extractPartOfSeq(partOfRead, startToExtract,endToExtract );
	}
	int startOnPreviousNode =  lengthPreviousNode - lengthToLook;
	size_t incl_end1 = startOnPreviousNode + lengthToLook -1;
	partOfPreviousNode = extractPartOfSeq(data.getSequence(previousNode.getVertex().getIdVertex()), startOnPreviousNode,incl_end1);
	std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfPreviousNode,optionNeedleman);
	size_t idx1 = previousNode.getVertex().getIdVertex();
	size_t idx2 = read_id;//the read is the last sequence in fasta file
	std::cout << "here! "<<std::endl;
	pw_alignment al(resultNeedleman.first,resultNeedleman.second, startOnPreviousNode,startToExtract, incl_end1,endToExtract, idx1, idx2);
	data.alignment_fits_ref(al);
	previousNode.getVertex().setStartOnRead(startToExtract);
	previousNode.getVertex().setEndOnRead(endToExtract);
	previousNode.getVertex().setStartOnNode(startOnPreviousNode);
	previousNode.getVertex().setEndOnNode(incl_end1+1);
	double score = countCostScore(model,al);
	previousNode.getVertex().setCostScore(idOrigin,score);
}

void Graph::lookNext(std::string partOfRead,Graph::OrientedVertex& onode, all_data data,size_t startOnRead,use_model model,int optionNeedleman,int idOrigin){
	if(VERBOSE) std::cout << "lookNext"<<std::endl;
	clock_t startTime = clock();
	size_t lengthToLook = partOfRead.size();
	std::string newPartOfRead, partOfNextNode;
	int startToExtract, endToExtract;
	size_t lengthNextNode = data.get_seq_size(onode.getVertex().getIdVertex());
	startToExtract = 0;
	if( lengthToLook > lengthNextNode ){
		lengthToLook = lengthNextNode;
		endToExtract = lengthToLook-1;
		partOfRead = extractPartOfSeq(partOfRead, startToExtract,endToExtract);
	}
	else
	endToExtract = lengthToLook -1;
	partOfNextNode = extractPartOfSeq(data.getSequence(onode.getVertex().getIdVertex()), startToExtract,endToExtract);
	std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfNextNode,optionNeedleman);
	size_t idx1 = onode.getVertex().getIdVertex();
	size_t idx2 = data.numSequences()-1; // the read is the last one in the fasta file
	size_t endOnRead =  startOnRead + lengthToLook-1;
	pw_alignment al(resultNeedleman.first,resultNeedleman.second, startToExtract,startOnRead,endToExtract,endOnRead, idx1, idx2);

	// Does the alignment is ok ?
	data.alignment_fits_ref(al);
	double score = countCostScore(model,al);
	onode.getVertex().setStartOnRead(startOnRead);
	onode.getVertex().setEndOnRead(endOnRead);
	onode.getVertex().setStartOnNode(startToExtract);
	onode.getVertex().setEndOnNode(endToExtract);
	onode.getVertex().setCostScore(idOrigin,score);
	clock_t endTime = clock()-startTime;
	if(VERBOSE)std::cout <<"lookNext time " << ((float)endTime)/CLOCKS_PER_SEC<< std::endl;
}

void Graph::lookGap(size_t & read_id , std::string partOfRead,Graph::OrientedVertex& onode, all_data data,size_t startOnRead,use_model model,int optionNeedleman,int idOrigin){
	//when we look at a Gap we always want the complete node for the alignment !
	std::cout << "look Gap " << std::endl;
	clock_t startTime = clock();
	size_t lengthToLook = partOfRead.size();
	std::string newPartOfRead, seqOfNextNode;
	size_t startToExtract, endToExtract;
	size_t lengthNextNode = data.get_seq_size(onode.getVertex().getIdVertex());
	startToExtract = 0;
	endToExtract = lengthToLook -1;
	if( lengthToLook > lengthNextNode ){
		lengthToLook = lengthNextNode;
		endToExtract = lengthToLook-1;
		partOfRead = extractPartOfSeq(partOfRead, startToExtract,endToExtract);
	}
	seqOfNextNode = convertDnastringToString(data.getSequence(onode.getVertex().getIdVertex()));
	std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,seqOfNextNode,optionNeedleman);
	size_t idx1 = onode.getVertex().getIdVertex();
	size_t idx2 = read_id; 
	size_t endOnRead =  startOnRead + lengthToLook-1;
	pw_alignment al(resultNeedleman.first,resultNeedleman.second, startToExtract,startOnRead,lengthNextNode-1,endOnRead, idx1, idx2);
	if(VERBOSE) std::cout << "look gap idx1 " << idx1 << " idx2 " << idx2 << " startToExtract " << startToExtract << " startOnRead " << startOnRead;
	if(VERBOSE) std::cout << " endToExtract " << lengthNextNode-1 << " endOnRead "<< endOnRead<< " " << onode.getVertex().getName() << std::endl;

	clock_t timeFits = clock();
	data.alignment_fits_ref(al);
	timeFits = clock() - timeFits;
	if(VERBOSE) std::cout <<"alignment_fits_ref time " <<  ((float)timeFits)/CLOCKS_PER_SEC<< std::endl;
	if(VERBOSE) std::cout << " start On node "<< startToExtract<< std::endl;
	double score = countCostScore(model,al);
	onode.getVertex().setStartOnRead(startOnRead);
	onode.getVertex().setEndOnRead(endOnRead);
	onode.getVertex().setStartOnNode(startToExtract);
	onode.getVertex().setEndOnNode(endToExtract);
	onode.getVertex().setCostScore(idOrigin,score);

	startTime = clock() - startTime;
	if(VERBOSE) std::cout <<"lookGap time " << ((float)startTime)/CLOCKS_PER_SEC << std::endl;
}
void Graph::set_first_and_last_node(size_t alsnumber, size_t start, size_t end){
	//Creation of virtual start and end node
		Vertex v1= Vertex(alsnumber+1);
		v1.setName("startNode");
		v1.setCostScore(alsnumber+1,0);
		v1.setEndOnRead(start);//End of a vertix is excluded.
		this->addVertex(v1);
		Vertex v2= Vertex(alsnumber+2);
		v2.setName("endNode");
		v2.setCostScore(alsnumber+2,0);
		v2.setStartOnRead(end+1);
		this->addVertex(v2);
}
/*bool Graph::checkVertex(std::string name, int start, int end,Graph g){ ///TODO Do it better
	for (std::map<int,Vertex>::const_iterator it = g.getVertices().begin(); it != g.getVertices().end(); ++it){
		std::string tmp = name;
		if((it->second.getName() == name || it->second.getName() == tmp) && it->second.getStartOnRead() == start && it->second.getEndOnRead() == end){
			return true;
		}
	}
	return false;
}
*/
void Graph::updateNewGraph(std::vector<Graph::OrientedVertex> path, all_data& data, size_t & read_id, std::string gap, use_model model, Graph& newGraph, int option, std::string whereToLook, std::vector<pw_alignment> & all_als){
	std::cout << "update New graph "<< std::endl;
	clock_t time = clock();
	size_t lengthOfGap = gap.size();
	dnastring read = data.getSequence(read_id); 
	OrientedVertex startOV = *(path.begin());
	path.erase(path.begin());
	OrientedVertex endOV = path.back();
	path.erase(path.end());
	std::string seqToAlign;
	std::string newName;
	if(!path.empty()){
	for(std::vector<OrientedVertex>::iterator it = path.begin(); it != path.end(); ++it){
		size_t idOV = it->getVertex().getIdVertex();

		 if(lengthOfGap > data.get_seq_size(idOV)){
			lengthOfGap = lengthOfGap - data.get_seq_size(idOV);
			seqToAlign += extractPartOfSeq(it->getDna(data, idOV),0, data.get_seq_size(idOV)-1);
                }
		else{
			if(whereToLook == "end"){
				seqToAlign += extractPartOfSeq(it->getDna(data, idOV),0, lengthOfGap-1);
			}
			else{ //gap
				seqToAlign += extractPartOfSeq(it->getDna(data, idOV),0, data.get_seq_size(idOV)-1);
			}
		}
		newName += it->getVertex().getName();
	}

	std::pair<std::string,std::string>resultNeedleman =runNeedleman(gap,seqToAlign,option);//TODO option for needleman !!! //XXX What does it do?
	size_t idx1 = 0; //generate new id ? //XXX Why??!!!
	size_t idx2 = read_id;
	std::cout << "read id "<< read_id <<std::endl;
	size_t incl_start2 = startOV.getVertex().getEndOnRead()+1;
	size_t excl_end2 =  endOV.getVertex().getStartOnRead();
	size_t incl_start, excl_end;
	double score;
	if(whereToLook == "end"){
		incl_start = 0 ;
		excl_end = gap.size()-1;
		pw_alignment al(resultNeedleman.first,resultNeedleman.second, incl_start,incl_start2, excl_end,excl_end2, idx1, idx2);
		model.cost_function(al);
		score = al.get_modify1();
		if(data.alignmentAlreadyExist(al) == false){
			std::cout << "alignemnt does not exist!"<<std::endl;
			data.add_pw_alignment(al);
		}

	}
	else{
		incl_start = 0;
		excl_end = seqToAlign.size()-1;
		pw_alignment al(resultNeedleman.first,resultNeedleman.second, incl_start,incl_start2, excl_end,excl_end2, idx1, idx2);
		score = countCostScore(model,al);
		model.cost_function(al);
		score = al.get_modify1();
		if(data.alignmentAlreadyExist(al) == false){
			data.add_pw_alignment(al);
			std::cout << " alignment does not exist"<< std::endl;
		}
	}
	//add link between : startPath -> newNode -> EndPAth
	size_t newIdStart, newIdEnd;
	std::cout << endOV.getVertex().getName() << " " << startOV.getVertex().getName() <<" " <<data.get_seq_name(idx1)<< std::endl;
	if(startOV.getVertex().getIdVertex() == all_als.size()+1){ // I cannot find virtual node in alignments 
		newIdStart = startOV.getVertex().getIdVertex();
	//	newIdEnd = data.findIdAlignment(endOV.getVertex().getName(),endOV.getVertex().getStartOnRead(), endOV.getVertex().getEndOnRead(),endOV.getVertex().getStartOnNode(),endOV.getVertex().getEndOnNode());//TODO
		newIdEnd = get_distance(all_als,endOV.getVertex().getIdVertex(),read_id,endOV.getVertex().getStartOnRead(), endOV.getVertex().getEndOnRead(),endOV.getVertex().getStartOnNode(),endOV.getVertex().getEndOnNode());
	}
	else if ( endOV.getVertex().getIdVertex() == all_als.size()+2){
	//	newIdStart = data.findIdAlignment(startOV.getVertex().getName(),startOV.getVertex().getStartOnRead(),startOV.getVertex().getEndOnRead(),startOV.getVertex().getStartOnNode(),startOV.getVertex().getEndOnNode());//TODO
		newIdStart = get_distance(all_als,startOV.getVertex().getIdVertex(),read_id,startOV.getVertex().getStartOnRead(),startOV.getVertex().getEndOnRead(),startOV.getVertex().getStartOnNode(),startOV.getVertex().getEndOnNode());
		newIdEnd = endOV.getVertex().getIdVertex();
	}
	else{

	//	newIdStart = data.findIdAlignment(startOV.getVertex().getName(),startOV.getVertex().getStartOnRead(),startOV.getVertex().getEndOnRead(),startOV.getVertex().getStartOnNode(),startOV.getVertex().getEndOnNode());//TODO
		newIdStart = get_distance(all_als,startOV.getVertex().getIdVertex(),read_id,startOV.getVertex().getStartOnRead(),startOV.getVertex().getEndOnRead(),startOV.getVertex().getStartOnNode(),startOV.getVertex().getEndOnNode());
		std::cout << endOV.getVertex().getName() << " " <<idx1<< " " << std::endl;
	//	newIdEnd = data.findIdAlignment(endOV.getVertex().getName(),endOV.getVertex().getStartOnRead(), endOV.getVertex().getEndOnRead(),endOV.getVertex().getStartOnNode(),endOV.getVertex().getEndOnNode());//TODO
		newIdEnd = get_distance(all_als,endOV.getVertex().getIdVertex(),read_id,endOV.getVertex().getStartOnNode(), endOV.getVertex().getEndOnNode(),endOV.getVertex().getStartOnRead(), endOV.getVertex().getEndOnRead());
	}
	size_t newIdMiddle = get_distance(all_als, idx1 ,read_id, incl_start, excl_end, incl_start2, excl_end2); //XXX Why did she set idx1 == 0 ???!!!

	int type;
	if(startOV.getIsForward() == true)
		type = 3;
	else
		type = 1;
	newGraph.addEdge(type,newIdStart,startOV.getVertex().getName(),newIdMiddle,newName);
	if(endOV.getIsForward() == true)
		type = 3;
	else
		type = 4;
	newGraph.addEdge(type,newIdMiddle,newName,newIdEnd,endOV.getVertex().getName());
	newGraph.setScore(newIdMiddle,score,newIdStart);
	newGraph.setStartOnRead(newIdMiddle,incl_start2);
	newGraph.setStartOnNode(newIdMiddle,incl_start);
	newGraph.setEndOnRead(newIdMiddle,excl_end2);
	newGraph.setEndOnNode(newIdMiddle,excl_end);
	}
	time = clock() - time;
	std::cout <<"updateNewGraph time " <<((float)time)/CLOCKS_PER_SEC<< std::endl;
}
void Graph::parseData(all_data & data, use_model & model,Graph& newGraph,std::vector<pw_alignment> & als, std::vector<pw_alignment> & all_als, dnastring seqOfRead){
	size_t maxNumberOfBases = 400; // use a distance that change regarding read length : take 1/4 of read length ? //XXX I really don't get the rationale behind this number!!!
	std::cout << "size of vectorAl is " <<als.size()<<std::endl;
	for(size_t itP =0; itP < als.size(); ++itP){//set of als that have less than 400b gap
		std::cout << "itp "<<itP<<std::endl;
		pw_alignment p = als.at(itP);
		std::map<size_t,Vertex>::iterator mapVertexItP = vertices.find(als.at(itP).getreference1()); //Walking over the reference graph. (From the dot file)
		if(mapVertexItP == vertices.end()){
			std::cerr << "error when find vertex " <<als.at(itP).getreference1() << std::endl;
			exit(1);
		}
		//update position on node base on alignment
		size_t left1,right1,left2,right2;
		als.at(itP).get_lr1(left1,right1);//on reference
		als.at(itP).get_lr2(left2,right2);//on read
		//first add all pw_alignment to the graph?
		size_t ref1 = p.getreference1();
		size_t ref2 = p.getreference2();
		size_t node_id = get_distance(all_als, ref1 , ref2 ,left1, right1,left2,right2); //check for all the alignments before the current one and returns its position in alignments vector as its node id.

		Vertex v1= Vertex(node_id);
		v1.setEndOnRead(right2+1);
		v1.setStartOnRead(left2);
		v1.setStartOnNode(left1);
		v1.setEndOnNode(right1+1);
		v1.setName(mapVertexItP->second.getName());
		v1.setCostScore(v1.getIdVertex(),als.at(itP).get_modify1());
		newGraph.addVertex(v1);

		OrientedVertex OV1 = OrientedVertex(mapVertexItP->second,als.at(itP));
		OV1.getVertex().setCostScore(OV1.getVertex().getIdVertex(),als.at(itP).get_modify1());
		OV1.getVertex().setStartOnRead(left2);
		OV1.getVertex().setEndOnRead(right2+1);
		OV1.getVertex().setStartOnNode(left1);
		OV1.getVertex().setEndOnNode(right1+1);
		// Link the current node to the startNode
		if(left2 == newGraph.getVertex(all_als.size()+1).getEndOnRead()){
			int type;
			if(OV1.getIsForward() == true)
				type = 3; 
			else
				type = 4;
			size_t start_node = all_als.size()+1;
			newGraph.addEdge(type, start_node,"startNode",node_id,mapVertexItP->second.getName());
		}
		// Link the current node to the endNode
		if(right2 == newGraph.getVertex(all_als.size()+2).getStartOnRead()-1){
			int type;
			if(OV1.getIsForward() == true)
				type = 3;
			else
				type = 1;
			size_t end_node = all_als.size()+2;
			newGraph.addEdge(type,node_id,mapVertexItP->second.getName(),end_node,"endNode");
		}
		size_t als_size = all_als.size();
		//Verify it in a way that the 5' end of the read fits to the graph 
		if(left2 < maxNumberOfBases && left2 !=0){
			std::string partOfPreviousNode;
			std::cout << "Look at the start of the read from position "<< left1 << " node "<<als.at(itP).getreference1() <<std::endl; //call needleman option 2 //XXX what is that for?
			std::string partOfRead = extractPartOfSeq(seqOfRead, 0, left2-1);
			std::vector<OrientedVertex> path;
			path = this->dijkstra(OV1,partOfRead,data,model,"previous",newGraph,ref2,als_size); //XXX what previous is for?
			if(path.size() > 1) //Artifact, I save the start node even if there is no path. then it is always bigger than one
				newGraph.updateNewGraph(path,data,ref2,partOfRead,model,newGraph,2,"end",all_als);
				std::cout<< "here!!!"<<std::endl;
		}
		//Verify it in a way that the 3' end of the read fits to the graph
		size_t tmp ;
		tmp = seqOfRead.length() - maxNumberOfBases;
		if (right2 > tmp && right2 != seqOfRead.length() - 1 ){ //&& vectorAl[itP].getend2() != seqOfRead.length() - 1 //XXX don't see why?
			std::string partOfRead, partOfNextNode;
			std::cout << "Look at the end of the read from position "<< right1 << " node "<<als.at(itP).getreference1() <<std::endl; 
			partOfRead = extractPartOfSeq(seqOfRead, right2+1, seqOfRead.length()-1);
			std::string endSeq = extractPartOfSeq(seqOfRead,right2+1, seqOfRead.length()-1);
			std::vector<OrientedVertex> path;
			path = this->dijkstra(OV1,endSeq,data,model,"next",newGraph,ref2,als_size);
			if(path.size() > 1)
				updateNewGraph(path,data,ref2,endSeq,model,newGraph,3,"end",all_als);
		}
  		for(size_t itPOther = itP+1 ; itPOther < als.size(); ++itPOther){
			//loop start at itP+1 because we dont need to look at node locate before and at the node[itP] with is the same
			std::cout << als.at(itPOther).getreference1()<<std::endl;
			std::map<size_t,Graph::Vertex>::iterator mapVertexItPOther = vertices.find(als.at(itPOther).getreference1());
			if(mapVertexItPOther == vertices.end()){
				std::cerr << "error cannot find the vertex! " << std::endl;
				als.at(itPOther).print();
				exit(1);
			}
			size_t left_1,right_1,left_2,right_2;
			als.at(itPOther).get_lr1(left_1,right_1);
			als.at(itPOther).get_lr2(left_2,right_2);
			//update position on node base on alignment
			vertices.find(als.at(itPOther).getreference1())->second.setStartOnRead(left_2);
			vertices.find(als.at(itPOther).getreference1())->second.setEndOnRead(right_2+1);
			vertices.find(als.at(itPOther).getreference1())->second.setStartOnNode(left_1);
			vertices.find(als.at(itPOther).getreference1())->second.setEndOnNode(right_1+1);
			//vertices.find(vectorAl[itPOther].getreference1())->second.set(rightNode2);
			std::map<size_t,Graph::Vertex>::iterator itForOV2 = vertices.find(als.at(itPOther).getreference1());
			if( itForOV2 == vertices.end()){
				std::cerr << "cannot find the vertex in parseData "<< als.at(itPOther).getreference1()<<std::endl;
				exit(1);
			}
			OrientedVertex OV2 = OrientedVertex(itForOV2->second, als.at(itPOther));
			size_t ref1 = als.at(itPOther).getreference1();
			size_t ref2 = als.at(itPOther).getreference2();
			size_t idOnNewGraph = get_distance(all_als, ref1,ref2, left_1, right_1, left_2, right_2);
			Graph::Vertex v1= Graph::Vertex(idOnNewGraph);
			v1.setEndOnRead(right_2+1);
			v1.setStartOnRead(left_2);
			v1.setStartOnNode(left_1);
			v1.setEndOnNode(right_1+1);
			v1.setCostScore(v1.getIdVertex(), als.at(itPOther).get_modify1());
			v1.setName(mapVertexItPOther->second.getName());
			newGraph.addVertex(v1);

			if((left2 == left_2 && right2 == right_2) || (right2 >= left_2)) // right >= leftNode2 == overlap
				continue;
			if(left_2 - right2 > maxNumberOfBases)
				break; // if we reach the max distance between two node, we stop looking at itP and go to the next one

			// Case 1 : No gap. when the next alignment has overlap with the current alignment or starts exatly where the previous one ends.
			if( right2 +1  >= left_2){

				std::vector<OrientedVertex> successors;
				OV1.getSuccessors(successors,*this);
				if( std::find(successors.begin(),successors.end(), OV2)!= successors.end()){//XXX???
					int type = OV1.typeOfEdge(OV2);
					newGraph.addEdge(type, node_id, OV1.getVertex().getName(), idOnNewGraph, OV2.getVertex().getName());
					newGraph.setStartOnRead(node_id,left2);//XXX Why did she add them twice?
					newGraph.setEndOnRead(node_id,right2+1);
					newGraph.setStartOnNode(node_id,left1); 
					newGraph.setEndOnNode(node_id,right1+1); 
					newGraph.setStartOnRead(idOnNewGraph,left_2);
					newGraph.setEndOnRead(idOnNewGraph,right_2+1);
					newGraph.setStartOnNode(idOnNewGraph,left_1);
					newGraph.setEndOnNode(idOnNewGraph,right_1+1);
					newGraph.setScore(node_id,als.at(itP).get_modify1(),node_id);
					newGraph.setScore(idOnNewGraph,als.at(itPOther).get_modify1(),node_id);
				}
			}
			//Case 2 : Gap
			else if( ((right2+1 < left_2) && (left_2 - right2 < maxNumberOfBases))){
				// take the gap sequence of the read, take the node after vectorAl[itP] that is different from vectorAl[itPOther]
				std::cout<< "read length "<< seqOfRead.length() <<" gap " << right2 << " - " << left_2 << std::endl;
				std::string gapSeq = extractPartOfSeq(seqOfRead,right2+1, left_2-1);
				std::vector<OrientedVertex> path;
				path = this->dijkstra(OV1,OV2,ref2,gapSeq,data, model);
				std::cout << "look gap"<<  path[0].getVertex().getStartOnNode()<< std::endl;
				if(path.size() > 1)
					updateNewGraph(path,data,ref2,gapSeq,model,newGraph,1,"gap", all_als);
			}
		}
	}

}

void Graph::prepare_read(Graph& newgraph, std::vector<pw_alignment> & als, const dnastring & current_read, use_model & m , all_data & data, std::ostream & output){
	//Remove part of read if gaps are too important, will find more than one path per read.
	std::sort(als.begin(),als.end(),sort_pw_alignment());
	std::cout << "als size " << als.size()<<std::endl;
	size_t maxDistance = 400;//XXX Why 400????!!!!
	//size_t posLastGap = 0 ;
	int nbRead = 1;
	size_t maxRightPos = 0;
	size_t firstLeft = 0;
	size_t firstRight;
	std::vector<pw_alignment> tmpVector;
	const pw_alignment & pal = als.at(0);//first alignment on the read
	std::string read_name = data.get_seq_name(pal.getreference2());//Read name
	pal.get_lr2(firstLeft,firstRight); 
	if(firstLeft  > maxDistance ){// A big part from the beginning of the read didn't map to any part of the reference //XXX what should be done in this case?
		std::cout << "big gap at the beginning"<< firstLeft << std::endl;
	}
	size_t begin = als.size()+1;
	size_t end = als.size()+2;
	for(size_t i =0; i < als.size()-1; i++){
		const pw_alignment & p = als.at(i);
		tmpVector.push_back(p);
		size_t left, right, left2,right2;
		p.get_lr2(left,right);
		const pw_alignment al = als.at(i+1);//XXX Wrong! Doesnt work for i = size-1 unless adding a break at i = size-2!!!
		al.get_lr2(left2,right2);
		std::cout << "left2 " << left2 <<std::endl;
		maxRightPos = std::max(maxRightPos,right);
		std::cout << "max " << maxRightPos << std::endl;
		if(maxRightPos + maxDistance < left2){//When two aligned regions on the the read are so far from each other 
			newgraph.print(); 
			std::cout << "big gap in the data from " << maxRightPos << " " << left2 << std::endl;

			std::cout << "first left "<< firstLeft <<std::endl;
			newgraph.set_first_and_last_node(als.size(),firstLeft,maxRightPos);
			parseData(data, m, newgraph,tmpVector, als ,current_read);//update the newgraph.
			OrientedVertex OV1 = OrientedVertex(newgraph.getVertex(begin),true);
			OrientedVertex OV2 = OrientedVertex(newgraph.getVertex(end),true);
			
			std::vector<OrientedVertex> path = newgraph.dijkstra(OV1,OV2);

			if(path.size()>1){
				alignmentOut(read_name, path,data,output);
			}

			firstLeft = left2;
			++nbRead;
			tmpVector.clear();
			newgraph = Graph();
		}

		if(i == als.size()-2){//XXX This one seems wrong too. Why the last alignment were ignored?
			std::cout << "i == als.size -2"<<std::endl;
			if(maxRightPos + maxDistance < current_read.length() ){//XXX ???
				std::cout << "A big gap at the end of the read" << maxRightPos<< std::endl;
				newgraph.set_first_and_last_node(als.size(),firstLeft,maxRightPos);
				parseData(data,m, newgraph,tmpVector, als,current_read);
				OrientedVertex OV1 = OrientedVertex(newgraph.getVertex(begin),true);
				OrientedVertex OV2 = OrientedVertex(newgraph.getVertex(end),true);
				std::vector<OrientedVertex> path = newgraph.dijkstra(OV1,OV2);
				if(path.size()>1){
					alignmentOut(read_name,path,data,output);
				}
				++nbRead;
			}
		}

	}
	if(nbRead == 1){
		newgraph.set_first_and_last_node(als.size(),0,current_read.length()-1);
		parseData(data,m,newgraph,als ,als ,current_read);
		OrientedVertex OV1 = OrientedVertex(newgraph.getVertex(begin),true);
		OrientedVertex OV2 = OrientedVertex(newgraph.getVertex(end),true);
		std::vector<OrientedVertex> path = newgraph.dijkstra(OV1,OV2);
		if(path.size()>1){
			alignmentOut(read_name,path,data,output);
		}
	}
	std::cout << "new nb of read " << nbRead -1<< std::endl;
}
std::string Graph::extractPartOfSeq(dnastring seq, int start, int end){
	std::cout << "read length is "<< seq.length() << std::endl;
	if( start > end){
		std::cerr << "Problem to extract part of seq start " << start << " end "<< end<< std::endl;
		exit(1);
	}
	std::string subSeq = "";
	for(int i= start; i <= end ; ++i){
	//	std::cout<< seq.at(i)<< " ";
		subSeq += seq.at(i);
	}
//	std::cout << ""<<std::endl;
	return subSeq;
}
std::string Graph::convertDnastringToString(dnastring seq){
	std::string newSeq;
	return newSeq = extractPartOfSeq(seq, 0,seq.length()-1);
}

void Graph::printFinalPath(std::vector<OrientedVertex> path){
	for(std::vector<OrientedVertex>::iterator it = path.begin(); it != path.end() ; ++it){
		std::cout << it->getVertex().getIdVertex()<<it->getIsForward() << "\t"<<it->getVertex().getName() <<"\t"<<it->getVertex().getStartOnRead()<<"-"<<it->getVertex().getEndOnRead()<<std::endl;
	}
}

void Graph::addVirtualEdges(use_model model,dnastring read){
	
	for(std::map<size_t,Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it){
		//check if node is not the first or last node
		if(it->second.getIdVertex() == -1 || it->second.getIdVertex() == -2){
			continue;
		}
		for(std::map<size_t,Vertex>::iterator it2 = vertices.begin(); it2 != vertices.end(); ++it2){
			if(it2->second.getIdVertex() == -1 || it2->second.getIdVertex() == -2)
						continue;
			if(it->second.getEndOnRead()+1 < it2->second.getStartOnRead()){
				addEdge(3,it->second.getIdVertex(),it->second.getName(),it2->second.getIdVertex(),it2->second.getName());
				size_t lengthSeq = it2->second.getStartOnRead() - it->second.getEndOnRead() -1;
				std::string seqN = std::string(lengthSeq, 'N');
				std::string seqRead = extractPartOfSeq(read,it->second.getEndOnRead()+1,it2->second.getStartOnRead()-1);
				pw_alignment al(seqN,seqRead,0,0,lengthSeq,lengthSeq,0,1);
				model.cost_function(al);
				double score = al.get_modify1();
				setScore(it2->second.getIdVertex(),score,it->second.getIdVertex()); //TODO score depending of the sequence
			}
		}
	}
}

const size_t Graph::get_distance(std::vector<pw_alignment> & all_als, const size_t & ref1, const size_t & ref2, const size_t & left1, const size_t & right1, const size_t & left2, const size_t &right2)const{
	for(size_t i = 0; i < all_als.size(); i++){
		pw_alignment p = all_als.at(i);
		size_t l1,r1,l2,r2;
		p.get_lr1(l1,r1);
		p.get_lr2(l2,r2);
		if(p.getreference1() == ref1 && p.getreference2() == ref2 && l1 == left1 && l2 == left2 && r1 ==right1 && r2 == right2){
			return i;
		}
	
	}



}
/*
bool sortOrientedVertex::operator()(Graph::OrientedVertex &a, Graph::OrientedVertex &b) const {
	return a.getVertex().getDistance() < b.getVertex().getDistance();
}
*/
int testSamFile(std::string file){
	std::ifstream inFile(file.c_str());
		int cmp = 0;
		if(!inFile){
				std::cerr << "Error : Cannot open " << file.c_str() << std::endl;
				exit(1);
			}
		else{
			std::string line;
			while(getline(inFile,line)){
				if(line[0] != '@')
					++cmp ;
				else
					continue;
			}
		}
	inFile.close();
	return cmp;
}

void outStats(std::string file, int id, int length, int nbAlignments, int nbNodes,int time){
	std::ofstream outFile(file.c_str(),std::ios_base::app);
	outFile << id << "\t" << length << "\t"<< nbAlignments << "\t"<<nbNodes <<"\t" << time << "\n";
	outFile.close();
}

void Graph::degreeNodes(std::string file){
	std::ofstream outFile(file.c_str());
	for( std::map<size_t,Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it){
		int degree = it->second.getFtoF().size();
		degree += it->second.getFtoR().size();
		degree += it->second.getRtoF().size();
		degree += it->second.getRtoR().size();
		degree += it->second.getPreviousFtoF().size();
		degree += it->second.getPreviousFtoR().size();
		degree += it->second.getPreviousRtoF().size();
		degree += it->second.getPreviousRtoR().size();
		outFile << it->first <<"\t" << it->second.getName() << "\t" << degree << "\n";
	}
	outFile.close();
}
