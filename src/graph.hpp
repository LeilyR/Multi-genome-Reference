/*
 * needleman.hpp
 *
 *  Created on: May 12, 2015
 *      Author: Marion Dubarry
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <typeinfo>
#include <iterator>
#include <functional>

#include "data.hpp"
#include "pw_alignment.hpp"
#include "dynamic_mc.hpp"
#include "needleman.hpp"

// Global var
const bool VERBOSE = false;
//typedef dynamic_mc_model use_model; TODO chenage it to this one
typedef model use_model;

//const all_data  data;

// ----------------------------
//		Graph
//-----------------------------
class Graph{

private:
	// ----------------------------
	//		Vertex
	//-----------------------------
	class Vertex{

	public:
		Vertex(size_t idVertex);
		Vertex(const Vertex& v);
		const size_t& getIdVertex()const;
		friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
		void printVertex(std::ostream &os)const;
		void printVertex()const;
		void setRtoF(const size_t& id);
		void setRtoR(const size_t& id);
		void setFtoF(const size_t& id);
		void setFtoR(const size_t& id);
		void setPreviousRtoF(const int& v);
		void setPreviousRtoR(const int& v);
		void setPreviousFtoF(const int& v);
		void setPreviousFtoR(const int& v);
		void setCostScore(const int& id,const double& score);
		void setDistance(const double& distance);
		void setPrevious(const int& idPrevious);
		void setStartOnRead(const int& start);
		void setEndOnRead(const int& end);
		void setStartOnNode(const int& start);
		void setEndOnNode(const int& end);
		void setName(const std::string& name);
		const std::vector<int>& getRtoF()const;
		const std::vector<int>& getRtoR()const;
		const std::vector<int>& getFtoF()const;
		const std::vector<int>& getFtoR()const;
		const std::vector<int>& getPreviousRtoF()const;
		const std::vector<int>& getPreviousRtoR()const;
		const std::vector<int>& getPreviousFtoF()const;
		const std::vector<int>& getPreviousFtoR()const;
		double getCostScore(const int& id)const;
		double getCostScore()const;
		const double& getDistance()const;
		const size_t& getStartOnRead()const;
		const size_t& getEndOnRead()const;
		const size_t& getStartOnNode()const;
		const size_t& getEndOnNode()const;
		const std::string& getName()const;
		bool operator < (const Graph::Vertex& v)const
		{
			return (distance < v.distance);
		}

	private:
		size_t idVertex;
		std::string name;
		std::map<int,double> costScore; //id from the edge where it come from and the score
		double distance;
		std::vector<int> RtoF;
		std::vector<int> RtoR;
		std::vector<int> FtoF;
		std::vector<int> FtoR;
		std::vector<int> previousRtoF;
		std::vector<int> previousRtoR;
		std::vector<int> previousFtoF;
		std::vector<int> previousFtoR;
		int startOnRead;
		int endOnRead;
		int startOnNode;
		int endOnNode;
		//TODO add length of the node, I can stop using data in a lot of function
	};
public:
	class OrientedVertex{

	public:
		// Constructors
		OrientedVertex(const Graph::Vertex& v, const bool& b);
		OrientedVertex(const Graph::Vertex& v, const pw_alignment& p);

		void getSuccessors(std::vector<Graph::OrientedVertex>& oriented,Graph g) const;
		void getPredecessors(std::vector<OrientedVertex>& oriented,Graph g) const;
		Graph::Vertex& getVertex();
		const bool& getIsForward()const;
		std::string getDna(all_data & , size_t &) const;

		bool compareOrientedVertex( Graph::OrientedVertex &v1, Graph::OrientedVertex &v2);

		// regarding 2 OrientedVertex gives the type of edge connecting them
		int typeOfEdge(OrientedVertex ov2);
	private:
		Vertex vertex;
		bool isForward;

	};

public:
	// get/set functions
	void setStartOnRead(const size_t&id, const size_t&start);
	void setEndOnRead(const size_t&id, const size_t&end);
	void setStartOnNode(const size_t&id, const size_t&start);
	void setEndOnNode(const size_t&id, const size_t&end);
	void setScore(const size_t& id, const double& score,const size_t& idOrigin);
	const double& getScore(const size_t& id);
	const Vertex& getVertex(const size_t& id);
	const std::map<size_t,Vertex>& getVertices()const;
	void set_vertex(const size_t & seq_id, const std::string & seq_name);

	// read/write files
	void read_dot_file(std::string & , std::map< std::string, size_t> &, std::string &);
	void writeDotFile(std::string file, std::vector<Graph::OrientedVertex> path);
	void writeDotFile(std::string file, Graph graph);
	void degreeNodes(std::string file);
	// output the final path :
	void alignmentOut(std::string& , std::vector<OrientedVertex> path, all_data data,std::ostream &);


	void addVertex(Vertex v);
	void addEdge(int typeOfEdge, const size_t & id1, const std::string & name1, const size_t & id2, const std::string & name2);
//	void add_edge()
	void addEdge(OrientedVertex v1,OrientedVertex v2);




	// print the object graph
	void printGraph(std::ostream &os);
	void print();

	// principal function of the program,
		//looks if the seed are everywhere on the long read, if no splits the read, and call parseData
	void prepare_read(Graph& , std::vector<pw_alignment> &, const dnastring&, use_model &, all_data & data, std::ostream &);
		// loop on the seeds, find link and update new graph
	void parseData(all_data & , use_model & ,Graph& newGraph, std::vector<pw_alignment> & als, std::vector<pw_alignment> & all_als , dnastring seq);

	// There is 3 different dijkstra algo :
		//To find a path between 2 seeds
	std::vector<Graph::OrientedVertex> dijkstra(OrientedVertex orientedStart,OrientedVertex orientedEnd,size_t & read_id ,std::string part,all_data data, use_model model);
		//To find a path at the ends of the read
	std::vector<Graph::OrientedVertex> dijkstra(OrientedVertex orientedStart,std::string partOfRead,all_data data ,use_model model, std::string option,Graph newGraph, size_t & , size_t &);
		// To find the final path in the new graph
	std::vector<Graph::OrientedVertex> dijkstra(OrientedVertex orientedStart,OrientedVertex orientedEnd);
	std::vector<Graph::OrientedVertex> findPath(std::map<int,std::vector<Graph::OrientedVertex> >previous,OrientedVertex start, OrientedVertex end);
	// There is 3 different function to take a new node and create the alignment,
	// those functions are call by the corresponding Dijkstra
		// If we are at the 5' end of the read and want to look at the predecessors
	void lookPrevious(std::string partOfRead,Graph::OrientedVertex& previousNode, all_data data, size_t startOnread, use_model model, int optionNeedleman,int idOrigin, size_t &);
		// If we are at the 3' end of the read and want to look at the successors
	void lookNext(std::string partOfRead,Graph::OrientedVertex& ov, all_data data,size_t startOnRead,use_model model,int optionNeedleman,int idOrigin);
		// If we are in between
	void lookGap(size_t & read_id,std::string partOfRead,Graph::OrientedVertex& onode, all_data data,size_t startOnRead,use_model model,int optionNeedleman,int idOrigin);

	// Update the distance for the dijkstra algo
	void updateDistance(Graph::OrientedVertex& ov1, Graph::OrientedVertex& ov2, std::map<int,std::vector<Graph::OrientedVertex> >& previous);

	//Create the strat(node -1) and end (node -2) nodes for the new graph
	void set_first_and_last_node(size_t alsnumber, size_t start, size_t end);
//	bool checkVertex(std::string name, int start, int end,Graph g);
	void printFinalPath(std::vector<OrientedVertex> path);

	// add new node/edge in the newGraph
	void updateNewGraph(std::vector<Graph::OrientedVertex> path,all_data& data, size_t & ,std::string gap,use_model model,Graph& newGraph,int option,std::string whereToLook, std::vector<pw_alignment> & all_als);

	// compute a score for a new node
	double countCostScore(use_model m,pw_alignment p);

	// TODO change those functions : just use dnastring format and not string
	std::string extractPartOfSeq(dnastring seq, int start, int end);
	std::string convertDnastringToString(dnastring seq);


	//add connection between all nodes
	void addVirtualEdges(use_model model,dnastring read);

	const size_t get_distance(std::vector<pw_alignment> & , const size_t & , const size_t & , const size_t & , const size_t & , const size_t & , const size_t &)const;


	friend std::ostream& operator<<(std::ostream &os, Graph & g);
	friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
	friend bool operator== ( OrientedVertex ov1, OrientedVertex ov2);
private:
	std::map<size_t,Vertex> vertices;
};

/*
class sortOrientedVertex {
	public:
	bool operator()(Graph::OrientedVertex &a,  Graph::OrientedVertex &b) const ;
};
*/
/*
class cmpOrientedVertex {
public:
	bool operator()(Graph::OrientedVertex& ov1, Graph::OrientedVertex& ov2) const;
};
*/
template<class InputIterator, class T>
  InputIterator findVertex (InputIterator first, InputIterator last, const T& val);
//int testSamFile(std::string file);


bool customComp(Graph::OrientedVertex ov1, Graph::OrientedVertex ov2);
struct special_compare : public std::unary_function<Graph::OrientedVertex, bool>
{
  explicit special_compare(const Graph::OrientedVertex &baseline) : baseline(baseline) {}
  bool operator() (const Graph::OrientedVertex &arg)
  { return customComp(arg, baseline); }
  Graph::OrientedVertex baseline;
};

void outStats(std::string file, int id, int length, int nbAlignments,int nbNodes, int time);
#endif // GRAPH_HPP


