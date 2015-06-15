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

#include "data.hpp"
#include "pw_alignment.hpp"
#include "needleman.hpp"

// Global var
const bool VERBOSE = true;
// ----------------------------
//		Graph
//-----------------------------
class Graph{

private:
//	class OrientedVertex;
	// ----------------------------
	//		Vertex
	//-----------------------------
	class Vertex{

	public:
		//Vertex();
		Vertex(int idVertex);
		//Vertex(Vertex const& v);
		//~Vertex();
		const int& getIdVertex()const;
		friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
		void printVertex(std::ostream &os)const;
		void setRtoF(const int& id);
		void setRtoR(const int& id);
		void setFtoF(const int& id);
		void setFtoR(const int& id);
		void setCostScore(const double& score);
		void setDistance(const double& distance);
		void setPrevious(const int& idPrevious);
		void setStartOnRead(const int& start);
		void setEndOnRead(const int& end);
		void setName(const std::string& name);
		const std::vector<int>& getRtoF()const;
		const std::vector<int>& getRtoR()const;
		const std::vector<int>& getFtoF()const;
		const std::vector<int>& getFtoR()const;
		const std::vector<int>& getPrevious();
		const double& getCostScore();
		const double& getDistance();
		const int& getStartOnRead();
		const int& getEndOnRead();
		const std::string& getName();
		bool operator < (const Graph::Vertex& v)const
				    {
				        return (distance < v.distance);
				    }

	private:
		int idVertex;
		std::string name;
		bool visited;
		std::vector<int> previous;
		double costScore;
		double distance;
		std::vector<int> RtoF;
		std::vector<int> RtoR;
		std::vector<int> FtoF;
		std::vector<int> FtoR;
		int startOnRead;
		int endOnRead;
		//TODO add length of the node, I can stop using data in a lot of function
		//friend void Graph::OrientedVertex::getSuccessors(std::vector<Graph::OrientedVertex>& oriented) const{
	};

	class OrientedVertex{
	public:
		//OrientedVertex();
		OrientedVertex(const Graph::Vertex& v, bool b);
		bool vertexIsForward(pw_alignment p);
		const std::vector<int>& getVertexFtoF(Graph::Vertex v)const;
		const std::vector<int>& getVertexFtoR(Graph::Vertex v)const;
		const std::vector<int>& getVertexRtoR(Graph::Vertex v)const;
		const std::vector<int>& getVertexRtoF(Graph::Vertex v)const;
		void getSuccessors(std::vector<Graph::OrientedVertex>& oriented) const;
		void getPredecessors(std::vector<OrientedVertex>) const;
		const Graph::Vertex& getVertex()const ;
		std::string getDna() const;
	
	private:
		Vertex vertex;
		bool isForward;

	};

public:
	//Graph();
	//Graph(Graph const& g);
	//~Graph();
	void setStart(const int&id, const int&start);
	void setEnd(const int&id, const int&end);
	void readDotFile(std::string dotFile,std::map<std::string, size_t> longname2seqidx); // Graph from a dot file -> Transform to constructor ?
	void addVertex(Vertex v);
	void addEdge(int typeOfEdge,int id1,std::string name1, int id2, std::string name2);
	const std::map<int,Vertex>& getVertices()const;
	void setScore(const int& id, const double& score);
	const double& getScore(const int& id);
	const Vertex& getVertex(const int& id);
	friend std::ostream& operator<<(std::ostream &os, Graph & g);
	friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
	void printGraph(std::ostream &os);
	void updateDistance(Graph::Vertex& v1, Graph::Vertex& v2, std::map<int,int>& previous);
	std::vector<int> dijkstra(int start,int end);
	void initAllCostScore(all_data data, Graph& newGraph);
	void lookPrevious(std::string partOfRead,Graph::Vertex it, all_data& data, Graph& newGraph);
	void lookNext(std::string partOfRead,Graph::Vertex it, all_data& data, Graph& newGraph);
	void DFS(Graph::Vertex v1, int target, int length, all_data data, std::vector<int> path);
	void lookGap(std::string gapSeq, Graph::Vertex v1, int idEnd, all_data data,Graph& newGraph);
	void parseData(all_data& data, Graph& newGraph);
	friend bool operator== ( const Vertex &v1, const Vertex &v2);
	std::string extractPartOfSeq(dnastring seq, int start, int end);
	void findFinalPath(size_t lengthRead);
	int findIdAlignment(Graph::Vertex v, all_data data);
//	std::vector<OrientedVertex>
private:
	std::map<int,Vertex> vertices;
	//std::vector<OrientedVertex> OV;
};
/*
class sortVertexByStart{
	public:
	bool operator () (const Vertex &v1, const Vertex &v2 )const;
};
*/
#endif // GRAPH_HPP


