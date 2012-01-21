// graph.h : Defines the classes defined in this application
//

#ifndef GRAPH_H
#define GRAPH_H
#include "graphLoader.h"
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cstring>
using namespace std;

class Complex;
class Arc
{
public:
	int m_iNodeTo;		//the ID of the node at one end
	int m_iNodeFrom;	//the ID of the other node of this arc
	Arc* m_pNextArc;	//the point of the next arc
	Arc* m_pImage;		//image arc: the first arc associated with the node at the other end of this arc
	double	m_fWeight;		//the weight of an graph-arc

public:
	Arc(int iNodeFrom, int iNodeTo, float fWeight, Arc* pNextArc)
	{
		m_iNodeFrom = iNodeFrom;
		m_iNodeTo=iNodeTo;
		m_fWeight=fWeight;
		m_pNextArc=pNextArc;
	}
	~Arc()
	{
		if(m_pNextArc)	delete m_pNextArc;
	}
};
class Node
{
public:
	int		m_iNode;		//the ID of the node
	char	m_szName[16];	//the name of the node
	int		m_iDegree;		//the degree of the node
	Arc*	m_pFirst;		//the point of the first arc associated with this node
	vector<int> subComplexes;	//array of complexes including this node

public:
	Node(int iNode,char szName[16])
	{
		strcpy(m_szName,szName);
		m_iNode=iNode;
		m_pFirst=0;
		m_iDegree=0;
	}

	~Node(){if(m_pFirst) delete m_pFirst;}
	void InsertArc(int iFrom, int iTo,float fWeight)
	{
		Arc* newArc=new Arc(iFrom,iTo,fWeight,m_pFirst);
		m_pFirst=newArc;
	}
};
class Clique	//the structure of the clique
{
public:
	vector<int>			m_CliqueNodes;		//the vector of the nodes in this clique
	int					m_iNumNodes;		//the number of nodes
	bool				m_bSubordinate;		//the flag showing if this clique is subordinate
	int					m_iCliqueID;		//the ID of this clique
public:
	Clique(int ID)
	{
		m_iNumNodes=0;
		m_iCliqueID=ID;
	}
	void sortNodes()
	{
		sort(m_CliqueNodes.begin(),m_CliqueNodes.end());
	}
};
class Complex		// the structure of the Complex
{
public:
	vector<int>			m_ComplexNodes;		//the vector of all the nodes in the complex
	vector<Arc*>		m_Arcs;				//the vector of the arc pointers
	int					m_iComplexNo;		//the ID of the complex
	int					m_iNumNodes;		//the number of the nodes
	int                 m_iInDegree;		//the indegree of the complex
	int                 m_iTotalDegree;		//the total degree of the complex
	bool				m_bMergeable;		//the flag showing if the complex is mergeable
	bool                m_bModule;			//the flag showing if this complex can be defined as a module
public:
	Complex(int iComplexNo)
	{	
		m_iComplexNo = iComplexNo;
		m_iNumNodes = 0;
		m_iInDegree = 0;
		m_iTotalDegree = 0;
		m_bMergeable = true;
		m_bModule = false;
	}
};

class Graph
{
public:
	char					m_szFileName[120];	//the name of the input file including the information of the graph
	string					m_strFileName;		//the name of the output file		
	char					m_eagleOut[60];		//the name of eagel output file
	vector<Node*>			m_NodeArray;		//the vector of the node pointers in this graph	
	vector<Clique*>			m_CliqueArray;		//the vector of the clique pointers	
	vector<Complex*>		m_ComplexArray;		//the vector of the complex pointers

	int						m_nNumNode;			//the number of Nodes	
	int						m_nNumEdge;			//the number of Edges
	int						m_nNumClique;		//the number of Cliques	
	int						m_nNumComplex;		//the number of Complexes
	bool					trace;				//record the merging process of not
	time_t					m_startT;
	time_t					m_endT;

	int						m_nCliThr;			//the clique threshold
	vector<int>				m_KeyArray;			//the key proteins
	int						m_nKeyNum;			//the number of keys
public:
	~Graph();
	void ClearNodes();
	void ClearCliques();
	void ClearComplex();
	void ClearComplex(vector<Complex *> &complexes);

	//functions to load file
	Graph(const char* szFileName);
	void LoadGraphFromRawGraph(GraphLoader&	rawGraph);
	int loadKeys(const char* szfileName);
	int locateNode(char m_szName[16]);
	//fuctions to generate cliques
	void generateCliques();
	vector<int> getNeighbors(int node);
	vector<int> getNeighbors(Clique* pc);
	bool searchInClique(Clique* pc, int id);
	bool exists(int index, vector<Clique*> cliques);
	void erase(vector<int> &nodes,int node);
	double calFitness(int node, vector<int> graph);
	double calModularity(vector<int> graph);
	bool searchNode(vector<int> graph, int node);

	int searchInComplexes(vector<Complex *> subgraph,int node);
	Complex* createComplex(int id, const string message);
	//called when use clustering algorithm
	void initializeComplex();
	void generateComplexes();
	vector<Complex*> mergeComplex(vector<Complex*> complexes);
	void mergeComplex(Complex* pc1, Complex* pc2);
	vector<Complex*> copyComplexSet(vector<Complex*> complexes);
	double calSimilarity(Complex* pc1,Complex* pc2);
	double calModularity();
	double arcWeights();
	double weightedDegree(int index);
	double isAdjacent(int a, int b);

	void getCliques();
	void getComplexes(const char* szFileName);
	void recordMerge(Complex* pc1, Complex* pc2);
	void dumpEagleComplex(vector<Complex *> optimal);
	void dump_EagleComplex(vector<Complex *> optimal);
	//functions to output informations
	void dumpNodeInfo(const char* szFileName);	//dump the node information
	void dumpCompleteComplexInfo(vector<Complex *> complexes,char fileName[]);
};
//print out of memory error
void OutMemoryError(const string message);

#endif
