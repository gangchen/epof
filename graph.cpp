// graph.cpp : Implementations of the classes used in the graph
//

#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <queue>
#include <algorithm>
#include <iostream>
#include <map>
#include <stack>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
using namespace std;

void OutMemoryError(const string message)
{
	ofstream errorFile;
	errorFile.open("error.txt");
	errorFile<<"Out Of Memory Error:"<<endl;
	errorFile<<message<<endl;
	exit(-1);
}

Graph::Graph(const char* szFileName)
{
	trace=false;
	strcpy(m_szFileName,szFileName);
	m_szFileName[strlen(m_szFileName)-4]='\0';
	m_strFileName=szFileName;	//variables of string type can get their values by simple equalation???
	string::size_type	pos=m_strFileName.find('.');	//the method: find()
	m_strFileName=m_strFileName.substr(0,pos);
	GraphLoader	rawGraph(szFileName);	//C++ criterion
	LoadGraphFromRawGraph(rawGraph);	
}
Graph::~Graph()
{
	ClearNodes();
	ClearComplex();
}

void Graph::ClearNodes()
{
	for(size_t i=0;i<m_NodeArray.size();++i)
	{
		if(m_NodeArray[i]!=0)	//  !=0??? not !=NULL???
			delete m_NodeArray[i];
		m_NodeArray[i]=0;
	}
	m_NodeArray.clear();
}
void Graph::ClearCliques()
{
	for(size_t i=0;i<m_CliqueArray.size();++i)
	{
		if(m_CliqueArray[i]!=0)
			delete m_CliqueArray[i];	//free memory
		m_CliqueArray[i]=0;
	}
	m_CliqueArray.clear();
}
void Graph::ClearComplex()
{
	for(size_t i=0;i<m_ComplexArray.size();i++)
	{
		if(m_ComplexArray[i]!=0)	//  !=0??? not !=NULL???
			delete m_ComplexArray[i];
		m_ComplexArray[i]=0;
	}
	m_ComplexArray.clear();
}
void Graph::ClearComplex(vector<Complex *> &complexes)
{
	for(size_t i=0;i<complexes.size();i++)
	{
		if(complexes[i]!=0)
			delete complexes[i];
		complexes[i]=0;
	}
	complexes.clear();
}

void Graph::LoadGraphFromRawGraph(GraphLoader&	rawGraph)	//load all the needed information from the raw-graph
{
	for(int i=0;i<rawGraph.m_nArcs;i++)		//read information from rawGraph
	{
		int j,k;
		bool bFoundFrom=false;		//if the start node exists
		bool bFoundTo=false;		//if the end node exists
		for(j=0;j<m_NodeArray.size();j++)	//check if the start node exists, j is used to store the position
		{
			if(!strcmp(m_NodeArray[j]->m_szName,rawGraph.m_RawGraph[i].m_szFrom))
			{
				bFoundFrom=true;
				break;
			}
			
		}		
		for(k=0;k<m_NodeArray.size();k++)	//check if the end node exists, k is used to store the position
		{
			if(!strcmp(m_NodeArray[k]->m_szName,rawGraph.m_RawGraph[i].m_szTo))
			{
				bFoundTo=true;
				break;
			}
		}
		if(bFoundFrom==false)		//if not create a node
		{
			Node* pNode=new Node(m_NodeArray.size(),rawGraph.m_RawGraph[i].m_szFrom);
			m_NodeArray.push_back(pNode);
			j=m_NodeArray.size()-1;						//j stores the position of the start node of the arc
		}
		if(bFoundTo==false)
		{
			Node* pNode=new Node(m_NodeArray.size(),rawGraph.m_RawGraph[i].m_szTo);
			m_NodeArray.push_back(pNode);
			k=m_NodeArray.size()-1;
		}
		m_NodeArray[j]->InsertArc(j,k,rawGraph.m_RawGraph[i].m_fWeight);    //add an arc in the link
		m_NodeArray[j]->m_iDegree++;
		m_NodeArray[k]->InsertArc(k,j,rawGraph.m_RawGraph[i].m_fWeight);   //add another arc
		m_NodeArray[k]->m_iDegree++;
		m_NodeArray[j]->m_pFirst->m_pImage = m_NodeArray[k]->m_pFirst;
		m_NodeArray[k]->m_pFirst->m_pImage = m_NodeArray[j]->m_pFirst;
	}
	m_nNumEdge=rawGraph.m_nArcs;
	m_nNumNode=m_NodeArray.size();
}

vector<string> Graph::strsplit(string str, string del){
  unsigned end=0;
  vector<string> fields;
  while(1){
    end = str.find(del);
    if(end > str.size()){
      fields.push_back(str);
      break;
    }
    fields.push_back(str.substr(0, end));
    str = str.substr(end+1, str.size());
  }
  return fields;
}


int Graph::loadKeyCliques(const char* szfileName){
  m_KeyCliqueArray.clear();
  ifstream InFile(szfileName);
  if(!InFile){
    cout << "Cannot open the key clique file" << endl;
    return -1;
  }
  int cliqueID = 0;
  string line;
  while(getline(InFile, line)){
    Clique *cq = new Clique(cliqueID++);
    vector<string> proteins = strsplit(line, " ");
    for(vector<string>::iterator it = proteins.begin();
	it != proteins.end();
	it++){
      char* proteinName = new char[strlen((*it).c_str())+1];
      strcpy(proteinName, (*it).c_str());
      cq->addNode(locateNode(proteinName));
    }
    m_KeyCliqueArray.push_back(cq);
  }
  InFile.close();
  return m_KeyCliqueArray.size();
}

int Graph::loadKeys(const char* szfileName)
{
	m_KeyArray.clear();
	ifstream InFile(szfileName);
	if(!InFile)
	{
		cout<<"Can't Open The key File 'key.txt'!"<<endl;
		return -1;
	}
	char	m_szName[16];	//the name of the node
	int index;
	GraphArc arc;
	for(;;)
	{
		if(!InFile)
			break;
		InFile>>m_szName;
		if(strlen(m_szName)>0){
		  //cout << m_szName << endl;
			index = locateNode(m_szName);
			if(index<0){
				cout<<"Error!! Can't load node ' "<<m_szName<<" '"<<endl;
				return -1;
			}
			m_KeyArray.insert(index);
		}
		InFile.ignore(100,'\n');
	}
	InFile.close();
	m_nKeyNum = (int)m_KeyArray.size();
/*
	for(size_t i=0;i<m_nKeyNum;i++)
	{
		cout<<m_KeyArray[i]<<endl;
	}
*/
	return 0;
}
int Graph::locateNode(char m_szName[16])
{
	for(int i=0;i<m_NodeArray.size();i++){
		if(!(strcmp(m_NodeArray[i]->m_szName,m_szName)))
			return i;
	}
	return -1;
}

vector<int> Graph::getNeighbors(int node)
{
	vector<int> ret;
	Node* pNode=m_NodeArray[node];
	Arc* pArc=pNode->m_pFirst;
	while(pArc!=NULL){
		ret.push_back(pArc->m_iNodeTo);
		pArc=pArc->m_pNextArc;
	}
	return ret;
}
bool Graph::searchInClique(Clique* pc, int node)
{
	for(int i=0;i<pc->m_CliqueNodes.size();++i)
		if(node == pc->m_CliqueNodes[i])
			return true;
	return false;
}
vector<int> Graph::getNeighbors(Clique* pc)
{
	vector<int> ret;
	int flag;

	for(int i=0;i<pc->m_CliqueNodes.size();i++){
		vector<int> temp = getNeighbors(pc->m_CliqueNodes[i]);
		for(size_t k=0;k<temp.size();k++){	//add non-duplicated nodes
			if(!searchInClique(pc,temp[k])){	//the node must not be contained in pc
				flag = 0;
				for(size_t l=0;l<ret.size();l++){
					if(temp[k] == ret[l]){
						flag = 1;
						break;
					}
				}
				if(flag == 0)	ret.push_back(temp[k]);
			}
		}
	}

	return ret;
}

int Graph::distanceNode2Node(Node* sNode, Node* dNode){
  vector<bool> visited;
  vector<int> depth;
  for(int i = 0; i < m_nNumNode; i++){
    visited.push_back(false); 
    depth.push_back(0);
  }
  queue<Node*> q;
  q.push(sNode);
  Node* curNode;
  while(q.size() != 0){
    curNode = q.front();
    q.pop();
    vector<int> neighbors = getNeighbors(curNode->getID());
    for(vector<int>::iterator it = neighbors.begin();
	it != neighbors.end();
	it++){
      if(visited[*it]) continue;
      depth[*it] = depth[curNode->getID()]+1;
      if(dNode->getID() == *it){
	return depth[*it];
      }
      q.push(m_NodeArray[*it]);
    }
  }
  return 100000;
}

bool Graph::within(int node, Clique* cq, int dist=2){
  int minDist = 0;
  for(vector<int>::iterator cqIt = (*cq).m_CliqueNodes.begin();
      cqIt != (*cq).m_CliqueNodes.end();
      cqIt++){
    Node* sNode = m_NodeArray[node];
    Node* dNode = m_NodeArray[*cqIt];
    if(distanceNode2Node(sNode, dNode) > dist){
      return false;
    }
  }
  return true;
}

void Graph::generateCliques(){

  cout << "OK" << endl;
  for(vector<Clique*>::iterator itClique = m_KeyCliqueArray.begin();
      itClique != m_KeyCliqueArray.end();
      itClique++){

    cout << (*itClique)->getID() << endl;//for progress monitor

    Clique* keyClique = *itClique;
    Clique* pc=new Clique(m_CliqueArray.size());
    pc->m_CliqueNodes = (*keyClique).m_CliqueNodes;

    vector<int> neighbors = getNeighbors(pc);
    while(neighbors.size()>0){
      //find the node with max fitness
      double maxFitness = -1, maxKeyFitness = -1, maxNonFitness = -1;
      int keyNode = -1, node = -1, nonNode = -1;
      for(vector<int>::iterator itNeighbor = neighbors.begin();
	  itNeighbor != neighbors.end();
	  itNeighbor++){
	pc->m_CliqueNodes.push_back(*itNeighbor);
	double curFitness = calFitness(*itNeighbor, pc->m_CliqueNodes);
	erase(pc->m_CliqueNodes, *itNeighbor);
	if(curFitness > maxFitness){
	  node = *itNeighbor;
	  maxFitness = curFitness;
	  if(m_KeyArray.find(node) != m_KeyArray.end()){
	    keyNode = node;
	    maxKeyFitness = maxFitness;
	  }else{
	    nonNode = node;
	    maxNonFitness = maxFitness;
	  }
	}
      }

      if(maxFitness < 0 ){
	break;
      }
      
      //key protein first
      if(keyNode != -1 && maxKeyFitness > 0){
      	node = keyNode;
      }else if(maxFitness < 0.015){
      	break;
      }

      // non-key protein first
      /*if(nonNode != -1 && maxNonFitness > 0){
      	node = nonNode;
      }else if(maxFitness < 0.015){
      	break;
	}*/
      
      // ignore node that far away
      //if(!within(node, pc, 4)){
      // 	erase(neighbors, node);
      //	continue;
      //}

      pc->m_CliqueNodes.push_back(node);
				
      if(calFitness(node,pc->m_CliqueNodes)<0){
	erase(pc->m_CliqueNodes,node);
	erase(neighbors,node);
	cout << "ERROR" << endl;
	cin.get();
      }else{
	neighbors = getNeighbors(pc);
      }
    }//while
    pc->m_iNumNodes = pc->m_CliqueNodes.size();
    m_CliqueArray.push_back(pc);
    

    for(vector<Clique*>::iterator itNext = itClique + 1;
	itNext != m_KeyCliqueArray.end();){
      bool flag = true;
      for(vector<int>::iterator itNode = (*(*itNext)).m_CliqueNodes.begin();
	  itNode != (*(*itNext)).m_CliqueNodes.end();
	  ){
	if(searchInClique(pc, *itNode)){
	  (*(*itNext)).m_CliqueNodes.erase(itNode);
	  if((*(*itNext)).m_CliqueNodes.size() == 1){
	    m_KeyCliqueArray.erase(itNext);
	    flag = false;
	  }
	}else{
	  itNode++;
	}
      }
      if(flag){
	itNext++;
      }

      /*if(flag){
	m_KeyCliqueArray.erase(itNext);
      }else{
	itNext++;
	}*/
    }

  }//for
}


double Graph::calFitness(int node, vector<int> graph)	//node is included in graph
{
	vector<int> temp;
	for(size_t i=0;i<graph.size();++i)
		if(graph[i]!=node)
			temp.push_back(graph[i]);
	double f1,f2;
	f1 = calModularity(graph);
	f2 = calModularity(temp);
	return (f1-f2);
}

void Graph::keyFilterByModularity(double threshold){
  for(vector<Clique*>::iterator itKeyClique = m_KeyCliqueArray.begin();
      itKeyClique != m_KeyCliqueArray.end();)
    {
    double modularity = calModularity( (*(*itKeyClique)).m_CliqueNodes );
    if(modularity < threshold){
      m_KeyCliqueArray.erase(itKeyClique);
    }else{
      itKeyClique++;
    }
  }
}


void Graph::sortKeyBySize(){
  sort(m_KeyCliqueArray.begin(), m_KeyCliqueArray.end());
}



double Graph::calModularity(vector<int> graph)
{
	double inWeight=0.0,outWeight=0.0;
	sort(graph.begin(), graph.end());// for faster search
	//for(size_t i=0;i<graph.size();++i){
	for(vector<int>::iterator it = graph.begin();
	    it != graph.end();
	    it++){
		Node* pNode=m_NodeArray[*it];
		Arc* pArc=pNode->m_pFirst;
		while(pArc!=NULL){//for all the adjacent nodes
		  if(binary_search(graph.begin(), graph.end(), pArc->m_iNodeTo))
			inWeight += pArc->m_fWeight;
		  else
			outWeight += pArc->m_fWeight;
		  pArc=pArc->m_pNextArc;
		}
	}
	double mod = ((inWeight)/(inWeight+outWeight));
	return mod;
}
bool Graph::searchNode(vector<int> graph, int node)
{
  return binary_search(graph.begin(), graph.end(), node);
}

bool Graph::exists(int index, vector<Clique*> cliques)
{
	for(size_t i=0;i<cliques.size();i++)
		for(size_t j=0;j<cliques[i]->m_CliqueNodes.size();j++)
			if(cliques[i]->m_CliqueNodes[j]==index)
				return true;
	return false;
}
void Graph::erase(vector<int> &nodes,int node)
{
	vector<int>::iterator it = nodes.begin( );
	for ( size_t i=0;i<nodes.size();i++ ){
		if(nodes[i] == node){
			nodes.erase(it);
			break;
		}
		it++;
	}
}

int Graph::searchInComplexes(vector<Complex *> subgraph,int node)
{
	int ret=0;
	int i,j;
	for(i=0;i<subgraph.size();++i)
		for(j=0;j<subgraph[i]->m_iNumNodes;++j)
			if(node==subgraph[i]->m_ComplexNodes[j])
				ret++;
	return ret;
}

//initialize all the maximal cliques larger than k as original complexes
void Graph::initializeComplex()
{
	Clique *pc;
	int i;
	//set the maximal cliques larger than setted threshold each as a complex
	for(i=0;i<m_CliqueArray.size();++i){
		pc=m_CliqueArray[i];
		if(pc->m_iNumNodes >= m_nCliThr){
			Complex *pcx=createComplex(m_ComplexArray.size(),
				"In initializeComplex():EAGLE, creating complex...");
			pc->m_bSubordinate=false;
			pcx->m_iNumNodes=pc->m_iNumNodes;
			for(int j=0;j<pc->m_iNumNodes;++j)
				pcx->m_ComplexNodes.push_back(pc->m_CliqueNodes[j]);
			m_ComplexArray.push_back(pcx);
		}
		else
			pc->m_bSubordinate=true;
	}
	//set each of the other sobordinate nodes as a complex
/*	int node,j;
	for(i=0;i<m_CliqueArray.size();++i){
		pc=m_CliqueArray[i];
		if(pc->m_bSubordinate)
			for(j=0;j<pc->m_iNumNodes;++j){
				node=pc->m_CliqueNodes[j];	//the node to be searched
				if(searchInComplexes(m_ComplexArray,node)==0){	//a subordinate node
					Complex *pcx=createComplex(m_ComplexArray.size(),
						"In initializeComplex():EAGLE, creating complex...");
					pcx->m_iNumNodes=1;
					pcx->m_ComplexNodes.push_back(node);
					m_ComplexArray.push_back(pcx);
				}
			}
	}*/
	ClearCliques();
	//dumpInfo(m_ComplexArray,"initialComplexes_Fitness.txt");
	dumpCompleteComplexInfo(m_ComplexArray,"initialComplexes_Fitness.txt");
}

void Graph::generateComplexes()
{
  //float mod,cur;
	vector<Complex *> optimal;
	vector<Complex *> initial;

	ClearComplex();
	m_startT = clock();
	initializeComplex();	//initialize the complex set
	initial = copyComplexSet(m_ComplexArray);	//get the intial set
	optimal = initial;	//set the initial set as optimal set
	/*mod = calModularity();

	//merge the set of complexes until they are merged as a whole
		do{
		m_ComplexArray = mergeComplex(m_ComplexArray);
		cur = calModularity();
		if(cur>mod){
			mod=cur;
			ClearComplex(optimal);	//delete the prior set of complexes	
			optimal = copyComplexSet(m_ComplexArray);	//get the new set
		}
	}while(m_ComplexArray.size()>1);
	m_endT = clock();
	//dump result
	dumpEagleComplex(optimal);
	//dump_EagleComplex(optimal);
	ClearComplex(optimal);*/

}

vector<Complex*> Graph::mergeComplex(vector<Complex*> complexes)//the number of complexes must be no less than two
{
	vector<Complex*> ret;
	double max=0,simval;
	int index1=0,index2=1;
	int i,j;
	//get the maximum similarity
	max=calSimilarity(complexes[0],complexes[1]);//what's the proper initial value for max??//
	for(i=0;i<complexes.size()-1;++i){
		for(j=i+1;j<complexes.size();++j){
			simval=calSimilarity(complexes[i],complexes[j]);
			if(simval>max){	//switch
				max=simval;
				index1=i;
				index2=j;
			}
		}
	}
	Complex *pc1=complexes[index1];
	Complex *pc2=complexes[index2];
	if(trace)	recordMerge(pc1, pc2);	//keep record of the merging process
	//merge the two selected complexes
	if(pc1->m_iNumNodes >= pc2->m_iNumNodes){
		mergeComplex(pc1,pc2);
		for(i=0;i<complexes.size();++i){
			if ( i!=index2 )
				ret.push_back(complexes[i]);
		}
	}
	else{
		mergeComplex(pc2,pc1);
		for(i=0;i<complexes.size();++i){
			if ( i!=index1 )
				ret.push_back(complexes[i]);
		}
	}
	return ret;
}

//merge complex pc2 into complex pc1
void Graph::mergeComplex(Complex* pc1, Complex* pc2){
	int i,j,flag;
	for(i=0;i<pc2->m_iNumNodes;++i){
		flag=0;
		for(j=0;j<pc1->m_iNumNodes;++j){
			if(pc2->m_ComplexNodes[i] == pc1->m_ComplexNodes[j]){
				flag=1;	//this node already exists in C1
				break;
			}
		}
		if(flag==0){
			pc1->m_iNumNodes++;
			pc1->m_ComplexNodes.push_back(pc2->m_ComplexNodes[i]);
		}
	}
	delete pc2;
}

void Graph::recordMerge(Complex* pc1, Complex* pc2){
	FILE * fp;
	fp=fopen("process_EALGE.txt", "a");
	if( fp == NULL) 
		cout<<"\t\tCan't open file process_EAGLE.txt\n"<<endl; 
	else{
		int i;
		fprintf(fp,"merge Complex%d with Complex%d:\n",pc1->m_iComplexNo,pc2->m_iComplexNo);
		for(i=0;i<pc1->m_ComplexNodes.size();++i)
			fprintf(fp,"%s ",m_NodeArray[pc1->m_ComplexNodes[i]]->m_szName);
		fprintf(fp,"   ");
		for(i=0;i<pc2->m_ComplexNodes.size();++i)
			fprintf(fp,"%s ",m_NodeArray[pc2->m_ComplexNodes[i]]->m_szName);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

vector<Complex*> Graph::copyComplexSet(vector<Complex*>complexes){
	vector<Complex *> temp;
	temp.reserve(complexes.size());
	int i,j;
	for(i=0;i<complexes.size();++i){
		Complex * comp=createComplex(i,"In copyComplexSet():EAGLE, creating complex...");
		comp->m_iNumNodes=complexes[i]->m_iNumNodes;
		for(j=0;j<complexes[i]->m_iNumNodes;++j)
			comp->m_ComplexNodes.push_back(complexes[i]->m_ComplexNodes[j]);
		temp.push_back(comp);
	}
	return temp;
}
double Graph::arcWeights()
{
	double w = 0;
	for(size_t i=0;i<m_NodeArray.size();i++){
		Node* pNode = m_NodeArray[i];
		Arc* pArc = pNode->m_pFirst;
		while(pArc!=NULL){
			if(pArc->m_iNodeFrom<pArc->m_iNodeTo)
				w += pArc->m_fWeight;
			pArc = pArc->m_pNextArc;
		}
	}
	return w;
}
double Graph::weightedDegree(int index)
{
	double wd = 0;
	Node* pNode = m_NodeArray[index];
	Arc* pArc = pNode->m_pFirst;
	while(pArc!=NULL){
		wd += pArc->m_fWeight;
		pArc = pArc->m_pNextArc;
	}
	return wd;
}

/**calculate the simiarities between all pairs of complexes which have been constructed before
 ****
 *S=1/2m*( sumof(Aij-ki*kj/2m) )	of which i IEO C1,j IEO C2, and i!=j
 */
double Graph::calSimilarity(Complex* pc1,Complex* pc2)
{
	double ret=0.0,temp;
	int n1,n2;
	double adj;
	double m = arcWeights();
	for(int i=0;i<pc1->m_ComplexNodes.size();++i){
		n1=pc1->m_ComplexNodes[i];
		for(int j=0;j<pc2->m_ComplexNodes.size();++j){
			n2=pc2->m_ComplexNodes[j];
			if(n1!=n2){//pn1 and pn2 is not the same node
				adj = isAdjacent(n1,n2);
				temp = weightedDegree(n1)*weightedDegree(n2);
				temp = -temp/(2*m);
				temp += adj;
				ret += temp;
			}
		}
	}
	//ret=ret/(2*m);
	//cout<<"\tThe similarity: "<<ret<<endl;
	return ret;
}
//evaluate the quality of the resulting complexes, by which we decide which division is optimal
/*
	EQ=1/2m*(sumof( (Aij-ki*kj/2m)/(Oi*Oj) ))	of which i,j IEO Cx
 */
double Graph::calModularity()
{
	double ret=0.0,temp;
	int i,j,k;
	double m = arcWeights();
	int n1,n2,c1,c2;
	double adj;
	for(i=0;i<m_ComplexArray.size();++i){
		for(j=0;j<m_ComplexArray[i]->m_iNumNodes-1;++j){
			n1=m_ComplexArray[i]->m_ComplexNodes[j];
			for(k=j+1;k<m_ComplexArray[i]->m_iNumNodes;++k){
				n2=m_ComplexArray[i]->m_ComplexNodes[k];
				/*cout<<"Complex: "<<i<<" n1: "<<n1<<" n2: "<<n2<<endl;*/
				adj=isAdjacent(n1,n2);
				c1=searchInComplexes(m_ComplexArray,n1);
				c2=searchInComplexes(m_ComplexArray,n2);
				temp = weightedDegree(n1)*weightedDegree(n2);
				temp=-temp/(2*m);
				temp+=adj;
				temp=temp/c1/c2;
				ret+=temp;
			}
		}
	}
	//ret=ret/(2*m);
	//cout<<"\tThe modularity: "<<ret<<endl;//
	return ret;
}
double Graph::isAdjacent(int a, int b)
{
	Node* pNode=m_NodeArray[a];
	Arc* pArc=pNode->m_pFirst;
	while(pArc!=NULL){
		if(pArc->m_iNodeTo==b){
			return pArc->m_fWeight;
		}
		pArc=pArc->m_pNextArc;
	}
	return 0;
}

Complex* Graph::createComplex(int id, const string message){
	Complex *pcx=new Complex(id);
	if(pcx==NULL)
		OutMemoryError("In initializeComplex():EAGLE, creating complex...");
	return pcx;
}

void Graph::dumpNodeInfo(const char* szFileName)
{
	ofstream OutFile;
	OutFile.open(szFileName);
	OutFile<<"There are "<<m_NodeArray.size()<<" nodes in total"<<endl;
	OutFile<<"Index\tName\tDegree"<<endl;
	size_t i;
	for (i=0; i<m_NodeArray.size(); ++i)
		OutFile<<"N"<<i+1<<":\t"<<m_NodeArray[i]->m_szName<<"\t"<<m_NodeArray[i]->m_iDegree<<" "<<endl;
	OutFile<<"================================================="<<endl;
	int arcnum=0;
	for (i=0; i<m_NodeArray.size(); ++i)
	{
		Node* pNode = m_NodeArray[i];
		Arc*  pArc  = pNode->m_pFirst;
		while(pArc!=NULL)
		{
			if(pNode->m_iNode < pArc->m_iNodeTo){
				arcnum++;
				OutFile<<"Arc"<<arcnum<<": "<<pNode->m_szName<<" "<<m_NodeArray[pArc->m_iNodeTo]->m_szName<<" "<<endl;
			}
			pArc = pArc->m_pNextArc;
		}
	}
	cout<<"NodeNum: "<<m_NodeArray.size()<<endl;	
	cout<<"ArcNum: "<<arcnum<<endl;
	OutFile.close();
}

void Graph::getCliques()
{
	ofstream outFile;
	outFile.open("cliques.txt");
	outFile<<"We've got "<<m_CliqueArray.size()<<" cliques from the graph:"<<endl<<endl;
	outFile<<"Index\tSize\tMembers"<<endl;
	for(int i=0;i<m_CliqueArray.size();++i){
		outFile<<i+1<<"\t"<<m_CliqueArray[i]->m_iNumNodes<<"\t";
		for(size_t j=0;j<m_CliqueArray[i]->m_iNumNodes;++j)
		{
			int index=m_CliqueArray[i]->m_CliqueNodes[j];
			outFile<<m_NodeArray[index]->m_szName<<" ";
		}
		outFile<<endl;
	}
	outFile.close();
}

void Graph::getComplexes(const char* szFileName)
{
	ofstream outFile;
	outFile.open(szFileName);
	int counter=0,i;
	for(i=0;i<m_ComplexArray.size();++i)
		 if(m_ComplexArray[i]->m_ComplexNodes.size()>0)
			 counter++;
	outFile<<"Complex amount: "<<counter<<endl<<endl;
	//outFile<<"Index Size Indegree TotalDegree isModule Members"<<endl;
	counter=0;
	vector<int> nodes;
	for(i=0;i<m_ComplexArray.size();++i){
		if(m_ComplexArray[i]->m_ComplexNodes.size()>0){
			counter++;
			nodes = m_ComplexArray[i]->m_ComplexNodes;
			//stable_sort(nodes.begin(),nodes.end(),LowerByName);
			/*outFile<<"Complex"<<counter<<"\tsize: "<<m_ComplexArray[i]->m_iNumNodes
				<<"\tInDegree: "<<m_ComplexArray[i]->m_iInDegree
				<<"  TotalDegree: "<<m_ComplexArray[i]->m_iTotalDegree
				<<"  isModule? "<<m_ComplexArray[i]->m_bModule<<endl;*/
			outFile<<"Complex "<<counter<<"\n";
			for(size_t j=0;j<nodes.size();++j)
			{
				int index = nodes[j];
				outFile<<m_NodeArray[index]->m_szName<<"\n";
			}
		}
	}
	m_endT=clock();
	if(strcmp(szFileName,"initial_complexes.txt"))
		outFile<<endl<<endl<<"Time used: "<<m_endT-m_startT<<" ms."<<endl;
	outFile.close();
}

void Graph::dumpCompleteComplexInfo(vector<Complex *> complexes,char fileName[])
{
	ofstream outFile;
	outFile.open(fileName);
	//outFile<<"Size of initial complexes: "<<complexes.size()<<endl<<endl;
	for(int i=0;i<complexes.size();++i){
	//	outFile<<"Complex"<<i<<"\tsize: "<<complexes[i]->m_iNumNodes
	//		<<"\tInDegree: "<<complexes[i]->m_iInDegree
	//		<<"  TotalDegree: "<<complexes[i]->m_iTotalDegree
	//		<<"  isModule? "<<complexes[i]->m_bModule<<endl;
		outFile<<"Complex "<<i<<"  "<<complexes[i]->m_iNumNodes<<endl;
		for(size_t j=0;j<complexes[i]->m_iNumNodes;++j)
		{
			int index=complexes[i]->m_ComplexNodes[j];
			outFile<<m_NodeArray[index]->m_szName<<endl;
		}
		//outFile<<endl;
	}
	outFile.close();
}

void Graph::dumpEagleComplex(vector<Complex *> optimal){
	ofstream outFile;
	outFile.open(m_eagleOut);
	outFile<<"Parameters"<<endl<<"\tclique threshold: "<<m_nCliThr<<endl<<endl;
	outFile<<"The optimal complexes:"<<endl;
	for(int i=0;i<optimal.size();++i){
		outFile<<"Complex "<<(i+1)<<"\n";
		for(int j=0;j<optimal[i]->m_iNumNodes;++j)
			outFile<<m_NodeArray[optimal[i]->m_ComplexNodes[j]]->m_szName<<"\n";
	}
	outFile<<"\n\nTime Used:  "<<m_endT-m_startT<<" ms.\n";
	outFile.close();
}

void Graph::dump_EagleComplex(vector<Complex *> optimal){
	char szFile[20];
	ofstream outFile;
	for(int i=0;i<optimal.size();i++){
		int num = i+1;
		sprintf(szFile,"%d",num); 
		strcat(szFile,".txt");
		outFile.open(szFile);
		for(int j=0;j<optimal[i]->m_iNumNodes;++j)
			outFile<<m_NodeArray[optimal[i]->m_ComplexNodes[j]]->m_szName<<"\n";
		outFile.close();
	}
	time_t te=clock();
	outFile.open("Time.txt");
	outFile<<"\nComplex amount: "<<optimal.size()<<"\nTime Used:  "<<m_endT-m_startT<<" ms.\n";
	outFile.close();
}

