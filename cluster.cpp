//cluster.cpp
#include "graph.h"
using namespace std;

int main(int argc, char* argv[])
{
	char netFile[80];
	char compFile[80];
	char keyFile[80];
	if(argc!=3)	//get the number of arguments
	{
		cout<<"Command line parameter error!\n";
		cout<<"Usage: "<<argv[0]<<" %InputGraph%"<<endl;
		exit(-1);
	}
	//parameters used for clustering
	strcpy(netFile,argv[1]);
/*	cout<<"Input file name: ";
	cin>>netFile;*/
	strcpy(keyFile,argv[2]);
	strcpy(compFile, netFile);
	compFile[strlen(compFile)-4] = '\0';
	strcat(compFile,"_Complex_Fitness.txt");
	
	fstream InFile(netFile);
	if(!InFile)
	{
		cout<<"Can't Open Network "<<netFile<<endl;
		return -1;
	}
	//load graph...
	Graph graph(netFile);
	graph.loadKeys("key.txt");
	if(graph.loadKeyCliques(keyFile)<0)
		return -1;
	graph.keyFilterByModularity(0.04);
	graph.sortKeyBySize();
	//generate maximal cliques...
	graph.generateCliques();
	cout << "Generate Cliques Done" << endl;

	//generate complexes...
	//int cliThr = atoi(argv[2]);
	int cliThr = 1;
	graph.m_nCliThr = cliThr;
	strcpy(graph.m_eagleOut,compFile);
	graph.generateComplexes();/**/

	cout<<endl<<"Succeed!!!Congratulations!"<<endl<<endl;


	return 0;
}

