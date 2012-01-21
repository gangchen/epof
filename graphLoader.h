//graphLoader.h head file

#ifndef GRALOADER_H
#define GRALOADER_H
#include <cstdlib>
#include <string>
#include <vector>
#include <cstring>


using namespace std;

class GraphArc
{
public:
	char	m_szFrom[16];	//name of the source protein
	char	m_szTo[16];		//name of the target protein
	float	m_fWeight;		//the weight of an graph-arc

public:
	GraphArc()
	{
		m_fWeight=1.0f;
	}

	GraphArc(char szFrom[16],char szTo[16],float fWeight)
	{
		strcpy(m_szFrom,szFrom);
		strcpy(m_szTo,szTo);
		m_fWeight=fWeight;
	}

};
class GraphLoader
{
public:
	vector<GraphArc>		m_RawGraph; 
	int						m_nArcs;
public:	
	GraphLoader();
	GraphLoader(const char* szFileName);
	void LoadGraphFromFile(const char* szFileName);
};

#endif
