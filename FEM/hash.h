#include <App.h>
#include <list>
#include <vector>
struct node
{
	int objID;
	int index;
	glm::vec3 position;
	node(int id, int i, glm::vec3 v)
	{
		objID = id;
		index = i;
		position = v;
	}
};

class HashCell
{
public:
	double T = 0;
	list<node> nodes;
};

class HashMap
{
public:
	double T = 0;
	vector<HashCell> cell;
};
