#include "common.h"
#include "Graph.h"


class Nodes{
public:
	int vertexID;
	int degree;
	
};

class Heu_return{				//用来存储结果，不同的虚拟边ij对应vector中不同的元素
public:
	int isSuccess;
	vector<int> grid_a;				//ij的a的开始位置
	vector<int> demand_tail;			//ij的头节点i在虚拟网络中的编号
	vector<int> demand_head;			//ij的尾节点j在虚拟网络中的编号
	vector<int> total_distance;		//总的代价
	vector<int> physical_distance;	//总的物理距离，用来判断是否超过mod的最大距离
	vector<vector<int> > edge_head;	//映射的所有物理边的头节点编号
	vector<vector<int> > edge_tail;	//映射的所有物理边的尾节点编号
	vector<int> bandwidth;
	vector<int> modtype;
	double proccessTime, resultWeight; 
};


class Heuristic{
private:
	Network de, re, *ori , *aux , *demand, *resource;
	Network auxGraph;
	Event result;
	vector<Edge> edges_before,edges_after, edges_middle;
	vector<Nodes> nodes;
	vector<vector<int>> FIk;//记录节点映射
	vector<int> last_vertex_ID;		//记录一次返回dijkstra 返回的网络，记得清零。
		

	static const int inf = 10000000;
	FILE *logFile, *optimalData;
	
	
	double solve();
	void sortnetwork();
	Network auxGraph_generate(int demand_head, int demand_tail, int a, int g, int maxdistance);
	int findpath();
	int findpath(Network *r,int source_vertexID, int dest_vertexID);
	
public:
	
	Heuristic();
	Heuristic(Network *d,Network *r,FILE *lf);
	Event solveByHeuristic(Network *d, Network *r, FILE *logFile);
	Heu_return work();
	~Heuristic();


};


