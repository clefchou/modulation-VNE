#include"Graph.h"
#include"Schedule.h"
#include"common.h"


void main(){
	//srand((unsigned)time(NULL));
	double lamda = 1;
	double mu = 1;
	//scanf("%lf", &lamda);
	
	
	//for(lamda = 10; lamda <= 100; lamda += 10)
	{
		char route[30] = "outputData_00.txt";
		route[11] = ((int)lamda)/10%10+'0';
		route[12] = ((int)lamda)%10+'0';
		//int nodeNum = 5;
		//int modNum = 4;
		//topoGenerate("inputData//whole_r.txt", nodeNum, nodeNum + 1, modNum, 10, 40,2000);

		Network *G = new Network("inputData\\whole_r_sixnodes.txt");
		Center ct(G, lamda, mu);
		
		//newDemandGenerate();
		
		ct.push(Event(0, 0, 0, new Network("inputData\\demand_r.txt")));
		while(!ct.empty()){
			ct.pop();
		}
		
		ct.show(route);
	}
}
