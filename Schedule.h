#include"Graph.h"
#include"ILPSolver.h"
#include"common.h"
#include "Heuristic.h"

#define WAIT 5
#define DEMANDNUM 30

int uniform(int a, int b){
	return a + rand()%(b-a+1);
}
void topoGenerate(char route[], int n, int mod, int maxWeight,int maxdistance){
/*
	生成随机拓扑
	route为拓扑文件存储路径
	n为点数，m为边数
	maxCapacity为带宽及节点计算能力
	maxWeight为点权边权最大值，具体值为[0,maxWeight]均匀分布
*/
	FILE *out = fopen(route, "w");
	int topo[4][4];
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			topo[i][j] = 0;
	int m;
	switch(n)
	{
	case 2:
		m = uniform(1, 2);
	    break;
	case 3:
		m = uniform(2, 5);
		break;
	case 4:
		m = uniform(3, 6);
		break;
	}
	fprintf(out, "%d %d %d\n", n, m, mod);
	for(int i = 0; i < n; i++)//每个节点的计算需求
		fprintf(out, "%d ", uniform(1, 4));
	fprintf(out, "\n");
	for(int i = 0; i < n; i++)
		fprintf(out, "%d ", maxWeight);
	fprintf(out, "\n");
	for(int i = 0; i < mod; i++)
		fprintf(out, "%d ", maxdistance);//只有物理图需要moddistance且为固定值
	fprintf(out, "\n");
	for(int i = 1; i < n; i++){
		int randomi = rand()%i;
		//fprintf(out, "%d %d %d %d %d\n", i, randomi, maxCapacity, maxWeight, maxdistance);
		topo[i][randomi] = 1;
	}
	for(int i = 0; i < m - n + 1; i++){
		int u = rand()%n, v ,ucount = 0;
		int uvalready = 0;
		do{
			uvalready = 0;
			v = rand()%n;
			ucount++;
			if(topo[u][v] == 1)
				uvalready = 1;
			if(ucount > n){		//这里容易产生死循环，如果v = 1而1的0,2,3全都满了，就会一直在里面转
				u++;
				break;
			}
		}while(u == v || uvalready);
	    topo[u][v] = 1;
	}
	for(int i = 0 ; i < 4; i++) 
		for(int j = 0; j < 4; j++)
			if(topo[i][j] == 1)
		       fprintf(out, "%d %d %d %d %d\n", i, j, uniform(1, 10), maxWeight, maxdistance);
	fclose(out);
	//system("type D:\\毕设\\xuni\\inputData\\demand_r.txt>>D:\\毕设\\xuni\\inputData\\alldemand.txt");
}

double U_Random()   /* 产生一个0~1之间的随机数 */
{
	double f;
	//srand( (unsigned)time( NULL ) );
	f = (float)(rand() % 100);
	/* printf("%f\n",f); */
	return f/100;
}
int possion(double lamda)  /* 产生一个泊松分布的随机数，Lamda为总体平均数*/
{
    int k = 0;
    long double p = 1.0;
    long double l=exp(-lamda);  /* 为了精度，才定义为long double的，exp(-Lambda)是接近0的小数*/
    //printf("%.15Lf\n",l);
    while (p>=l)
    {
        double u = U_Random();
        p *= u;
        k++;
    }
    return (k-1)/2;
}
double Exponential( double Beta)
{
    double Random_Number;
    double Return_Value;
    
    Random_Number = U_Random();
    while(Random_Number<0.0000001)
    Random_Number = U_Random();
    Return_Value = -1.0 * Beta * log( Random_Number);
    return Return_Value ;

}

void newDemandGenerate(){
	int vn = uniform(4, 4);
	topoGenerate("inputData\\demand_r.txt", vn, 4, 1, 300);
}

class Center{
    private:
		int eventCount;
		int numILPSuccess, numILPBlock, numHeuSuccess;
		double totalILPTime, totalILPWeight, totalHeuTime, totalHeuWeight, HeuWeight;
		

		Network *resource;
		Network resource_original;
		Network *resource_Hr;
		double lamda, mu;
        priority_queue<Event> schedule;
		FILE *logFile;

        void arrange(Event cur){
			fprintf(logFile, "Event No.%d Eventtype:%d process\n", cur.ID, cur.type);
			//if(eventCount == 4 && cur.type == 0){
			    /*for(int k = 0; k < resource->n; k++)
			        cout<<resource->vertexWeight[k];
				cout<<"\n";
				int m = resource->edges.size();
			    for(int ij = 0; ij < m; ij++){
				        int bandwidth = resource->edges[ij].wave.size();
				        for(int b = 0; b < bandwidth; b++)
					        cout<<resource->edges[ij].wave[b];
						cout<<"\n";
			}*///用于检查
			int checkError;
			//**********************************************************************

			Heuristic Hr(cur.demand, resource, logFile);


			//**********************************************************************
			//利用ILP求解并检查
			 printf("***************************");

			 Event retILP = Hr.solveByHeuristic(cur.demand, resource, logFile);
			//Event retILP = solveBySimpleILP(cur.demand, resource, logFile);
			if(retILP.isSuccess){
				numILPSuccess += 1;
				totalILPTime += retILP.proccessTime;
				totalILPWeight += retILP.resultWeight;
			}
			else
				numILPBlock += 1;
			checkError = 0;
			if(retILP.isSuccess){
				int m = resource->edges.size();
				for(int ij = 0; ij < m; ij++){
					int bandwidth = resource->edges[ij].wave.size();
					for(int b = 0; b < bandwidth; b++)
						if(retILP.demand->edges[ij].wave[b] > resource->edges[ij].wave[b])
							checkError = 1;
				}
			}
			if(checkError)
				fprintf(logFile, "error, ILPSolver use unavaiable wave\n");
			else
				fprintf(logFile, "ILPSolver correct after check\n");
			/*if(retILP.isSuccess)
				fprintf(logFile, "ILP success, get object of %f, using Time %f\n ", retILP.resultWeight, retILP.proccessTime);*/
				if(retILP.isSuccess){
					//if(eventCount == 4){
				    fprintf(logFile, "ILP success, get object of %f, using Time %f\n ", retILP.resultWeight, retILP.proccessTime);
				    for(int i = 0; i < cur.demand->n; i++)
					   for(int k = 0; k < resource->n; k++)
					       if(retILP.a[i][k] == 1){
								fprintf(logFile, "A: virtualNod; %d,--> phsicalNod: %d\n ", i, k);}
					for(int i = 0; i < cur.demand->m; i++)
		                for(int j = 0; j < resource->m ; j++)
				            for (int a =0; a < resource->edges[j].bandwidth - cur.demand->edges[i].bandwidth + 1; a++)
					            for(int mod =0; mod < resource->mod; mod++)
									if(retILP.b[i][j][a][mod] == 1){
										fprintf(logFile, "Z: virtualedge; %d,--> phsicaledge: %d, startgrid: %d, modulation: %d\n", i, j, a, mod);}
					/*for(int i = 0; i < cur.demand->m; i++)
		                for(int j = 0; j < resource->m ; j++)
				            for (int a =0; a < resource->edges[j].bandwidth - cur.demand->edges[i].bandwidth + 1; a++)
								if(retILP.c[i][j][a] == 1){
									fprintf(logFile, "M: virtualedge; %d,--> phsicaledge: %d, startgrid: %d\n ", i, j, a);
							}
					for(int i = 0; i < cur.demand->m; i++)
			            for (int a =0; a < resource->maxBandwidth - cur.demand->edges[i].bandwidth + 1; a++)
							if(retILP.d[i][a] == 1){
								fprintf(logFile, "T: virtualedge; %d,--> startgrid: %d\n ", i, a);
							}
					 for(int i = 0; i < cur.demand->m; i++)
		                 for(int mod =0; mod < resource->mod; mod++)
							 if(retILP.e[i][mod] == 1){
								 fprintf(logFile, "S: virtualedge; %d,--> modulation: %d\n ", i, mod);
							 }*/

				}
				else
					fprintf(logFile, "ILP failed\n");
				

			int thisTime = cur.arrivalTime;
			if(retILP.isSuccess){
				//分配资源
				int m = resource->edges.size();
				for(int ij = 0; ij < m; ij++){
					int bandwidth = resource->edges[ij].wave.size();
					for(int b = 0; b < bandwidth; b++)
						if(retILP.demand->edges[ij].wave[b] == 1)
							resource->edges[ij].wave[b] = 0;
				}
				for(int k = 0; k < resource->n; k++)
				    resource->vertexWeight[k] -= retILP.demand->vertexWeight[k];
				//if(eventCount == 4)
				
			       /* for(int k = 0; k < resource->n; k++)
			            cout<<resource->vertexWeight[k];
					cout<<"\n";*/
			       /* for(int ij = 0; ij < m; ij++){
				        int bandwidth = resource->edges[ij].wave.size();
				        for(int b = 0; b < bandwidth; b++)
					        cout<<resource->edges[ij].wave[b];
						cout<<"\n";
			  
				}*///用于检查
				//增加离去事件
				*(cur.demand) = *(retILP.demand);
				cur.arrivalTime +=  Exponential(mu);
				cur.type = 1;
				push(cur);
			}
			delete retILP.demand;


			

			//下一个到达事件
			if(eventCount < DEMANDNUM){
				newDemandGenerate();
				push(Event(0, thisTime + 1 /*Exponential(1/lamda)*/, eventCount, new Network("inputData\\demand_r.txt")));
			}
        }
		//释放资源
        void release(Event cur){
			fprintf(logFile, "Event No.%d Eventtype:%d process\n", cur.ID, cur.type);
			int m = resource->edges.size();
			for(int ij = 0; ij < m; ij++){
				int bandwidth = resource->edges[ij].wave.size();
				for(int b = 0; b < bandwidth; b++)
					if(cur.demand->edges[ij].wave[b] == 1)
						resource->edges[ij].wave[b] = 1;
			}
			for(int k = 0; k < resource->n; k++)
			    resource->vertexWeight[k] += cur.demand->vertexWeight[k];
			delete cur.demand;
        }
    public:
        Center() {
			logFile = fopen("outputData\\logFile.txt", "w");
		}
        Center(Network *t, double lam, double muu){
            resource = t;
			

			lamda = lam;
			mu = muu;
			eventCount = 0;
			numILPSuccess = numILPBlock = 0;
			totalILPTime = totalILPWeight = 0;
			logFile = fopen("outputData\\logFile.txt", "w");
        }

        void push(Event cur){
            schedule.push(cur);
			if(cur.type == 0)
			    eventCount++;
		}
        void pop(){
            Event cur = schedule.top();
            schedule.pop();

			if(!cur.type){
                arrange(cur);
				show("outputData\\processInfo.txt", cur.ID);
			}
			else{
				resource_original = *resource;
                release(cur);								//释放资源影响到resource变量的么，resource明明不是cur的变量啊。
				*resource = resource_original;
			}
        }
		bool empty(){
			return schedule.empty();
		}
		void show(char *outputFile, int curID = -1){
			FILE *info = fopen(outputFile, "w");
			if(curID != -1)
				fprintf(info, "lamda %.2f: demand NO.%d proccess\n\n", lamda, curID);
			else
				fprintf(info, "lamda %.2f:", lamda);
			fprintf(info, "ILP:\n\tblock: %d\n\tsuccess %d\n\ttime %f\n\tweight %f\n\n", numILPBlock, numILPSuccess, totalILPTime, totalILPWeight);
			fclose(info);
		}

		~Center(){
			fclose(logFile);
		}
};
