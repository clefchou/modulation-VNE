#ifndef GraphH
#define GraphH

#include"common.h"

class Edge{
    public:
        int head, tail, bandwidth, distance;
		double cmn;
        int headU, tailU;
        vector<bool> wave;
		Edge() {;}
        Edge(int a, int b, int c, double d ,int e){
            head = b;
            tail = a;
			bandwidth = c;
            cmn = d;
            distance=e;
            for(int i = 0; i < c; i++)
                wave.push_back(true);
        }
};

class Network{
    public:
        vector<Edge> edges;
        vector<int> vertexWeight;
		vector<double> ck;
		vector<int> degree;
		vector<int> moddistance;
		vector<vector<int> >neighbor;
        int n, m, mod, maxBandwidth;
		Network() {;}
        Network(char f[]){
            FILE *input = fopen(f, "r");
            int tn, tm, tmod;
            fscanf(input, "%d%d%d", &tn, &tm, &tmod);
            n = tn, m = tm, mod = tmod, maxBandwidth = 0;
            for(int i = 0; i < tn; i++){
                int t;
                fscanf(input, "%d", &t);
                vertexWeight.push_back(t);
				degree.push_back(0);
            }
			for(int i = 0; i < tn; i++){
				double t;
				fscanf(input, "%lf", &t);
				ck.push_back(t);
			}
			for(int i =0; i < tmod; i++){
				int t;
				fscanf(input, "%d", &t);
				moddistance.push_back(t);
			}
            for(int i = 0; i < tm; i++){
                int a, b, c, e; 
				double d;
                fscanf(input, "%d%d%d%lf%d", &a, &b, &c, &d, &e);
                edges.push_back(Edge(a, b, c, d, e));
				maxBandwidth = max(maxBandwidth, c);
				degree[a]++;
				degree[b]++;
            }
            fclose(input);
        }
		void show(){
			for(int i = 0; i < m; i++){
				printf("Edge_%d: ", i);
				for(int j = 0; j < edges[i].bandwidth; j++)
					printf("%d ", edges[i].wave[j]);
				printf("\n");
			}
		}
};

class Event{
    public:
        bool type;
        int arrivalTime, ID;
		int isSuccess;
		double proccessTime, resultWeight; 
		int proccessTreeNode, proccessCircle;
		//vector<vector<int>> a;
		int a[10][10];
		int b[10][20][30][4];
		//int c[10][20][30];
		//int d[10][30];
		//int e[10][4];
        Network *demand;
        Event() {
			demand = NULL;
		}
        Event(bool t, int curTime, int id, Network *d){
            type = t;
            arrivalTime = curTime;
			ID = id;
			demand = new Network(*d);
			proccessTime = resultWeight = proccessTreeNode = proccessCircle = isSuccess = 0;
			//memset(a,0,sizof(a));
			/*for(int i = 0; i < 6; i++)
			      for(int j = 0; j < 6; j++){
				    a[i][j] = 0;
			}*/
			/*a.resize(demand -> n);
			for(int i = 0; i < demand-> n; i++)
			    a[i].push_back(0);*/
			
        }
        bool operator<(const Event &b) const{
            return this->arrivalTime > b.arrivalTime;
		}
};

#endif
