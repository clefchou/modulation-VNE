#include "Heuristic.h"
#include "UsTime.h"
#include "Graph.h"



static int nodeID = -1;
static vector<Nodes> nodes_temp;
bool sort_nodes(const Nodes &n1, const Nodes &n2){				//排序函数，用来对sort进行排序
	return n1.degree>n2.degree;
}


void Heuristic::sortnetwork(){
	
	bool sort_edges(const Edge &e1, const Edge &e2);

	for(int vertexID = 0; vertexID < demand->n; vertexID++){
		Nodes tempnode;
		tempnode.vertexID = vertexID, tempnode.degree = demand->degree[vertexID];
		nodes.push_back(tempnode);
	}
	
	sort(nodes.begin(),nodes.end(),sort_nodes);
	for(int k = 0; k < demand->n; k++){
		int edge_reserve = demand->edges.size();	//随着demand里edges被选入，vector大小会变小，所以判断条件要跟着变小
		vector<Edge>::iterator it = demand->edges.begin();
		for(int ij = 0 ; ij < demand->edges.size(); ij++){
			if(demand->edges[ij].head == nodes[k].vertexID || demand->edges[ij].tail == nodes[k].vertexID){
				edges_middle.push_back(demand->edges[ij]);
				demand->edges.erase(it);
				ij--;			//每次删掉内容后 ij对应的值也要回退一格，因为整个vector都被压缩了一格
			}
			else
				it++;
		}
		nodeID = k;
		nodes_temp = nodes;
		sort(edges_middle.begin(),edges_middle.end(),sort_edges);			//创建中间变量，对边的相邻节点度进行排序，再把它加进一总的vector里面
		for (int i = 0; i < edges_middle.size(); i++){
			edges_after.push_back(edges_middle[i]);
		}
		edges_middle.clear();
	}	
}
	//sort(nodes.begin(),nodes.end(),demand->degree(n));


 bool sort_edges(const Edge &e1, const Edge &e2){

	int degree_e1, degree_e2;
	if (e1.head == nodeID)
		degree_e1 = nodes_temp[e1.tail].degree;
	if (e1.tail == nodeID)
		degree_e1 = nodes_temp[e1.head].degree;
	if (e2.head == nodeID)
		degree_e2 = nodes_temp[e2.tail].degree;
	if (e2.tail == nodeID)
		degree_e2 = nodes_temp[e2.head].degree;
	return degree_e1>degree_e2;
}

Network Heuristic::auxGraph_generate(int demand_head, int demand_tail,int a,int g){
	//int mod = demand->mod;
	vector<int> r_in, r_mj, d_in, d_mj;
	int lanmda_v_i = 0; 
	for(int k = 0; k < resource->m; k++){
		r_in.push_back(100000);
		r_mj.push_back(100000);	
		d_in.push_back(0);
		d_mj.push_back(0);//公式17的r_in,r_jn 未完成
	}


//	int F_g_1 = resource->maxBandwidth - mod +1;
//	for (int a = 0; a < F_g_1; a++){
		{auxGraph.edges.clear();
		auxGraph.ck.clear();
		auxGraph.degree.clear();
		auxGraph.vertexWeight.clear();
		
		}
		 auxGraph = *resource ;



		int nshift = auxGraph.n, mshift = auxGraph.m;
		int n_number= auxGraph.n, m_number = auxGraph.m;
			for(int i = 0; i < n_number; i++){
				
				auxGraph.vertexWeight.push_back(auxGraph.vertexWeight[i]);		
				auxGraph.degree.push_back(auxGraph.degree[i]);
				auxGraph.ck.push_back(auxGraph.ck[i]);
				auxGraph.n += 1;
			}
			for(int ij = 0; ij < m_number; ij++){	
				vector<Edge>::iterator it = auxGraph.edges.begin();
				int if_cut_the_edge = 0;
				
				int dest_a = a + g;
				for( int ak = a; ak < dest_a; ak++){
					if(auxGraph.edges[ij].wave[ak] != true)
						if_cut_the_edge = 1;
				}
				if(if_cut_the_edge){
					auxGraph.edges.erase(it);
					ij--;
					m_number--;
					auxGraph.m--;
				}	
				else
					it++;
			}													//删除了容量不够的边，接下来把所有边翻倍形成新图
			for(int ij = 0; ij < m_number; ij++){
				Edge te = auxGraph.edges[ij];
				te.head += nshift;
				te.tail += nshift;
				te.cmn *= auxGraph.edges[ij].distance;	//更新cmn 乘上距离
				auxGraph.edges[ij].cmn = te.cmn;
				auxGraph.edges.push_back(te);
				auxGraph.maxBandwidth = max(auxGraph.maxBandwidth, te.bandwidth);
				auxGraph.m += 1;

				auxGraph.edges[ij].head += nshift;				//让初始的点连到新的图的新点上去

			}
		
			//接下来要根据 i,j是否被映射选择增加虚边的数量，如果被映射了只增加一条，没有的话就全加上。
			
			auxGraph.n +=2;	//加上了两个虚拟节点
			int FIk_i = -1, FIk_j = -1, if_node_embeded = 0;	//初始化映射flag 假定没有被映射到任何一点上
			for (int k = 0; k < resource->n; k++){
					if(FIk[demand_head][k] == 1)
						FIk_j = k;
					if(FIk[demand_tail][k] == 1)
						FIk_i = k;
			}
			if(FIk_i >= 0){
				Edge *te = new Edge(auxGraph.n-2,FIk_i,auxGraph.maxBandwidth ,r_in[FIk_i],0);	//为了让i,j编号不与实际节点k编号搞混，另j为最后一个节点，i为倒数第二个
				auxGraph.edges.push_back(*te);
				auxGraph.m +=1;
			}			//如果i点被映射过了，那就只加一条被映射的边↑；如果没有被映射过，那就把其他所有边都加上去↓
			else{
				for(int k = 0; k < resource->n; k++){
					for (int demand_k = 0; demand_k < demand->n; demand_k ++){
						if(FIk[demand_k][k] == 1 )
							if_node_embeded = 1;
					}
					int dmn,rin = 0;
					if(if_node_embeded == 0){
					Edge *te = new Edge(auxGraph.n-2,k,auxGraph.maxBandwidth ,r_in[k],0);
/************************************************************************
					int count=0;
					for (int edgeID = 0; edgeID < auxGraph.m; edgeID++){
						
			            if(auxGraph.edges[edgeID].head == k||auxGraph.edges[edgeID].tail == k){
						    count++;
							for(int waveID = 0; waveID+te->bandwidth <= resource->maxBandwidth; waveID++)
							rin ++;//怎么表达?
						}
					}
					if (count != 0)
						rin=rin/count;		   
					rin = pow(float(rin),resource->degree[k]);
					
					
					int Dis = 0,aa=0,mid=0; //aa记录deta(v,n)的个数
					for(int ij = 0 ; ij < edges_after.size(); ij++){
						
			           if(edges_after[ij].head == nodes[k].vertexID ){
						     for(int kk = 0; kk < resource->n; kk++){
								 if(FIk[edges_after[ij].tail][kk] == 1 ){
							         if_node_embeded = 1;
								     mid=kk;//用mid记录节点v映射的物理节点
								     break;
						         }
						      }	
							 if(if_node_embeded = 1){
								int dis = findpath(resource,k,mid);//调用一次find找到Dis(M,n),
								Dis+=dis;
								aa++;
							 }
					    }

			           if(edges_after[ij].tail == nodes[k].vertexID ){
						     for(int kk = 0; kk < resource->n; kk++){
						           if(FIk[edges_after[ij].head][kk] == 1 ){
							           if_node_embeded = 1;
								       mid=kk;//用mid记录节点v映射的物理节点
								       break;
						           	}
						     }
							if( if_node_embeded = 1){
								int dis = findpath(resource,k,mid);//调用一次find找到Dis(M,n),
								Dis+=dis;
								aa++;
							}
						}  							  	
					if( Dis != 0)
					 dmn=Dis/aa;	
                    int cmn= rin * dmn;
					te->cmn = cmn;

					}
					*///************************************************************************
					auxGraph.edges.push_back(*te);
					auxGraph.m +=1;
					}
					if_node_embeded = 0;
				}
			}	


			if(FIk_j >= 0){
				Edge *te = new Edge(FIk_j+resource->n,auxGraph.n-1,auxGraph.maxBandwidth ,r_mj[FIk_j],0);	//复制的图的终点，每个节点编号都要shift。
				auxGraph.edges.push_back(*te);
				auxGraph.m +=1;
				
			}			//如果i点被映射过了，那就只加一条被映射的边↑；如果没有被映射过，那就把其他所有边都加上去↓
			else{
				for(int k = 0; k < resource->n; k++){
					for (int demand_k = 0; demand_k < demand->n; demand_k ++){
						if(FIk[demand_k][k] == 1 )
							if_node_embeded = 1;
					}
				
					int dmn,rmj = 0;
					if(if_node_embeded == 0){
					Edge *te = new Edge(k+resource->n,auxGraph.n-1,auxGraph.maxBandwidth ,r_mj[k],0);
/************************************************************************
					int count=0;
					for (int edgeID = 0; edgeID < auxGraph.m; edgeID++){
						
			            if(auxGraph.edges[edgeID].head == k||auxGraph.edges[edgeID].tail == k){
						    count++;
							for(int waveID = 0; waveID+te->bandwidth <= resource->maxBandwidth; waveID++)
							rmj ++;//怎么表达?
						}
					}
					if (count != 0)
						rmj=rmj/count;		   
					rmj = pow(float(rmj),resource->degree[k]);
					
					
					int Dis = 0,aa=0,mid=0; //aa记录deta(v,n)的个数
					for(int ij = 0 ; ij < edges_after.size(); ij++){
						
			           if(edges_after[ij].head == nodes[k].vertexID ){
						     for(int kk = 0; kk < resource->n; kk++){
								 if(FIk[edges_after[ij].tail][kk] == 1 ){
							         if_node_embeded = 1;
								     mid=kk;//用mid记录节点v映射的物理节点
								     break;
						         }
						      }	
							 if(if_node_embeded = 1){
								int dis = findpath(resource,k,mid);//调用一次find找到Dis(M,n),
								Dis+=dis;
								aa++;
							 }
					    }

			           if(edges_after[ij].tail == nodes[k].vertexID ){
						     for(int kk = 0; kk < resource->n; kk++){
						           if(FIk[edges_after[ij].head][kk] == 1 ){
							           if_node_embeded = 1;
								       mid=kk;//用mid记录节点v映射的物理节点
								       break;
						           	}
						     }
							if( if_node_embeded = 1){
								int dis = findpath(resource,k,mid);//调用一次find找到Dis(M,n),
								Dis+=dis;
								aa++;
							}
						}  							  	
					if( Dis != 0)
					 dmn=Dis/aa;	
                    int cmn= rmj * dmn;
					te->cmn = cmn;

					}
*///**********************************************************
					auxGraph.edges.push_back(*te);
					auxGraph.m +=1;
					}
					if_node_embeded = 0;
				}
			}
		
		return auxGraph;
			
}

int Heuristic::findpath(){
	int temp_dist=100000000,temp_node_id ,findcounts = 0;
	vector<double>dist(auxGraph.n,1000000000);
	//vector<int> last_vertex_ID;							//记录上一个节点的ID

	last_vertex_ID.clear();
	for(int i = 0; i < auxGraph.n; i++)
		last_vertex_ID.push_back(-1);
	
	vector<int> used_vertex_ID;							//记录已经到达的点ID
	used_vertex_ID.push_back(auxGraph.n-2);				//auxGraph.n-1代表的就是源节点 auxGraph.n代表目的节点
	dist[auxGraph.n-2] = 0;						//初始化距离数组
	do{
	for (int i = 0; i < used_vertex_ID.size(); i++){
		for (int edgeID = 0; edgeID < auxGraph.m; edgeID++){
			if(auxGraph.edges[edgeID].tail == used_vertex_ID[i]  )
				if(dist[auxGraph.edges[edgeID].tail]+auxGraph.edges[edgeID].cmn < dist[auxGraph.edges[edgeID].head]){
					dist[auxGraph.edges[edgeID].head] = dist[auxGraph.edges[edgeID].tail]+auxGraph.edges[edgeID].cmn;
					last_vertex_ID[auxGraph.edges[edgeID].head] = auxGraph.edges[edgeID].tail;		//记录下上一个节点位置。
				}
		}
	}
	for (int i = 0; i < auxGraph.n; i++){
		int if_i_is_used = 0;
		for (int k = 0; k < used_vertex_ID.size(); k++){
			if( i == used_vertex_ID[k])
				if_i_is_used = 1;
		}
		if (dist[i] < temp_dist && if_i_is_used == 0){
			temp_dist = dist[i];
			temp_node_id = i;
		}
	}								//最终求得了一个节点使距离最小，把节点加入重新计算直到到达目的节点
	used_vertex_ID.push_back(temp_node_id);
	temp_dist=100000000;
	findcounts ++;
	if( findcounts >= auxGraph.m){
		return 0;						//	如果寻找次数超过边数了还没到达终点，则寻找失败返回0 ？？？？不知道判定是否正确
	}
	}while(temp_node_id != auxGraph.n-1);

	return dist[auxGraph.n-1];			//这个距离相当于cmn之和，物理距离要另算

}

int Heuristic::findpath(Network *r,int source_vertexID, int dest_vertexID){
	int temp_dist=100000000,temp_node_id ,findcounts = 0;
	Network *tempnet = r;
	vector<double>dist(tempnet->n,1000000000);
	//vector<int> last_vertex_ID;							//记录上一个节点的ID

	last_vertex_ID.clear();
	for(int i = 0; i < tempnet->n; i++)
		last_vertex_ID.push_back(-1);
	
	vector<int> used_vertex_ID;							//记录已经到达的点ID
	used_vertex_ID.push_back(source_vertexID);				
	dist[source_vertexID] = 0;						//初始化距离数组
	do{
	for (int i = 0; i < used_vertex_ID.size(); i++){
		for (int edgeID = 0; edgeID < tempnet->m; edgeID++){
			if(tempnet->edges[edgeID].tail == used_vertex_ID[i])
				if(dist[tempnet->edges[edgeID].tail]+tempnet->edges[edgeID].distance < dist[tempnet->edges[edgeID].head]){//这里改成了distance距离而不是cmn
					dist[tempnet->edges[edgeID].head] = dist[tempnet->edges[edgeID].tail]+tempnet->edges[edgeID].distance;
					last_vertex_ID[tempnet->edges[edgeID].head] = tempnet->edges[edgeID].tail;		//记录下上一个节点位置。
				}
		}
	}
	for (int i = 0; i < tempnet->n; i++){
		int if_i_is_used = 0;
		for (int k = 0; k < used_vertex_ID.size(); k++){
			if( i == used_vertex_ID[k])
				if_i_is_used = 1;
		}
		if (dist[i] < temp_dist && if_i_is_used == 0){
			temp_dist = dist[i];
			temp_node_id = i;
		}
	}								//最终求得了一个节点使距离最小，把节点加入重新计算直到到达目的节点
	used_vertex_ID.push_back(temp_node_id);
	temp_dist=100000000;
	findcounts ++;
	if( findcounts >= tempnet->m){
		return 0;						//	如果寻找次数超过边数了还没到达终点，则寻找失败返回0 ？？？？不知道判定是否正确
	}
	}while(temp_node_id != dest_vertexID);

	return dist[dest_vertexID];			

}




Heu_return Heuristic::work(){
	UsClock clk;
	clk.setStartTime();

	Event ret = Event(0, 0, 0, resource);
	Heu_return heu_ret;
	for(int ij = 0; ij < edges_after.size(); ij++){
		
		int e_has_been_embeded = 0;
		for(int k = 0; k < resource->n; k++){
			if(FIk[edges_after[ij].head][k] == 1 || FIk[edges_after[ij].tail][k] == 1)
				e_has_been_embeded++;
		}
				{vector<int> t;
		heu_ret.demand_tail.push_back(edges_after[ij].tail);
		heu_ret.demand_head.push_back(edges_after[ij].head);
		heu_ret.edge_head.push_back(t);
		heu_ret.edge_tail.push_back(t);
		heu_ret.total_distance.push_back(100000000);
		heu_ret.grid_a.push_back(-1);	
		heu_ret.physical_distance.push_back(0);				//初始化结果类
		heu_ret.bandwidth.push_back(-1);
		heu_ret.modtype.push_back(-1);
		}

		//if(e_has_been_embeded == 2)							//这一段判断是否已经被映射，如果两个点都被映射过了，则这条边已被映射，continue		
			//continue;
		

		int flag_FIk_same_point = 0 , same_point_ID;		//判断是否映射到同一点
		int DijisSuccess;
		for (int mod = 16; mod > 1; mod = mod/2){ 
			int g = ceil(double(edges_after[ij].bandwidth)/(mod));
			int F_g_1 = resource->maxBandwidth - g + 1;	//这里擅自取了F=maxBandwidth 循环多了应该也不影响，能找到会在之前就找到，找不到循环到底也找不到
			for (int a = 1; a < F_g_1; a++){
				heu_ret.grid_a[ij] = a;

				auxGraph_generate(edges_after[ij].head,edges_after[ij].tail,a,g);	//生成辅助图

				if(flag_FIk_same_point == 1){
					flag_FIk_same_point = 0;
					int m_number = auxGraph.m;
					for(int mn = 0; mn < m_number; mn++){	
						vector<Edge>::iterator it = auxGraph.edges.begin();
						if(mn == same_point_ID){
							auxGraph.edges.erase(it);
							mn--;
						}	
						else
							it++;
					}				
				}							//删除其中一条虚边，保证不会映射到同一个点

				 DijisSuccess = findpath();		//根据生成辅助图寻找最短路
				if (!DijisSuccess)
					continue;
				else{
					if(DijisSuccess < heu_ret.total_distance[ij]){
						heu_ret.total_distance[ij] = DijisSuccess;
						heu_ret.edge_head[ij].clear();
						heu_ret.edge_tail[ij].clear();			//结果更小，删除上一次保留的结果

							int j_m ,i_n ,temp = auxGraph.n-1;
							j_m = last_vertex_ID[auxGraph.n-1];
							heu_ret.edge_head[ij].push_back(auxGraph.n-1);
							heu_ret.edge_tail[ij].push_back(j_m);
					
							do{
								i_n = temp;
								temp = last_vertex_ID[temp];
								heu_ret.edge_head[ij].push_back(temp);
								heu_ret.edge_tail[ij].push_back(last_vertex_ID[temp]);


							}while(last_vertex_ID[temp] != auxGraph.n-2);
							i_n = temp;

							if(i_n + resource->n == j_m){			//两个点I,J映射到同一个物理节点的情况
								flag_FIk_same_point = 1;
								same_point_ID = j_m;
								a--;
								continue;
							}

							for (int demand_k = 0; demand_k<resource->n;demand_k++){
								if(demand_k == j_m || demand_k == j_m-resource->n)
								FIk[edges_after[ij].head][demand_k]=1;
							}
							for (int demand_k = 0; demand_k<resource->n;demand_k++){
								if(demand_k == i_n || demand_k == i_n-resource->n)
								FIk[edges_after[ij].tail][demand_k]=1;
							}														//将节点映射记录下来
						}
					break;
					}
			}
			int physical_distance = 0;
			if(DijisSuccess){				//如果找到了一条路径，break 并将资源情况更新
				for(int mn = 0; mn < resource->m; mn++){
					for(int k = 0; k < heu_ret.edge_head[ij].size(); k++){
						if(resource->edges[mn].head == heu_ret.edge_head[ij][k] ||resource->edges[mn].head == heu_ret.edge_head[ij][k]-resource->n)
							if(resource->edges[mn].tail == heu_ret.edge_tail[ij][k] ||resource->edges[mn].tail == heu_ret.edge_tail[ij][k]-resource->n){
								physical_distance += resource->edges[mn].distance;
								for(int a = heu_ret.grid_a[ij]; a < heu_ret.grid_a[ij]+g; a++){
									resource->edges[mn].wave[a] = false;
								}
							}
					}					
				}
				int mod_num = 0-2,tempmod = mod;
				do{tempmod = tempmod/ 2;
					mod_num++ ;
				}while(tempmod != 0);
				if( physical_distance >= resource->moddistance[mod_num])
					continue;
				else{
					heu_ret.physical_distance[ij] = physical_distance;
					heu_ret.modtype[ij] = mod;
					heu_ret.bandwidth[ij] = g;
					printf("******total distance = %d\n",heu_ret.total_distance[ij]);
					break;
				}
			}
			if(mod == 2){
				heu_ret.isSuccess = 0;
				return heu_ret;
			}
		}
	}
	clk.setEndTime();
	heu_ret.proccessTime = clk.getTime();
	heu_ret.isSuccess = 1;
	return heu_ret;
}




Heuristic::Heuristic(Network *d, Network *r, FILE *lf){
	de = *d, re = *r ,ori = r;
	resource = new Network (*r);
	demand = new Network (*d);
	resource->mod++;
	logFile = lf;
	result.isSuccess = result.proccessCircle = result.proccessTime = result.proccessTreeNode = result.resultWeight = 0;
	
	for(int nodei=0;nodei<demand->n;nodei++){
		vector<int> tempFIk;
		FIk.push_back(tempFIk);
	}
	for(int nodei=0;nodei<demand->n;nodei++){
		for(int k=0;k<resource->n;k++){
			FIk[nodei].push_back(0);
		}
	}
	
	sortnetwork();
	for (int i = 0; i < resource->n; i++){
		vector<int> tempneighbor;
		resource->neighbor.push_back(tempneighbor);
	}

/*	for (int i = 0; i < resource->n; i++)				//定义每个节点的邻居节点
		for( int k = 0; k < resource->n; k++)
			resource->neighbor[i].push_back(0);*/
	for(int k = 0; k < resource->m; k++){
		int head = resource->edges[k].head;
		int tail = resource->edges[k].tail;
		resource->neighbor[head].push_back(tail);
		resource->neighbor[tail].push_back(head);
	}
	for(int i = 0; i < resource->n; i++){
		sort( resource->neighbor[i].begin(),resource->neighbor[i].end());
		resource->neighbor[i].erase( unique(resource->neighbor[i].begin(),resource->neighbor[i].end()),resource->neighbor[i].end());
		
	}
}

Heuristic::~Heuristic(){

	
	delete resource;
	
//	fclose(optimalData);
}

Event Heuristic::solveByHeuristic(Network *d, Network *r, FILE *logFile){
	Event ret = Event(0, 0, 0, r);
	int numILPSuccess, numILPBlock, numHeuSuccess;
		double totalILPTime, totalILPWeight, totalHeuTime, totalHeuWeight, HeuWeight;
		Network original_resource = *r;
		Network *resource_Hr = &original_resource;
			 Heu_return heu_ret;
			 heu_ret = work();
			 HeuWeight = 0;
			 if(heu_ret.isSuccess){
				 int n = d->edges.size();
				for(int ij = 0; ij < n; ij++){
					for (int mn = 0; mn < resource_Hr->edges.size(); mn++){
						for (int k = 0; k < heu_ret.edge_head[ij].size(); k++){
							if(resource_Hr->edges[mn].head == heu_ret.edge_head[ij][k] || resource_Hr->edges[mn].head == heu_ret.edge_head[ij][k]-resource_Hr->n)
								if(resource_Hr->edges[mn].tail == heu_ret.edge_tail[ij][k] || resource_Hr->edges[mn].tail == heu_ret.edge_tail[ij][k]-resource_Hr->n)
									for(int ak = 0; ak < heu_ret.bandwidth[ij]; ak++){
										HeuWeight += resource_Hr->edges[mn].cmn*heu_ret.bandwidth[ij];
										int wave_ak = heu_ret.grid_a[ij] + ak;
										resource_Hr->edges[mn].wave[wave_ak] = false;
									}
						}
					}
				}
			 }
			 

			 if(heu_ret.isSuccess){
				numHeuSuccess += 1;
				totalILPTime += heu_ret.proccessTime;
				totalILPWeight += HeuWeight;
				fprintf(logFile, "Heu success, get object of %f, using Time %f\n ", HeuWeight, heu_ret.proccessTime);
				ret.proccessTime = heu_ret.proccessTime;
				for(int ij = 0; ij < heu_ret.demand_tail.size(); ij++){
					for(int k = 0; k < heu_ret.edge_head[ij].size(); k++){
						if(heu_ret.edge_head[ij][k] >= resource_Hr->n )
							heu_ret.edge_head[ij][k] -= resource_Hr->n;
						if(heu_ret.edge_tail[ij][k] >= resource_Hr->n )
							heu_ret.edge_tail[ij][k] -= resource_Hr->n;
					}
				}


				for(int ij = 0; ij < heu_ret.demand_tail.size(); ij++){
					printf("Virtual node: %d, %d ==> Physical node: %d, %d\n",heu_ret.demand_tail[ij],heu_ret.demand_head[ij],heu_ret.edge_head[ij][heu_ret.edge_head[ij].size()-1],heu_ret.edge_tail[ij][0]);


				}
			}
			 else
					fprintf(logFile, "Heu failed\n");
			 ret.resultWeight = HeuWeight;
			for(int ij = 0; ij < ret.demand->edges.size(); ij++)
				for(int b = 0; b < ret.demand->edges[ij].wave.size(); b++)
					ret.demand->edges[ij].wave[b] = false;

			for(int k = 0; k < heu_ret.modtype.size(); k++){
				if(heu_ret.modtype[k] == 2)
					heu_ret.modtype[k] = 0;
				if(heu_ret.modtype[k] == 4)
					heu_ret.modtype[k] = 1;
				if(heu_ret.modtype[k] == 8)
					heu_ret.modtype[k] = 2;
				if(heu_ret.modtype[k] == 16)
					heu_ret.modtype[k] = 3;
			}
			for(int i = 0; i < 10; i++)
		      for(int j = 0; j < 20; j++)
				  for (int a =0; a < 30; a++)
					  for(int mod =0; mod < 4; mod++){
						  ret.b[i][j][a][mod] = 0;
					  }
			 for(int i = 0; i < 10; i++)
		      for(int j = 0; j < 10; j++)
			 {
			      ret.a[i][j] = 0;
			}

		for(int mn = 0; mn < r->m; mn++){
			for(int ij = 0; ij < d->m; ij++){
				for( int k = 0; k < heu_ret.edge_head[ij].size(); k++){
					if(r->edges[mn].head == heu_ret.edge_head[ij][k] )
						if(r->edges[mn].tail == heu_ret.edge_tail[ij][k]){
							int a_start = heu_ret.grid_a[ij];
							for(int a = heu_ret.grid_a[ij]; a < heu_ret.grid_a[ij] + heu_ret.bandwidth[ij]; a++){
								ret.b[ij][mn][a_start][heu_ret.modtype[ij]] = 1;
							    ret.demand->edges[mn].wave[a] = true;
							}
						}
				}
			}
		}
		for (int i = 0; i < d->n ; i ++){
			for(int k = 0; k < r->n; k++){
				if(FIk[i][k] == 1)
					ret.a[i][k] = 1;
			}
		}
		
		for(int k = 0; k < r->n; k++)
			ret.demand->vertexWeight[k] = 0;
		for(int k = 0; k < r->n; k++)
			for(int i = 0; i < d->n; i++)
				if(ret.a[i][k] == 1)
					ret.demand->vertexWeight[k] += d->vertexWeight[i];
		
		ret.isSuccess = heu_ret.isSuccess;
		return ret;

					  
					  
					  


}