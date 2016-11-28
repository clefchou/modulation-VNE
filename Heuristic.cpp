#include "Heuristic.h"
#include "UsTime.h"
#include "Graph.h"

#define Ucmn 10	//控制rin rmj在0.5-1.5
#define Cr 0.05	//0.5当前结果最好
static int nodeID = -1;
static vector<Nodes> nodes_temp;

/*排序函数，用来对节点的度进行排序*/
bool sort_nodes(const Nodes &n1, const Nodes &n2){				
	return n1.degree>n2.degree;
}
/*排序函数，用来对边的度进行排序*/
bool sort_edges(const Edge &e1, const Edge &e2) {

	int degree_e1, degree_e2;
	if (e1.head == nodeID) {
		for(int i = 0; i < nodes_temp.size();i++)
			if(nodes_temp[i].vertexID==e1.tail)
				degree_e1 = nodes_temp[i].degree;		
	}
	if (e1.tail == nodeID)
		for (int i = 0; i < nodes_temp.size(); i++)
			if (nodes_temp[i].vertexID == e1.head)
				degree_e1 = nodes_temp[i].degree;
		
	if (e2.head == nodeID)
		for (int i = 0; i < nodes_temp.size(); i++)
			if (nodes_temp[i].vertexID == e2.tail)
				degree_e2 = nodes_temp[i].degree;		
	if (e2.tail == nodeID)
		for (int i = 0; i < nodes_temp.size(); i++)
			if (nodes_temp[i].vertexID == e2.head)
				degree_e2 = nodes_temp[i].degree;		

	return degree_e1>degree_e2;
}
/*排序函数，用来对带宽进行排序*/
bool sort_edges_bandwidth(const Edge &e1, const Edge &e2){
	
	return e1.bandwidth>e2.bandwidth;
}


void Heuristic::sortnetwork(int type){
	for(int vertexID = 0; vertexID < demand->n; vertexID++){
		Nodes tempnode;
		tempnode.vertexID = vertexID, tempnode.degree = demand->degree[vertexID];
		nodes.push_back(tempnode);
	}	
	sort(nodes.begin(),nodes.end(),sort_nodes);
	/* type == 0时对带宽进行排序,带宽占用大的边排在前
	 * type == 1时对度进行排序,先按节点的度对节点排序,再按节点的边的另一个点的度进行排序
	 * type == 2时先对节点的度进行排序,再根据节点边的带宽进行排序
	 */
	if(type == 0){//根据带宽排序
		int edge_reserve = demand->edges.size();	//随着demand里edges被选入，vector大小会变小，所以判断条件要跟着变小
		vector<Edge>::iterator it = demand->edges.begin();
		for(int ij = 0 ; ij < demand->edges.size(); ij++){
			if(demand->edges.size()){
				edges_middle.push_back(demand->edges[ij]);
				demand->edges.erase(it);
				ij--;			//每次删掉内容后 ij对应的值也要回退一格，因为整个vector都被压缩了一格
			}
			else
				it++;
		}
		sort(edges_middle.begin(),edges_middle.end(),sort_edges_bandwidth);			//创建中间变量，对边的占用带宽进行排序，再把它加进一总的vector里面
		for (int i = 0; i < edges_middle.size(); i++){
			edges_after.push_back(edges_middle[i]);
		}
	}
	if(type){
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
		if(type == 1)//按照度排序
			sort(edges_middle.begin(),edges_middle.end(),sort_edges);
		if(type == 2)//先按照度排序，再按照带宽排序
			sort(edges_middle.begin(),edges_middle.end(),sort_edges_bandwidth);			//创建中间变量，对边的相邻节点度进行排序，再把它加进一总的vector里面
		for (int i = 0; i < edges_middle.size(); i++){
			edges_after.push_back(edges_middle[i]);
		}
		edges_middle.clear();		
	}
	}
}


/* 生成辅助图的函数*/
Network Heuristic::auxGraph_generate(int demand_head, int demand_tail,int a,int g,int maxdistance){
	vector<int> r_in, r_mj;
	int lanmda_v_i = 0; 
	for(int k = 0; k < resource->m; k++){
		r_in.push_back(10000000);
		r_mj.push_back(10000000);	
	}
		auxGraph.edges.clear();
		auxGraph.ck.clear();
		auxGraph.degree.clear();
		auxGraph.vertexWeight.clear();
		/*首先将原图的全部信息复制到辅助图中*/
		auxGraph = *resource ;
		int nshift = auxGraph.n, mshift = auxGraph.m;
		int n_number= auxGraph.n, m_number = auxGraph.m;
			for(int i = 0; i < n_number; i++){				
				auxGraph.vertexWeight.push_back(auxGraph.vertexWeight[i]);
				auxGraph.degree.push_back(auxGraph.degree[i]);
				auxGraph.ck.push_back(auxGraph.ck[i]);
				auxGraph.n += 1;
			}
		/*删除容量不足的边,如果一段边的物理距离已经超过了mod限制,则直接将其从辅助图中删除*/
			vector<Edge>::iterator it = auxGraph.edges.begin();
			for(int ij = 0; ij < m_number; ij++){					
				int if_cut_the_edge = 0;				
				int dest_a = a + g;
				if(auxGraph.edges[ij].distance > maxdistance)			//判断边的物理距离是否足够；
					if_cut_the_edge = 1;
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
			}													
			/*接下来把所有边翻倍形成新图,生成和原图一样的一组点和边,将序号更新并将原图的点和边连接到新图的点和边上*/
			m_number = auxGraph.edges.size();
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
		
			/*接下来要根据 i,j是否被映射选择增加虚边的数量，如果被映射了只增加一条，没有的话就全加上。*/

			auxGraph.n +=2;	//加上了两个虚拟节点
			int FIk_i = -1, FIk_j = -1, if_node_embeded = 0;	//初始化映射flag 假定没有被映射到任何一点上
			for (int k = 0; k < resource->n; k++){
					if(FIk[demand_head][k] == 1)
						FIk_j = k;
					if(FIk[demand_tail][k] == 1)
						FIk_i = k;
			}
			if(FIk_i >= 0){
				if(resource->vertexWeight[FIk_i]>demand->vertexWeight[demand_tail]){
				Edge *te = new Edge(auxGraph.n-2,FIk_i,auxGraph.maxBandwidth ,r_in[FIk_i],0);	//为了让i,j编号不与实际节点k编号搞混，另j为最后一个节点，i为倒数第二个
				auxGraph.edges.push_back(*te);
				auxGraph.m +=1;
				delete te;
				}
				
			}			//如果i点被映射过了，那就只加一条被映射的边↑；如果没有被映射过，那就把其他所有边都加上去↓
			else{
				for(int k = 0; k < resource->n; k++){
					for (int demand_k = 0; demand_k < demand->n; demand_k ++){
						if(FIk[demand_k][k] == 1 )
							if_node_embeded = 1;
					}
					int dmn = 0;
					double rin = 0;
					if(if_node_embeded == 0){
						if(resource->vertexWeight[k]>demand->vertexWeight[demand_tail]){
					Edge *te = new Edge(auxGraph.n-2,k,auxGraph.maxBandwidth ,r_in[k],0);	//添加新的虚边
/************************************************************************/
					int count=0;
					/*计算∑ae/|Eo(n)|*/
					for (int edgeID = 0; edgeID < resource->m; edgeID++){	
			            if(resource->edges[edgeID].tail == k){
						    count++;
							for(int waveID = 0; waveID < resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									rin ++;
							
						}
					}
					if (count != 0)
						rin=rin/g/count;

					/*p是物理资源上有的容量数,r是虚边需要的容量数*/
					int pi_in,pi_out,ri_in,ri_out ,new_rin;
					pi_in = pi_out = ri_in = ri_out = 0;
					int out_degree = 0;
					for (int edgeID = 0; edgeID < resource->m; edgeID++)
						if(resource->edges[edgeID].tail == k )
							for(int waveID = 0; waveID < resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									pi_out ++;
					for (int edgeID = 0; edgeID < resource->m; edgeID++)
			            if(resource->edges[edgeID].head == k )
							for(int waveID = 0; waveID < resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									pi_in ++;
					for (int edgeID = 0; edgeID < edges_after.size(); edgeID++)
						if(edges_after[edgeID].head == demand_tail)
							for(int waveID = 0; waveID < edges_after[edgeID].bandwidth; waveID++)
									ri_in ++;	
					for (int edgeID = 0; edgeID < edges_after.size(); edgeID++)
						if(edges_after[edgeID].tail == demand_tail){
							out_degree++;
							for(int waveID = 0; waveID < edges_after[edgeID].bandwidth; waveID++)
									ri_out ++;
						}
					
					int pin_rin = pi_in - ri_in;
					int pout_rout = pi_out - ri_out;
					/*如果p-r<-0.75r,即为p<0.25r,说明该物理边无法提供最低限度的资源,将其删去*/
					if(pin_rin <= -(double(ri_in)*0.75) || pout_rout <= -(double(ri_out)*0.75)){
					//if(1){
						te->cmn = 10000000;
						auxGraph.edges.push_back(*te);
						auxGraph.m +=1;
						delete te;
						if_node_embeded = 0;
						continue;
					}
					else{
						//rin = abs(rin+pin_rin+pout_rout);
						/*这块公式比较乱,改了很多版本,要配合最后的paper看*/
						//double temp = pout_rout/out_degree;
						double temp = (double(pi_out)/count-resource->maxBandwidth/2)/out_degree;
						 new_rin = double(pi_out)/count ;//+ rin;	2016/11/21 删了后面那个 看看是不是差别不大
						//printf("%lf\n",temp);
						rin = 1+Cr*(0.5-1/(1+pow(2.718,-temp)));
					}
					

					//rin = pow(double(rin),-1)*0.1;
					
					/*计算∑dis(M,n)·δ/∑δ*/
					int Dis = 0,aa=0,mid=-1; //aa记录deta(v,n)的个数	
					for(int edgeID = 0;edgeID < resource->n;edgeID++){
						for(int kk = 0; kk < demand->n;kk++){
							if(FIk[kk][edgeID] == 1)
								//if(kk == demand_head || kk == demand_tail)	//2016/11/21	为了保证δ(v,i)直连
								if( int dis = findpath(resource,k,edgeID)){
									Dis+= dis;
									aa++;
								} 
						}
					}/*
					for (int ij = 0; ij < edges_after.size(); ij++) {

						if (edges_after[ij].tail == k) {                      //？？？？？这里一直有问题 什么意思？？？ 2016/9/4
							for (int kk = 0; kk < resource->n; kk++) {
								if (FIk[edges_after[ij].head][kk] == 1) {
									if_node_embeded = 1;
									mid = kk;//用mid记录节点v映射的物理节点
									break;
								}	
							}
							if (if_node_embeded = 1) {
								int dis = findpath(resource, k, mid);//调用一次find找到Dis(M,n),
								if (dis != 0) {
									Dis += dis;
									aa++;
								}
							}
						}
					}*/

					/*第一次没有任何点被映射上去时,设置默认距离为1000*/
					if( Dis != 0)
						dmn=Dis/aa;
					else
						dmn = 1000;
					//rin =1;
                   // double cmn= rin * dmn*Ucmn;
					double cmn = Ucmn * dmn + new_rin;
					if (cmn == 0)
						cmn = 10000000;
					te->cmn = cmn;

					///************************************************************************
					auxGraph.edges.push_back(*te);
					auxGraph.m +=1;
					delete te;
					}
						
					}
					if_node_embeded = 0;
				}
			}	


			if(FIk_j >= 0){
				if(resource->vertexWeight[FIk_j]>demand->vertexWeight[demand_head]){
				Edge *te = new Edge(FIk_j+resource->n,auxGraph.n-1,auxGraph.maxBandwidth ,r_mj[FIk_j],0);	//复制的图的终点，每个节点编号都要shift。
				auxGraph.edges.push_back(*te);
				auxGraph.m +=1;
				delete te;
				}
				else
					printf("????????????????");
			}			//如果i点被映射过了，那就只加一条被映射的边↑；如果没有被映射过，那就把其他所有边都加上去↓
			else{
				for(int k = 0; k < resource->n; k++){
					for (int demand_k = 0; demand_k < demand->n; demand_k ++){
						if(FIk[demand_k][k] == 1 )
							if_node_embeded = 1;
					}//判断某个点在这轮映射里有没有被用过，避免重复映射到同一点
				
					int dmn = 0;
					double rmj = 0;
					if(if_node_embeded == 0){
						if(resource->vertexWeight[k]>demand->vertexWeight[demand_head]){
					Edge *te = new Edge(k+resource->n,auxGraph.n-1,auxGraph.maxBandwidth ,r_mj[k],0);
/************************************************************************/
					int count=0;
					for (int edgeID = 0; edgeID < resource->m; edgeID++){
						
			            if(resource->edges[edgeID].head == k){
						    count++;
							for(int waveID = 0; waveID <= resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									rmj ++;//怎么表达?
						}
					}
					if (count != 0)
						rmj=rmj/g/count;	
					/************2016/10/21****/
					int pj_in,pj_out,rj_in,rj_out,new_rmj;
					pj_in = pj_out = rj_in = rj_out = 0;
					int in_degree = 0;
					for (int edgeID = 0; edgeID < resource->m; edgeID++)
						if(resource->edges[edgeID].tail == k )
							for(int waveID = 0; waveID < resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									pj_out ++;
					for (int edgeID = 0; edgeID < resource->m; edgeID++)
						if(resource->edges[edgeID].head == k )
							for(int waveID = 0; waveID < resource->edges[edgeID].bandwidth; waveID++)
								if(resource->edges[edgeID].wave[waveID] == true)
									pj_in ++;
					for (int edgeID = 0; edgeID < edges_after.size(); edgeID++)
						if(edges_after[edgeID].head == demand_head){
							in_degree++;
							for(int waveID = 0; waveID < edges_after[edgeID].bandwidth; waveID++)
									rj_in ++;			
						}
					for (int edgeID = 0; edgeID < edges_after.size(); edgeID++)
						if(edges_after[edgeID].tail == demand_head)
							for(int waveID = 0; waveID < edges_after[edgeID].bandwidth; waveID++)
									rj_out ++;
					//pj_in /= 2; //由于辅助图的关系，入边被算了两次
					int pin_rin = pj_in - rj_in;
					int pout_rout = pj_out - rj_out;
					if(pin_rin <= -(double(rj_in)*0.75) || pout_rout <= -(double(rj_out)*0.75)){
					//if(1){
						te->cmn = 10000000;
						auxGraph.edges.push_back(*te);
						auxGraph.m +=1;
						delete te;
						if_node_embeded = 0;
						continue;
					}
					else{
						//rmj = abs(rmj+pin_rin+pout_rout);
						
						//double temp = pin_rin/in_degree;
						double temp = (double(pj_in)/count-resource->maxBandwidth/2)/in_degree;
						 new_rmj = double(pj_in)/count ;//+ rmj;	2016/11/21 删了后面那个 看看是不是差别不大
						rmj = 1+Cr*(0.5-1/(1+pow(2.718,-temp)));
					}					
					//rmj = pow(double(rmj),-1)*0.1;
					
					
					int Dis = 0,aa=0,mid=-1; //aa记录deta(v,n)的个数
					
					for(int edgeID = 0;edgeID < resource->n;edgeID++){//16/11/3 原先是edgeID <= resource->n 内存报错后纠正,不知道为什么之前没错.
						for(int kk = 0; kk < demand->n;kk++){
							if(FIk[kk][edgeID] == 1)
								//if(kk == demand_head || kk == demand_tail)	//2016/11/21	为了保证δ(v,i)直连
								if( int dis = findpath(resource,k,edgeID)){
									Dis+= dis;
									aa++;
								} 
						}
					}/*
					for (int ij = 0; ij < edges_after.size(); ij++) {
						if (edges_after[ij].head == k ) {
							for (int kk = 0; kk < resource->n; kk++) {
								if (FIk[edges_after[ij].tail][kk] == 1) {
									if_node_embeded = 1;
									mid = kk;//用mid记录节点v映射的物理节点
									break;
								}
							}
							if (if_node_embeded = 1) {
								int dis = findpath(resource, k, mid);//调用一次find找到Dis(M,n),
								if (dis != 0) {
									Dis += dis;
									aa++;
								}
							}
						}

					}*/
					
					if( Dis != 0)
					 dmn=Dis/aa;	
					else
						dmn = 1000;
					//rmj = 1;
                   // double cmn= rmj * dmn*Ucmn;
					double cmn = Ucmn * dmn + new_rmj;
 					if (cmn == 0)
						cmn = 10000000;
					te->cmn = cmn;

					
///**********************************************************
					auxGraph.edges.push_back(*te);
					auxGraph.m +=1;
					delete te;
						}
					}
					if_node_embeded = 0;
				}
			}
			r_in.swap(vector<int>());
			r_mj.swap(vector<int>());
		return auxGraph;
			
}

/*dijkstra算法的实现*/
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
/*带参数的dijkstra,用来在生成图中求dis(M,n)*/
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


/*work函数,根据生成图和dijkstra找到一条路并进行判断,如果可行将其存入结果中*/

Heu_return Heuristic::work(){
	UsClock clk;
	clk.setStartTime();

	Heu_return heu_ret;
	heu_ret.isSuccess = 1;//判断整个映射是否成功
	heu_ret.resultWeight = 0;
	/*初始化heu_ret,具体定义在Heuristic.h中*/
	for (int ij = 0; ij < edges_after.size(); ij++) {
		if(heu_ret.isSuccess){
		vector<int> t;
		heu_ret.demand_tail.push_back(edges_after[ij].tail);
		heu_ret.demand_head.push_back(edges_after[ij].head);
		heu_ret.edge_head.push_back(t);
		heu_ret.edge_tail.push_back(t);
		heu_ret.total_distance.push_back(100000000);
		heu_ret.grid_a.push_back(-1);
		heu_ret.physical_distance.push_back(0);				//初始化结果类
		heu_ret.bandwidth.push_back(-1);
		heu_ret.modtype.push_back(-2);
		

		int flag_FIk_same_point = 0, same_point_ID;		//判断是否映射到同一点
		int DijisSuccess;
		int modtype = 4;
		int argmod,argi;
		
		/*argvmin mod 使算法不至于太贪心*/
		argmod = 10000;
		for(int i = 1; i<= 4; i++){
			int p = ceil(float(edges_after[ij].bandwidth)/i)*10+i;
			if( p < argmod ){
				argi = i;
				argmod = p;
			}
		}
		printf("%d %d",edges_after[ij].bandwidth,argi);
		//modtype = argi;
		//if(modtype < 4) modtype += 1;
		for (int mod = 3; mod >= 0; mod--) {
			modtype--;
			int g = ceil(float(float(edges_after[ij].bandwidth) / (mod + 1)));//g为在当前mod下需要使用的带宽
			int F_g_1 = resource->maxBandwidth - g + 1;	//这里擅自取了F=maxBandwidth 循环多了应该也不影响，能找到会在之前就找到，找不到循环到底也找不到
			for (int a = 0; a < F_g_1; a++) {
				heu_ret.grid_a[ij] = a;
				int maxdistance = resource->moddistance[modtype];
				
				if (flag_FIk_same_point == 0)
					auxGraph_generate(edges_after[ij].head, edges_after[ij].tail, a, g, maxdistance);	//生成辅助图
				else
				{
					/*如果在上一次选路中起点和终点被选为同一点,则在这里将其中一条虚边删除确保此次选路不会选到同一点*/
					flag_FIk_same_point = 0;
					int m_number = auxGraph.m;
					vector<Edge>::iterator it = auxGraph.edges.begin();
					for (int mn = 0; mn < m_number; mn++) {
						if (auxGraph.edges[mn].tail == same_point_ID &&auxGraph.edges[mn].head == auxGraph.n - 1) {
							printf("delete same point: %d", same_point_ID);
							auxGraph.edges.erase(it);
							auxGraph.m--;
							break;
							//mn--;
						}
						else
							it++;
					}
				}							//删除其中一条虚边，保证不会映射到同一个点

				DijisSuccess = findpath();		//根据生成辅助图寻找最短路
				if (!DijisSuccess) {

					continue;
				}
				else {
					if (DijisSuccess <= heu_ret.total_distance[ij]) {//为啥要小于上一次距离来着..忘记了 16/11/5
						heu_ret.total_distance[ij] = DijisSuccess;
						heu_ret.edge_head[ij].clear();
						heu_ret.edge_tail[ij].clear();			//结果更小，删除上一次保留的结果

						int j_m, i_n, temp = auxGraph.n - 1;
						j_m = last_vertex_ID[auxGraph.n - 1];
						heu_ret.edge_head[ij].push_back(auxGraph.n - 1);
						heu_ret.edge_tail[ij].push_back(j_m);
						/*将dijkstra中的结果存储到heu_ret结构中*/
						do {
							i_n = temp;
							temp = last_vertex_ID[temp];
							heu_ret.edge_head[ij].push_back(temp);
							heu_ret.edge_tail[ij].push_back(last_vertex_ID[temp]);
						} while (last_vertex_ID[temp] != auxGraph.n - 2);


						i_n = temp;
						/*两个点映射到同一点了,重新选路并在生成辅助图时将其中一条虚边删除*/
						if (i_n + resource->n == j_m) {			//两个点I,J映射到同一个物理节点的情况
							DijisSuccess = 0;
							printf("same_point = %d,a = %d\n!!!!!!!!!!!!!!!!!!", i_n, a);
							flag_FIk_same_point = 1;
							same_point_ID = j_m;
							heu_ret.total_distance[ij] = 100000000;
							if (a == F_g_1 - 1) {
								if(mod > 1){	//之前没有加这个，找完a没有就直接判定失败了，其实应该把mod降低继续找的 10/21
									break;
								}
								heu_ret.isSuccess = 0;
								printf("no enough bandwidth!\n");
								break;
							}
							a--;
							continue;
						}
						
						
					}
					
					break;
				}
			}


			int physical_distance = 0;
			if (DijisSuccess) {				//如果找到了一条路径，break 并将资源情况更新
				for (int mn = 0; mn < resource->m; mn++)
					for (int k = 1; k < heu_ret.edge_head[ij].size() - 1; k++)		//不能从0开始 第0个点是虚边起点，最后一个点也不能算是虚边终点,只算物理边节点
						if (resource->edges[mn].head == heu_ret.edge_head[ij][k] || resource->edges[mn].head == heu_ret.edge_head[ij][k] - resource->n)
							if (resource->edges[mn].tail == heu_ret.edge_tail[ij][k] || resource->edges[mn].tail == heu_ret.edge_tail[ij][k] - resource->n) {
								physical_distance += resource->edges[mn].distance;
								for (int a = heu_ret.grid_a[ij]; a < heu_ret.grid_a[ij] + g; a++)
									resource->edges[mn].wave[a] = false;
							}
				/*判断总的物理长度是否大于当前mod下的极限值*/
				if (physical_distance >= resource->moddistance[mod]) {
					if (mod == 0) {
						heu_ret.isSuccess *= 0;
						printf("physical failed!!\n");
						break;//物理映射失败 跳出
					}
					else
						continue;
				}
				else {
					/*物理距离和资源都满足,此次映射成功,将节点映射的情况存到FIk数组中*/
					for (int demand_k = 0; demand_k < resource->n; demand_k++) {
						if (demand_k == heu_ret.edge_head[ij][1] || demand_k == heu_ret.edge_head[ij][1] - resource->n)
								FIk[edges_after[ij].head][demand_k] = 1;
						}
					int heusize = heu_ret.edge_tail[ij].size();
						for (int demand_k = 0; demand_k < resource->n; demand_k++) {
							if (demand_k == heu_ret.edge_tail[ij][heusize-2] || demand_k == heu_ret.edge_tail[ij][heusize-2] - resource->n)
								FIk[edges_after[ij].tail][demand_k] = 1;
						}														//将节点映射记录下来

					heu_ret.physical_distance[ij] = physical_distance;
					heu_ret.modtype[ij] = mod;
					heu_ret.bandwidth[ij] = g;

					int demand_t = ij;
							for(int t = 1; t < heu_ret.edge_head[demand_t].size()-1;t++){
								int e_head = heu_ret.edge_head[demand_t][t],e_tail = heu_ret.edge_tail[demand_t][t];
								if(e_head > resource->n)
									e_head -= resource->n;
								if(e_tail > resource->n)
									e_tail -= resource->n;
								for(int mn = 0; mn < resource->m; mn++)
									if(resource->edges[mn].head == e_head && resource->edges[mn].tail == e_tail)
										heu_ret.resultWeight += heu_ret.bandwidth[demand_t] * resource->edges[mn].cmn;
								
							}
						
					printf("******total distance = %d\n", heu_ret.total_distance[ij]);
					heu_ret.isSuccess *= 1;
					break;
				}
			}

		}
		}
	}
	for (int is = 0; is < heu_ret.edge_head.size(); is++)			//判断是否findpath全部成功，如果有一个不成功则映射失败。
		if (heu_ret.edge_head[is].size() == 0)
			heu_ret.isSuccess = 0;
	clk.setEndTime();
	heu_ret.proccessTime = clk.getTime();

	return heu_ret;
}




Heuristic::Heuristic(Network *d, Network *r, FILE *lf){
	ori = r;
	resource = new Network (*r);
	demand = new Network (*d);	
	resource->mod++;
	logFile = lf;
	/*初始化FIk数组*/
	for(int nodei=0;nodei<demand->n;nodei++){
		vector<int> tempFIk;
		FIk.push_back(tempFIk);
	}
	for(int nodei=0;nodei<demand->n;nodei++){
		for(int k=0;k<resource->n;k++){
			FIk[nodei].push_back(0);
		}
	}
	
	sortnetwork(1); //0->带宽排序	1->度排序	2->先度后带宽排序
}

Heuristic::~Heuristic(){

	edges_before.swap(vector<Edge>());
	edges_after.swap(vector<Edge>());
	edges_middle.swap(vector<Edge>());
//	fclose(optimalData);
}
/*这里将heu_ret的结果转存到Event ret中,为了在schedule.h中保持接口统一*/
Event Heuristic::solveByHeuristic(Network *d, Network *r, FILE *logFile){
	Event ret = Event(0, 0, 0, r);
	int numILPSuccess, numILPBlock, numHeuSuccess=0;
		double totalILPTime=0, totalILPWeight=0, totalHeuTime, totalHeuWeight, HeuWeight;
		Network original_resource = *r;
		Network *resource_Hr = &original_resource;
			 Heu_return heu_ret;
			 heu_ret = work();
			 HeuWeight = 0;

			  if(heu_ret.isSuccess){
				numHeuSuccess += 1;
				totalILPTime += heu_ret.proccessTime;
				totalILPWeight += HeuWeight;
				fprintf(logFile, "Heu success, get object of %f, using Time %f\n ", HeuWeight, heu_ret.proccessTime);
				ret.proccessTime = heu_ret.proccessTime;
				/*将辅助图的节点编号还原成实际物理图上的节点编号*/
				for(int ij = 0; ij < heu_ret.demand_tail.size(); ij++){
					for(int k = 0; k < heu_ret.edge_head[ij].size(); k++){
						if(heu_ret.edge_head[ij][k] >= resource_Hr->n )
							heu_ret.edge_head[ij][k] -= resource_Hr->n;
						if(heu_ret.edge_tail[ij][k] >= resource_Hr->n )
							heu_ret.edge_tail[ij][k] -= resource_Hr->n;
					}
				}

			 if(heu_ret.isSuccess){

				for(int ij = 0; ij < demand->n; ij++)					
					for(int k = 0; k < resource->n; k++)				
						if(FIk[ij][k] == 1){
							HeuWeight += demand->vertexWeight[ij]*resource->ck[k];
							printf("**** HeuWeight = %f, vertexID: %d, vertexWeight : %d, ck = : %f\n",HeuWeight,ij, demand->vertexWeight[ij],resource->ck[k]);
						}

						HeuWeight += heu_ret.resultWeight;
			 }
			 

				for(int ij = 0; ij < heu_ret.demand_tail.size(); ij++){
					printf("Virtual node: %d, %d ==> Physical node: %d, %d",heu_ret.demand_tail[ij],heu_ret.demand_head[ij],heu_ret.edge_head[ij][heu_ret.edge_head[ij].size()-1],heu_ret.edge_tail[ij][0]);
					printf("	modtype: %d , start grid : %d ,bandwidth: %d\n",heu_ret.modtype[ij],heu_ret.grid_a[ij],heu_ret.bandwidth[ij]);
					for(int pk = 0; pk < heu_ret.edge_tail[ij].size(); pk++)
						printf("--%d",heu_ret.edge_tail[ij][pk]);
					printf("\n");
				}
			}
			 else{
					fprintf(logFile, "Heu failed\n");
					printf("\n********************************\n embedding failed!\n");
					
			 }
			 ret.resultWeight = HeuWeight;

if(heu_ret.isSuccess){
			for(int ij = 0; ij < ret.demand->edges.size(); ij++)
				for(int b = 0; b < ret.demand->edges[ij].wave.size(); b++)
					ret.demand->edges[ij].wave[b] = false;

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
					if(r->edges[mn].head == heu_ret.edge_head[ij][k])
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
}
		ret.isSuccess = heu_ret.isSuccess;

		return ret;
		
}