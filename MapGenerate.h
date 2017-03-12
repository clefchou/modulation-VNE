#pragma once
#include"common.h"


#define d_same 300		//定义相同区域平均距离
#define d_diff 1000		//定义不同区域平均距离
#define max_slot 20
//#define m_degree 3		//平均每个节点的度

int uniform(int a, int b){
	return a + rand()%(b-a+1);
}

void mapGenerate(char route[], int n, int m, int area, int maxWeight,int maxdistance){
/*
	生成随机拓扑
	route为拓扑文件存储路径
	n为点数，m为边数
	area为需要划分的区域数
*/
	FILE *out = fopen(route, "w");
	
	/*在一个vector中存area个数，其中保存同区域最大节点序号
	 *默认所有节点按顺序排列，例如1-10在一个区域，11-15在一个区域以此类推
	 *vector中则存[10,15,..]
	 */
	vector<int> area_n;
	for(int i = 0; i < area; i++){
		area_n.push_back(rand()%n-1);
	}
	sort(area_n.begin(),area_n.end());
	vector<int>::iterator k = area_n.begin();
    area_n.erase(k);
	area_n.push_back(n);//最后添加一个点作为终点
	//随机计算平均节点度
	int n_degree = uniform(1,2*m/n);
	/* n_ 和当前点连接的点, n_index 从属于哪个区域, flag 是否被映射第一条边
	 * distance 随机生成的距离, m_ 现在已经映射的边数
	 */
	int n_,n_index,flag,distance,m_;
	n_index = 0;
	flag = 0;
	m_ = 0;
	for(int i = 0; i < n; i++)
		fprintf(out, "%d ", 1000); //节点的资源
	fprintf(out, "\n");
	for(int i = 0; i < m; i++)
		fprintf(out, "%d ", 2); //???
	fprintf(out, "\n");

	for(int i = 0; i < n; i++){
		//判断在哪个区域
		if(i >= area_n[n_index])
			n_index++;
		for(int j = 0; j < uniform(0,n_degree);j++){
			//如果至少有一条边那么必须和区域内的点相连
			if(flag==0){
				if(n_index == 0){
					do{
						n_ = rand()%(area_n[n_index]);//防止溢出
					}while(n_ == i);
				}
				else{
					do{
						n_ = rand()%(area_n[n_index]-area_n[n_index-1])+area_n[n_index-1];
					}while(n_ == i);
				}
				distance = uniform(0,2*d_same);
				fprintf(out, "%d %d %d %d %d\n", i, n_,max_slot, 1, distance);
				fprintf(out, "%d %d %d %d %d\n", n_, i,max_slot, 1, distance);
				flag = 1;
				if(++m_ == m)
					break;
			}
			do{
				n_ = rand()%n;
			}while(n_ == i);
			distance = uniform(0,2*d_diff);
			fprintf(out, "%d %d %d %d %d\n", i, n_,max_slot, 1, distance);
			fprintf(out, "%d %d %d %d %d\n", n_, i,max_slot, 1, distance);
			if(++m_ == m)
				break;
		}
		flag = 0;
		if(m_ == m)
			break;
	}
	//如果边数不足的情况
	if(m_< m){
		int n__;
		for(int i = m_; i < m; i++){
			do{
				n_ = rand()%n;
				n__ = rand()%n;
			}while(n_ == n__);
			distance = uniform(0,2*d_diff);
			fprintf(out, "%d %d %d %d %d\n", n__, n_,max_slot, 1, distance);
			fprintf(out, "%d %d %d %d %d\n", n_, n__,max_slot, 1, distance);
		}
	}

	fclose(out);
	system("type \inputData\\demand_r.txt>>\inputData\\alldemand.txt");
}




int findpath(Network *r,int source_vertexID, int dest_vertexID){
	int temp_dist=100000000,temp_node_id ,findcounts = 0;
	Network *tempnet = r;
	vector<double>dist(tempnet->n,1000000000);
	vector<int> last_vertex_ID;							//记录上一个节点的ID

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


void chooseLi(Network *r,int boundary){
	Network *tempr = r;
	vector<vector<int>>node;
	//调用最短路算法算出每个点到其他所有点距离
	for(int i = 0; i < tempr->n; i++){
		vector<int>ndist;
		for(int j = 0; j < tempr->n; j++){
			int dist = findpath(tempr,i,j);
			//根据boundary限制得到每个点的相邻区域
			if(dist < boundary)
				ndist.push_back(j);	  
		}
		node.push_back(ndist);
	}
	r->neighbor = node;


}
//随机给demand节点安排所要映射的区域
void randomGenerate(Network *d,Network *r){

	vector<int> startArray;
	for(int i = 0; i <r->n; i++){
		startArray.push_back(i);
	}	    
	int N = d->n;//随机数个数    
	vector<int> resultArray;  
	for(int i = 0; i < N; i++)    
	{    
		int seed = rand()%(startArray.size()-i);    
		vector<int> pp;
		d->neighbor.push_back(pp);
		d->neighbor[i].push_back(startArray[seed]);
		//resultArray.push_back(startArray[seed]);//赋值给结果数组    
		startArray[seed] = startArray[startArray.size() - i - 1];//把随机数产生过的位置替换为未被选中的值。    
	}    

}
