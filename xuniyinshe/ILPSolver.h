#include "common.h"
#include "Graph.h"
#include <ilcplex/ilocplex.h>
#include "UsTime.h"




extern double ILPT;

Event solveBySimpleILP(Network *d, Network *r, FILE *logFile){
	Event ret = Event(0, 0, 0, r);

	IloEnv environment;
    IloModel model(environment);
	//IloCplex cplex(environment);
	//cplex.setParam(IloCplex::NodeFileInd,2);

	
		IloArray<IloArray<IloIntVarArray>> M(environment, d->m);
	    for(int ij = 0; ij < d->m; ij++){
		    M[ij] = IloArray<IloIntVarArray>(environment, r->m);
		    for(int mn = 0; mn < r->m; mn++)
			    M[ij][mn] = IloIntVarArray(environment, r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1, 0, 1);
	    }
		IloArray<IloArray<IloArray<IloIntVarArray>>>Z(environment, d->m);
	    for(int ij = 0; ij < d->m; ij++){
		    Z[ij] = IloArray<IloArray<IloIntVarArray>>(environment, r->m);
		    for(int mn = 0; mn < r->m; mn++){
				Z[ij][mn]=IloArray<IloIntVarArray>(environment, r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1);
				for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++){
					Z[ij][mn][a]=IloIntVarArray(environment, r->mod, 0, 1);}}}
	    IloArray<IloIntVarArray> A(environment, d->n);
	    for(int i = 0; i < d->n; i++)
		    A[i] = IloIntVarArray(environment, r->n, 0, 1);
	    IloArray<IloIntVarArray> T(environment, d->m);
	    for(int ij = 0; ij < d->m; ij++)
		    T[ij] = IloIntVarArray(environment, r->maxBandwidth - d->edges[ij].bandwidth + 1, 0, 1);
		IloArray<IloIntVarArray> S(environment, d->m);
		for(int ij = 0; ij < d->m; ij++){
			S[ij]=IloIntVarArray(environment, r->mod, 0, 1);}


		/**************************************************优化目标*********************************************************************/
	
		IloExpr BP(environment);
	    for(int ij = 0; ij < d->m; ij++)
		    for(int mn = 0; mn < r->m; mn++)
			    for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					for(int mod = 0; mod < r->mod; mod++){
						//double p1=0;
						//p1 += d->edges[ij].bandwidth/( mod+1 );
				        BP += r->edges[mn].cmn * ceil(float (float(d->edges[ij].bandwidth)/( mod+1 ))) * Z[ij][mn][a][mod];
					}
	    for(int i = 0; i < d->n; i++)
		    for(int k = 0; k < r->n; k++)
			   ;// BP += r->ck[k] * d->vertexWeight[i] * A[i][k];
	    model.add(IloMinimize(environment, BP));
		

		/****************************************以下是约束段***************************************************/
		
		//4
		for(int i = 0; i < d->n; i++){
		    IloExpr constraint(environment);
		    for(int k = 0; k < r->n; k++)
			    constraint += A[i][k];
		    model.add(constraint == 1);
	    }
		//5
		for(int k = 0; k < r->n; k++){
		    IloExpr constraint(environment);
		    for(int i = 0; i < d->n; i++)
			    constraint += A[i][k];
		    model.add(constraint <= 1);
	    }
		//6
		for(int ij = 0; ij < d->m; ij++)
		    for(int k = 0; k < r->n; k++){
			    IloExpr constraint(environment);
			    for(int mn = 0; mn < r->m; mn++)
				    if(r->edges[mn].head == k)
					    for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
						    constraint += M[ij][mn][a];
			    for(int mn = 0; mn < r->m; mn++)
				    if(r->edges[mn].tail == k)
					    for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
						    constraint -= M[ij][mn][a];
			    constraint -= A[d->edges[ij].head][k] - A[d->edges[ij].tail][k];
			    model.add(constraint == 0);

		}
		//7
		for(int ij = 0; ij < d->m; ij++){
		IloExpr constraint(environment);
		for(int a = 0; a < r->maxBandwidth - d->edges[ij].bandwidth + 1; a++)
			constraint += T[ij][a];
		model.add(constraint == 1);
	    }
		//8
		for(int ij = 0; ij < d->m; ij++)
		    for(int a = 0; a < r->maxBandwidth - d->edges[ij].bandwidth + 1; a++){
			    IloExpr constraint(environment);
			    for(int mn = 0; mn < r->m; mn++)
				    if(a + d->edges[ij].bandwidth <= r->edges[mn].bandwidth)
				 	    constraint += M[ij][mn][a];
			    constraint -= T[ij][a];
			    model.add(constraint >= 0);
		}
		//9
		for(int ij = 0; ij < d->m; ij++)
		    for(int a = 0; a < r->maxBandwidth - d->edges[ij].bandwidth + 1; a++){
			    IloExpr constraint(environment);
			    for(int mn = 0; mn < r->m; mn++)
				    if(a + d->edges[ij].bandwidth <= r->edges[mn].bandwidth)
					    constraint += M[ij][mn][a];
			    constraint -= 1000 * T[ij][a];
			    model.add(constraint <= 0);
		}
		//10
		for(int ij = 0; ij < d->m; ij++){
			IloExpr constraint(environment);
		    for(int mod = 0; mod < r->mod; mod++)
				constraint += S[ij][mod];
			model.add(constraint == 1);
		}
		//11
		for(int ij = 0; ij < d->m; ij++){
			IloExpr constraint(environment);
			for(int mn = 0; mn < r->m; mn++)
			    for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					constraint += M[ij][mn][a] * r->edges[mn].distance;
			for(int mod = 0; mod < r->mod; mod++)
				constraint -= S[ij][mod] * r->moddistance[mod];
			model.add(constraint <= 0);
		}
		//12
		for(int k = 0; k < r->n; k++){
		IloExpr constraint(environment);
		for(int i = 0; i < d->n; i++)
			constraint += d->vertexWeight[i] * A[i][k];
		model.add(constraint <= r->vertexWeight[k]);
	    }
		//13
		for(int ij = 0; ij < d->m; ij++)
			for(int mn = 0; mn < r->m; mn++)
				for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					for(int mod = 0; mod < r->mod; mod++){
						IloExpr constraint(environment);
						//double p1=0;
						//p1 += d->edges[ij].bandwidth/( mod+1 );
						constraint += Z[ij][mn][a][mod] * ceil(float (float(d->edges[ij].bandwidth)/( mod+1 )));
						for(int t = a; t < a + ceil(float (float(d->edges[ij].bandwidth)/( mod+1 ))); t++)
						    constraint -= r->edges[mn].wave[t];
						model.add(constraint <= 0);
					}
		//14
					
        for(int mn = 0; mn < r->m; mn++)
		    for(int a = 0; a < r->edges[mn].bandwidth ; a++)
			    for(int mod = 0; mod < r->mod; mod++){
			        IloExpr constraint(environment);
					constraint -= r->edges[mn].wave[a];
				    for(int ij = 0; ij < d->m; ij++)
						for(int t = 0; t< r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1 ; t++){
							//double p1=0;
							//p1 += d->edges[ij].bandwidth/( mod+1 );
						    if(t <= a && t + ceil(float (float(d->edges[ij].bandwidth)/( mod+1 ))) -1>= a)
								 constraint += Z[ij][mn][t][mod];}
					    model.add(constraint <= 0);
				}
		//15
	    for(int ij = 0; ij < d->m; ij++)
			for(int mn = 0; mn < r->m; mn++)
				for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					for(int mod = 0; mod < r->mod; mod++){
						IloExpr constraint(environment);
						constraint += M[ij][mn][a] + S[ij][mod];
						constraint -= Z[ij][mn][a][mod];
						model.add(constraint <= 1);
					}
		//16
        for(int ij = 0; ij < d->m; ij++)
			for(int mn = 0; mn < r->m; mn++)
				for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					for(int mod = 0; mod < r->mod; mod++){
						IloExpr constraint(environment);
						constraint += M[ij][mn][a] + S[ij][mod];
						constraint -= 2 * Z[ij][mn][a][mod];
						model.add(constraint >= 0);
					}
	
	IloCplex solver(model);
	//environment.out()<<"   "<<solver.getObject();
	solver.setOut(environment.getNullStream());
	//FILE *file;
	//file = fopen("D:\\毕设\\xuni\\outputData\\file.lp", "w");
	//solver.exportModel("D:\\毕设\\xuni\\outputData\\file.lp");
	//fclose(file);

	UsClock clk;
	clk.setStartTime();
	ret.isSuccess = solver.solve();
	clk.setEndTime();
	ret.proccessTime = clk.getTime();

	if(ret.isSuccess){
		ret.resultWeight = solver.getObjValue();
		for(int ij = 0; ij < ret.demand->edges.size(); ij++)
			for(int b = 0; b < ret.demand->edges[ij].wave.size(); b++)
				ret.demand->edges[ij].wave[b] = false;
		for(int mn = 0; mn < r->m; mn++)
			for(int ij = 0; ij < d->m; ij++)
				for(int a = 0; a < r->edges[mn].bandwidth - d->edges[ij].bandwidth + 1; a++)
					for(int mod = 0; mod < r->mod; mod++)
					    if(solver.getValue(Z[ij][mn][a][mod]) + 0.5 > 1)
						    for(int t = a; t <= a + ceil(float (float(d->edges[ij].bandwidth)/( mod+1 ))) - 1; t++)
							    ret.demand->edges[mn].wave[t] = true;
		for(int i = 0; i < 10; i++)
		      for(int j = 0; j < 20; j++)
				  for (int a =0; a < 30; a++)
					  for(int mod =0; mod < 4; mod++){
						  ret.b[i][j][a][mod] = 0;
					  }
		for(int i = 0; i < d->m; i++)
		      for(int j = 0; j < r->m; j++)
				  for (int a =0; a < r->edges[j].bandwidth - d->edges[i].bandwidth + 1; a++)
					  for(int mod =0; mod < r->mod; mod++)
					  {
					      ret.b[i][j][a][mod] = solver.getValue(Z[i][j][a][mod]);
					  }//Z
		/*for(int i = 0; i < 10; i++)
		      for(int j = 0; j < 20; j++)
				  for (int a =0; a < 30; a++){
					  ret.c[i][j][a] = 0;
				  }
		for(int i = 0; i < d->m; i++)
		      for(int j = 0; j < r->m; j++)
				  for (int a =0; a < r->edges[j].bandwidth - d->edges[i].bandwidth + 1; a++){
					  ret.c[i][j][a] = solver.getValue(M[i][j][a]);
				  }//M*/
		for(int i = 0; i < 10; i++)
		      for(int j = 0; j < 10; j++)
			 {
			      ret.a[i][j] = 0;
			}
		for(int i = 0; i < d->n; i++)
			for(int k = 0; k < r->n; k++)
			{
			    ret.a[i][k] = solver.getValue(A[i][k]);
					}//A
		/*for(int i = 0; i < 10; i++)
			for (int a =0; a < 30; a++){
				ret.d[i][a] = 0;
			}
		for(int i = 0; i < d->m; i++)
			for (int a =0; a < r->maxBandwidth - d->edges[i].bandwidth + 1; a++){
				ret.d[i][a] = solver.getValue(T[i][a]);
			}//T
	  for(int i = 0; i < 10; i++)
		  for(int mod =0; mod < 4; mod++){
		      ret.e[i][mod] = 0;
		  }
      for(int i = 0; i < d->m; i++)
		  for(int mod =0; mod < r->mod; mod++){
			  ret.e[i][mod] = solver.getValue(S[i][mod]);
		  }//S */
		for(int k = 0; k < r->n; k++)
			ret.demand->vertexWeight[k] = 0;
		for(int k = 0; k < r->n; k++)
			for(int i = 0; i < d->n; i++)
				if(solver.getValue(A[i][k]) + 0.5 > 1)
					ret.demand->vertexWeight[k] += d->vertexWeight[i];
	}
	else
		;

	A.end();
	for(int i = 0; i < d->m; i++)
		T[i].end();
	T.end();
	for(int i = 0; i < d->m; i++)
		S[i].end();
	S.end();
	for(int i = 0; i < d->m; i++){
		for(int j = 0; j < r->m; j++)
			M[i][j].end();
		M[i].end();
	}
	M.end();
	for(int i = 0; i < d->m; i++){
		for(int j = 0; j < r->m; j++)
			Z[i][j].end();
		Z[i].end();
	}
	Z.end();
	solver.end();
    model.end();
    environment.end();
	BP.end();

	return ret;
}


