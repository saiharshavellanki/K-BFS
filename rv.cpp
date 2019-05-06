#include<bits/stdc++.h>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

unordered_map<int,int> bfs(int v,vector<vector<int> > &ind,int &maxi)
{
	unordered_map<int,int> frontier;
	queue<int> q;
	int par,i;
	frontier[v]=0;
	q.push(v);
	while(!q.empty())
	{
		par = q.front();
		q.pop();
		if(q.size()==0)
			maxi=frontier[par];
		for(i=0;i<ind[par].size();i++)
		{
			if(frontier.find(ind[par][i])==frontier.end())
			{
				frontier[ind[par][i]] = frontier[par]+1;
				q.push(ind[par][i]);
			}
		}
	}
	return frontier;
}

unordered_map<int,int> nearestneighbour(int v,vector<vector<int> > &ind,unordered_set<int> &NSw)
{
	unordered_map<int,int> frontier,nearneigh;
	for(auto it=NSw.begin();it!=NSw.end();++it)
		nearneigh[*it]=*it;
	queue<int> q;
	int par,i;
	frontier[v]=0;
	q.push(v);
	while(!q.empty())
	{
		par = q.front();
		q.pop();
		for(i=0;i<ind[par].size();i++)
		{
			if(frontier.find(ind[par][i])==frontier.end())
			{
				frontier[ind[par][i]] = 1;
				q.push(ind[par][i]);
				if(NSw.find(ind[par][i])==NSw.end())
				{
					nearneigh[ind[par][i]] = nearneigh[par];
				}
			}
		}
	}
	return nearneigh;
}

bool compare(pair<int,int> p1, pair<int,int> p2)
{
	return p1.second>p2.second;
}
vector<vector<int> > ind;
vector<int> ecc,eprime;

int main()
{
	int max_threads = omp_get_max_threads();
	cout<<"No of threads : "<<max_threads<<endl;
	omp_set_num_threads(max_threads);

	int n,m,i,j,k,a,b,q,t,tid,s,w;
	scanf("%d %d",&n,&m);
	ind.resize(n+1);

	cout<<"No of vertices : "<<n<<" and No of edges : "<<m<<endl;

	for(i=0;i<m;i++)
	{
		scanf("%d	%d",&a,&b);
		ind[a].push_back(b);
		ind[b].push_back(a);
	}
	cout<<"scanned\n";
	ecc.resize(n+1, 0);
	eprime.resize(n+1, 0);

	auto start = high_resolution_clock::now();
	for(i=1;i<=n;i++)
	{
		if(ecc[i]==0)
		{
//			printf("Eccentricity of %d started\n",i);
			vector<int> comp;
			unordered_map<int,int> frontier,nearneigh;
			unordered_map<int,int>::iterator it;
			int l;
			frontier = bfs(i,ind,l);

			for(it=frontier.begin();it!=frontier.end();++it)
			{
				comp.push_back(it->first);
			}

			int num = comp.size();
			if(num==1)
				continue;
			if(num==2)
			{
				ecc[comp[0]]=1;
				ecc[comp[1]]=1;
				continue;
			}
			s = (int)(sqrt((double)num*(double)log(num)))/50;

			for(j=0;j<s;j++)
			{
				q=rand()%num;
				t = comp[j];
				comp[j]=comp[q];
				comp[q] = t;
			}

			vector<int> S,NSw;
			for(j=0;j<s;j++)
				S.push_back(comp[j]);

			unordered_map<int,unordered_map<int,int>> dists, distw;
			vector<int> mini(num, INT_MAX);
			
			int minecc_ins=INT_MAX,tasks,subtasks;
			tasks = s/max_threads+(int)(s%max_threads!=0);
			subtasks = num/max_threads+(int)(num%max_threads!=0);

//			printf("bfs of S started with %d tasks and %d subtasks %d vertices\n",tasks,subtasks,s);
			for(int j=0;j<tasks;j++){
				#pragma omp parallel private(tid)
				{
					tid = omp_get_thread_num();
					int y = j*max_threads+tid;
					if(y<s)
						dists[S[y]] = bfs(S[y], ind, ecc[S[y]]);
				}
				int garb = max_threads;
				if(j+1==tasks)
					garb = s%max_threads;
				for(int y=0;y<subtasks;y++)
				{
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num();
                                                int z = y*max_threads+tid;
						if(z<num)
						{
							for(int p=j*max_threads;p<j*max_threads+garb && p<s;p++)
							{
								mini[z]=min(mini[z],dists[S[p]][comp[z]]);
							}
						}
					}
//					printf("task %d subtask %d complete\n",j,y);
				}
			}
			for(j=0;j<s;j++)
				minecc_ins = min(minecc_ins,ecc[S[j]]);
			w=0;
			for(j=1;j<num;j++)
				if(mini[w]<mini[j])
					w=j;
			w=comp[w];

			frontier=bfs(w, ind,l);
			vector<pair<int,int>> temp(frontier.begin(), frontier.end());
			sort(temp.begin(), temp.end(), compare);

			for(j=0;j<s;j++)
				NSw.push_back(temp[j].first);

			tasks = s/max_threads+(int)(s%max_threads!=0);
//			printf("bfs of NSw started\n");
			for(j=0;j<tasks;j++){
				#pragma omp parallel private(tid)
				{
					tid = omp_get_thread_num();
					int y = j*max_threads+tid;
					if(y<s)
						distw[NSw[y]] = bfs(NSw[y], ind,ecc[NSw[y]]);
				}
			}

			unordered_set<int> others_set(comp.begin(), comp.end());
			for(j=0;j<s;j++)
			{
				others_set.erase(S[j]);
				others_set.erase(NSw[j]);
			}
			vector<int> others(others_set.begin(), others_set.end());

//			printf("eprime calculation started\n");
			#pragma omp prallel for
			{
			for(j=0;j<others.size();j++){
				eprime[others[j]]=distw[w][others[j]];
				for(k=0;k<s;k++)
					eprime[others[j]]=max(eprime[others[j]], dists[S[k]][others[j]]);
			}
			}
//			printf("eprime calculation ended\n");
			unordered_set<int> NSW(NSw.begin(),NSw.end());
			nearneigh = nearestneighbour(w,ind,NSW);
			
			for(j=0;j<others.size();j++)
			{
				int v = others[j],vt=nearneigh[v];
				if(distw[vt][v]<=distw[vt][w])
				{
					ecc[v]=max(eprime[v],ecc[vt]);
				}
				else
				{
					ecc[v]=max(eprime[v],minecc_ins);
				}
			}
		}
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << "Time taken for rv is " << duration.count() << endl;
	return 0;
}
