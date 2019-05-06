#include<bits/stdc++.h>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

unordered_map<int,int> bfs(int v,vector<vector<int> > &ind)
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

void kbfs(vector<int> frontier,vector<vector<int> > &ind,int num,unordered_map<int,int> &d)
{
	int i,tid;

	while(frontier.size()>0){
		int tasks = frontier.size()/num+(frontier.size()%num!=0);
		vector<int> new_frontier;
		for(i=0;i<tasks;i++)
		{
			#pragma omp parallel private(tid)
			{
				tid = omp_get_thread_num();
				int z = i*num+tid;
				if(z<frontier.size())
				{
					int par = frontier[z];
					int level = d[par]+1;
					for(int j=0;j<ind[par].size();j++)
					{
						if(d.find(ind[par][j])==d.end())
						{
							#pragma omp critical
							{
								d[ind[par][j]] = level;
								new_frontier.push_back(ind[par][j]);
							}
						}
					}
				}
			}
		}
		int sz = new_frontier.size();
		frontier = new_frontier;
		frontier.resize(sz);
	}
}

bool cmpfunc(pair<int,int> p1,pair<int,int> p2)
{
	return p1.second<p2.second;
}
vector<vector<int> > ind;
vector<int> ecc;

int main()
{
	int max_threads = omp_get_max_threads();
	cout<<"No of threads : "<<max_threads<<endl;
	omp_set_num_threads(max_threads);

	int n,m,i,j,k,r,a,b,q,t,tid,K;
//	n=1157827;
//	m=2987624;
	scanf("%d %d", &n, &m);
	cout<<"No of vertices : "<<n<<" and No of edges : "<<m<<endl;
	ind.resize(n+1);
	for(i=0;i<m;i++)
	{
		scanf("%d	%d",&a,&b);
		ind[a].push_back(b);
		ind[b].push_back(a);
	}
	for(int val = 6 ; val <= 13 ; val++)
	{
		ecc.resize(n+1,0);
		K = 1<<val;
		auto start = high_resolution_clock::now();
		for(i=1;i<=n;i++)
		{
			if(ecc[i]==0)
			{
				k = K;
				vector<int> comp;
				unordered_map<int,int> frontier;
				unordered_map<int,int>::iterator it;
				unordered_map<int,int> d,d1;

				frontier = bfs(i,ind);

				for(it=frontier.begin();it!=frontier.end();++it)
				{
					comp.push_back(it->first);
					d[it->first]=0;
				}
				if(comp.size()<k)
					k=comp.size();
				for(j=0;j<k;j++)
				{
					q=rand()%comp.size();
					t = comp[j];
					comp[j]=comp[q];
					comp[q] = t;
				}

				vector<int> S,temp,S1;
				for(j=0;j<k;j++)
					S.push_back(comp[j]);

				kbfs(S,ind,max_threads,d);


				for(it=d.begin();it!=d.end();++it)
				{
					ecc[it->first]=it->second;
				}

				vector<pair<int,int> >newd(d.begin(),d.end());
				sort(newd.begin(),newd.end(),cmpfunc);
				for(int l=0;l<k;l++)
				{
					S1.push_back(newd[l].first);
					d1[newd[l].first]=0;
				}
				kbfs(S1,ind,max_threads,d1);

				for(it=d1.begin();it!=d1.end();++it)
				{
					ecc[it->first]=max(ecc[it->first],it->second);
				}
			}
		}
		auto stop = high_resolution_clock::now();

		auto duration = duration_cast<seconds>(stop - start);
		cout<<" K is : "<< K << endl;
		cout << "Time taken by k-bfs on cit-Patents graph : " << duration.count() << " seconds" << endl;
	}
	return 0;
}
