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
	scanf("%d %d",&n,&m);
	ind.resize(n+1);

	//	n=1157827;
	//	m=2987624;

	cout<<"No of vertices : "<<n<<" and No of edges : "<<m<<endl;

	for(i=0;i<m;i++)
	{
		scanf("%d	%d",&a,&b);
		ind[a].push_back(b);
		ind[b].push_back(a);
	}
	cout<<"scanned\n";
	ecc.clear();
	ecc.resize(n+1);
	for(i=1 ; i<=n ; i++)
		ecc[i] = 0;
	K = 256;
	cout<<"started\n";
	auto start = high_resolution_clock::now();
	for(i=1;i<=n;i++)
	{
		if(ecc[i]==0)
		{
			k = K;
			vector<int> comp;
			unordered_map<int,int> frontier;
			unordered_map<int,int>::iterator it;
			unordered_map<int,int> d;

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

			int num = frontier.size();

			int tasks = k/max_threads+(int)(k%max_threads!=0);
			int subtasks = (num)/max_threads+int(num%max_threads!=0);
			vector<unordered_map<int,int> > temps(max_threads);

			cout << "first " <<tasks<< endl;
			cout << "first " <<subtasks<< endl;
			for(int x=0;x<tasks;x++){
				cout<<x<<endl;
				#pragma omp parallel private(tid)
				{
					tid = omp_get_thread_num();
					int y = x*max_threads+tid;
					if(y<k)
					{
						temps[tid] = bfs(S[y],ind);
					}
				}
				#pragma omp barrier 
				int garb = max_threads;
				if(x+1 == tasks)
					garb = k%max_threads;
//				cout << "first_sub" <<endl;
				for(int y=0;y<subtasks;y++){
					cout<<"subtask "<<y<<endl;
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num();
						int z = y*max_threads+tid;
						if(z<num){
							for(int p=0;p<garb;p++)
							{
//								if(comp[z]>=d.size())
//									cout<<"1."<<d.size()<<" "<<comp[z]<<endl;
//								if(comp[z]>=temps[p].size())
//									cout<<"2."<<temps[p].size()<<" "<<comp[z]<<endl;
								d[comp[z]] = max(temps[p][comp[z]],d[comp[z]]);
							}
						}
					}
					#pragma omp barrier 
				}
//				cout << "sub_end" << endl;
			}
			vector<pair<int,int> >newd(d.begin(),d.end());
			sort(newd.begin(),newd.end(),cmpfunc);
//			cout<<newd.size()<<endl;
			for(int l=0;l<k;l++)
			{
				S1.push_back(newd[l].first);
			}
		//	cout << "second" << endl;

			cout << "second " <<tasks<< endl;
			cout << "second " <<subtasks<< endl;
			for(int x=0;x<tasks;x++){
				cout<<x<<endl;
				#pragma omp parallel private(tid)
				{
					tid = omp_get_thread_num();
					int y = x*max_threads+tid;
					if(y<k)
					{
						temps[tid] = bfs(S1[y],ind);
					}
				}
				int garb = max_threads;
				if(x+1 == tasks)
					garb = k%max_threads;
				for(int y=0;y<subtasks;y++){
					cout<<"subtask "<<y<<endl;
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num();
						int z = y*max_threads+tid;
						if(z<num){
							for(int p=0;p<garb;p++)
							{
								d[comp[z]] = max(temps[p][comp[z]],d[comp[z]]);
							}
						}
					}
				}
			}


			for(it=d.begin();it!=d.end();++it)
			{
				ecc[it->first]=it->second;
			}
		}
	}
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout<<" K is :"<<K<<endl;
		cout << "Time taken by k-bfs " << duration.count() << " seconds" << endl;
	return 0;
}
