#include<cuda.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

using namespace std;

//For nodes in adjacency list
typedef struct node {
	int val;
	struct node* next;
}node;

//Stores visit array's old and new values
typedef struct node1
{
	int oldval,newval;
}node1;


//Compare function to sort based on decreasing order of oldvalue
int cmpfunc(const void* a,const void* b)
{
	node1 x,y;
	x = *(node1*)a;
	y = *(node1*)b;
	if(x.oldval<y.oldval)
		return 1;
	else if(x.oldval==y.oldval)
		return 0;
	else
		return -1;
}


//Function to update depth of nodes in next level of k-bfs
__global__
void Updatenextlevel(int *d_g_edges,int *d_g_edgepos,int *d_depth,node1 *d_visit,int *d_n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	int i, j;
	for(i=index;i<d_n[0];i+=stride)
	{
		if(d_visit[i].oldval==d_depth[0])
		{
			for(j = d_g_edgepos[i];j<d_g_edgepos[i+1];j++)
			{
				if(d_visit[d_g_edges[j]].oldval==-1)
					d_visit[d_g_edges[j]].newval = d_depth[0]+1;
			}
		}
	}
}

//Function to update visit values after completion of iteration of k-bfs
__global__
void UpdateVisit(node1 *visit,int *d_n,int *d_depth,int *d_sz)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for(int j=index;j<d_n[0];j+=stride)
	{
		if(visit[j].newval==d_depth[0]+1)
			d_sz[0]+=1;
		visit[j].oldval = visit[j].newval;
    visit[j].newval = 0;
	}
}

//Array of vectors to store input graph in host memory
node *head[1157828];

//Inserting edge a->b in graph
void insert(int a,int b)
{
	node* temp;
	temp=(node*)malloc(sizeof(node));
	temp->val = b;
	temp->next = head[a];
	head[a] = temp;
}

int main()
{
	int n,m,i,a,b,start,end,j,k,q,t,K,num_blocks=1,num_threads=64;
	scanf("%d %d",&n,&m);

  //Initialising all the edge lists as NULL
	for(i=0;i<n+1;i++)
		head[i]=NULL;

  //Scanning Input graph
	for(i=0;i<m;i++)
	{
		scanf("%d %d",&a,&b);
		insert(a,b);
		insert(b,a);
	}

  //Value of K in K-bfs algorithm
	K = 32;

  /*G_edges is used to store edge end points in CSR format (stored in host memory)
  G_edgepos is used to store number of edges from vertex (stored in host memory)
  G_ecc to store our approximated eccentricity (stored in host memory)
  visit to store whether a vertex is visited (stored in host memory)*/
	int *G_edges,*G_edgepos,*G_ecc;
	node1 *visit;
	G_ecc = (int*)malloc((n+1)*sizeof(int));
	visit = (node1*)malloc((n+1)*sizeof(node1));
	G_edges = (int*)malloc((2*m)*sizeof(int));
	G_edgepos=(int*) malloc((n+2)*sizeof(int));
	G_edgepos[0] = 0 ;
	G_edgepos[1] = 0 ;

  //Converting graph to CSR format
	j = 0;
	for(i=1;i<=n;i++)
	{
		G_ecc[i] = 0;
		visit[i].oldval = -1;
		visit[i].newval = -1;
		node *temp = head[i];
		start = j;
		while(temp!=NULL)
		{
			G_edges[j] = temp->val;
			j++;
			temp=temp->next;
		}
		end = j-1;
		G_edgepos[i+1] = G_edgepos[i]+(end-start+1);
	}


  /*d_g_edges is used to store edge end points in CSR format (stored in device memory)
  d_g_edgepos is used to store number of edges from vertex (stored in device memory)
  d_g_ecc to store our approximated eccentricity (stored in device memory)
  d_visit to store whether a vertex is visited (stored in device memory)*/
	int *d_g_edges,*d_g_edgepos,*d_sz,*d_n;
	node1 *d_visit;

  /*Memory Allocation for variables in device memory and copying corresponding variables from
    host memory to device memory  */
	cudaMalloc( (void**) &d_g_edges, sizeof(int)*(2*m)) ;
	cudaMemcpy( d_g_edges,G_edges, sizeof(int)*(2*m), cudaMemcpyHostToDevice) ;


	cudaMalloc( (void**) &d_g_edgepos, sizeof(int)*(n+2)) ;
	cudaMemcpy( d_g_edgepos,G_edgepos, sizeof(int)*(n+2), cudaMemcpyHostToDevice) ;

	cudaMalloc( (void**) &d_visit, sizeof(node1)*(n+1)) ;

	int *d_depth;
	cudaMalloc( (void**) &d_depth, sizeof(int)) ;

	cudaMalloc( (void**) &d_sz, sizeof(int)) ;

	cudaMalloc( (void**) &d_n, sizeof(int)) ;
	cudaMemcpy( d_n, &n, sizeof(int), cudaMemcpyHostToDevice) ;


	for(i=1;i<=n;i++)
	{
    /*If vertex is not visited yet, then approximate the eccentricity
     of all vertices to which this vertex belongs.This vertex is taken as source vertex */
		if(G_ecc[i]==0)
		{
      //If a vertex is isolated vertex
			if(head[i]==NULL)
				continue;

      //Initialising value of K and visit array (host memory) before bfs
			k = K;
			for(j=1;j<=n;j++)
			{
				visit[j].oldval=-1;
				visit[j].newval=-1;
			}
      //Mark visit of first vertex found with zero approximated eccentricity as zero
			visit[i].oldval=0;

      //Copying host visit array to device visit array
			cudaMemcpy(d_visit,visit, sizeof(node1)*(n+1), cudaMemcpyHostToDevice) ;
			int sz=1,depth=0,comp_size=0;
			cudaMemcpy(d_depth, &depth, sizeof(int), cudaMemcpyHostToDevice) ;


      //Loop runs bfs on source vertex
      //The condition in while loop means we will stop when there are no nodes in current level of bfs
			while(sz>0)
			{
				comp_size+=sz;

        //Update next level in k-bfs
				Updatenextlevel<<<num_blocks,num_threads>>>(d_g_edges,d_g_edgepos,d_depth,d_visit,d_n);
				sz=0;

        //Copy size variable (which is 0) from host to device
				cudaMemcpy(d_sz, &sz, sizeof(int), cudaMemcpyHostToDevice) ;

        //Update visit array in device memory
				UpdateVisit<<<num_blocks,num_threads>>>(d_visit,d_n,d_depth,d_sz);

        //Increase the depth (host memory) in bfs
				depth++;

        //Copy depth variable in host memory to device memory
				cudaMemcpy(d_depth, &depth, sizeof(int), cudaMemcpyHostToDevice) ;

        //Copy the number of nodes in current level from device memory to host memory
				cudaMemcpy(&sz, d_sz, sizeof(int), cudaMemcpyDeviceToHost) ;
			}

      //Copying the visit array which has distances from the source back to host memory
			int *comp;
			cudaMemcpy(visit,d_visit, (n+1)*sizeof(node1), cudaMemcpyDeviceToHost) ;
			int l=0;

      //Getting number of nodes in current component
			for(j=1;j<=n;j++)
			{
				if(visit[j].oldval!=-1)
				{
					l++;
				}
			}

      //Adding values of nodes in current component to comp array (host memory)
			comp_size = l;
			comp = (int*)malloc(l*sizeof(int));
			l = 0;
			for(j=1;j<=n;j++)
			{
				if(visit[j].oldval!=-1)
				{
					comp[l]=j;
					l++;
				}
			}

      //If component size is less than k then k is changed to component size
			if(comp_size<k);
			k=comp_size;

      //Selecting k random nodes from component array
			for(j=0;j<k;j++)
			{
				q=rand()%comp_size;
				t = comp[j];
				comp[j]=comp[q];
				comp[q] = t;
			}

      //Initialise visit for all vertices as -1
			for(j=1;j<=n;j++)
			{
				visit[j].oldval = -1;
				visit[j].newval = -1;
			}

      //Mark visit for all vertices in currentcomponent as 0
			for(j=0;j<k;j++)
			{
				visit[comp[j]].oldval = 0;
				visit[comp[j]].newval = 0;
			}

      //Copy visit array from device memory to host memory
			cudaMemcpy(d_visit,visit, sizeof(node1)*(n+1), cudaMemcpyHostToDevice) ;

      //Initialise number of nodes in 1st level as k and their depth as 0
			sz=k,depth=0;

      //Copy depth variable from host memory to device memory
			cudaMemcpy(d_depth,&depth, sizeof(int), cudaMemcpyHostToDevice) ;

      //Running bfs with above selected k nodes in first level
			while(sz>0)
			{
        //Update next level in k-bfs
				Updatenextlevel<<<num_blocks,num_threads>>>(d_g_edges,d_g_edgepos,d_depth,d_visit,d_n);
				sz=0;

        //Copy size variable (which is 0) from host to device
				cudaMemcpy(d_sz,&sz, sizeof(int), cudaMemcpyHostToDevice) ;

        //Update visit array in device memory
				UpdateVisit<<<num_blocks,num_threads>>>(d_visit,d_n,d_depth,d_sz);

        //Increase the depth (host memory) in bfs
				depth++;

        //Copy depth variable in host memory to device memory
				cudaMemcpy(d_depth,&depth, sizeof(int), cudaMemcpyHostToDevice) ;

        //Copy the number of nodes in current level from device memory to host memory
				cudaMemcpy(&sz,d_sz, sizeof(int), cudaMemcpyDeviceToHost) ;
			}

      //Copying the visit array which has distances from the source back to host memory
			cudaMemcpy(visit,d_visit, (n+1)*sizeof(node1), cudaMemcpyDeviceToHost) ;

      //Update the eccentricities of nodes in this component based on visit array
			for(j=0;j<comp_size;j++)
			{
				G_ecc[comp[j]] = visit[comp[j]].oldval;
			}

      //newd array of type struct node1 which stores depth along with node value
			node1* newd;
			newd = (node1*)malloc(sizeof(node1));
			for(j=0;j<comp_size;j++)
			{
				newd[j].oldval = visit[comp[j]].oldval;
				newd[j].newval = comp[j];
			}

      //Sort newd array based on decreasing order of depth
			qsort(newd,comp_size,sizeof(node1),cmpfunc);

      //Initialise visit for all vertices as -1
			for(j=1;j<=n;j++)
			{
				visit[j].oldval=-1;
				visit[j].newval=-1;
			}

      //Mark visit for all vertices in currentcomponent as 0
			for(j=0;j<k;j++)
			{
				visit[newd[j].newval].oldval=0;
				visit[newd[j].newval].newval=0;
			}
      //Copy visit array from device memory to host memory
      cudaMemcpy(d_visit,visit, sizeof(node1)*(n+1), cudaMemcpyHostToDevice) ;

      //Initialise number of nodes in 1st level as k and their depth as 0
      sz=k,depth=0;

      //Copy depth variable from host memory to device memory
      cudaMemcpy(d_depth,&depth, sizeof(int), cudaMemcpyHostToDevice) ;

      //Running bfs with above selected k nodes in first level
      while(sz>0)
      {
        //Update next level in k-bfs
        Updatenextlevel<<<num_blocks,num_threads>>>(d_g_edges,d_g_edgepos,d_depth,d_visit,d_n);
        sz=0;

        //Copy size variable (which is 0) from host to device
        cudaMemcpy(d_sz,&sz, sizeof(int), cudaMemcpyHostToDevice) ;

        //Update visit array in device memory
        UpdateVisit<<<num_blocks,num_threads>>>(d_visit,d_n,d_depth,d_sz);

        //Increase the depth (host memory) in bfs
        depth++;

        //Copy depth variable in host memory to device memory
        cudaMemcpy(d_depth,&depth, sizeof(int), cudaMemcpyHostToDevice) ;

        //Copy the number of nodes in current level from device memory to host memory
        cudaMemcpy(&sz,d_sz, sizeof(int), cudaMemcpyDeviceToHost) ;
      }

      //Copying the visit array which has distances from the source back to host memory
      cudaMemcpy(visit,d_visit, (n+1)*sizeof(node1), cudaMemcpyDeviceToHost) ;

      /*Compare the value in visit array(depth)
      with previous approximated eccentricity value and update it if it is more*/
			for(j=0;j<comp_size;j++)
			{
				if(visit[comp[j]].oldval>G_ecc[comp[j]])
					G_ecc[comp[j]] = visit[comp[j]].oldval;
			}

		}

	}
	return 0;
}
