
#include<thrust/scan.h>
#include<thrust/sort.h>
#include <thrust/unique.h>
#include<thrust/copy.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<vector>
#include<set>
#include<algorithm>
#include<numeric>
#include<set>
#include<iostream>
#include<stdio.h>
#include<fstream>
#include <sys/time.h>
#include<cuda_runtime.h>
using namespace std;


//structure which stores source vertex,target vertex and weight....

struct edgeingraph{
    int u;
    int v;
    int w;
};
bool comp(edgeingraph o1, edgeingraph o2)
    {
  
        if(o1.u==o2.u){
          return (o1.w<o2.w);
        }
        return (o1.u<o2.u);
    }
//set is just for finding number of vertices....
set<int> SET;
int verticescount_n=0;
int edgescount_m=0;
vector<struct edgeingraph> edgeStream;


//this function will read the input from file
//one edge is added twice
//example if 0 1 2 is taken 
//then 1 0 2 is also added
//making it easier to generate csr representation

void readEdgeStream(){
    ifstream fin;
    fin.open("edgesinput.txt"); //can give your input file name here
    struct edgeingraph edge;
    struct edgeingraph edge2;
        fin>>edge.u;
        SET.insert(edge.u);
        fin>>edge.v;
        SET.insert(edge.v);
        fin>>edge.w;
        edgescount_m++;
        edge2.u=edge.v;
        edge2.v=edge.u;
        edge2.w=edge.w;
        if(edge.u!=edge.v){
        edgeStream.push_back(edge);
        edgeStream.push_back(edge2);
        }
    while (fin) {
      struct edgeingraph edge1;
      struct edgeingraph edge3;
        fin>>edge1.u;
        if( fin.eof() ) break;
        SET.insert(edge1.u);
        fin>>edge1.v;
        SET.insert(edge1.v);
        fin>>edge1.w;
        edgescount_m++;
        edge3.u=edge1.v;
        edge3.v=edge1.u;
        edge3.w=edge1.w;
        if(edge1.u!=edge1.v){
        edgeStream.push_back(edge1);
        edgeStream.push_back(edge3);
        }
    }
    fin.close();
}

//will create and sorts edgestreamarray
//from edgeStream vector

void sortEdgeStream(struct edgeingraph edgestreamarray[]){
    int *keys=new int[edgescount_m*2];
    for(int i=0;i<2*edgescount_m;i++){
      keys[i]=static_cast<struct edgeingraph>(edgeStream[i]).u;
      edgestreamarray[i].u=edgeStream[i].u;
      edgestreamarray[i].v=edgeStream[i].v;
      edgestreamarray[i].w=edgeStream[i].w;
    }
    thrust::sort_by_key(keys,keys+(2*edgescount_m),edgestreamarray);
    delete[] keys;
}



void edgeStreamToCSR(int edgesarray[],int verticesarray[], int weightsarray[],struct edgeingraph edgestreamarray[]){
  for(int i=0;i<2*edgescount_m;i++){
        edgesarray[i]=edgestreamarray[i].v;
  }

  for(int i=0;i<2*edgescount_m;i++){
        weightsarray[i]=edgestreamarray[i].w;
  }
  
  int it=0;
  for(int i=0;i<2*edgescount_m;i++){
   if(edgestreamarray[i].u==it){
     verticesarray[it]=i;
     it++;
   }     
  }
  verticesarray[verticescount_n]=2*edgescount_m;

}



__global__ void createSusscessorArray(int *su,int *w,int *ed,int *v,int n,int m)
{
    int tid=(blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid<n&& tid>=0){
    int minindex=v[tid];
    for(int i=v[tid]+1;i<v[tid+1];i++){
        if(w[i]<w[minindex]){
          minindex=i;
        }
    }
    su[tid]=ed[minindex];
    }
}

__global__ void createSuperVertices(int *su,int n)
{
    int tid=(blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid<n&& tid>=0){
    while(su[tid]!=su[su[tid]]){
      su[tid]=su[su[tid]];
    }
    }
}
__global__ void createnewedgestream(int *s,edgeingraph *edsa, int n,int m,int *ts,int *ind)
{
    int tid=(blockIdx.x * blockDim.x) + threadIdx.x;
    if(tid<2*m && tid>=0){
      (edsa[tid]).u=ind[s[(edsa[tid]).u]];
      (edsa[tid]).v=ind[s[(edsa[tid]).v]];
    }
}

int main(){

float totaltime=0;
  cout<<"EDGES OF MINIMUM SPANNING TREE: "<<endl;
  int cost=0;
  int count=0;
  readEdgeStream();


  //cleans multigraph as simple graph by deleting larger weights and keeping small weight between two vertices.

  for (auto i = edgeStream.begin(); i != edgeStream.end(); ++i) {
    for(auto j=i+1;j != edgeStream.end();++j){
      if((*j).u==(*i).u && (*j).v == (*i).v){
        if((*j).w>=(*i).w){
          edgeStream.erase(j);
          j--;
        }
        else{
             edgeStream.erase(i);
             i--;break;
        }
      }
    }  
}
edgescount_m=edgeStream.size()/2;

  int f=SET.size();
  int flag[f];
  for(int i=0;i<f;i++){
    flag[i]=i;
  }
  vector<struct edgeingraph> backup;
  backup=edgeStream;



while(count<f-1){
  verticescount_n=SET.size();
  struct edgeingraph edgestreamarray[2*edgescount_m];
  sortEdgeStream(edgestreamarray);
  sort(edgestreamarray,edgestreamarray+2*edgescount_m,comp);
  int *edgesarray=new int[2*edgescount_m];
  int *weightsarray=new int[2*edgescount_m];
  int *verticesarray=new int[verticescount_n+1];
  edgeStreamToCSR(edgesarray,verticesarray,weightsarray,edgestreamarray);

  int *successorarray=new int[verticescount_n];
    int *v;
    cudaMalloc((void **)&v, (verticescount_n+1)*sizeof(int));
    int *ed;
    cudaMalloc((void **)&ed, (edgescount_m*2)*sizeof(int));
    int *w;
    cudaMalloc((void **)&w, (edgescount_m*2)*sizeof(int));
    int *su;
    cudaMalloc((void **)&su, (verticescount_n+1)*sizeof(int));
    

  cudaMemcpy(v, verticesarray, (verticescount_n+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(ed, edgesarray, (edgescount_m*2)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(w, weightsarray, (edgescount_m*2)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(su, successorarray, (verticescount_n)*sizeof(int), cudaMemcpyHostToDevice);
  int th=(verticescount_n+1023)/1024;
  float time=0;
  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
  createSusscessorArray<<<th, 1024>>>(su,w,ed,v,verticescount_n,edgescount_m);
cudaEventRecord(stop,0);
  cudaDeviceSynchronize();
  cudaEventElapsedTime(&time,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  totaltime+=time;
  cout<<endl<<endl<<time<<endl<<endl;
  cudaMemcpy(successorarray,su,(verticescount_n)*sizeof(int), cudaMemcpyDeviceToHost);
cudaFree(v);
cudaFree(ed);
cudaFree(w);


 for(int i=0;i<verticescount_n;i++){
   if(successorarray[successorarray[i]]==i){
     successorarray[successorarray[i]]=successorarray[i];
   }
 }
 for(int i=0;i<verticescount_n;i++){
   if(successorarray[i]!=i){
     int weight=0;
     for(int j=verticesarray[i];j<verticesarray[i+1];j++){
       if(edgesarray[j]==successorarray[i]){
         weight=weightsarray[j];
       }
     }
    for(int k=0;k<backup.size();k++){
        if(backup[k].w==weight){
          if((flag[backup[k].u]==i && flag[backup[k].v]==successorarray[i]) || (flag[backup[k].v]==i && flag[backup[k].u]==successorarray[i])){
            cout<<backup[k].u<<" "<<backup[k].v<<" "<<weight<<endl;
            break;
          }
        }
    }
    cost+=weight;
    count++;
   }
 }

  cudaMemcpy(su, successorarray, (verticescount_n)*sizeof(int), cudaMemcpyHostToDevice);
  time=0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
  createSuperVertices<<<th, 1024>>>(su,verticescount_n);
  cudaEventRecord(stop,0);
  cudaDeviceSynchronize();
  cudaEventElapsedTime(&time,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  totaltime+=time;
  cout<<endl<<endl<<time<<endl<<endl;
  cudaMemcpy(successorarray,su,(verticescount_n)*sizeof(int), cudaMemcpyDeviceToHost);
cudaFree(su);
//relabeling

int *tempsuc=new int[verticescount_n];
thrust::copy(successorarray,successorarray+verticescount_n,tempsuc);
thrust::sort(tempsuc,tempsuc+verticescount_n);
thrust::unique(tempsuc,tempsuc+verticescount_n);

int index[verticescount_n]={-1};
fill(index,index+verticescount_n,-1);
for(int i=0;i<verticescount_n;i++){
  if(index[tempsuc[i]]==-1){
  index[tempsuc[i]]=i;
  }
}
for(int i=0;i<f;i++){

  int temp=flag[i];
   
  flag[i]=index[successorarray[temp]];

}
int *ind;
cudaMalloc((void **)&ind, verticescount_n*sizeof(int));
int *ts;
cudaMalloc((void **)&ts, sizeof(tempsuc));
int *s;
cudaMalloc((void **)&s, (verticescount_n)*sizeof(int));
struct edgeingraph *edsa;
cudaMalloc((void **)&edsa, (2*edgescount_m)*sizeof(edgeingraph));
th=((2*edgescount_m)+1023)/1024;
cudaMemcpy(ind,index,verticescount_n*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(ts,tempsuc,sizeof(tempsuc),cudaMemcpyHostToDevice);
cudaMemcpy(edsa, edgestreamarray, (edgescount_m*2)*sizeof(edgeingraph), cudaMemcpyHostToDevice);
cudaMemcpy(s, successorarray, (verticescount_n)*sizeof(int), cudaMemcpyHostToDevice);
 time=0; 
 /*cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);*/
createnewedgestream<<<th, 1024>>>(s,edsa,verticescount_n,edgescount_m,ts,ind);
/*cudaEventRecord(stop,0);
  cudaDeviceSynchronize();
  cudaEventElapsedTime(&time,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  totaltime+=time;*/
  //cout<<endl<<endl<<time<<endl<<endl;
cudaMemcpy(edgestreamarray, edsa,(edgescount_m*2)*sizeof(edgeingraph), cudaMemcpyDeviceToHost);
cudaFree(s);
cudaFree(edsa);
cudaFree(ts);
cudaFree(ts);
edgeStream.clear();
SET.clear();
edgeStream.insert(edgeStream.begin(),edgestreamarray,edgestreamarray+edgescount_m*2);
 for (auto i = edgeStream.begin(); i != edgeStream.end(); ++i) {
     SET.insert((*i).u);
     SET.insert((*i).v);
        if ((*i).u == (*i).v) {
            edgeStream.erase(i);
            i--;
        }
}
 for (auto i = edgeStream.begin(); i != edgeStream.end(); ++i) {
    for(auto j=i+1;j != edgeStream.end();++j){
      if((*j).u==(*i).u && (*j).v == (*i).v){
        if((*j).w>=(*i).w){
          edgeStream.erase(j);
          j--;
        }
        else{
             edgeStream.erase(i);
             i--;break;
        }
      }
    }  
}
edgescount_m=edgeStream.size()/2;
}
cout<<endl<<"TOTAL COST = "<<cost<<endl;
cout<<totaltime<<endl;
}
