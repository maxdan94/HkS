/*
Maximilien Danisch
April 2016
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

to compile:
gcc improvegreedy.c -o improvegreedy -O9

to execute:
./improvegreedy net.txt init.txt res.txt
net.txt should contain the adjacency list: "source targe weight"
init.txt should contain the initial subgraph in the format: "weight k nodeID1 nodeID2 nodeID3 ... nodeIDk"
Will print in res.txt "weight k nodeID1 nodeID2 nodeID3 ... nodeIDk"
Where weight is the sum of the weights of the subgraph induced on the k nodes.

Info:
This local search algorithm improves a solution by switching a node inside the subgraph of size k and a node outside the subgraph such that the weight of the subgraph increases. It stops when no the solution becomes a stable local optimum with respect to this switch operation. It picks the best switch everytime (the one leading the the largest increase in the weight of the subgraph).

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define NLINKS 1000000000 //maximum number of links, will increase if needed

// graph datastructure:

typedef struct {
	unsigned s;//source node
	unsigned t;//target node
	double w;//edge weight
} edge;

//sparse graph structure
typedef struct {
	unsigned n;	//dimensions of the squared matrix
	unsigned e;	//number of edges (nonzero entries simetric)
	unsigned *cd;	//cumulative degree
	unsigned *a; //list of neighbors
	double *w; //weights
	edge *el;//edge list
	unsigned *new;//maping: new[oldID] = newID;
	unsigned *old;//maping: old[newID] = oldID;
} sparse;

typedef struct {
	unsigned n;//number of nodes in the subgraph
	unsigned *l;//list of the nodes in the subgraph
	bool *in;//boolean table in[i]=1 iff i is the subgraph

	unsigned nb;//number of nodes at the border of the subgraph
	unsigned *lb;//list of the nodes at the border of the subgraph
	unsigned *bor;////bor[i]=number of neighbors in the subgraph

	double wt;//sum of the weights in the subgraph
	double *wd;//wd[i]=induced weighted degree of node l[i]

} subgraph;

//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the weighted edgelist
sparse* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	g->el=malloc(e1*sizeof(edge));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	while (fscanf(file,"%u %u %le", &(g->el[g->e].s), &(g->el[g->e].t), &(g->el[g->e].w))==3) {
		g->n=max3(g->n,g->el[g->e].s,g->el[g->e].t);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->el=realloc(g->el,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->el=realloc(g->el,g->e*sizeof(edge));

	return g;
}

//relabel nodes from 0 to n
void relabel(sparse *g) {
	unsigned i,j;

	g->new=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++) {
		g->new[i]=g->n;
	}
	g->old=malloc(g->n*sizeof(unsigned));
	j=0;
	for (i=0;i<g->e;i++) {
		if (g->new[g->el[i].s]==g->n){
			g->new[g->el[i].s]=j;
			g->old[j++]=g->el[i].s;
		}
		if (g->new[g->el[i].t]==g->n){
			g->new[g->el[i].t]=j;
			g->old[j++]=g->el[i].t;
		}
		g->el[i].s=g->new[g->el[i].s];
		g->el[i].t=g->new[g->el[i].t];
	}
	g->n=j;
	g->old=realloc(g->old,g->n*sizeof(unsigned));
}

//build the weighted graph structure
void mkgraph(sparse *g){
	unsigned *tmp;
	unsigned i;

	tmp=calloc(g->n,sizeof(unsigned));
	for (i=0;i<g->e;i++) {
		tmp[g->el[i].s]++;
		tmp[g->el[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+tmp[i-1];
		tmp[i-1]=0;
	}

	g->a=calloc(2*g->e,sizeof(unsigned));
	g->w=calloc(2*g->e,sizeof(double));
	for (i=0;i<g->e;i++) {
		g->a[ g->cd[g->el[i].s] + tmp[g->el[i].s] ]=g->el[i].t;
		g->a[ g->cd[g->el[i].t] + tmp[g->el[i].t] ]=g->el[i].s;
		g->w[ g->cd[g->el[i].s] + tmp[g->el[i].s]++ ]=g->el[i].w;
		g->w[ g->cd[g->el[i].t] + tmp[g->el[i].t]++ ]=g->el[i].w;
	}
	free(tmp);
}


//initialise the subgraph with theinput subgraph
subgraph *initsubgraph(sparse *g,char *inputsg){
	unsigned i,j,node,neigh;
	subgraph *sg=malloc(sizeof(subgraph));
	FILE* file=fopen(inputsg,"r");

	fscanf(file,"%le",&(sg->wt));
	fscanf(file," %u",&(sg->n));
	sg->l=malloc((sg->n+1)*sizeof(unsigned));//one more empty spot used in addnode()
	sg->in=calloc(g->n,sizeof(bool));
	for (i=0;i<sg->n;i++){
		fscanf(file," %u",sg->l+i);
		sg->l[i]=g->new[sg->l[i]];
		sg->in[sg->l[i]]=1;
	}
	fclose(file);

	sg->nb=0;
	sg->lb=malloc((g->n-sg->n)*sizeof(unsigned));
	sg->bor=calloc(g->n,sizeof(unsigned));
	sg->wd=calloc(g->n,sizeof(double));
	for (i=0;i<sg->n;i++){
		node=sg->l[i];
		for (j=g->cd[node];j<g->cd[node+1];j++){
			neigh=g->a[j];
			sg->wd[neigh]+=g->w[j];
			if ((sg->bor[neigh]++==0) && (sg->in[neigh]==0)){
				sg->lb[sg->nb++]=neigh;
			}
		}
	}
	return sg;
}

double bestswich(sparse* g, subgraph* sg, unsigned *bi, unsigned *bj){
	unsigned i,j;
	unsigned n1,n2;
	static double *delta=NULL;
	double max=0,tmp;

	if (delta==NULL){
		delta=calloc(g->n,sizeof(double));
	}

	for (i=0;i<sg->n;i++){
		n1=sg->l[i];
		for (j=g->cd[n1];j<g->cd[n1+1];j++){
			delta[g->a[j]]=g->w[j];
		}
		for (j=0;j<sg->nb;j++){
			n2=sg->lb[j];
			tmp=sg->wd[n2]-delta[n2]-sg->wd[n1];
			if (tmp>max){
				max=tmp;
				*bi=i;
				*bj=j;
			}
		}
		for (j=g->cd[n1];j<g->cd[n1+1];j++){
			delta[g->a[j]]=0;
		}
	}

	return max;
}

void addnode(sparse* g, subgraph* sg, unsigned j){
	unsigned node,neigh;

	node=sg->lb[j];
	sg->l[sg->n++]=node;
	sg->in[node]=1;
	for (j=g->cd[node];j<g->cd[node+1];j++){
		neigh=g->a[j];
		sg->wd[neigh]+=g->w[j];
		if (sg->in[neigh]){
			sg->wt+=g->w[j];
		}
		if ((sg->bor[neigh]++==0) && (sg->in[neigh]==0)){
			sg->lb[sg->nb++]=neigh;
		}
	}
}

void rmnode(sparse* g, subgraph* sg, unsigned i){
	unsigned node,neigh;

	node=sg->l[i];
	sg->l[i]=sg->l[--sg->n];
	sg->lb[sg->nb++]=node;
	sg->in[node]=0;
	for (i=g->cd[node];i<g->cd[node+1];i++){
		neigh=g->a[i];
		sg->wd[neigh]-=g->w[i];
		if (sg->in[neigh]){
			sg->wt-=g->w[i];
		}
		sg->bor[neigh]--;
	}

	for (i=0;i<sg->nb;i++){
		if ((sg->bor[sg->lb[i]]==0) || (sg->in[sg->lb[i]])){
			sg->lb[i--]=sg->lb[--sg->nb];
		}
	}
}

void optimize(sparse *g,subgraph *sg){
	double delta;
	unsigned i,j,k=0;
	while (1){
		delta=bestswich(g,sg,&i,&j);
		printf("it delta w = %u %le %le\n",++k,delta,sg->wt);
		if (delta==0){
			break;
		}
		addnode(g,sg,j);
		rmnode(g,sg,i);
	}
}

//printing the result in file output: "size node density" on each line
void printres(sparse *g, subgraph *sg, char* output){
	unsigned i;
	FILE *file=fopen(output,"w");
	fprintf(file,"%.10le %u",sg->wt,sg->n);
	for (i=0;i<sg->n;i++){
		fprintf(file," %u",g->old[sg->l[i]]);
	}
	fprintf(file,"\n");
	fclose(file);
}

//free the graph stucture
void freegraph(sparse *g){
	free(g->cd);
	free(g->a);
	free(g->w);
	free(g->el);
	free(g->old);
	free(g->new);
	free(g);
}

//free the graph stucture
void freesg(subgraph *sg){
	free(sg->in);
	free(sg->l);
	free(sg->lb);
	free(sg->bor);
	free(sg->wd);
	free(sg);
}

int main(int argc,char** argv){
	sparse* g;
	subgraph* sg;

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Building the datastructure\n");
	printf("Number of nodes = %u\n",g->n);

	relabel(g);
	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);
	mkgraph(g);
	printf("Initialising the dense subgraph with %s\n",argv[2]);
	sg=initsubgraph(g,argv[2]);
	printf("Optimizing\n");
	optimize(g,sg);
	printf("Printing the result to file %s\n",argv[3]);
	printres(g,sg,argv[3]);
	freegraph(g);
	freesg(sg);
	return 0;
}
