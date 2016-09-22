/*
Maximilien Danisch
April 2016
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

to compile:
gcc wcore.c -o wcore -O9

to execute:
./wcore net.txt res.txt
Will print in res.txt "size nodeID weight" on each line
Where weight is the sum of the weights of the subgraph induced on all previous nodes.

Info:
Compute an ordering of the vertices according to the weighted kcore (or weighted degeneracy ordering) and compute the weight of the subgraphs induced by all prefixes.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define NLINKS 1000000000 //maximum number of links, will increase if needed


//heap data structure:

typedef struct {
	unsigned key;
	double val;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max){
	unsigned i;
	bheap *heap=malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=malloc(n_max*sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap,unsigned i, unsigned j) {
	keyvalue kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].val>heap->kv[i].val) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].val<heap->kv[j1].val) ) ? j2 : j1 ;
		if (heap->kv[j].val < heap->kv[i].val) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

void update(bheap *heap,unsigned key,double value){
	unsigned i=heap->pt[key];
	if (i!=-1){
		(heap->kv[i]).val-=value;
		bubble_up(heap,i);
	}
}

keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}

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
	double *wd;//weighted degree of each node
	double tw;//total weight
	unsigned *map;//maping: map[newID] = oldID;
} sparse;

typedef struct {
	unsigned id;
	unsigned idr;
	double val;
} idval;


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
	unsigned *newlabel;

	newlabel=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++) {
		newlabel[i]=g->n;
	}
	g->map=malloc(g->n*sizeof(unsigned));
	j=0;
	for (i=0;i<g->e;i++) {
		if (newlabel[g->el[i].s]==g->n){
			newlabel[g->el[i].s]=j;
			g->map[j++]=g->el[i].s;
		}
		if (newlabel[g->el[i].t]==g->n){
			newlabel[g->el[i].t]=j;
			g->map[j++]=g->el[i].t;
		}
		g->el[i].s=newlabel[g->el[i].s];
		g->el[i].t=newlabel[g->el[i].t];
	}
	g->n=j;
	free(newlabel);
	g->map=realloc(g->map,g->n*sizeof(unsigned));
}

//build the weighted graph structure
void mkgraph(sparse *g){
	unsigned *tmp;
	unsigned i;

	tmp=calloc(g->n,sizeof(unsigned));
	g->wd=calloc(g->n,sizeof(double));
	g->tw=0;
	for (i=0;i<g->e;i++) {
		tmp[g->el[i].s]++;
		tmp[g->el[i].t]++;
		g->wd[g->el[i].s]+=g->el[i].w;
		g->wd[g->el[i].t]+=g->el[i].w;
		g->tw+=g->el[i].w;
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

//Building the heap structure
bheap* mkheap(sparse *g){
	unsigned i;
	keyvalue kv;
	bheap* heap=construct(g->n);
	for (i=0;i<g->n;i++){
		kv.key=i;
		kv.val=g->wd[i];
		insert(heap,kv);
	}
	return heap;
}

//compute the ranking and edges weights
keyvalue* mkdense(sparse* g,bheap* heap){
	unsigned i,j,k;
	keyvalue kv;
	unsigned e=g->e;
	double ew=0;
	keyvalue *rank=malloc(g->n*sizeof(keyvalue));

	ew=g->tw;
	rank[g->n-1].val=ew;
	for (i=g->n-1;i>0;i--){
		kv=popmin(heap);
		ew-=kv.val;
		rank[i].key=kv.key;
		rank[i-1].val=ew;
		for (j=g->cd[kv.key];j<g->cd[kv.key+1];j++){
			k=g->a[j];
			update(heap,k,g->w[j]);
		}
	}
	kv=popmin(heap);
	rank[0].key=kv.key;
	return rank;
}

//printing the result in file output: "size node density" on each line
void printres(keyvalue* rank, sparse *g, char* output){
	unsigned i;
	FILE *file=fopen(output,"w");
	for (i=0;i<g->n;i++){
		fprintf(file,"%u %u %.10le\n",i+1,g->map[rank[i].key],rank[i].val);
	}
	fclose(file);
}

void freeheap(bheap *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//free the graph stucture
void freegraph(sparse *g){
	free(g->cd);
	free(g->a);
	free(g->w);
	free(g->el);
	free(g->map);
	free(g);
}

int main(int argc,char** argv){
	sparse* g;
	double *vect;
	keyvalue* rank;
	bheap* heap;

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Building the datastructure\n");
	relabel(g);
	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);
	mkgraph(g);
	printf("Computing subgraphs of high weight using the weighted core ordering\n");
	heap=mkheap(g);
	rank=mkdense(g,heap);
	free(heap);
	printf("Printing result in file %s\n",argv[2]);
	printres(rank,g,argv[2]);
	free(rank);
	freegraph(g);
	return 0;
}
