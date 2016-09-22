/*
Maximilien Danisch
Adapted from the java version of Manthos Letsios
Technical consultants: Oana Balalau, Emmanuel Orsini and Mauro Sozio
September 2016
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish. This program finds an alpha approximation of the densest subgraph of k nodes. The algorithm is based on a smart branch and bound: the branching phase is based on deciding whether to add a specific edge (and thus its endpoints) in the corresponding solution or not.

To compile:
gcc bb_dks.cpp -o bb_dks -O9

To execute:
./bb_dks k alpha net.txt
net.txt should contain the weighted graph "nodeID1 nodeID2 weight". It is better if the nodeIDs are from 0 to n-1.

Will print:
- some info
- total_weight
- the k nodes
*/

#define NLINKS 10000000 //Maximum number of links, will increase if needed.
#define NHEAP 10000000 //Maximum number of sets of l edges, will NOT increase.

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
	unsigned s;//source node
	unsigned t;//target node
	double w;//edge weight
} edge;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of nodes
	edge* edges;//edge list
	unsigned* map;//maping: map[newID] = oldID;
} edgelist;

typedef struct {
	unsigned n;//node
	unsigned e;//edge id
	double v;//edge weight
} nodeidval;

//The following structure is used to check adjencency between two nodes
typedef struct {
	unsigned n;//number of nodes
	unsigned *cd;//cd[u+1]=cummulative degree of node u. cd[0]=0
	nodeidval *adj;//adj[cd[u]..cd[u+1]] = list of neighbors of u, weights on edges, and id of edges
} graph;

typedef struct {
	unsigned e;//edge id
	double v;//edge weight
} idval;

typedef struct {
	unsigned n;//node id
	unsigned d;//deg
} nodedeg;

//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//Reads the weighted edgelist from file named input
edgelist* readedgelist(char* input){
	unsigned s=NLINKS;
	edgelist *el=malloc(sizeof(edgelist));
	el->edges=malloc(s*sizeof(edge));
	FILE *file;
	el->n=0;
	el->e=0;
	file=fopen(input,"r");
	while (fscanf(file,"%u %u %le", &(el->edges[el->e].s), &(el->edges[el->e].t), &(el->edges[el->e].w))==3) {
		el->n=max3(el->n,el->edges[el->e].s,el->edges[el->e].t);
		if (el->e++==s) {
			s+=NLINKS;
			el->edges=realloc(el->edges,s*sizeof(edge));
		}
	}
	fclose(file);
	el->n++;
	el->edges=realloc(el->edges,el->e*sizeof(edge));
	return el;
}

void free_edgelist(edgelist* el){
	free(el->edges);
	free(el->map);
	free(el);
}

//For futur use in qsort.
static int compare_edges (void const *a, void const *b){
	edge const *pa = a;
	edge const *pb = b;
	if ((*pa).w<=(*pb).w)
		return 1;
	return -1;
}

//For futur use in qsort.
static int compare_nodedeg (void const *a, void const *b){
	nodedeg const *pa = a;
	nodedeg const *pb = b;
	if ((*pa).d<=(*pb).d)
		return 1;
	return -1;
}

//Sorts edges in non-increasing order of weight, relabels nodes from 0 to n-1 (in non-increasing order of degree) and transforms edges such that sourceID > targetID.
void relabel(edgelist *el) {
	unsigned i,j,u,v;
	unsigned *newlabel=malloc(el->n*sizeof(unsigned));
	nodedeg *nd=malloc(el->n*sizeof(nodedeg));

	//Ordering nodes in non-increasing order of degree
	for (i=0;i<el->n;i++) {
		nd[i].n=i;
		nd[i].d=0;
	}
	for (i=0;i<el->e;i++) {
		nd[el->edges[i].s].d++;
		nd[el->edges[i].t].d++;
	}
	qsort(nd,el->n,sizeof(nodedeg),compare_nodedeg);
	el->map=malloc(el->n*sizeof(unsigned));
	for (i=0;i<el->n;i++) {
		if (nd[i].d==0)
			break;
		el->map[i]=nd[i].n;
		newlabel[nd[i].n]=i;
	}
	j=i;
	for (i=j;i<el->n;i++) {
		newlabel[nd[i].n]=el->n;
	}
	free(nd);
	el->n=j;
	el->map=realloc(el->map,el->n*sizeof(unsigned));

	//Sorting edges in non-increasing order of weight
	qsort(el->edges,el->e,sizeof(edge),compare_edges);

	//Relabeling nodes in edgelist
	for (i=0;i<el->e;i++) {
		u=newlabel[el->edges[i].s];
		v=newlabel[el->edges[i].t];
		if (u>v){
			el->edges[i].s=u;
			el->edges[i].t=v;
		}
		else {
			el->edges[i].s=v;
			el->edges[i].t=u;
		}
	}
	free(newlabel);
}

//For futur use in qsort.
static int compare_nodeidval (void const *a, void const *b){
	nodeidval const *pa = a;
	nodeidval const *pb = b;
	if ((*pa).n>(*pb).n)
		return 1;
	return -1;
}

//Builds the graph structure (note that we have sourceID > tardetID).
graph* mkgraph(edgelist *el) {
	unsigned i,j;
	graph *g=malloc(sizeof(graph));
	g->n=el->n;
	unsigned *d=calloc(g->n,sizeof(unsigned));
	unsigned core=0;

	for (i=0;i<el->e;i++) {
		d[el->edges[i].s]++;
	}
	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		core=(d[i-1]>core)?d[i-1]:core;
		d[i-1]=0;
	}
	printf("Core value of the graph <= %u\n",core);
	g->adj=malloc(el->e*sizeof(nodeidval));
	for (i=0;i<el->e;i++) {
		j=g->cd[el->edges[i].s]+d[el->edges[i].s]++;
		g->adj[j].n=el->edges[i].t;
		g->adj[j].e=i;
		g->adj[j].v=el->edges[i].w;
	}

	for (i=0;i<g->n;i++) {
		qsort(&(g->adj[g->cd[i]]),d[i],sizeof(nodeidval),compare_nodeidval);
	}
	free(d);
	return g;
}

void free_graph(graph* g){
	free(g->cd);
	free(g->adj);
	free(g);
}

//Checks if u and v are adjacent. Returns the ID and the weight of edge (u,v) if they are adjacent and {-1,-1} if they are not. Do binary search in the neighbor list of the node with smaller index. Time O(log(c)) where c is the core value of the graph (a degree ranking is implemented instead of a core ranking, this does not seem to change much the runing time in practice in real world graphs).
inline idval isadj(unsigned u, unsigned v, graph *g){
	unsigned first,middle,last;
	idval res;
	if (v>u){
		unsigned w=u;
		u=v;
		v=w;
	}
	first=g->cd[u];
	last=g->cd[u+1]-1;
	middle=(first+last)/2;
	while (first<=last) {
		if (g->adj[middle].n<v)
			first=middle+1;
		else if (g->adj[middle].n>v)
			last=middle-1;
		else {
			res.e=g->adj[middle].e;
			res.v=g->adj[middle].v;
			return res;
		}
		middle=(first+last)/2;
	}
	res.e=-1;
	res.v=-1.;
	return res;
}

//A node in the branch and bound explored tree
typedef struct {
	unsigned size;//size of the subgraph induced on the selected edges
	double val;//sum of weights of the edges induced on the nodes
	double ub;//upper bound on the weights the subgraph can achieve
	unsigned eid;//number of edges already considered.
	unsigned *nodes;//id of the nodes in the subgraph
} treenode;

//Special heap data structure
typedef struct {
	unsigned n_max;//max number of elements.
	unsigned n;//number of elements.
	treenode **trns;//heap of pointers to elements
	treenode **avail;//list of available pointers
	treenode *elems;//array of elements
	unsigned *nodes;//nodes associated to elements
} bheap;

bheap *construct(unsigned k){
	unsigned i;
	bheap *heap=malloc(sizeof(bheap));
	heap->n_max=NHEAP;
	heap->n=0;
	heap->trns=malloc(NHEAP*sizeof(treenode*));
	heap->avail=malloc(NHEAP*sizeof(treenode*));
	heap->elems=malloc((NHEAP+3)*sizeof(treenode));//3 more treenodes are used: trn, best and child.
	heap->nodes=malloc(k*(NHEAP+3)*sizeof(unsigned));
	for (i=0;i<NHEAP;i++){
		heap->avail[i]=&(heap->elems[i+3]);
		(heap->elems[i+3]).nodes=heap->nodes+k*(i+3);
	}
	(heap->elems[0]).nodes=heap->nodes;
	(heap->elems[1]).nodes=heap->nodes+k;
	(heap->elems[2]).nodes=heap->nodes+2*k;
	return heap;
}

void free_bheap(bheap* heap){
	free(heap->trns);
	free(heap->avail);
	free(heap->elems);
	free(heap->nodes);
	free(heap);
}

inline void swap(bheap *heap,unsigned i, unsigned j) {
	treenode* trn=heap->trns[i];
	heap->trns[i]=heap->trns[j];
	heap->trns[j]=trn;
}

inline void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->trns[j]->val<heap->trns[i]->val) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

inline void bubble_down(bheap *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->trns[j2]->val>heap->trns[j1]->val) ) ? j2 : j1 ;
		if (heap->trns[j]->val>heap->trns[i]->val) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

//Inserts trn in heap and returns an available pointers to a treenode
inline treenode* insert(bheap *heap,treenode *trn){
	if (heap->n==NHEAP){
		printf("Error: NHEAP is too small\n");
		exit(1);
	}
	heap->trns[heap->n]=trn;
	trn=heap->avail[heap->n];
	bubble_up(heap,(heap->n)++);
	return trn;
}

//Returns maximum in heap, should provide an available pointer
inline treenode* popmax(bheap *heap,treenode *trn){
	heap->avail[--(heap->n)]=trn;
	trn=heap->trns[0];
	heap->trns[0]=heap->trns[heap->n];
	bubble_down(heap);
	return trn;
}

//Checks if the 2 nodes are already in the set
inline char edgeinset(edge ed,treenode* trn){
	char c=0;
	unsigned i;
	for (i=0;i<trn->size;i++){//simple linear search, improve to binary search later maybe. trn.s<k and k is small anyway
		if (trn->nodes[i]==ed.s){
			if (c==2){
				return 3;
			}
			c=1;
		}
		else if (trn->nodes[i]==ed.t){
			if (c==1){
				return 3;
			}
			c=2;
		}
	}
	return c;
}

void mkchild1(treenode* trn,treenode* child,edgelist *el,graph *g,unsigned k){
	edge ed=el->edges[trn->eid];
	unsigned i,end,beg;
	double add;
	idval ev;
	char c=edgeinset(ed,trn);

	if (c==3){//both nodes in set
		for (i=0;i<trn->size;i++){
			child->nodes[i]=trn->nodes[i];
		}
		child->size=trn->size;
		child->eid=trn->eid+1;
		child->val=trn->val;
		//computing new upper bound
		child->ub=trn->ub-ed.w;
		end=trn->eid+(k*(k-1))/2-(trn->size*(trn->size-1))/2;
		if (end<el->e){
			child->ub+=el->edges[end].w;
		}
		return;
	}

	if (c==2){//ed.t in set, ed.s not in set
		add=0;
		child->size=trn->size+1;
		child->eid=trn->eid+1;
		for (i=0;i<trn->size;i++){
			child->nodes[i]=trn->nodes[i];
			ev=isadj(trn->nodes[i],ed.s,g);
			if (ev.v>0){
				if (ev.e<trn->eid){
					child->ub=-1;
					return;
				}
				add+=ev.v;
			}
		}
		child->nodes[i]=ed.s;
	}

	else if (c==1){//ed.s in set, ed.t not in set
		add=0.;
		child->size=trn->size+1;
		child->eid=trn->eid+1;
		for (i=0;i<trn->size;i++){
			child->nodes[i]=trn->nodes[i];
			ev=isadj(trn->nodes[i],ed.t,g);
			if (ev.v > 0.){
				if (ev.e<trn->eid){
					child->ub=-1;
					return;
				}
				add+=ev.v;
			}
		}
		child->nodes[i]=ed.t;
	}

	else {//c==0 nodes not in set
		add=0.;
		child->size=trn->size+2;
		if (child->size>k){
			child->ub=-1;
			return;
		}
		child->eid=trn->eid+1;
		for (i=0;i<trn->size;i++){
			child->nodes[i]=trn->nodes[i];
			ev=isadj(trn->nodes[i],ed.s,g);
			if (ev.v>0){
				if (ev.e<trn->eid){
					child->ub=-1;
					return;
				}
				add+=ev.v;
			}
			ev=isadj(trn->nodes[i],ed.t,g);
			if (ev.v>0){
				if (ev.e<trn->eid){
					child->ub=-1;
					return;
				}
				add+=ev.v;
			}
		}
		child->nodes[i]=ed.s;
		child->nodes[i+1]=ed.t;
		add+=ed.w;
	}

	child->val=trn->val+add;
	//Computing new upper bound: add the weight of the edges that were added to the subgraph, then remove the weight of the last edges used to compute the upper bound
	child->ub=trn->ub+add-ed.w;
  beg=trn->eid+(k*(k-1))/2-(child->size*(child->size-1))/2+1;
  beg=(beg<el->e)?beg:el->e;
  end=trn->eid+(k*(k-1))/2-(trn->size*(trn->size-1))/2;
  end=(end<el->e)?end:el->e;
  for (i=beg;i<end;i++){
    child->ub-=el->edges[i].w;
  }

}

void mkchild2(treenode* trn,treenode* child,edgelist *el,graph *g,unsigned k){
	edge ed=el->edges[trn->eid];
	unsigned i,end;

	if (edgeinset(ed,trn)==3){//both nodes in set
		child->ub=-1;
		return;
	}

	for (i=0;i<trn->size;i++){
		child->nodes[i]=trn->nodes[i];
	}
	child->val=trn->val;
	child->size=trn->size;
	child->eid=trn->eid+1;

	//computing new upper bound: remove the weight of the examined edge and add the weight of the one that comes next after the edges used to compute the previous upper bound
	child->ub=trn->ub-ed.w;
	end=trn->eid+(k*(k-1))/2-(trn->size*(trn->size-1))/2;
	if (end<el->e){
		child->ub+=el->edges[end].w;
	}

}

//copy trn_in in trn_out
inline void copytrn(treenode* trn_in,treenode* trn_out){
	unsigned i;
	trn_out->size=trn_in->size;
	trn_out->val=trn_in->val;
	trn_out->ub=trn_in->ub;
	trn_out->eid=trn_in->eid;
	for (i=0;i<trn_in->size;i++){
		trn_out->nodes[i]=trn_in->nodes[i];
	}
}

treenode* find_dks(bheap* heap,edgelist* el, graph* g, unsigned k, double alpha){
	unsigned i,j,l=(k*(k-1))/2;
	treenode *trn,*child,*best;
	trn=heap->elems;
	child=trn+1;
	best=child+1;

	trn->size=2;
	trn->val=el->edges[0].w;
	trn->ub=0;
	for (i=0;i<l;i++){
		trn->ub+=el->edges[i].w;
	}
	trn->eid=1;
	trn->nodes[0]=el->edges[0].s;
	trn->nodes[1]=el->edges[0].t;
	trn=insert(heap,trn);

	trn->size=0;
	trn->val=0.;
	trn->ub=0;
	for (i=1;i<l+1;i++){
		trn->ub+=el->edges[i].w;
	}
	trn->eid=1;
	trn=insert(heap,trn);

	best->val=0.;

	while (heap->n>0){
		//printf("%u ",heap->n);//checking the size of the heap. Interestingly it is small.
		trn=popmax(heap,trn);
		if (trn->val>best->val){
			copytrn(trn,best);
		}
		if (trn->ub>best->val*alpha){
			if ((trn->size<k) && (trn->eid<el->e)){//subgraph of less than k nodes and not all edges considered yet
				mkchild1(trn,child,el,g,k);
				if (child->ub>best->val*alpha){
					child=insert(heap,child);
				}
				mkchild2(trn,child,el,g,k);
				if (child->ub>best->val*alpha){
					child=insert(heap,child);
				}
			}
		}
	}

	return best;
}

//Printing results
void printres(treenode* trn,edgelist *el,double alpha) {
	unsigned i=0;
	printf("size: %u\n",trn->size);
	printf("approximation: %lf\n",alpha);
	printf("val: %.10le\n",trn->val);
	printf("nodes:");
	for (i=0;i<trn->size;i++){
		printf(" %u",el->map[trn->nodes[i]]);
	}
	printf("\n");
}

int main(int argc,char** argv){
	unsigned k=atoi(argv[1]);
	double alpha=atof(argv[2]);
	edgelist *el;
	graph *g;
	treenode* trn;
	printf("Reading weighted edgelist from file %s\n",argv[3]);
	el=readedgelist(argv[3]);
	printf("Relabeling edgelist\n");
	relabel(el);
	printf("Number of nodes = %u\n",el->n);
	printf("Number of edges = %u\n",el->e);
	printf("Building the graph structure\n");
	g=mkgraph(el);
	printf("Building the heap structure\n");
	bheap *heap=construct(k);
	printf("Computing DkS\n");
	trn=find_dks(heap,el,g,k,alpha);
	printf("Printing results\n");
	printres(trn,el,alpha);
	printf("Releasing memory\n");
	free_edgelist(el);
	free_graph(g);
	free_bheap(heap);
	return 0;
}
