#include <hre/config.h>

//#include "Ddlts.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <mpi.h>
#include <hre/runtime.h>
#include <assert.h>
#include "Dtaudlts.h"
#include "bufs.h"
#include "sortcount.h"

#define INVERSE_TAG 10
#define COLLAPSE_TAG 15
#define SOME_TAG 20
#define REACH_TAG 30
#define WEIGHT_TAG 40

#define MAX_MANAGER_GRAPH_SIZE 100000

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wuninitialized"

//#define DEBUG
//#define VERBOSE

// in all types and functions
// i made the hypothesis that 
// the number of segments = the number of workers in _comm_



#ifdef DEBUG

int* totsize;

#endif

int Nme, Nnodes;

// the number of vertices and arcs ..
// final
extern int Nfinal, Mfinal;
// hooked, i.e. on a cycle
extern int Nhooked, Mhooked;
// not yet classified as "final" or "hooked"
extern int Nleft, Mleft;





/****************************************************


taudlts_t taudlts_create(MPI_Comm communicator)


*****************************************************/

taudlts_t taudlts_create(MPI_Comm communicator){
	taudlts_t t;
	t=(taudlts_t)malloc(sizeof(struct taudlts));
	if (!t) Fatal(1,error,"out of memory in taudlts_create");
  t->comm=communicator;
	t->N = 0;
  t->M = 0;
	t->begin=NULL;
	t->w=NULL;
	t->o=NULL;	
	return t;
}










/****************************************************


taudlts_t taudlts_free(taudlts_t t)


*****************************************************/

void taudlts_free(taudlts_t t){
 free(t->w);
 free(t->o);
 free(t->begin);
 free(t);
}


void taudlts_reset(taudlts_t t){
 free(t->w);t->w=NULL;
 free(t->o);t->o=NULL;
 free(t->begin);t->begin=NULL;
 t->M=t->N=0;
}



void taudlts_simple_join(taudlts_t t1, taudlts_t t2){
 // t1 and t2 are in the aux form, i.e. begin is 
 // array of transitions' sources and not array of start indeces!
 // effect: t1 := t1+t2, t2 gets reset
 int i, newM;
 
 MPI_Comm_rank(t1->comm, &Nme);
#ifdef DEBUG
 assert(t1->N == t2->N);
#endif
 newM = t1->M + t2->M + 1;
 
 // Warning(info,"%d a %d (%d+%d),%p,%p  a   ",
 //				 Nme,newM, t1->M, t2->M, t1, t1->begin); 
 if ((t1->begin = (int*)realloc(t1->begin, newM*sizeof(int))) == NULL)
	Fatal(1,error,"out of memory in simple_join1");

 // Warning(info,"%d b ",Nme); fflush(stdout);
 if ((t1->w = (int*)realloc(t1->w, sizeof(int)*newM)) == NULL)
	Fatal(1,error,"out of memory in simple_join2"); 
 // Warning(info,"%d c ",Nme); 
 if ((t1->o = (int*)realloc(t1->o, sizeof(int)*newM)) == NULL)
	Fatal(1,error,"out of memory in simple_join3");

 //Warning(info,"\n\n%d:M1=%d, M2=%d, newM=%d  t1 %p t1->begin %p",Nme,
 //	 t1->M, t2->M, newM, t1, t1->begin);

 for(i = 0; i < t2->M; i++){
	t1->begin[(t1->M)+i] = t2->begin[i];
	t1->w[(t1->M)+i] = t2->w[i];
	t1->o[(t1->M)+i] = t2->o[i];
 }

 t1->M = newM-1; 
 taudlts_reset(t2);
}


// transform the format of t between
// normal (as defined in Dtaudlts.h) and
// aux (t->begin is t->M long and 
//          holds sources instead of begin indeces)

void taudlts_aux2normal(taudlts_t t){
 int *osrc, *tmp;
 int i;
 MPI_Comm_rank(t->comm, &Nme);
 osrc = t->begin;
 if ((t->begin = (int*)calloc(t->N + 1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in aux2normal");
 if (t->M > 0){
	/*	for(i=0;i<t->M;i++){
	 if (osrc[i] >= t->N)
		Warning(info,"%d: i=%d, osrc=%d, N=%d  ",Nme,i,osrc[i],t->N);
	 assert(osrc[i] < t->N);
	 }
*/
	ComputeCount(osrc, t->M, t->begin); 
	Count2BeginIndex(t->begin, t->N);
	tmp=t->w; SortArray(&tmp, t->M, osrc, t->begin, t->N); t->w=tmp;
	tmp=t->o; SortArray(&tmp, t->M, osrc, t->begin, t->N); t->o=tmp;
	free(osrc);
 }
 // Warning(info,"\n%d:aux2normal: w=%p, o=%p, begin=%p, N=%d, M=%d\n",
 //				 Nme,t->w, t->o, t->begin, t->N, t->M);
}


void taudlts_normal2aux(taudlts_t t){
 int i,j;
 int* b;
 b = t->begin;
 if ((t->begin = (int*)calloc(t->M + 1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in normal2aux");
 for(i = 0; i < t->N; i++)
	for(j = b[i] ;j < b[i+1]; j++)
	 t->begin[j] = i;
 if (b!=NULL) free(b);
}



/****************************************************


void taudlts_copy(taudlts_t tfrom, taudlts_t tto)


*****************************************************/


void taudlts_copy(taudlts_t tfrom, taudlts_t tto){
 int i;
 tto->comm = tfrom->comm;
 tto->N = tfrom->N;
 tto->M = tfrom->M;
 tto->w = (int*)calloc(tto->M, sizeof(int));
 tto->o = (int*)calloc(tto->M, sizeof(int)); 
 tto->begin = (int*)calloc(tto->N+1, sizeof(int));
 for(i=0;i<tto->M;i++){
	tto->w[i] = tfrom->w[i];
	tto->o[i] = tfrom->o[i];
 }
 for(i=0;i<=tto->N;i++)
	tto->begin[i] = tfrom->begin[i];
}














/****************************************************


void taudlts_extract_from_dlts(taudlts_t t, dlts_t lts)


*****************************************************/

void taudlts_extract_from_dlts(taudlts_t t, dlts_t lts){
// extracts the FORWARD tau transitions from a given dlts
// and orders them 
// The resulted taults has ALL the states of the dlts.
// ANOTHER WAY (?): renumber the states and give a map dlts->taudlts
// need the following fields in dlts:
//    state_count[Nme]
//    transition_count[Nme][i], forall i
//    src,label,dest [Nme][i], forall i
//    tau

 int i, j, k, l, index, Nme, Nnodes, Mtot, Dtot, Stot;

 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme);
 t->comm = lts->comm;

 t->N=lts->state_count[Nme]; 
 t->begin=(int*)calloc(t->N+1, sizeof(int));
 if (!(t->begin))
	Fatal(1,error,"out of memory in extract_from_dlts"); 

 Mtot=0;
 for(i=0;i<Nnodes;i++){
	Mtot += lts->transition_count[Nme][i];
	for(j=0;j<lts->transition_count[Nme][i];j++)
	 if (lts->label[Nme][i][j] == lts->tau)
		t->begin[lts->src[Nme][i][j]]++;
 }
 Count2BeginIndex(t->begin, t->N);

#ifdef VERBOSE
 MPI_Allreduce(&(t->N), &Dtot, 1, MPI_INT, MPI_SUM, t->comm); 
 j=Mtot; Mtot=0; MPI_Allreduce(&j, &Mtot, 1, MPI_INT, MPI_SUM, t->comm); 
 if(Nme==0) Warning(info,"%12d     transitions and %12d states", Mtot, Dtot);
#endif

 t->M = t->begin[t->N];

 t->w=(int*)calloc(t->M, sizeof(int));
 t->o=(int*)calloc(t->M, sizeof(int));
 if ((!(t->w)) || (!(t->o)))
	Fatal(1,error,"out of memory in extract_from_dlts");

 for(i = Nnodes-1; i >= 0; i--){
	index = 0;
	for(j=0; j<lts->transition_count[Nme][i]; j++)
	 if (lts->label[Nme][i][j] == lts->tau){
		t->w[t->begin[lts->src[Nme][i][j]]]=i;
		t->o[t->begin[lts->src[Nme][i][j]]]=lts->dest[Nme][i][j];
		t->begin[lts->src[Nme][i][j]]++;
	 }
	 else{
		lts->src[Nme][i][index] = lts->src[Nme][i][j];
		lts->label[Nme][i][index] = lts->label[Nme][i][j];
		lts->dest[Nme][i][index] = lts->dest[Nme][i][j];
		index++;
	 }
	//	Warning(info,"%3d: %12d non-tau transitions to %3d .. %12d taus OUT (%12d)", 
	//					me, index, i, lts->transition_count[Nme][i] - index,t->M);
	lts->transition_count[Nme][i] = index;
	lts->src[Nme][i] = realloc(lts->src[Nme][i], index*sizeof(int));
	lts->label[Nme][i] = realloc(lts->label[Nme][i], index*sizeof(int));
	lts->dest[Nme][i] = realloc(lts->dest[Nme][i], index*sizeof(int));
 }

 EndIndex2BeginIndex(t->begin, t->N);
 j = l = 0;
 for(i=0;i<t->N;i++)
	if (t->begin[i+1]==t->begin[i]) j++;
	else for(k = t->begin[i]; k < t->begin[i+1]; k++)
	 if ((t->w[k]==Nme)&&(t->o[k]==i))
		l++;
 // Warning(info,"%3d has %12d tau transitions and %12d deadlocks",Nme, t->M, j);
#ifdef VERBOSE
 i=Mtot;
 MPI_Allreduce(&j, &Dtot, 1, MPI_INT, MPI_SUM, t->comm); 
 MPI_Allreduce(&(t->M), &Mtot, 1, MPI_INT, MPI_SUM, t->comm); 
 MPI_Allreduce(&l, &Stot, 1, MPI_INT, MPI_SUM, t->comm); 
 if(Nme==0) Warning(info,"%12d tau transitions (%3.3f pct.)\nand %12d deadlocks and %12d self-cycles", Mtot, ((float)Mtot/(float)i) * 100, Dtot, Stot);
#endif
}



/**************************
 ***
 ***      taudlts_is_dlts_without_labels(taudlts_t t, dlts_t lts)
 ***
 **************************/
/* DON'T KNOW IF REALLY NEEDED
taudlts_is_dlts_without_labels(taudlts_t t, dlts_t lts){
todo
}
*/











void taudlts_insert_to_dlts(taudlts_t t, dlts_t lts){

 int i, j, index;
 int* count_to_w;

 MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme);
 if (Nme==0) Warning(info,"INSERT TAU TRANSITIONS");
#ifdef VERBOSE
 Warning(info,"%d: %12d taus  ",Nme,t->M);
#endif
 if ((count_to_w = (int*)calloc(Nnodes, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in insert_to_dlts");
 for(i=0; i<t->M; i++){
	//	assert(t->w[i]<Nnodes);
	count_to_w[t->w[i]]++;
 }
 for(i = 0; i < Nnodes; i++){
	index = lts->transition_count[Nme][i] + count_to_w[i]+1;
	//	Warning(info,"%3d: %12d non-tau transitions to %3d .. %d new taus (%d) index=%d ", 
	//					me, lts->transition_count[Nme][i], i, count_to_w[i], t->M, index);
	lts->src[Nme][i] = realloc(lts->src[Nme][i], index*sizeof(int));
	//	Warning(info,"%d i=%d  ",Nme,i);
	lts->label[Nme][i] = realloc(lts->label[Nme][i], index*sizeof(int));
	lts->dest[Nme][i] = realloc(lts->dest[Nme][i], index*sizeof(int));
	if ((!lts->src[Nme][i])||(!lts->label[Nme][i])||(!lts->dest[Nme][i]))
	 Fatal(1,error,"out of memory in insert_to_dlts");
	//Warning(info,"%d oooo ",Nme);
 }
#ifdef DEBUG
 Warning(info,"%3d: lts reallocated  %d = %d",Nme, t->begin[t->N], t->M);
#endif
 for(i = 0; i < t->N; i++){
	//assert(t->begin[i]<=t->begin[i+1]);
	for(j = t->begin[i]; j < t->begin[i+1]; j++){
	 //assert(j < t->M);
	 index = lts->transition_count[Nme][t->w[j]];
	 lts->src[Nme][t->w[j]][index] = i;
	 lts->label[Nme][t->w[j]][index] = lts->tau;
	 lts->dest[Nme][t->w[j]][index] = t->o[j];
	 lts->transition_count[Nme][t->w[j]]++;
	}
 }
 free(count_to_w); 
 MPI_Barrier(lts->comm);
 if (Nme==0) Warning(info,"END INSERT TAU TRANSITIONS");
}







/****************************************************


void taudlts_load(taudlts_t t, char* filename, int type)


*****************************************************/

void taudlts_load(taudlts_t t, char* filename, int type){}
// reads actually a dlts, but only remembers the tau transitions..















/****************************************************


void taudlts_write(taudlts_t t, char* filename)


*****************************************************/

void taudlts_write(taudlts_t t, char* filename){
// dumps t as .aut

 FILE* output;
 MPI_Status status;
 char buffer[1024];
 int root, states, transitions, count, mytrans;
 int i,j,k,s,l,d;
 int *auxdata, *datapair, *first;
 
 ///// MPI_Barrier(t->comm);
 MPI_Comm_rank(t->comm, &Nme);
 MPI_Comm_size(t->comm, &Nnodes);

 // send all-to-all the number of states and outgoing transitions
 // store everything in the array _auxdata_
 auxdata=(int*)calloc(2*Nnodes, sizeof(int));
 first=(int*)calloc(Nnodes, sizeof(int));
 datapair=(int*)calloc(2, sizeof(int));
 datapair[0] = t->N;
 datapair[1] = t->M;
 ///// MPI_Barrier(t->comm); 
 MPI_Allgather(datapair, 2, MPI_INT, auxdata, 2, MPI_INT, t->comm);
 // compute all (i.e. for every worker) first global indexes 
 count=0;
 for(i=0;i<Nnodes;i++){
	first[i]=count;
	count+=auxdata[2*i];
 }

 for(i=0;i<Nnodes;i++) {
	if (Nme==i) { 
	 // open the file to write
	 if (Nme==0) {
		output=fopen(filename,"w");
		// 0: compute and write the first row: root, no. transitions, no. states 
		states=transitions=0; 
		for(k=0;k<Nnodes;k++) {
		 states+=auxdata[2*k];
		 transitions+=auxdata[2*k+1];
		}
		root=0;
		Warning(info,"Root: %d  States: %d  Transitions: %d", root, states, transitions); fflush(stdout);	
		fprintf(output,"des(%d,%d,%d)\n", root, transitions,states);
	 } 
	 else
		output=fopen(filename,"a");
	 Warning(info,"%d starts writing %d-%d",Nme,t->N,t->M); 
	 fflush(stdout);
	 // all: dump outgoing transitions
	 if (t->M > 0){
		Warning(info,"%d: actual M is %d (%d .. )",
						Nme,t->begin[t->N], t->begin[t->N-1]);
		for(k = 0; k < t->N; k++)
		 for(j = t->begin[k]; j < t->begin[k+1]; j++){
			s = first[Nme] + k;
			d = first[t->w[j]] + t->o[j];
			fprintf(output,"(%d,%s,%d)\n", s,"i",d);
		 }
	 }
	 fclose(output);
	} // end if Nme==i	
	//	Warning(info,"%d %d happy",Nme,i);
	///// MPI_Barrier(t->comm);
 }
 free(auxdata);free(first);free(datapair);
 // Warning(info,"\n>>>>>> %d finished!   N %d  M %d <<<<<<<",
 //		 Nme, t->N, t->M);
}



















/****************************************************


void taudlts_fwd2back(taudlts_t t)


*****************************************************/

void taudlts_fwd2back(taudlts_t t){
// transform the fwd representation of the LTS
// (i.e.: for all states, know the outgoing transitions)
// to the backward representation

 int i, j, x, Nnodes, Nme, src, aux;
 int *count_to_w, *begin_to_w, *count_from_w, *begin_from_w;
 int *to_w, *from_w;
 MPI_Request* request_array;
 MPI_Status* status_array;

 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 /*
 aux=x=0; 
 for(i=0;i<t->N;i++)
	for(j=t->begin[i];j<t->begin[i+1];j++)
	 if (t->w[j]==Nme){
		aux++;
		if (t->o[j]==i) x++;
	 }
 Warning(info,"%d BEFORE ...... %d self transitions, %d self loops",Nme, aux, x);
 */
 // alloc..
 if ((begin_to_w = (int*)calloc(Nnodes+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in fwd2back");
 if ((count_to_w = (int*)calloc(Nnodes+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in fwd2back");
 if ((begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in fwd2back"); 
 if ((count_from_w = (int*)calloc(Nnodes+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in fwd2back");
 if ((request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request)))==NULL)
	Fatal(1,error,"out of memory in fwd2back");
 if ((status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status)))==NULL)
	Fatal(1,error,"out of memory in fwd2back");
 // if ((!begin_to_w) || (!begin_from_w) || (!count_to_w) || (!count_from_w) || (!request_array) || (!status_array))
 //	Fatal(1,error,"out of memory in fwd2back"); 

 // count outgoing transitions
 for(i=0;i<t->M;i++)
	count_to_w[t->w[i]]++;
 begin_to_w[0]=0;
 for(i=1;i<=Nnodes;i++)
	begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];
	
 // exchange number of outgoing transitions
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm);
 // if (Nme==0) Warning(info,"EXCHANGED count_from_w");

 // prepare receive buffers 
 begin_from_w[0]=0;
 for(i=1;i<=Nnodes;i++)
  begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 from_w=(int*)calloc(2*begin_from_w[Nnodes], sizeof(int));
 if (!(from_w))
	Fatal(1,error,"out of memory in fwd2back");

 // prepare send buffers
 aux=count_to_w[Nme];
 to_w = (int*)calloc(2*(begin_to_w[Nnodes]-aux), sizeof(int));
 if (!(to_w))
	Fatal(1,error,"out of memory in fwd2back");
 for(i=Nme+1;i<=Nnodes;i++)
	begin_to_w[i]-=aux;

 /*
Warning(info,"%d: %d local transitions, %d expected and %d to be sent",
 				 Nme, aux, begin_from_w[Nnodes]-aux, begin_to_w[Nnodes]); 
Warning(info,"%d,before: local trans begin from %d and end in %d",
				me,begin_from_w[Nme], begin_from_w[Nme+1]);
 */

 for(src=0;src<t->N;src++)	
	for(i=t->begin[src];i<t->begin[src+1];i++)
	 if (t->w[i]==Nme){
		aux=begin_from_w[Nme]++;
		from_w[2*aux]=src;
		from_w[2*aux+1]=t->o[i];
	 }
	 else {
		aux=begin_to_w[t->w[i]]++;
		to_w[2*aux]=src;
		to_w[2*aux+1]=t->o[i];
	 }
 begin_from_w[Nme]-=count_from_w[Nme];

 // Warning(info,"%d,after: local trans begin from %d and end in %d begin: %p  ",
 // 		me,begin_from_w[Nme], begin_from_w[Nme+1],t->begin);

 for(i=Nnodes;i>0;i--)
	begin_to_w[i]=begin_to_w[i-1];
 begin_to_w[0]=0;
 free(t->w); t->w=NULL;
 free(t->o); t->o = NULL;
 // if (Nme==0) Warning(info,"BEEN HERE3!");
 free(t->begin); t->begin = NULL;

 // send and receive
 ///// MPI_Barrier(t->comm);
 // if (Nme==0) Warning(info,"EXCHANGING TRANSITIONS");
 aux=0;
 for(i=0;i<Nnodes;i++)
	if (i!=Nme){
	 //	 Warning(info,"%d->%d: %d   ",Nme,i,2 * count_to_w[i]);
	 MPI_Isend(to_w + 2 * begin_to_w[i],
						 2 * count_to_w[i], MPI_INT, 
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	 MPI_Irecv(from_w + 2 * begin_from_w[i],
						 2 * count_from_w[i], MPI_INT,
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	}
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(t->comm);
 // if (Nme==0) Warning(info,"EXCHANGED new transitions");

 // reorganize  t
 free(to_w);
 t->M = begin_from_w[Nnodes];
 t->w=(int*)calloc(t->M, sizeof(int));
 t->o=(int*)calloc(t->M, sizeof(int));
 t->begin=(int*)calloc(t->N+1, sizeof(int));
 if ((!(t->w)) || (!(t->o)) || (!(t->begin)))
	Fatal(1,error,"out of memory in fwd2back");
 
 for(i=0;i<t->M;i++)
	t->begin[from_w[2*i+1]]++;
 for(i=1;i <= t->N;i++)
	t->begin[i] += t->begin[i-1];
 for(i=Nnodes-1;i>=0;i--)
	for(j=begin_from_w[i+1]-1;j>=begin_from_w[i];j--){
	 t->begin[from_w[2*j+1]]--;
	 t->w[t->begin[from_w[2*j+1]]]=i;
	 t->o[t->begin[from_w[2*j+1]]]=from_w[2*j];
 }

 /*
 aux=x=0; 
 for(i=0;i<t->N;i++){
	for(j=t->begin[i];j<t->begin[i+1];j++)
	 if (t->w[j]==Nme){
		aux++;
		if (t->o[j]==i) x++;
	 }
 }
 Warning(info,"%d AFTER ...... %d self transitions, %d self loops",Nme, aux, x);
 */
 
 free(from_w);
 free(count_to_w); free(count_from_w);
 free(begin_to_w); free(begin_from_w);
 free(request_array); free(status_array);
}













/****************************************************


void taudlts_elim_trivial(taudlts_t t, int* oscc)


*****************************************************/

void taudlts_elim_trivial(taudlts_t t,  taudlts_t ta, int* oscc){
// topological sort
// input: t, oscc
//        same restrictions as for the input of reduce_some
// output: t without the transitions not belonging to a cycle
//         oscc is not changed!!
// All states are kept.

 int Nzeros, Nzeros_all, Mviz;
 int i,j, src, aux;
 int* degree;
 int *count_to_w, *begin_to_w, *count_from_w, *begin_from_w;
 int *to_w, *from_w;
 intbuf_t buflocal;
 MPI_Request* request_array;
 MPI_Status* status_array;
 
 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 
 if (Nme==0) Warning(info,"\nELIM TRIVIAL");

#ifdef DEBUG
	for(i=0; i<t->N; i++)
	 if (t->begin[i] < t->begin[i+1])
		assert(oscc[i]==i);
#endif

 // alloc..
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 count_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
 degree=(int*)calloc(t->N, sizeof(int));
 if ((!begin_to_w) || (!begin_from_w) || (!count_to_w) || (!count_from_w) 
		 || (!request_array) || (!status_array) || (!degree))
	Fatal(1,error,"out of memory in elim_trivial"); 
 buflocal=newBuffer(0); 

 // compute degree of all states
 //                  count trans./dest. worker
 Mviz=0;
 ComputeCount(t->w, t->M, count_to_w);
 count_to_w[Nme] = 0;
 begin_to_w[0] = 0; 
 for(i=1;i<=Nnodes;i++)
	begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];

 /*for(i=0;i<t->M;i++)
	count_to_w[t->w[i]]++;
 count_to_w[Nme]=0;
 begin_to_w[0]=0;
 for(i=1;i<=Nnodes;i++)
	begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];
 */

 //                  exchange number of outgoing transitions
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm);
 //                  prepare receive buffers 
 begin_from_w[0]=0;
 for(i=1;i<=Nnodes;i++)
	begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 from_w=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 if (!(from_w))
	Fatal(1,error,"out of memory in elim_trivial");
 //                  prepare send buffers
 to_w = (int*)calloc(begin_to_w[Nnodes], sizeof(int));
 if (!(to_w))
	Fatal(1,error,"out of memory in elim_trivial");
 for(src=0;src<t->N;src++)	
	for(i=t->begin[src];i<t->begin[src+1];i++)
	 if (t->w[i]==Nme)
		degree[t->o[i]]++;
	 else {
		aux=begin_to_w[t->w[i]]++;
		to_w[aux]=t->o[i];
	 }
 for(i=Nnodes;i>0;i--)
	begin_to_w[i]=begin_to_w[i-1];
 begin_to_w[0]=0;
 //                     send/receive

 //// ///// MPI_Barrier(t->comm);

 aux=0;
 for(i=0;i<Nnodes;i++)
	if (i!=Nme){
	 MPI_Isend(to_w + begin_to_w[i],
						 count_to_w[i], MPI_INT, 
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	 MPI_Irecv(from_w + begin_from_w[i],
						 count_from_w[i], MPI_INT,
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	}
 MPI_Waitall(aux, request_array, status_array);
 //// ///// MPI_Barrier(t->comm);
 //                      really compute degrees
 free(to_w);
 for(i=Nnodes-1;i>=0;i--)
	for(j=begin_from_w[i+1]-1;j>=begin_from_w[i];j--)
	 degree[from_w[j]]++;
 free(from_w);
 // Warning(info,">>>>>>>>>>>>>>>%d: computed the degrees",Nme); 
 //// ///// MPI_Barrier(t->comm);

 // iterations
 Nzeros_all=1;
 for(;;) {
	//        init
	resetBuffer(buflocal);
	for(i=0;i<Nnodes;i++)
	 count_to_w[i]=count_from_w[i]=0;
	//        count 
	Nzeros=0;aux=0;Nzeros_all=0;
	for(i=0;i<t->N;i++)
	 if (degree[i]==0) {
#ifdef DEBUG
		assert((oscc[i] == i) || (t->begin[i] == t->begin[i+1]));
#endif
		Nzeros++;
		aux+=t->begin[i+1]-t->begin[i];
		for(j=t->begin[i];j<t->begin[i+1];j++)
		 count_to_w[t->w[j]]++;
		if (t->begin[i] < t->begin[i+1]){
		 Nfinal++;
		 Nleft--;
		}
	 };	 
	//        test termination
		 //	Warning(info,">>>>>>>>>>>>>>>%d: %d zeros %d all %d eliminated transitions (%d total, %d local)",
		 //	me,Nzeros,Nzeros_all, aux,t->begin[t->N],count_to_w[Nme]); 
		 //  for(aux=0;aux<Nnodes;aux++) printf("COUNT%d-%d-%d\n",Nme,aux,count_to_w[aux]);	

	////	///// MPI_Barrier(t->comm);

	MPI_Allreduce(&Nzeros, &Nzeros_all, 1, MPI_INT, MPI_SUM, t->comm );
	//	if (Nme==0) Warning(info,"\n**************************\n%d total zeros",Nzeros_all);
	if (Nzeros_all==0) break;
	count_to_w[Nme]=0; begin_to_w[0]=0;
	for(i=1;i<=Nnodes;i++)
	 begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];
	MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm); 
	begin_from_w[0]=0;
	for(i=1;i<=Nnodes;i++)
	 begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
	from_w=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
	if (!(from_w))
	 Fatal(1,error,"out of memory in elim_trivial");
	//	 	Warning(info,">>>>>>>>>>>>>>>%d: before marking...",Nme); 
	//        mark transitions and prepare send buffers
	to_w = (int*)calloc(begin_to_w[Nnodes], sizeof(int));
	if (!(to_w))
	 Fatal(1,error,"out of memory in elim_trivial");
	//	for(aux=0;aux<Nnodes;aux++) printf("B%d-%d-%d\n",Nme,aux,begin_to_w[aux]);
	for(i=0;i<t->N;i++)
	 if (degree[i] == 0) {
		for(j = t->begin[i]; j < t->begin[i+1]; j++){
		 //		 Warning(info,"%d: i=%d, j=%d",Nme,i,j);
		 if (t->w[j]==Nme) 
			add1(buflocal,t->o[j]);
		 else{		
			//			if ((t->w[j] < 0)||(t->w[j] >= Nnodes))
			//			 Fatal(1,error,"%d: i=%d, wdest=%d outdegree=%d",Nme,i,t->w[j], t->begin[i+1]-t->begin[i]);
			to_w[begin_to_w[t->w[j]]]=t->o[j];
			begin_to_w[t->w[j]]++;
			 // if (begin_to_w[t->w[j]]>begin_to_w[Nnodes]){
			 //			 for(aux=0;aux<Nnodes;aux++) printf("B%d-%d-%d\n",Nme,aux,begin_to_w[aux]);
			 // Fatal(1,error,"%d: i=%d, j=%d, t->w[j]=%d, count=%d, index=%d, totalM=%d",
			 //			 Nme,i,j,t->w[j], count_to_w[t->w[j]],begin_to_w[t->w[j]],begin_to_w[Nnodes]);
			 //}
		 }
		 t->w[j] = -(t->w[j]) -1;    // to mark elimination
		 Mviz++;
		 Mfinal++;
		 Mleft--;
		};
		degree[i]=-1;
	 }; // end if degree==0		
	//Warning(info,">>>>>>>>>>>>>>>%d: before exchange...",Nme); 
	//        exchange
	for(i=Nnodes;i>0;i--)
	 begin_to_w[i]=begin_to_w[i-1];
	begin_to_w[0]=0;

	////	///// MPI_Barrier(t->comm);

	aux=0;
	for(i=0;i<Nnodes;i++)
	 if (i!=Nme){
		MPI_Isend(to_w + begin_to_w[i],
							count_to_w[i], MPI_INT, 
							i, 
							INVERSE_TAG,
							t->comm,
							request_array + aux);
		aux++;
		MPI_Irecv(from_w + begin_from_w[i],
							count_from_w[i], MPI_INT,
							i, 
							INVERSE_TAG,
							t->comm,
							request_array + aux);
		aux++;
	 }
	MPI_Waitall(aux, request_array, status_array);

	////	///// MPI_Barrier(t->comm);	

	//        end exchange
	//Warning(info,">>>>>>>>>>>>>>>%d: before decrease...",Nme); 
	//        decrease the degrees
	free(to_w);
	for(i=0;i<buflocal->index;i++)
	 degree[buflocal->b[i]]--;
	for(i=Nnodes-1;i>=0;i--)
	 for(j=begin_from_w[i+1]-1;j>=begin_from_w[i];j--)
		degree[from_w[j]]--;

	//        free
	free(from_w);
	//	freeBuffer(buflocal);
 } 
 
 // free..
 free(degree);
 freeBuffer(buflocal);
 free(count_to_w); free(count_from_w);
 free(begin_to_w); free(begin_from_w);
 free(request_array); free(status_array);

 // really eliminate the "negative"(i.e. vizible) transitions
#ifdef VERBOSE
 Warning(info,"%3d: %12d vizible transitions",Nme, Mviz);
#endif
 ta->begin = (int*)calloc(Mviz, sizeof(int));
 ta->w = (int*)calloc(Mviz, sizeof(int));
 ta->o = (int*)calloc(Mviz, sizeof(int));
 if ((!(ta->begin)) || (!(ta->w)) || (!(ta->o)))
	Fatal(1,error,"out of memory in elim_trivial");
 ta->N = t->N;
 ta->M = Mviz;
 Mviz=0;
 i = 0;  // i is "first free"
 for(src = 0; src < t->N; src++){
	aux = t->begin[src];
	t->begin[src] = i;
	for(j = aux; j < t->begin[src+1]; j++)
	 if(t->w[j] >= 0){
		t->o[i] = t->o[j];
		t->w[i] = t->w[j];
		i++;
	 }
	 else{
		ta->begin[Mviz] = src;
		ta->w[Mviz] = - (t->w[j]) -1;
		ta->o[Mviz++] = t->o[j];
	 }
 }

 j=aux=0;
 MPI_Reduce(&i, &aux, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&(t->M), &j, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if (Nme==0) Warning(info,"old M: %12d     new M: %12d",j,aux);
 t->begin[t->N] = t->M = i;
 t->w=(int*)realloc(t->w, (t->M) * sizeof(int));
 t->o=(int*)realloc(t->o, (t->M) * sizeof(int));
}



/****************************************************


void taudlts_print_info(taudlts_t t, int* oscc)


*****************************************************/
// compute how many: 
// - entry, exit and entry+exit states, cross transitions, local transitions
// - indegree 1, outdegree 1, indegree 1 + outdegree 1
// - states reachable from 0.. and the number and size 
//                      of the other reachability oasis


void taudlts_printinfo(taudlts_t t, int* oscc){


 char *entry, *exit;
 int *in, *out;
 int i, j, cross, Nentry, Nexit, Nee, Nin0, Nout0, Nin1, Nout1, Nio, max;
 taudlts_t ta;

 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 if (Nme==0) Warning(info,"\nSTATISTICS");

 ta=taudlts_create(t->comm);
 taudlts_elim_trivial(t, ta, oscc);

 entry=(char*)calloc(t->N, sizeof(char));
 exit=(char*)calloc(t->N, sizeof(char));
 in=(int*)calloc(t->N, sizeof(int));
 out=(int*)calloc(t->N, sizeof(int));
 
 cross=0;
 for(i=0;i<t->N;i++){
	out[i] = t->begin[i+1] - t->begin[i];
	for (j=t->begin[i];j<t->begin[i+1];j++)
	 if (t->w[j] != Nme){
		exit[i] = 1;
		cross++;
	 }
 } 
 MPI_Reduce(&cross, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"cross transitions: %12d",i);

 taudlts_fwd2back(t);
 
 
 taudlts_elim_trivial(t, ta, oscc);
 taudlts_free(ta);

 for(i=0;i<t->N;i++){
	in[i] = t->begin[i+1] - t->begin[i];
	for (j=t->begin[i];j<t->begin[i+1];j++)
	 if (t->w[j] != Nme)
		entry[i] = 1;
 }  
 taudlts_fwd2back(t);

 Nentry=Nexit=Nee=Nin0=Nout0=Nin1=Nout1=Nio=0;
 for(i=0;i<t->N;i++){
	if (entry[i]) Nentry++;
	if (exit[i]) Nexit++;
	if ((entry[i])&&(exit[i])) Nee++;
	if (in[i]==0) Nin0++; else if (in[i]==1) Nin1++;
	if (out[i]==0) Nout0++; else if (out[i]==1) Nout1++;
	if ((in[i]==1)&&(out[i]==1)) Nio++;
 }

 MPI_Reduce(&Nentry, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"entries       : %12d",i); 
 MPI_Reduce(&Nexit, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"exits         : %12d",i); 
 MPI_Reduce(&Nee, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"both          : %12d",i); 
 MPI_Reduce(&Nin0, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"in-degree 0   : %12d",i); 
 MPI_Reduce(&Nout0, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"out-degree 0  : %12d",i); 
 MPI_Reduce(&Nin1, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"in-degree 1   : %12d",i); 
 MPI_Reduce(&Nout1, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"out-degree 1  : %12d",i); 
 MPI_Reduce(&Nio, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if(Nme==0) Warning(info,"both degrees 1: %12d",i);
 max=0;for(i=0;i<t->N;i++) if (out[i]>max) max = out[i];
 for(j=2;j<max;j++){
	Nio=0; for(i=0;i<t->N;i++)
	 if (out[i]==j) Nio++;
	MPI_Reduce(&Nio, &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
	if(Nme==0) Warning(info,"out degree %3d: %12d",j,i);
 }

 free(in);free(out);free(entry);free(exit);
}




/****************************************************


void taudlts_elim_trivial(taudlts_t t, int* oscc)


*****************************************************/

void taudlts_elim_trivial_A_BIT_OLDER(taudlts_t t, int* oscc){
// topological sort
// input: t
// output: t without the transitions not belonging to a cycle
// output: oscc[x] = x, for all x with degree 0
//        (oscc must be initialized outside of this function)  
// All states are kept.

 int Nzeros, Nzeros_all;
 int i,j, src, aux;
 int* degree;
 int *count_to_w, *begin_to_w, *count_from_w, *begin_from_w;
 int *to_w, *from_w;
 intbuf_t buflocal;
 MPI_Request* request_array;
 MPI_Status* status_array;
 
 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 
 if (Nme==0) Warning(info,"\nELIM TRIVIAL");
 // alloc..
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 count_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
 degree=(int*)calloc(t->N, sizeof(int));
 if ((!begin_to_w) || (!begin_from_w) || (!count_to_w) || (!count_from_w) 
		 || (!request_array) || (!status_array) || (!degree))
	Fatal(1,error,"out of memory in elim_trivial"); 
 buflocal=newBuffer(0); 

 // compute degree of all states
 //                  count trans./dest. worker
 for(i=0;i<t->M;i++)
	count_to_w[t->w[i]]++;
 count_to_w[Nme]=0;
 begin_to_w[0]=0;
 for(i=1;i<=Nnodes;i++)
	begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];
 //                  exchange number of outgoing transitions
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm);
 //                  prepare receive buffers 
 begin_from_w[0]=0;
 for(i=1;i<=Nnodes;i++)
	begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 from_w=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 if (!(from_w))
	Fatal(1,error,"out of memory in elim_trivial");
 //                  prepare send buffers
 to_w = (int*)calloc(begin_to_w[Nnodes], sizeof(int));
 if (!(to_w))
	Fatal(1,error,"out of memory in elim_trivial");
 for(src=0;src<t->N;src++)	
	for(i=t->begin[src];i<t->begin[src+1];i++)
	 if (t->w[i]==Nme)
		degree[t->o[i]]++;
	 else {
		aux=begin_to_w[t->w[i]]++;
		to_w[aux]=t->o[i];
	 }
 for(i=Nnodes;i>0;i--)
	begin_to_w[i]=begin_to_w[i-1];
 begin_to_w[0]=0;
 //                     send/receive
 ///// MPI_Barrier(t->comm);
 aux=0;
 for(i=0;i<Nnodes;i++)
	if (i!=Nme){
	 MPI_Isend(to_w + begin_to_w[i],
						 count_to_w[i], MPI_INT, 
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	 MPI_Irecv(from_w + begin_from_w[i],
						 count_from_w[i], MPI_INT,
						 i, 
						 INVERSE_TAG,
						 t->comm,
						 request_array + aux);
	 aux++;
	}
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(t->comm);
 //                      really compute degrees
 free(to_w);
 for(i=Nnodes-1;i>=0;i--)
	for(j=begin_from_w[i+1]-1;j>=begin_from_w[i];j--)
	 degree[from_w[j]]++;
 free(from_w);
 //	 Warning(info,">>>>>>>>>>>>>>>%d: computed the degrees",Nme); 
 ///// MPI_Barrier(t->comm);

 // iterations
 Nzeros_all=1;
 for(;;) {
	//        init
	resetBuffer(buflocal);
	for(i=0;i<Nnodes;i++)
	 count_to_w[i]=count_from_w[i]=0;
	//        count 
	Nzeros=0;aux=0;Nzeros_all=0;
	for(i=0;i<t->N;i++)
	 if (degree[i]==0) {
		oscc[i] = i;
		Nzeros++;
		aux+=t->begin[i+1]-t->begin[i];
		for(j=t->begin[i];j<t->begin[i+1];j++)
		 count_to_w[t->w[j]]++;
	 };	 
	//        test termination
		 //	Warning(info,">>>>>>>>>>>>>>>%d: %d zeros %d all %d eliminated transitions (%d total, %d local)",
		 //	me,Nzeros,Nzeros_all, aux,t->begin[t->N],count_to_w[Nme]); 
		 //  for(aux=0;aux<Nnodes;aux++) printf("COUNT%d-%d-%d\n",Nme,aux,count_to_w[aux]);	
	///// MPI_Barrier(t->comm);
	MPI_Allreduce(&Nzeros, &Nzeros_all, 1, MPI_INT, MPI_SUM, t->comm );
	//	if (Nme==0) Warning(info,"\n**************************\n%d total zeros",Nzeros_all);
	if (Nzeros_all==0) break;
	count_to_w[Nme]=0; begin_to_w[0]=0;
	for(i=1;i<=Nnodes;i++)
	 begin_to_w[i]=begin_to_w[i-1]+count_to_w[i-1];
	MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm); 
	begin_from_w[0]=0;
	for(i=1;i<=Nnodes;i++)
	 begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
	from_w=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
	if (!(from_w))
	 Fatal(1,error,"out of memory in elim_trivial");
	 //	Warning(info,">>>>>>>>>>>>>>>%d: before marking...",Nme); 
	//        mark transitions and prepare send buffers
	to_w = (int*)calloc(begin_to_w[Nnodes], sizeof(int));
	if (!(to_w))
	 Fatal(1,error,"out of memory in elim_trivial");
	//	for(aux=0;aux<Nnodes;aux++) printf("B%d-%d-%d\n",Nme,aux,begin_to_w[aux]);
	for(i=0;i<t->N;i++)
	 if (degree[i] == 0) {
		for(j = t->begin[i]; j < t->begin[i+1]; j++){
		 //		 Warning(info,"%d: i=%d, j=%d",Nme,i,j);
		 if (t->w[j]==Nme) 
			add1(buflocal,t->o[j]);
		 else{		
			//			if ((t->w[j] < 0)||(t->w[j] >= Nnodes))
			//			 Fatal(1,error,"%d: i=%d, wdest=%d outdegree=%d",Nme,i,t->w[j], t->begin[i+1]-t->begin[i]);
			to_w[begin_to_w[t->w[j]]]=t->o[j];
			begin_to_w[t->w[j]]++;
			 // if (begin_to_w[t->w[j]]>begin_to_w[Nnodes]){
			 //			 for(aux=0;aux<Nnodes;aux++) printf("B%d-%d-%d\n",Nme,aux,begin_to_w[aux]);
			 // Fatal(1,error,"%d: i=%d, j=%d, t->w[j]=%d, count=%d, index=%d, totalM=%d",
			 //			 Nme,i,j,t->w[j], count_to_w[t->w[j]],begin_to_w[t->w[j]],begin_to_w[Nnodes]);
			 //}
		 }
		 t->w[j] = -(t->w[j])-1;    // to mark elimination
		};
		degree[i]=-1;
	 }; // end if degree==0		
	//	Warning(info,">>>>>>>>>>>>>>>%d: before exchange...",Nme); 
	//        exchange
	for(i=Nnodes;i>0;i--)
	 begin_to_w[i]=begin_to_w[i-1];
	begin_to_w[0]=0;
	///// MPI_Barrier(t->comm);
	aux=0;
	for(i=0;i<Nnodes;i++)
	 if (i!=Nme){
		MPI_Isend(to_w + begin_to_w[i],
							count_to_w[i], MPI_INT, 
							i, 
							INVERSE_TAG,
							t->comm,
							request_array + aux);
		aux++;
		MPI_Irecv(from_w + begin_from_w[i],
							count_from_w[i], MPI_INT,
							i, 
							INVERSE_TAG,
							t->comm,
							request_array + aux);
		aux++;
	 }
	MPI_Waitall(aux, request_array, status_array);
	///// MPI_Barrier(t->comm);	
	//        end exchange
	//	Warning(info,">>>>>>>>>>>>>>>%d: before decrease...",Nme); 
	//        decrease the degrees
	free(to_w);
	for(i=0;i<buflocal->index;i++)
	 degree[buflocal->b[i]]--;
	for(i=Nnodes-1;i>=0;i--)
	 for(j=begin_from_w[i+1]-1;j>=begin_from_w[i];j--)
		degree[from_w[j]]--;

	//        free
	free(from_w);
	//	freeBuffer(buflocal);
 } 
 
 // free..
 free(degree);
 freeBuffer(buflocal);
 free(count_to_w); free(count_from_w);
 free(begin_to_w); free(begin_from_w);
 free(request_array); free(status_array);

 // really eliminate the "negative" transitions
 i = 0;  // i is "first free"
 for(src = 0; src < t->N; src++){
	aux = t->begin[src];
	t->begin[src] = i;
	for(j = aux; j < t->begin[src+1]; j++)
	 if(t->w[j] >= 0){
		t->o[i] = t->o[j];
		t->w[i] = t->w[j];
		i++;
	 }
 }

 j=aux=0;
 MPI_Reduce(&i, &aux, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&(t->M), &j, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if (Nme==0) Warning(info,"old M: %12d     new M: %12d",j,aux);
 t->begin[t->N] = t->M = i;
}












/****************************************************


SCC LOCAL


*****************************************************/

int* dfsindex;
int *s,*b,*w,*o;
int* stack;
char* decided;
int dfsn = 0;
int stackn = 0;
int INF = 0xffffff;
int nnn=0;
int aaa=0;


int dfsTarjan(int x){
	int i,y;

	s[x] = x;
	dfsindex[x] = dfsn++;
	stack[stackn++] = x;
	decided[x] = 0;
	
	for (i=b[x];i<b[x+1];i++)
	 if (w[i] == Nme) {
		y=o[i];
		if (dfsindex[y] == INF)
		 dfsTarjan(y);
		if (!decided[y])
		 if (dfsindex[s[y]] < dfsindex[s[x]]){
			s[x] = s[y];		
			aaa++;
		 }
	 }
	
	if (s[x] == x){
	 stackn--;		
	 //	 assert(stackn>=0);
	 y = stack[stackn];
	 for (;x != y;){
		//		printf("\n**%d",y);
		s[y] = x;    // added just for try..
		decided[y] = 1;
		stackn--; 
		//		assert(stackn>=0);
		y = stack[stackn];
	 }
	 decided[x] = 1;
	}
	
	return 0;
}






void taudlts_scc_local(taudlts_t t, int* oscc){
// input: t
// output: ( t is not modified )
// output: map states -> their LOCAL scc head
// in the end, the "global" Nnodes are those for which oscc[x] = x
//                                                    and t->begin[x] < t->begin[x+1]
// i.e. the heads of local components
// and the "global" transitions are ...
// (to investigate whether "virtual" transitions would be useful)
// in the beginning: oscc[x]=x iff t->begin[x] = t->begin[x+1]; in rest, oscc[x] = -1
 

 int i,n,m,j, transout;

 // ///// MPI_Barrier(t->comm);
 MPI_Comm_rank(t->comm, &Nme); 

 if ((dfsindex=(int*)calloc(t->N, sizeof(int))) == NULL)
	Fatal(1,0,"cannot allocate dfsindex");
 if ((stack=(int*)calloc(t->N, sizeof(int))) == NULL)
	Fatal(1,0,"cannot allocate stack");
 if ((decided=(char*)calloc(t->N, sizeof(char))) == NULL)
	 Fatal(1,0,"cannot allocate decided");	
 
 s = oscc;
 b = (t->begin);
 w = (t->w);
 o = (t->o);

 for(i = 0; i < t->N; i++)
	 if (s[i]==-1) {
		dfsindex[i] = INF; 
		s[i] = INF;
	 }
	 else 
		dfsindex[i]=-1;

 for(i = 0; i < t->N; i++)
	 if (dfsindex[i] == INF)
		dfsTarjan(i);
 
 free(dfsindex);
 free(decided);
 free(stack);

 n=0;m=0;transout=0;
 for(i = 0; i < t->N; i++){
	if (oscc[i]==i)
	 n++;
	for(j=t->begin[i];j<t->begin[i+1];j++){
	 //	 assert(t->w[j]==Nme);
	 if (t->w[j]==Nme){
		if (oscc[i] == oscc[t->o[j]])
		 m++;
	 }	 
	 else transout++;
	}
 }

 // ///// MPI_Barrier(t->comm);
 Warning(info,"%3d: %d local components, %d transitions on local cycles (AND %d global, %d out)\n", 
				 Nme,n, m, t->M - m, transout);
#ifdef DEBUG
 for(i=0;i<t->N;i++){
	assert(oscc[i] >= 0);
	assert(oscc[i] < t->N);
	if (oscc[i] != oscc[oscc[i]])
	 Warning(info,"%d -> %d (%d) -> %d (%d) -> %d -> ..",
					 i,oscc[i], b[oscc[i]+1] - b[oscc[i]], 
					 oscc[oscc[i]], b[oscc[oscc[i]]+1] - b[oscc[oscc[i]]], 
					 oscc[oscc[oscc[i]]]);	
	 assert(oscc[oscc[i]] == oscc[i]);
 }
#endif
 /*
 m=0;
 for(i=0;i<t->N;i++){
	for(j=t->begin[i];((j < t->begin[i+1]) && (t->w[j]==Nme));) j++;
	if (j<t->begin[i+1]) m++;
 }
 Warning(info,"%d: %d exit states", Nme, m);
 */
}


//              END         scc    local                                    
//**********************************************************






/****************************************************




*****************************************************/

void taudlts_delete_transitions(taudlts_t t, char* deleted){
 // remove the transitions marked in _deleted_
 int i,j, index, index_old;

 index = 0; 
 for (i = 0; i < t->N; i++){
	index_old = t->begin[i];
	t->begin[i] = index;
	for(j = index_old; j < t->begin[i+1]; j++)
	 if (!deleted[j]){
		t->w[index] = t->w[j];
		t->o[index] = t->o[j];
		index++;
	 }
 }
 // Warning(info,"reduced from %d to %d transitions!", t->M, index);
 t->begin[t->N] = index;
 t->M = index;
 t->w = realloc(t->w, t->M * sizeof(int));
 t->o = realloc(t->o, t->M * sizeof(int));

}


void taudlts_elim_small(taudlts_t t, int* wscc, int* oscc){
 // eliminate very small cycles (<= 2)

 int i,j,k,x,y,Ndel,aux,Nme;
 taudlts_t tinverse;
 char* deleted;

 MPI_Comm_rank(t->comm, &Nme); 


 tinverse=taudlts_create(t->comm);
 taudlts_copy(t,tinverse);
 taudlts_fwd2back(tinverse);
 
 deleted=(char*)calloc(t->M, sizeof(char));
 Ndel=0;
 for(i=0;i<t->N;i++)
	for(j=t->begin[i]; (j<t->begin[i+1]) && (wscc[i]==Nme) && (oscc[i]==i);j++){
	 x=t->w[j]; y=t->o[j]; aux=0;
	 for (k = tinverse->begin[i]; ((k < tinverse->begin[i+1])&&(!aux)); k++) 
		if ((tinverse->w[k]=x)&&(tinverse->o[k]==y)){		 
		 deleted[j]=1;
		 Ndel++;
		 if (x < Nme) {wscc[i]=x; oscc[i]=y;}
		 else if ((x==Nme)&&(y<i)) {wscc[i]=x; oscc[i]=y;}
		 aux=1;
		}
	}

 taudlts_free(tinverse);

 taudlts_delete_transitions(t, deleted);
 free(deleted);
 
 j=0;
 for(i=0;i<t->N;i++)
	if ((wscc[i]==Nme)&&(oscc[i]==i))
	 j++;
 MPI_Reduce(&Ndel, &aux, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&(t->M), &i, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&j, &k, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if (Nme==0) Warning(info,"%12d transitions deleted, %12d left\n%12d head states left",aux, i, k); 
}




/****************************************************


void taudlts_reduce_some
   (taudlts_t* t, char* workers, 0, int* wscc, int* oscc, int* weight);


*****************************************************/


void BuffersToManager(char* workers, int Nnodes, int Nme, int manager, 
											int* bufout, int size, int* bufin, int* s, 
											MPI_Comm comm){

 int aux,i;
 MPI_Request* request_array;
 MPI_Status* status_array;
 /*
 Warning(info,"BuffersToManager. Nme %d (some %d) manager %d",
				 Nme, workers[Nme], manager);
 if (Nme != manager)
	Warning(info,"%d : send %d to M",Nme,size);
 else 
	for(i=0;i<Nnodes;i++)
	 if ((workers[i]) && ( i != manager))
		printf("expect from %d: %d",i, s[i+1]-s[i]);
  */

 aux = 0;
 if (Nme == manager){
	request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
	status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
	for(i = 0; i < Nnodes; i++)
	 if ((workers[i]) && (i != manager)){
		MPI_Irecv(bufin + s[i], s[i+1] - s[i], 
							MPI_INT, i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }
	//	Warning(info,"%d waiting",Nme);
	MPI_Waitall(aux, request_array, status_array);
 }
 else{
	MPI_Send(bufout, size, MPI_INT, manager, SOME_TAG, comm);
	//	Warning(info,"%d SENT",Nme);
 }
}



void BuffersFromManager(char* workers, int Nnodes, int Nme, int manager, 
											int* bufin, int size, int* bufout, int* s, 
											MPI_Comm comm){

 int aux,i;
 MPI_Request* request_array;
 MPI_Status* status_array;

 /*
 if (Nme != manager)
	Warning(info,"%d : send %d to M",Nme,size);
 else 
	for(i=0;i<Nnodes;i++)
	 if ((workers[i]) && ( i != manager))
		printf("expect from %d: %d",i, s[i+1]-s[i]);
 */

 aux = 0;
 if (Nme == manager){
	request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
	status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
	for(i = 0; i < Nnodes; i++)
	 if ((workers[i]) && ( i != manager)){
		MPI_Isend(bufout + s[i], s[i+1] - s[i], 
							MPI_INT, i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }
	MPI_Waitall(aux, request_array, status_array);
 }
 else
	MPI_Recv(bufin, size, MPI_INT, manager, SOME_TAG, comm,NULL);
}


void ExchangeBuffers(char* workers_from, char* workers_to, 
										 int Nnodes, int Nme, MPI_Comm comm,  
										 int* bufout, int* begin_to_w, int* size_to_w,
										 int** bufin, int* begin_from_w, int* size_from_w){
 // *bufin must be NULL if it's not the right size
 // if it is the right size, begin_from_w and size_from_w 
 // must also be filled in
 int aux,i,j;
 MPI_Request* request_array;
 MPI_Status* status_array;
 // if (Nme==0)
 //		Warning(info,"\n\n\n\n\n\n===========\n\n");
 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
 aux = 0;
 if (workers_from[Nme])
	for(i=0;i<Nnodes;i++)
	 if (workers_to[i]){
		MPI_Isend(size_to_w+i, 1, MPI_INT, 
							i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }	 
 if (workers_to[Nme])
	for(i=0;i<Nnodes;i++)
	 if (workers_from[i]){
		MPI_Irecv(size_from_w+i, 1, MPI_INT, 
							i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }
	 else size_from_w[i]=0;
 //	 else {begin_from_w[i] = 0; size_from_w[i]=0;}
 MPI_Waitall(aux, request_array, status_array);
 if (workers_to[Nme]){
	if (*bufin == NULL){
	 begin_from_w[0]=0;
	 for(i=1;i<Nnodes;i++)
		if (workers_from[i-1]) 
		 begin_from_w[i] = begin_from_w[i-1] + size_from_w[i-1];
		else begin_from_w[i] = begin_from_w[i-1];
	 begin_from_w[Nnodes] = begin_from_w[Nnodes-1] + size_from_w[Nnodes-1];
	 *bufin = (int*)calloc(begin_from_w[Nnodes]+1, sizeof(int));
	}
	//	else{ 
	 //	 Warning(info,"%d already set:  size %12d %12d %12d",
	 //					 Nme, size_from_w[0], size_from_w[1], size_from_w[2]);
	 //	 Warning(info,"%d already set:  begin %12d %12d %12d %12d",
	 // Nme, begin_from_w[0], begin_from_w[1], begin_from_w[2], begin_from_w[Nnodes]);
	 //	}
 }
 //  Warning(info,"%d counts exchanged",Nme);

	/*
 if (workers_from[Nme])
	Warning(info,"%d sizes TO  :  %12d %12d %12d",
					me, size_to_w[0], size_to_w[1], size_to_w[2]);
 if (workers_to[Nme])
	Warning(info,"%d sizes FROM:  %12d %12d %12d",
					me, size_from_w[0], size_from_w[1], size_from_w[2]); 
	*/

 aux = 0;
 if (workers_from[Nme])
	for(i = 0; i < Nnodes; i++)
	 if (workers_to[i]){
		MPI_Isend(bufout + begin_to_w[i], size_to_w[i], MPI_INT, 
							i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }
 if (workers_to[Nme])
	for(i=0;i<Nnodes;i++)
	 if (workers_from[i]){ 
		//		assert(begin_from_w[i]+size_from_w[i] <= begin_from_w[Nnodes]);
		MPI_Irecv(*bufin + begin_from_w[i], size_from_w[i], MPI_INT, 
							i, SOME_TAG, comm, request_array + aux);
		aux++;
	 }
 MPI_Waitall(aux, request_array, status_array); 
 // Warning(info,"%3d : received %d ",Nme, begin_from_w[Nnodes]);
 // Warning(info,"%d data exchanged",Nme);
}

























/****************************************************************/
void to_Manager(taudlts_t t, char* workers, int manager, 
								int* wscc, int* oscc, int* pns, taudlts_t tcopy, 
								int** back, 
								int* begintrans,
								int* beginstates){
/****************************************************************/
 
 int x, i, j, k, Nnodes, Nme, n, ns, first;
 int *osrc, *wdest, *odest;
 int *tmp, *begin_to_w, *begin_from_w, *size_to_w, *size_from_w, *gindex;
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 // Warning(info,"%d to Manager N: %d M: %d",Nme,t->N,t->M);
 // count SOME transitions
 n = ns = k= x = 0;
 for(i=0; i<t->N; i++){
	if ((oscc[i]==i)&&(wscc[i]==Nme)) ns++;
	if ((oscc[i]==i)&&(wscc[i]==Nme)&&(t->begin[i+1]>t->begin[i])) x++;
	//	if (t->begin[i+1] > t->begin[i])
	//		 assert((oscc[i]==i)&&(wscc[i]==Nme));
	for(j=t->begin[i]; j<t->begin[i+1]; j++)
	 if (workers[t->w[j]]) n++;
 }
 Warning(info, "%3d: %12d SOME (i.e. head) states (%12d all, %12d nonempty)\nand  %12d SOME transitions (%12d all)\n", 
				 Nme, ns, t->N, x, n, t->M);
 // send the states and transition counts
 if (Nme==manager){
	beginstates[Nme] = ns;
	begintrans[Nme] = n;
	for (i=0;i<Nnodes;i++)
	 if ((workers[i]) && (i!=manager)){
		MPI_Recv(beginstates+i, 1, MPI_INT, i, SOME_TAG, t->comm, NULL);
		MPI_Recv(begintrans+i, 1, MPI_INT, i, SOME_TAG, t->comm, NULL);
		//					Warning(info,"%d RECV",Nme);
	 }
	Count2BeginIndex(beginstates, Nnodes);
	Count2BeginIndex(begintrans, Nnodes);
 }
 else {	
	MPI_Send(&ns, 1, MPI_INT, manager, SOME_TAG, t->comm);
	MPI_Send(&n, 1, MPI_INT, manager, SOME_TAG, t->comm);
	//		Warning(info,"%d SSSSSSENT",Nme);
 } 
 // M->W  set first global index
 if (Nme==manager){
	first=beginstates[Nme];
	for (i=0;i<Nnodes;i++)
	 if ((workers[i]) && (i!=manager))
		MPI_Send(beginstates+i, 1, MPI_INT, i, SOME_TAG, t->comm);
 }
 else MPI_Recv(&first, 1, MPI_INT, manager, SOME_TAG, t->comm, NULL);


#ifdef DEBUG 
 MPI_Barrier(t->comm);
 if (Nme==manager)
	Warning(info, "MMM: %12d total SOME transitions", begintrans[Nnodes]);
#endif

 // prepare the buffers with SOME transitions  
 gindex = (int*)calloc(t->N, sizeof(int));
 osrc = (int*)calloc(n, sizeof(int));
 wdest = (int*)calloc(n, sizeof(int));
 odest = (int*)calloc(n, sizeof(int));	
 if (Nme != manager)
	tmp = (int*)calloc(ns, sizeof(int));
 else {
	*back = (int*)calloc(beginstates[Nnodes], sizeof(int));
	tmp = (*back) + beginstates[Nme];
 }
 x = 0; k = 0;
 for(i = 0; i < t->N; i++){
	if ((oscc[i] == i)&&(wscc[i]==Nme))
	 { gindex[i] = k; tmp[k] = i; k++; }
	//	else 
	//	 assert(t->begin[i+1] == t->begin[i]);
	for(j = t->begin[i]; j < t->begin[i+1]; j++)
	 if (workers[t->w[j]]){
		//		assert(x < n);
		osrc[x] = gindex[i]; 
		wdest[x] = t->w[j]; 
		odest[x] = t->o[j];
		x++;
	 }
 }	
 
 // build *back 
 if (Nme == manager){
	for(i = 0; i < Nnodes; i++)
	 if ((workers[i]) && (i!=manager)) 
		MPI_Recv((*back) + beginstates[i], 
						 beginstates[i+1]-beginstates[i], MPI_INT, i, SOME_TAG, t->comm,NULL);
		}
 else{ 
	MPI_Send(tmp, ns, MPI_INT, manager, SOME_TAG, t->comm);
	free(tmp); 
 }
 tmp = NULL;
#ifdef DEBUG 
  Warning(info,"%d transforminggggg! first=%d, have %d head states",Nme,first,ns);
#endif
 // transform local offsets to global
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 ComputeCount(wdest, n, begin_to_w);
 Count2BeginIndex(begin_to_w, Nnodes);
 //Warning(info,"%d q ",Nme);
#ifdef DEBUG
/* for(i = 0; i < Nnodes; i++){
	assert((workers[i]) || (begin_to_w[i] == begin_to_w[i+1]));
	assert((workers[i]) || (beginstates[i] == beginstates[i+1]));
	assert((workers[i]) || (begintrans[i] == begintrans[i+1]));
 }
*/
#endif
 SortArray(&odest, n, wdest, begin_to_w, Nnodes);
 SortArray(&osrc, n, wdest, begin_to_w, Nnodes); 
 free(wdest); 
#ifdef DEBUG
 Warning(info,"%d sort DONE",Nme);
#endif
 for(i=0;i<Nnodes;i++)
	size_to_w[i]=begin_to_w[i+1]-begin_to_w[i];
 ExchangeBuffers(workers, workers, Nnodes, Nme, t->comm, 
								 odest, begin_to_w, size_to_w, 
								 &tmp, begin_from_w, size_from_w);
 //   Warning(info,"%d exchanged destinations",Nme);
 for (i=0;i<begin_from_w[Nnodes];i++)
	tmp[i] = gindex[tmp[i]]+first;
 free(gindex);
 ExchangeBuffers(workers, workers, Nnodes, Nme, t->comm, 
								 tmp, begin_from_w, size_from_w,
								 &odest, begin_to_w, size_to_w);
 // now osrc has global names - first 
 // and odest has global names

 ////    Warning(info,"%d exchanged NEW GLOBAL destinations",Nme);
 // build tcopy->begin
 if (Nme==manager){
	tcopy->N = beginstates[Nnodes];
	tcopy->M = begintrans[Nnodes];
	tcopy->begin = (int*)calloc(tcopy->N+1, sizeof(int));
	tmp = tcopy->begin + beginstates[Nme];
 } 
 else tmp = (int*)calloc(ns+1, sizeof(int));
 ComputeCount(osrc, n, tmp);
 Count2BeginIndex(tmp, ns);
 SortArray(&odest, n, osrc, tmp, ns);
 free(osrc);
 //  Warning(info,"%d bla",Nme);
 if (Nme == manager){ 
	//	assert(tcopy->N == beginstates[Nnodes]);
	BuffersToManager(workers, Nnodes, Nme, manager, 
									 tmp, ns, tcopy->begin, beginstates, t->comm);
	for(i=0;i<Nnodes;i++)	 
	 for(j=beginstates[i];j<beginstates[i+1];j++)
		tcopy->begin[j] += begintrans[i];
	tcopy->begin[beginstates[Nnodes]] = begintrans[Nnodes];
#ifdef DEBUG
	for(i=0;i<Nnodes;i++){	 
	 //	 if (tcopy->begin[beginstates[i]] != begintrans[i])
	 //		Warning(info,"%d : %d != %d !!",
	 //						i, tcopy->begin[beginstates[i]],begintrans[i]);
	 assert((workers[i]) || (beginstates[i] == beginstates[i+1]));
	 assert((workers[i]) || (begintrans[i] == begintrans[i+1])); 
	 assert(tcopy->begin[beginstates[i]] == begintrans[i]);
	}
#endif
 }
 else { 
	BuffersToManager(workers, Nnodes, Nme, manager, 
									 tmp, ns, NULL, NULL, t->comm);
	free(tmp); tmp=NULL;
 }
#ifdef DEBUG
 Warning(info,"%d ver gekomen",Nme);
#endif
 // build tcopy->o
 if (Nme == manager){
	tcopy->o = (int*)calloc(tcopy->M, sizeof(int));
	for(i=0; i<n; i++) tcopy->o[begintrans[Nme]+i] = odest[i];
	BuffersToManager(workers, Nnodes, Nme, manager, 
									 odest, n, tcopy->o, begintrans, t->comm);
 }
 else BuffersToManager(workers, Nnodes, Nme, manager, 
											 odest, n, NULL, begintrans, t->comm);
#ifdef DEBUG
 if (Nme==manager){
	//		 	Warning(info, "MANAGER (%d): tcopy filled!! N=%d, M=%d", 
	//					me, tcopy->N, tcopy->M);
	assert(tcopy->M == tcopy->begin[tcopy->N]);
 }
#endif
 *pns = ns;
 // free unused memory
 free(odest);
}


/****************************************************************/
void to_Workers(taudlts_t t, char* workers, int manager, 
								int ns, int* wscc, int* oscc, 
								int* gscc, int* omapback, int* beginstates){
 /****************************************************************/
 
 int i,j,Nnodes,Nme;
 int *wmapback, *buf, *wsccold; 
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 

 if (Nme==manager) Warning(info,"%d's crowd : UPDATE wscc, oscc",Nme);
 // send/receive the new wscc values


 if (Nme == manager){
	wmapback = (int*)calloc(beginstates[Nnodes], sizeof(int));
	buf = (int*)calloc(beginstates[Nnodes], sizeof(int));
	for(i = 0; i < Nnodes; i++)
	 for(j = beginstates[i]; j < beginstates[i+1]; j++){
		wmapback[j] = i;
#ifdef DEBUG
		if(omapback[j]>=totsize[wmapback[j]])
		 Warning(info,"!!!! i=%d j=%d omapback[j]=%d wmapback[j]=%d totsize[..]=%d",
						 i,j,omapback[j],wmapback[j], totsize[wmapback[j]]);
		assert(omapback[j]<totsize[wmapback[j]]);
#endif
	 }

	for(i = 0; i < beginstates[Nnodes]; i++)
	 buf[i] = wmapback[gscc[i]];
	free(wmapback);
#ifdef DEBUG
	Warning(info,"%d wmapback DONE",Nme);
#endif
	BuffersFromManager(workers, Nnodes, Nme, manager, 
										 NULL, 0, buf, beginstates, t->comm);
	buf = buf + beginstates[Nme];
#ifdef DEBUG
	Warning(info,"%d wmapback SENT",Nme);
#endif
 }
 else{
	buf = (int*)calloc(ns, sizeof(int));
	BuffersFromManager(workers, Nnodes, Nme, manager, 
										 buf, ns, NULL, NULL, t->comm);
#ifdef DEBUG
	Warning(info,"%d wmapback RECV",Nme);
#endif
 }

 // save wscc
 wsccold = (int*)calloc(t->N, sizeof(int));
 for(i = 0; i < t->N; i++)
	wsccold[i]=wscc[i];
 // update wscc
 j = 0;
 for(i = 0; i < t->N; i++)
	if ((oscc[i] == i)&&(wscc[i] == Nme)){    // used to be a global state..
	 //	 if (j>=ns)
		//		Warning(info,"%d: i=%d, N=%d, j=%d, ns=%d",Nme,i,t->N, j, ns);
	 //	 assert(j < ns);
	 //	 assert(buf[j]<Nnodes);
	 wscc[i] = buf[j++];
	}



 // send/receive the new oscc values
 if (Nme==manager){
	free(buf);

#ifdef DEBUG
	for(i=0;i<Nnodes;i++)
	 for(j=beginstates[i];j<beginstates[i+1];j++)
		assert(omapback[j] < totsize[i]);
#endif

	for(i = 0; i < beginstates[Nnodes]; i++)
	 gscc[i] = omapback[gscc[i]];
	free(omapback);
	BuffersFromManager(workers, Nnodes, Nme, manager, 
										 NULL, 0, gscc, beginstates, t->comm);
	buf = gscc + beginstates[Nme];
	//	Warning(info,"omapback");
 }
 else
	BuffersFromManager(workers, Nnodes, Nme, manager, 
										 buf, ns, NULL, NULL, t->comm);
 //  Warning(info,"%d   oscc",Nme);
 j = 0;
 for(i = 0; i < t->N; i++)
	if ((oscc[i] == i)&&(wsccold[i] == Nme)){    // used to be a global state..
	 //	 assert(j < ns);
	 oscc[i] = buf[j++];
	}
#ifdef DEBUG	
 Warning(info,"%d   oscc updated",Nme);
#endif
 //free(buf); PROBLEM at worker 15, when it is a manager (WHY ?????)
 
 // statistics..
 j = 0;
 for(i = 0; i < t->N; i++)
	if ((wscc[i] == Nme) && (oscc[i] == i))
	 j++;
 Warning(info,"%3d: Had %12d global states, now I have %12d",
				 Nme,ns,j);

 // check that wscc and oscc are right
#ifdef DEBUG
 for(i=0;i<t->N;i++){
	if (oscc[i]>=totsize[wscc[i]])
	 Warning(info,"%d: i=%d, oscc=%d, wscc=%d, totsize=%d",
					 Nme,i,oscc[i],wscc[i],totsize[wscc[i]]);
	assert(oscc[i]<totsize[wscc[i]]);
 }
#endif

}






















/****************************************************************/
void update_destinations(taudlts_t t, int* wscc, int* oscc, char* workers){
/****************************************************************/

 int i,j,Nnodes,Nme;
 int *osrc, *odest, *tmp, *wtmp, *aux;
 int *begin_to_w, *begin_from_w, *size_to_w, *size_from_w;
 char *all;

 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 // Warning(info,"%d UPDATE DESTINATIONS ",Nme);
#ifdef DEBUG
 assert(t->begin[t->N] == t->M);
#endif
 // Warning(info,"%d blaa %d",Nme,t->M);

 // transform t->begin representation into osrc..
 if ((osrc=(int*)calloc(t->M, sizeof(int))) == NULL){
	Warning(info,"%3d: out of memory in update_destinations",Nme);
	exit(HRE_EXIT_FAILURE);
 }
	//	Fatal(1,error,"%3d: out of memory in update_destinations",Nme);
 // Warning(info,"%d blllllllllllllaa ",Nme);
 for(i = 0; i < t->N; i++){
	//	assert(t->begin[i]>=0);
	//	assert(t->begin[i]<=t->M);
	for(j = t->begin[i]; j < t->begin[i+1]; j++)
	 osrc[j] = i;
 }
 //Warning(info,"%d bbbbbbbbbbbbbbbblllllllllllllaa ",Nme);
 /*
 aux=(int*)calloc(t->N+1, sizeof(int));
 for(i=0;i<t->M;i++){
	assert(osrc[i]<t->N);
	aux[osrc[i]]++;
 }
 Count2BeginIndex(aux,t->N);
 for(i=0;i<t->N;i++)
	assert(aux[i]==t->begin[i]);
 */

 // Warning(info,"%d sorts ",Nme);
 ///// MPI_Barrier(t->comm);

 // sort transitions on t->w 
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 ComputeCount(t->w, t->M, begin_to_w);
 for(i = 0; i < Nnodes; i++) size_to_w[i] = begin_to_w[i];
 Count2BeginIndex(begin_to_w, Nnodes); 
 tmp = t->o; SortArray(&tmp, t->M, t->w, begin_to_w, Nnodes); t->o = tmp;
 SortArray(&osrc, t->M, t->w, begin_to_w, Nnodes); 
 for(i=0;i<Nnodes;i++)
	for(j=begin_to_w[i]; j<begin_to_w[i+1]; j++)
	 t->w[j] = i;

 ///// MPI_Barrier(t->comm);
 // Warning(info,"%d sorted ",Nme);
 // send the destinations to their owners
 all = (char*)calloc(Nnodes, sizeof(char));
 for(i = 0; i < Nnodes; i++) all[i] = 1; 
 tmp = NULL;
 ExchangeBuffers(all, workers, Nnodes, Nme, t->comm,
								 t->o, begin_to_w, size_to_w,
								 &tmp, begin_from_w, size_from_w);
 // Warning(info,"%d exchanged ",Nme);
 ///// MPI_Barrier(t->comm);

 //  for (i = 0 ; i < Nnodes ; i++)
 //	Warning (1,"%d -> %d: beginto %12d sizeto %12d",
 //					 Nme, i, begin_to_w[i], size_to_w[i]);
 wtmp = NULL;
 if (workers[Nme]){
	//	Warning(info,"%d : received %d %d %d",
	//					me, begin_from_w[0], begin_from_w[1], begin_from_w[Nnodes]);

	wtmp = (int*)calloc(begin_from_w[Nnodes], sizeof(int));
	for(i = 0; i < begin_from_w[Nnodes]; i++){ 
	 //	 assert(tmp[i] < t->N);
	 wtmp[i] = wscc[tmp[i]];
	}
	//	Warning(info,"zzzzzzzzz",Nme);
 }
 // Warning(info,"%d exchangeddddd ",Nme);

 ExchangeBuffers(workers, all, Nnodes, Nme, t->comm, 
								 wtmp, begin_from_w, size_from_w, 
								 &(t->w), begin_to_w, size_to_w);
 if (workers[Nme]){
	free(wtmp);
	for(i=0;i<begin_from_w[Nnodes];i++) 
	 tmp[i] = oscc[tmp[i]];
 }
 ExchangeBuffers(workers, all, Nnodes, Nme, t->comm,
								 tmp, begin_from_w, size_from_w,
								 &(t->o), begin_to_w, size_to_w);
 /*
 for (i = 0 ; i < Nnodes ; i++)
	Warning (1,"%d -> %d: beginfrom %d   beginto %d",
					 Nme, i, begin_from_w[i], begin_to_w[i]);
 */
 // Warning(info,"%d exchaaaaangedddddddd ",Nme);
 free(tmp); free(begin_to_w); free(begin_from_w); free(all);

 // normalize t
 // Warning(info,"%d fixes t ",Nme);
 /*  
 tmp = (int*)calloc(t->N, sizeof(int));
 for(i=0;i<t->M;i++){
	assert(osrc[i]<t->N);
	tmp[osrc[i]]++;
 }
 Count2BeginIndex(tmp,t->N);
 for(i=0;i<t->N;i++){
	assert(tmp[i]==t->begin[i]);
	assert(t->begin[i]<=t->begin[i+1]);
 } 
 */
#ifdef DEBUG
 assert(t->M == t->begin[t->N]);
#endif
 tmp = t->w; SortArray(&tmp, t->M, osrc, t->begin, t->N); t->w = tmp;
 tmp = t->o; SortArray(&tmp, t->M, osrc, t->begin, t->N); t->o = tmp;
 // PROBLEMMM  // free(osrc); osrc=NULL; 
 // Warning(info,"%d t fixed ",Nme);
}


/****************************************************************/
void	migrate_transitions(taudlts_t t, int* wscc, int* oscc, char* workers){
/****************************************************************/

 int i, j, Nnodes, Nme, back, firstnew;
 int *wsrc, *osrc, *tmp;
 int *size_to_w, *size_from_w, *begin_to_w, *begin_from_w;
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 

 /*
#ifdef DEBUG
 for(i=0;i<t->N;i++){
	if (oscc[i]>=totsize[wscc[i]])
	 Warning(info,"%d: i=%d, oscc=%d, wscc=%d, totsize=%d",
					 Nme,i,oscc[i],wscc[i],totsize[wscc[i]]);
	assert(oscc[i]<totsize[wscc[i]]);
 }
#endif DEBUG
 */


 //  Warning(info,"%3d MIGRATES",Nme);
 // fill in wsrc, osrc AND sort on wsrc
 wsrc=(int*)calloc(t->M, sizeof(int));
 osrc=(int*)calloc(t->M, sizeof(int));
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
 size_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 for(i = 0; i < t->N; i++)
	for(j = t->begin[i]; j < t->begin[i+1]; j++){
	 wsrc[j] = wscc[i];
	 osrc[j] = oscc[i];
	 begin_to_w[wsrc[j]]++;
	}
 free(t->begin);
 for(i=0;i<Nnodes;i++) size_to_w[i] = begin_to_w[i];
 Count2BeginIndex(begin_to_w, Nnodes);
 SortArray(&osrc, t->M, wsrc, begin_to_w, Nnodes);
 tmp=t->w;SortArray(&tmp, t->M, wsrc, begin_to_w, Nnodes);t->w=tmp;
 tmp=t->o;SortArray(&tmp, t->M, wsrc, begin_to_w, Nnodes);t->o=tmp;
 free(wsrc);
#ifdef DEBUG
 for(i=begin_to_w[Nme]; i<begin_to_w[Nme+1];i++)
	assert(osrc[i]<t->N);
#endif
 // exchange osrc
 // AND realloc and remove "left" transitions (only osrc)
 // AND fill in "arrived" transitions (osrc)
 tmp=NULL;
 ExchangeBuffers(workers, workers, Nnodes, Nme, t->comm, 
								 osrc, begin_to_w, size_to_w,
								 &tmp, begin_from_w,size_from_w);
#ifdef DEBUG
 for(i=0;i<begin_from_w[Nnodes];i++)
	assert(tmp[i] < t->N);
#endif
 back=0; 
 for(i=0;i<Nnodes;i++)
	if (workers[i])
	 back += (begin_to_w[i+1]-begin_to_w[i]);
	else
	 for(j = begin_to_w[i]; j < begin_to_w[i+1]; j++){
#ifdef DEBUG
		assert(j>=back);
		assert(j<t->M);
#endif
		osrc[j-back] = osrc[j];
	 }
 i=t->M; t->M = t->M - back + begin_from_w[Nnodes];
 osrc = realloc(osrc, t->M * sizeof(int));
 // Warning(info,"%3d: before migration %d transitions, after migration %d",
 //				 Nme, i, t->M);
 firstnew=i-back; 
 for(i=firstnew,j=0; j<begin_from_w[Nnodes]; j++,i++)
	osrc[i]=tmp[j];
#ifdef DEBUG
 for(i=0;i<t->M;i++)
	assert(osrc[i] < t->N);
#endif
 // exchange wdest
 // AND realloc and remove "left" transitions (only wdest)
 // AND fill in "arrived" transitions (wdest)
 ExchangeBuffers(workers, workers, Nnodes, Nme, t->comm, 
								 t->w, begin_to_w, size_to_w,
								 &tmp, begin_from_w,size_from_w);
 back=0; 
 for(i=0;i<Nnodes;i++)
	if (workers[i])
	 back += (begin_to_w[i+1]-begin_to_w[i]);
	else
	 for(j = begin_to_w[i]; j < begin_to_w[i+1]; j++)
		t->w[j-back] = t->w[j];
 t->w = realloc(t->w, t->M * sizeof(int));
 for(i=firstnew,j=0;j<begin_from_w[Nnodes];j++,i++)
	t->w[i]=tmp[j];
 // exchange odest
 // AND realloc and remove "left" transitions (only odest)
 // AND fill in "arrived" transitions (odest)
 ExchangeBuffers(workers, workers, Nnodes, Nme, t->comm, 
								 t->o, begin_to_w, size_to_w,
								 &tmp, begin_from_w,size_from_w);
 back=0; 
 for(i=0;i<Nnodes;i++)
	if (workers[i])
	 back += (begin_to_w[i+1]-begin_to_w[i]);
	else
	 for(j = begin_to_w[i]; j < begin_to_w[i+1]; j++)
		t->o[j-back] = t->o[j];
 t->o = realloc(t->o, t->M  * sizeof(int));
 for(i=firstnew,j=0;j<begin_from_w[Nnodes];j++,i++)
	t->o[i]=tmp[j];
 // sort on osrc and rebuild t->begin
 // Warning(info,"%3d now sorting M=%d",Nme,t->M);
 free(tmp);
 t->begin=(int*)calloc(t->N+1, sizeof(int));
#ifdef DEBUG
 for(i=0;i<t->M;i++)
	assert(osrc[i] < t->N);
#endif
 ComputeCount(osrc, t->M, t->begin);
 Count2BeginIndex(t->begin, t->N);
 // Warning(info,"%3d now sorting M=%d begin[%d]=%d",Nme,t->M, t->N,t->begin[t->N]);
 tmp=t->w;SortArray(&tmp, t->M, osrc, t->begin, t->N);t->w=tmp;
 // Warning(info,"%3d llll",Nme);
 tmp=t->o;SortArray(&tmp, t->M, osrc, t->begin, t->N);t->o=tmp;
 free(osrc);
#ifdef DEBUG
 for(i=0; i<t->N; i++)
	if (t->begin[i+1] > t->begin[i])
	 assert((oscc[i]==i)&&(wscc[i]==Nme));
#endif
	 // Warning(info,"%3d: migration finished",Nme);
}


/****************************************************************/
void taudlts_clear_useless_transitions1(taudlts_t t, char* workers, char verbose){
/****************************************************************/

 int i, j, k, Nnodes, Nme, max, size;
 int *begin_to_w, *tmp, *aux;
 char* deleted;
 if (t->M == 0) return;
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 

 begin_to_w=(int*)calloc(Nnodes+1, sizeof(int));
 // sort on source, t->w, t->o
 // max = 0; for(i=0;i<t->M;i++) if (t->o[i]>max) max = t->o[i]; max++;
 // Warning(info,"%3d: max is %d",Nme,max);
 for(i = 0; i < t->N; i++)
	if (t->begin[i+1] > t->begin[i]+1){
	 tmp = t->o + t->begin[i];
	 for(j = 0; j <= Nnodes; j++) begin_to_w[j]=0;
	 ComputeCount(t->w+t->begin[i],  t->begin[i+1] - t->begin[i], begin_to_w);
	 Count2BeginIndex(begin_to_w, Nnodes);
	 SortArray_copy(tmp, t->begin[i+1] - t->begin[i], 
									t->w + t->begin[i], begin_to_w, Nnodes);
	 if ( t->begin[i+1]-t->begin[i] > 1000)
		Warning(info,"%3d: sorted %3d (%d) on wdest", Nme, i, t->begin[i+1]-t->begin[i]);
	 
	 for(j=0;j<Nnodes;j++){
		for(k = t->begin[i] + begin_to_w[j] ; k < t->begin[i] + begin_to_w[j+1]; k++ )
		 t->w[k]=j;
		size = begin_to_w[j+1]-begin_to_w[j];
	 
	 //	 quicksort(t->o + t->begin[i]+begin_to_w[j], 
	 //						 0, begin_to_w[j+1]-begin_to_w[j]-1);
	 /*
	 if (size>1){
		//		Warning(info,"%d quick-sorting %d..%d",
		//			me, t->begin[i]+begin_to_w[j],t->begin[i]+begin_to_w[j+1]-1);
		//quicksort(t->o,  t->begin[i]+begin_to_w[j],  t->begin[i]+begin_to_w[j+1]-1);
		for(i=0;i<=max;i++) aux[i]=0;
		ComputeCount(tmp+ begin_to_w[j], size, aux);
		index=0;
		for(k=0;k<max;k++)
		for(t=0;t<aux[k];t++){
		(tmp+begin_to_w[j])[index+k]=i;
		index++;
		}
		}
	 */

	 if ((workers[j])&&(size>1))
		bucketsort(tmp+begin_to_w[j], size);
	 }
	}
 // Warning(info,"%d been so far",Nme);
 // mark what to delete
 deleted = (char*)calloc(t->M, sizeof(int));
 for(i=0;i<t->N;i++)
	if (t->begin[i+1] > t->begin[i]){
	 j=t->begin[i];if ((t->o[j] == i)&&(t->w[j]==Nme)) deleted[j]=1; j++;
	 for( ; j < t->begin[i+1]; j++){
		if ((t->o[j] == i)&&(t->w[j]==Nme)) deleted[j]=1;
		else if ((t->w[j]==t->w[j-1])&&(t->o[j]==t->o[j-1]))
		 deleted[j]=1;
	 }
	}
 // delete
 i=t->M;
 taudlts_delete_transitions(t, deleted);
 if (verbose)
	Warning(info,"%d useless transitions CLEARED. Old M: %12d New M: %12d", 
					Nme, i, t->M);
}



/****************************************************************/
void taudlts_clear_useless_transitions(taudlts_t t){
/****************************************************************/
 char* all;
 int i;
 MPI_Comm_size(t->comm, &Nnodes);
 all=(char*)calloc(Nnodes, sizeof(char));
 for(i=0;i<Nnodes;i++) all[i]=1;
 taudlts_clear_useless_transitions1(t,all,0);
}



/****************************************************************/
void taudlts_cleanup(taudlts_t t, int* wscc, int* oscc){
/****************************************************************/
 int i,j,n;
 char* deleted;
 MPI_Barrier(t->comm);
 MPI_Comm_rank(t->comm, &Nme); 
 if (Nme==0) Warning(info,"CLEANUP");
 taudlts_global_collapse(t, wscc, oscc);
 MPI_Barrier(t->comm);
 if (Nme==0) Warning(info,"after global collapse");
 MPI_Barrier(t->comm);
 // count self transitions
 if ((deleted = (char*)calloc(t->M, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in cleanup");
 n=0;
 for (i=0;i<t->N;i++)
	for (j=t->begin[i]; j < t->begin[i+1]; j++)
	 if ((t->w[j]==Nme)&&(t->o[j]==i)){
		Mhooked++;
		Mleft--;
		deleted[j]=1;
		n++;
	 }
 // taudlts_clear_useless_transitions(t);
 taudlts_delete_transitions(t, deleted);
 free(deleted);
 i = n; MPI_Reduce(&i, &n, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if (Nme==0) Warning(info,"%d intra transitions deleted",n);
}





/****************************************************************/
void taudlts_clear_useless_transitions_old(taudlts_t t){
/****************************************************************/

 int i, j, k, x, n, Nnodes, Nme, max, size;
 int *begin_to_w, *aux, *aw, *ao;
 char* deleted;
 if (t->M == 0) return;
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 
 Warning(info,"%d CLEARS USELESS TRANSITIONS",Nme);
 // sort transitions of every source ON odest!
 // and then on wdest
 for(x=0;x<t->N;x++)
	if (t->begin[x]+1 < t->begin[x+1]){
	 aw = t->w + t->begin[x];
	 ao = t->o + t->begin[x];
	 n = t->begin[x+1] - t->begin[x];
	 Warning(info,"%d %d.. %d %d %d .................",Nme,x, t->begin[x], t->begin[x+1], n);
	 sortpiece(aw, ao, 0, n-4);
	 Warning(info,"%d %d ************************",Nme,x);
	}
 /*
 for(x=0;x<t->N;x++)
	if (t->begin[x]+1 < t->begin[x+1]){
	 aw = t->w + t->begin[x];
	 ao = t->o + t->begin[x];
	 n = t->begin[x+1] - t->begin[x];
	 max=0; for(i=0;i<t->M;i++) if (t->o[i]>max) max = t->o[i]; max++;
	 aux=(int*)calloc(max+1, sizeof(int));
	 ComputeCount(ao, n, aux);
	 Count2BeginIndex(aux, max);
	 Warning(info,"%d %d .................",Nme,x);
	 SortArray_copy(aw, n, ao, aux, max);
	 Warning(info,"%d %d ************************",Nme,x);
	 for(i=0;i<max;i++)
		for(j=aux[i];j<aux[i+1];j++)
		 ao[j]=i;
	}
 */
 // mark what to delete
 deleted = (char*)calloc(t->M, sizeof(char));
 for(i=0;i<t->N;i++)
	if (t->begin[i+1] > t->begin[i]){
	 j=t->begin[i];if ((t->o[j] == i)&&(t->w[j]==Nme)) deleted[j]=1; j++;
	 for( ; j < t->begin[i+1]; j++){
		if ((t->o[j] == i)&&(t->w[j]==Nme)) deleted[j]=1;
		else if ((t->w[j]==t->w[j-1])&&(t->o[j]==t->o[j-1]))
		 deleted[j]=1;
	 }
	}
 // delete
 taudlts_delete_transitions(t, deleted);
 Warning(info,"%3d: new M: %d",Nme, t->M);
}

/****************************************************************/
void taudlts_reduce_some 
  (taudlts_t t, char* workers, int manager, 
	 int* wscc, int* oscc){
/****************************************************************/
// input: all
//    - only "head" states X with wscc[X]=Nme , oscc[X]=X  will be considered
//    - (I): the transitions in t must be only between head  states     
// output: 
//    - t modified (transitions found on cycles and doubles deleted)
//    - wscc, oscc modified to reflect the new SCC's found 
//    - t->w, t->o modified such that (I) still holds

 int i,j,n,ns;
 taudlts_t tcopy;
 int *omapback, *gscc, *begintrans, *beginstates;
 char* all;
 MPI_Request* request_array;
 MPI_Status* status_array;
 
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 

 /////  Warning(info,"%d starting REDUCE_SOME by manager %d",Nme, manager);
 ///// MPI_Barrier(t->comm);

 taudlts_scc_stabilize(t, wscc, oscc);
 MPI_Barrier(t->comm);


#ifdef DEBUG	
 i = t->N;
 totsize = (int*)calloc(Nnodes, sizeof(int));
 for(j=0;j<Nnodes;j++)
	MPI_Gather(&i, 1, MPI_INT, totsize, 1, MPI_INT, j, t->comm);
 for(i=0;i<t->N;i++){
	assert(wscc[i] >= 0);
	assert(wscc[i] < Nnodes);
	assert(oscc[i] >= 0);
	if (oscc[i] >= totsize[wscc[i]])
	 Warning(info,"%3d: i=%d, oscc=%d, wscc=%d, totsize=%d", Nme, i, oscc[i], wscc[i], totsize[wscc[i]]);
	assert(oscc[i] < totsize[wscc[i]]);
	if (t->begin[i]<t->begin[i+1]) assert((wscc[i]==Nme)&&(oscc[i]==i));
 }
 Warning(info,"%d: iiiiii ",Nme);
#endif

 //1. SOME: send to the manager all "global" transitions 
 //      with both ends in the SOME set
 //   SOME: Manager: create a taudlts structure for the "small" graph

 // if (Nme==manager){
	tcopy=taudlts_create(t->comm);
	// }
	// else tcopy=NULL;

 begintrans = beginstates = NULL;
 if (Nme == manager){
	begintrans=(int*)calloc(Nnodes+1, sizeof(int));
	beginstates=(int*)calloc(Nnodes+1, sizeof(int));
 }

 if (workers[Nme])
	to_Manager(t, workers, manager, wscc, oscc, &ns, 
						 tcopy, &omapback, begintrans, beginstates); 
 ///// MPI_Barrier(t->comm);

#ifdef DEBUG
 if (Nme==manager)
	for(i=0;i<Nnodes;i++)
	 for(j=beginstates[i];j<beginstates[i+1];j++){
		if (omapback[j] >= totsize[i])
		 Warning(info,"!!!! worker=%d, j=%d(%d), omapback=%d, totsize=%d",
						 i,j,j-beginstates[i],omapback[j],totsize[i]);
		assert(omapback[j] < totsize[i]);
	 }
#endif


 gscc = NULL;
 //2. SOME: Manager: find the local components
 if (Nme==manager){
	Warning(info,"\nSMALL GRAPH by manager %d: %12d states, %12d transitions \n", 
					Nme, tcopy->N, tcopy->M);
	gscc=(int*)calloc(tcopy->N, sizeof(int));
	for(i=0;i<tcopy->N;i++)
	 //	 if (tcopy->begin[i] == tcopy->begin[i+1])
	 //gscc[i]=i; else 
	 gscc[i] = -1;
	tcopy->w = (int*)calloc(tcopy->M, sizeof(int));
	for(i=0;i<tcopy->M;i++) tcopy->w[i] = Nme;
 }
	
 // taudlts_write(tcopy,"tcopy.aut");

	 if (Nme==manager){
		taudlts_scc_local(tcopy,gscc);		
#ifdef DEBUG
		for(i=0;i<tcopy->N;i++)
		 assert(gscc[gscc[i]] == gscc[i]);	
 Warning(info,"%d: ooooooooooooo ",Nme);
#endif
		taudlts_free(tcopy);
	 }

 ///// MPI_Barrier(t->comm);
 
 //3. SOME: receive the new scc heads and update the local wscc, oscc
 if (workers[Nme])
	to_Workers(t, workers, manager, ns, wscc, oscc, gscc, omapback, beginstates);
 omapback = NULL;
 ///// MPI_Barrier(t->comm);

 //4. ALL: update the destinations of transitions TO workers in SOME
 all=(char*)calloc(Nnodes, sizeof(char));
 for(i=0;i<Nnodes;i++) all[i]=1;
 update_destinations(t, wscc, oscc, all);
 ///// MPI_Barrier(t->comm);
 
 //5. SOME: migrate transitions to their new source : 
 //   the scc head of their source state
 if (workers[Nme])
	migrate_transitions(t, wscc, oscc, workers);
 ///// MPI_Barrier(t->comm);
 
 //6. SOME: delete useless transitions
 if (workers[Nme])
	taudlts_clear_useless_transitions1(t, workers,0);
 ///// MPI_Barrier(t->comm);

 // Warning(info,"%d finished REDUCE_SOME. New M: %d",Nme, t->M);
 //     taudlts_scc_stabilize(t, wscc, oscc);
}































void taudlts_scc_stabilize(taudlts_t t, int* wscc, int* oscc){

 int i, j, k, aux, total, total_all;
 int *tmp, *wdesttmp, *odesttmp, *wscctmp, *oscctmp;
 char *all, *newstable;
 int *begin_to_w, *begin_from_w, *size_to_w, *size_from_w;
 taudlts_t tscc;

 // initialize 
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme);  

 if (Nme==0)
	Warning(info,"\nSTABILIZE"); 
 ///// MPI_Barrier(t->comm);
 all =  (char*)calloc(Nnodes, sizeof(char));
 for(i=0;i<Nnodes;i++) all[i]=1;
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 size_from_w = (int*)calloc(Nnodes+1, sizeof(int));

 // build a graph wehere x->y if scc(x) = y
 tscc=taudlts_create(t->comm);
 tscc->N = tscc->M = t->N;
 tscc->w=(int*)calloc(t->N, sizeof(int));
 tscc->o=(int*)calloc(t->N, sizeof(int));
 tscc->begin=(int*)calloc(t->N+1, sizeof(int));
 for(i=0;i<t->N;i++) {
	tscc->w[i]=wscc[i]; 
	tscc->o[i]=oscc[i];
	tscc->begin[i]=i;
 };
 tscc->begin[tscc->N]=tscc->N;

 taudlts_clear_useless_transitions1(tscc, all,0);

 // then reverse it. 
 // Now arrows point to the "owned" states
 taudlts_fwd2back(tscc);

 //Warning(info,"\n%d reversed \n",Nme); 

 newstable=(char*)calloc(t->N, sizeof(char));
 j=0;
 for(i=0;i<t->N;i++) 
	if ((wscc[i]==Nme)&&(oscc[i]==i)) {
	 newstable[i]=1;
	 j++;
	}
 MPI_Allreduce(&j, &total, 1, MPI_INT, MPI_SUM, t->comm);
 total_all = 0;
 // Warning(info,"%d: initially %d stable (out of %d)", Nme, j,t->N);

 // check for cycles of 2: x->y->x
 for(i=0;i<t->N;i++)
	for(k=tscc->begin[i]; k<tscc->begin[i+1]; k++)
	 if ((tscc->w[k]==wscc[i])&&(tscc->o[k]==oscc[i]))
		Warning(info,"%d: %d is WRONG! ( %d , %d )\n",Nme,i,wscc[i],oscc[i]);

 // while
 while (total > 0) {
	if (Nme==0){
	 total_all += total;
	 //	 Warning(info,"Total new stable: %d total stable: %d",total, total_all);
	}

	aux=0;
	for(i=0;i<t->N;i++)
	 if (newstable[i])
		aux += tscc->begin[i+1] - tscc->begin[i];
	//		Warning(info,"%d: %d to send",Nme,aux);

	wdesttmp=(int*)calloc(aux, sizeof(int));
	odesttmp=(int*)calloc(aux, sizeof(int));
	wscctmp=(int*)calloc(aux, sizeof(int));
	oscctmp=(int*)calloc(aux, sizeof(int));
	
	j=0;
	for(i=0;i<t->N;i++)
	 if (newstable[i]){
		for(k=tscc->begin[i]; k<tscc->begin[i+1]; k++){
		 wdesttmp[j] = tscc->w[k];
		 odesttmp[j] = tscc->o[k];
		 wscctmp[j] = wscc[i];
		 oscctmp[j] = oscc[i];
		 j++;
		}
		newstable[i]=0;
	 }
	
	for(i=0;i<Nnodes;i++) begin_to_w[i]=0;
	ComputeCount(wdesttmp, aux, begin_to_w);
	for(i=0;i<Nnodes;i++) size_to_w[i]=begin_to_w[i];
	Count2BeginIndex(begin_to_w, Nnodes);
	SortArray(&odesttmp, aux, wdesttmp, begin_to_w, Nnodes);
	SortArray(&wscctmp, aux, wdesttmp, begin_to_w, Nnodes);
	SortArray(&oscctmp, aux, wdesttmp, begin_to_w, Nnodes);

	tmp=NULL;
	ExchangeBuffers(all, all, Nnodes, Nme, tscc->comm, 
									odesttmp, begin_to_w, size_to_w,
									&tmp, begin_from_w, size_from_w);
	odesttmp=(int*)realloc(odesttmp, begin_from_w[Nnodes] * sizeof(int));
	for(i=0;i<begin_from_w[Nnodes];i++) odesttmp[i] = tmp[i];

	///// MPI_Barrier(tscc->comm);
	free(tmp);tmp=NULL;
	ExchangeBuffers(all, all, Nnodes, Nme, tscc->comm, 
									wscctmp, begin_to_w, size_to_w,
									&tmp, begin_from_w, size_from_w);
	wscctmp=(int*)realloc(wscctmp, begin_from_w[Nnodes] * sizeof(int));
	for(i=0;i<begin_from_w[Nnodes];i++) wscctmp[i] = tmp[i];

	///// MPI_Barrier(tscc->comm);
	free(tmp);tmp=NULL;
	ExchangeBuffers(all, all, Nnodes, Nme, tscc->comm, 
									oscctmp, begin_to_w, size_to_w,
									&tmp, begin_from_w, size_from_w);
	oscctmp=(int*)realloc(oscctmp, begin_from_w[Nnodes] * sizeof(int));
	for(i=0;i<begin_from_w[Nnodes];i++) oscctmp[i] = tmp[i];
	free(tmp);tmp=NULL;

	//		Warning(info,"%d: %d received",Nme,begin_from_w[Nnodes]);

	for(i=0;i<begin_from_w[Nnodes];i++){
	 newstable[odesttmp[i]]=1;
	 wscc[odesttmp[i]]=wscctmp[i];
	 oscc[odesttmp[i]]=oscctmp[i];
	}

	free(wdesttmp); free(odesttmp); free(wscctmp); free(oscctmp);

	j=0;
	for(i=0;i<t->N;i++) 
	 if (newstable[i]) j++;
	//	Warning(info,"%d: %d stable", Nme, j);
	MPI_Allreduce(&j, &total, 1, MPI_INT, MPI_SUM, t->comm);
 }

 taudlts_free(tscc);
  MPI_Barrier(t->comm); 
 
  if (Nme==0)
	 Warning(info,"\nSTABILIZE ended\n\n");
} 






/****************************************************


void taudlts_local_collapse(taudlts_t t, int* map)


*****************************************************/

void taudlts_local_collapse(taudlts_t t, int* map){
// input: t, map
// output: t where the transitions x-> y become map[x] -> map[y]
//         and transitions map[x] -> map[x] are deleted
// input: domain(map) must be = codomain(map)
 
 int *b, *w, *o;
 int i,j;

 MPI_Comm_rank(t->comm, &Nme); 

 b = (int*)calloc(t->N + 1, sizeof(int));

 for (i=0;i<t->N;i++)
	for(j=t->begin[i];j<t->begin[i+1];j++)
	 if ((t->w[j]!=Nme) || ((t->w[j]==Nme) && (map[i] != map[t->o[j]])))
		b[map[i]]++;
 for (i=1;i<=t->N;i++)
	b[i] += b[i-1];

 // assert(b[t->N] <= t->begin[t->N]);
 t->M = b[t->N];
 w = t->w; o = t->o;
 t->w = (int*)calloc(t->M, sizeof(int));
 t->o = (int*)calloc(t->M, sizeof(int));
 for (i = 0;i <= t->N; i++)
	for (j = t->begin[i]; j < t->begin[i+1]; j++)
	 if ((w[j] != Nme) || ((w[j] == Nme) && (map[i] != map[o[j]]))){
		b[map[i]]--;
		t->w[b[map[i]]] = w[j];
		if (w[j] == Nme) t->o[b[map[i]]] = map[o[j]];
		else t->o[b[map[i]]] = o[j];
	 } 
 // Warning(info,"%d: had %d transitions.. %d deleted, %d left",
 //				 Nme, t->begin[t->N], t->begin[t->N] - t->M, t->M);
 if (t->begin != NULL) {free(t->begin); t->begin = b;}
 if (w!=NULL) free(w); if (o!=NULL) free(o);
}














/****************************************************


void taudlts_set_entry_exit(taudlts_t t, char* ee)


*****************************************************/

void taudlts_set_entry_exit(taudlts_t t, char* ee){}
// input: t
// output: ee[x] = 3 if x is both entry and exit;
//                 2, only exit, 1 only entry, 0 nothing












/****************************************************


void taudlts_to_virtual(taudlts_t t, char* ee, taudlts_t tabs)


*****************************************************/

void taudlts_to_virtual(taudlts_t t, char* ee, taudlts_t tabs){}
// input: t, ee
// output: tabs 
// a transition x->y in tabs Nmeans:
// x is entry, y is exit and there's a local path from x to y






























/****************************************************


void taudlts_scc_global(taudlts_t t, int* wscc, int* oscc)


*****************************************************/

void taudlts_scc_global(taudlts_t t, int* wscc, int* oscc){
// input: t
// input: oscc (= local map to component heads)
// output: map states -> their global scc head 
// SIDE EFFECT: it destroys t! (can be avoided by copying..)

 int NG;
 int i, j, first;
 int *count, *begin, *gscc, *wmapback_all, *omapback_all, *wscc_all, *oscc_all, *wmap, *omap, *omapback;
 int *aaux;
 int oldN;

 MPI_Comm_rank(t->comm, &Nme);
 // compute the number of owned "global" states
 NG=0;
 for(i=0;i<t->N;i++)
	if ((oscc[i]==i) && (t->begin[i+1]>t->begin[i]))
	 NG++;
 // compute "first"
 count = begin = NULL;
 if ( Nme == 0) { 
	MPI_Comm_size(t->comm, &Nnodes); 
	count = (int*)calloc(Nnodes,sizeof(int)); 
	begin = (int*)calloc(Nnodes+1,sizeof(int)); 
 } 
 MPI_Gather(&NG, 1, MPI_INT, count, 1, MPI_INT, 0, t->comm); 
 if (Nme==0){
	begin[0]=0;
	for(i=1;i<=Nnodes;i++) begin[i]=begin[i-1]+count[i-1];
	//	Warning(info,"begin: %d %d %d !!",begin[0],begin[1],begin[2]);
 }
 MPI_Scatter(begin, 1, MPI_INT, &first, 1, MPI_INT, 0, t->comm); 
 // "create" wmap, omap, wmapback, omapback
 wmap = (int*)calloc(t->N, sizeof(int));
 omap = (int*)calloc(t->N, sizeof(int));
 omapback = (int*)calloc(NG, sizeof(int));
 j = first;
 for(i=0;i<t->N;i++)
	if ((oscc[i]==i) && (t->begin[i+1]>t->begin[i])){
	 omap[i] = first++;
	 omapback[omap[i]-j]=i;
	}
 for(i=0;i<t->N;i++){
	wmap[i] = 0;
	omap[i] = omap[oscc[i]];
 }
 // move the "global" "small" graph to 0 and reduce it there
 oldN = t->N;
 taudlts_global_collapse(t, wmap, omap); // side effect at 0: t->N := begin[Nnodes]
 free(wmap); free(omap);
 ///// MPI_Barrier(t->comm);
 gscc = NULL;
 if (Nme==0){
	Warning(info,"\nSMALL GRAPH: %d states, %d transitions\n", t->N, t->M);
	gscc=(int*)calloc(t->N, sizeof(int));
	taudlts_scc_local(t,gscc);	
 }
 // compute wscc, oscc for the "global" graph
 wmapback_all = omapback_all = wscc_all = oscc_all = NULL;
 if (Nme==0){
	wmapback_all = (int*)calloc(t->N, sizeof(int));
	omapback_all = (int*)calloc(t->N, sizeof(int));
	wscc_all = (int*)calloc(t->N, sizeof(int));
	oscc_all = (int*)calloc(t->N, sizeof(int));
 }
 MPI_Gatherv(omapback, NG, MPI_INT, omapback_all, count, begin, MPI_INT, 0, t->comm);

 begin = NULL;
 if (Nme==0){
	// Warning(info,"begin: %d %d %d !!",begin[0],begin[1],begin[2]);
	aaux = (int*)calloc(Nnodes, sizeof(int));
	for(i=0;i<Nnodes;i++)
	 for(j=begin[i];j<begin[i+1];j++)
		wmapback_all[j]=i;
	for(i=0;i<t->N;i++){
	 wscc_all[i] = wmapback_all[gscc[i]];
	 aaux[wscc_all[i]]++;
	 oscc_all[i] = omapback_all[gscc[i]];
	}
	free(wmapback_all); free(omapback_all); free(gscc);
	Warning(info,"DISTRIBUTION OF THE SMALL GRAPH: %d in 0, %d in 1, %d in 2 !!",aaux[0],aaux[1],aaux[2]);
	free(aaux);
 }
 aaux = (int*)calloc(NG, sizeof(int));
 MPI_Scatterv(wscc_all, count, begin, MPI_INT, aaux, NG, MPI_INT, 0, t->comm);
 for(i=0;i<oldN;i++)
	wscc[i] = Nme;
 for(i=0;i<NG;i++)
	wscc[omapback[i]] = aaux[i];
 if (Nme==0) free(wscc_all);
 MPI_Scatterv(oscc_all, count, begin, MPI_INT, aaux, NG, MPI_INT, 0, t->comm);
 if (Nme==0) {free(count); free(begin);}
 for(i=0;i<NG;i++)
	oscc[omapback[i]] = aaux[i];
 if (Nme==0)free(oscc_all);
 free(aaux);

 // extend wscc, oscc to all states. 
 // omapback still contains the local "global" states

 t->N=oldN; j=0;
 /*
 for(i=0;i<t->N;i++)
	if (i < omapback[j]) {
	 wscc[i] = wscc[oscc[i]];
	 	 if ((oscc[i]==i)&&(wscc[i]!=Nme))
	 		Fatal(1,error,"bad scc of non-global state: %d,%d,wscc=%d",Nme,i,wscc[i]);
	 oscc[i] = oscc[oscc[i]];
	}
	else {
	 if (i>omapback[j])
		Fatal(1,error,"RRRRRRRR%d: i=%d(out of %d) j=%d(out of %d) omapback=%d",
					me,i,t->N,j,NG,omapback[j]);
	 j++;
	}  
 */
 for(i=0;i<t->N;i++)
	if (wscc[i]==Nme){
	 wscc[i] = wscc[oscc[i]];
	 oscc[i] = oscc[oscc[i]];
	}
 free(omapback);
}













/****************************************************


void taudlts_to_real(taudlts_t tabs, char* ee, int* wscc, int* oscc, taudlts_t t)


*****************************************************/

void taudlts_to_real(taudlts_t tabs, char* ee, int* wscc, int* oscc, taudlts_t t){
// input: tabs, ee, t 
// input: wscc, oscc defined for the states in ee
// output: wscc, oscc defined for all states of t

}











void taudlts_global_collapse(taudlts_t t, int* wscc, int* oscc){

 char* all;
 int i;

 MPI_Comm_size(t->comm, &Nnodes);
 all=(char*)calloc(Nnodes, sizeof(char));
 for(i=0;i<Nnodes;i++) all[i]=1;

 fflush(stdout);fflush(stderr);
 ///// MPI_Barrier(t->comm);
 update_destinations(t, wscc, oscc, all);
 ///// MPI_Barrier(t->comm);
 migrate_transitions(t, wscc, oscc, all);
 ///// MPI_Barrier(t->comm);
}


/****************************************************


void taudlts_global_collapse(taudlts_t t, int* wmap, int* omap)


*****************************************************/

void taudlts_global_collapse_wrong(taudlts_t t, int* wmap, int* omap){

// input: t, wmap, omap
// output: t where the transitions 
//    x -> y become wmap[x],omap[x] -> wmap[y], omap[y]

 int *wsrc, *osrc, *wdest, *odest;
 int *aux1, *aux2, *aux3;
 int *count_to_w, *begin_to_w, *count_from_w, *begin_from_w; 
 int i,j,aux;
 MPI_Request* request_array;
 MPI_Status* status_array;
 
 ///// MPI_Barrier(t->comm);
 MPI_Comm_size(t->comm, &Nnodes);
 MPI_Comm_rank(t->comm, &Nme); 

 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
 wsrc=(int*)calloc(t->M, sizeof(int));
 osrc=(int*)calloc(t->M, sizeof(int));

 // transform t-wmap-omap to wsrc-osrc-wdest-odest
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 for(i=0;i<t->N;i++)
	for(j = t->begin[i]; j < t->begin[i+1]; j++) {
	 wsrc[j] = wmap[i];
	 osrc[j] = omap[i];
	 begin_to_w[wsrc[j]]++;
	}
 free(t->begin); 
 wdest=t->w; odest=t->o;

 count_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 
 // sort on wsrc
 //     sort4(wsrc, osrc, wdest, odest, t->M, Nnodes, begin_to_w);
 ComputeCount(wsrc,t->M,begin_to_w);
 Count2BeginIndex(begin_to_w, Nnodes);
 SortArray(&osrc,t->M,wsrc,begin_to_w,Nnodes);
 SortArray(&wdest,t->M,wsrc,begin_to_w,Nnodes);
 SortArray(&odest,t->M,wsrc,begin_to_w,Nnodes);

 // exchange osrc-wdest-odest 
 for(i=0;i<Nnodes;i++) count_to_w[i] = begin_to_w[i+1]-begin_to_w[i];
 //  Warning(info,"%d: tosend %d to_myself %d",Nme, begin_to_w[Nnodes], count_to_w[Nme]);
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm); 
 begin_from_w[0]=0; 
 for(i=1;i<=Nnodes;i++) begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 // Warning(info,"%d: toreceive:%d from_myself %d",Nme,begin_from_w[Nnodes], count_from_w[Nme]);
 t->M = begin_from_w[Nnodes];
 aux1=(int*)calloc(t->M, sizeof(int));
 ///// MPI_Barrier(t->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(osrc + begin_to_w[i],
						count_to_w[i], MPI_INT, 
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(aux1 + begin_from_w[i],
						count_from_w[i], MPI_INT,
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(t->comm);
 free(osrc);osrc=aux1;
 aux2=(int*)calloc(t->M, sizeof(int));
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(wdest + begin_to_w[i],
						count_to_w[i], MPI_INT, 
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(aux2 + begin_from_w[i],
						count_from_w[i], MPI_INT,
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(t->comm);
 free(wdest);wdest=aux2;
 // for (i=0;i<t->M;i++)
 //	if (aux2[i]>=Nnodes)
 //	 Fatal (1,1,"%d: i=%d, wdest=%d M=%d",Nme,i,aux2[i],t->M);
 aux3=(int*)calloc(t->M, sizeof(int));
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(odest + begin_to_w[i],
						count_to_w[i], MPI_INT, 
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(aux3 + begin_from_w[i],
						count_from_w[i], MPI_INT,
						i, 
						COLLAPSE_TAG,
						t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(t->comm);
 free(odest);odest=aux3;

 // sort on wdest
 for(i=0;i<=Nnodes;i++) count_to_w[i]=count_from_w[i]=begin_to_w[i] = 0;
 //             sort3(wdest, osrc, odest, t->M, Nnodes, begin_to_w);
 // sort3(&wdest, &osrc, &odest, t->M, Nnodes, begin_to_w);
 ComputeCount(wdest,t->M,begin_to_w);
 Count2BeginIndex(begin_to_w, Nnodes);
 SortArray(&osrc,t->M,wdest,begin_to_w,Nnodes);
 SortArray(&odest,t->M,wdest,begin_to_w,Nnodes);


 // exchange mappings of destinations 
 //       echange counts
 for(i=0;i<Nnodes;i++) count_to_w[i] = begin_to_w[i+1]-begin_to_w[i];
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,t->comm);
 begin_from_w[0]=0; 
 for(i=1;i<=Nnodes;i++) begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 // Warning(info,"%d: i'm the new destination of %d transitions, %d local",
 //		 Nme, begin_from_w[Nnodes], begin_from_w[Nme]);

 aux1=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 aux2=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 //       send odest / receive aux1  
 ///// MPI_Barrier(t->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(odest + begin_to_w[i],
						count_to_w[i], MPI_INT, i, COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(aux1 + begin_from_w[i],
						count_from_w[i], MPI_INT,i,COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 //      aux2 := wmap[aux1]
 for(i=0;i<begin_from_w[Nnodes];i++)
	aux2[i] = wmap[aux1[i]];
 //      send aux2 / receive wdest
 ///// MPI_Barrier(t->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(aux2 + begin_from_w[i],
						count_from_w[i], MPI_INT, i, COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(wdest + begin_to_w[i],
						count_to_w[i], MPI_INT,i,COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array); 
 //      aux2 := omap[aux1]
 for(i=0;i<begin_from_w[Nnodes];i++)
	aux2[i] = omap[aux1[i]];
 //      send aux2 / receive odest
 ///// MPI_Barrier(t->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(aux2 + begin_from_w[i],
						count_from_w[i], MPI_INT, i, COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(odest + begin_to_w[i],
						count_to_w[i], MPI_INT,i,COLLAPSE_TAG, t->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);   
 ///// MPI_Barrier(t->comm);
 free(aux1); free(aux2);
 // re-create t
 t->N=-1;
 for(i=0;i<t->M;i++)
	if (osrc[i] > t->N) t->N = osrc[i];
 t->N++;
 t->begin=(int*)calloc(t->N+1, sizeof(int));

 //        sort3(osrc, wdest, odest, t->M, t->N, t->begin);
 //sort3(&osrc, &wdest, &odest, t->M, t->N, t->begin);
 
 ComputeCount(osrc, t->N, t->begin);
 Count2BeginIndex(t->begin, t->N);
 SortArray(&wdest,t->M,osrc,t->begin,t->N);
 SortArray(&odest,t->M,osrc,t->begin,t->N); 
 free(osrc);

 t->w = wdest;
 t->o = odest; 
 ///// MPI_Barrier(t->comm);
 Warning(info,"%d after shuffle: %d states, %d transitions",Nme, t->N, t->M);
 ///// MPI_Barrier(t->comm);
}


















/****************************************************


void dlts_shuffle(dlts_t lts, int* wmap, int* omap)


*****************************************************/

void dlts_shuffle(dlts_t lts, int* wmap, int* omap){

// input: lts, wmap, omap
// output: lts where the transitions 
//    x -a-> y become wmap[x],omap[x] -a-> wmap[y], omap[y]
// for every a =/= tau or  wmap[x],omap[x] =/= wmap[y], omap[y]
// input: wmap, omap should be lts->state_count[Nme] long
// SIDE EFFECT: 
//      destroys state_count and transition_count
//         - i.e. only state_count[Nme] is correct and transition_count[Nme][i]

 int *wsrc, *osrc, *lab, *wdest, *odest;
 int *aux1, *aux2, *aux3, *aux4;
 int *count_to_w, *begin_to_w, *count_from_w, *begin_from_w; 
 int i,j,aux;
 int M, N;
 MPI_Request* request_array;
 MPI_Status* status_array;
 
 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme); 

 N = lts->state_count[Nme];
 M = 0;
 for(i=0;i<Nnodes;i++) M += lts->transition_count[Nme][i];

#ifdef DEBUG
 Warning(info,"\n%d: N=%12d, M=%12d",N,M);
#endif
 (void) N;
 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));
 wsrc=(int*)calloc(M+1, sizeof(int));
 osrc=(int*)calloc(M+1, sizeof(int));
 lab=(int*)calloc(M+1, sizeof(int));
 wdest=(int*)calloc(M+1, sizeof(int));
 odest=(int*)calloc(M+1, sizeof(int));
 begin_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 if ((!request_array) || (!status_array) || 
		 (!wsrc) || (!osrc) || (!wdest) || (!odest) || (!begin_to_w))
	Fatal(1,error,"out of memory in dlts_shuffle");

 // fill in wsrc-osrc-lab-wdest-odest
 M = 0;
 for(i=0;i<Nnodes;i++)
	for(j = 0; j < lts->transition_count[Nme][i]; j++,M++) {
	 wsrc[M] = wmap[lts->src[Nme][i][j]];
	 begin_to_w[wsrc[M]]++;
	 osrc[M] = omap[lts->src[Nme][i][j]];
	 lab[M] = lts->label[Nme][i][j];
	 wdest[M] = i;
	 odest[M] = lts->dest[Nme][i][j];
	}
 for(i=0;i<Nnodes;i++) {
	free(lts->src[Nme][i]); lts->src[Nme][i]=NULL;
	free(lts->label[Nme][i]); lts->label[Nme][i]=NULL;
	free(lts->dest[Nme][i]); lts->dest[Nme][i]=NULL;
 }

 // sort on wsrc
 Count2BeginIndex(begin_to_w, Nnodes);
 SortArray(&osrc,M,wsrc,begin_to_w,Nnodes);
 SortArray(&lab,M,wsrc,begin_to_w,Nnodes);
 SortArray(&wdest,M,wsrc,begin_to_w,Nnodes);
 SortArray(&odest,M,wsrc,begin_to_w,Nnodes);
 free(wsrc);

 count_to_w = (int*)calloc(Nnodes+1, sizeof(int));
 begin_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 if ((!count_to_w)||(!begin_from_w)||(!count_from_w))
	Fatal(1,error,"out of memory in dlts_shuffle");

 // exchange osrc-wdest-odest 
 for(i=0;i<Nnodes;i++) count_to_w[i] = begin_to_w[i+1]-begin_to_w[i];
 // Warning(info,"%d: tosend %d to_myself %d",Nme, begin_to_w[Nnodes], count_to_w[Nme]);
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,lts->comm); 
 begin_from_w[0]=0; 
 for(i=1;i<=Nnodes;i++) begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
 
 // Warning(info,"%d: toreceive:%d from_myself %d",Nme,begin_from_w[Nnodes], count_from_w[Nme]);

 M = begin_from_w[Nnodes];
#ifdef DEBUG
 Warning(info,"%d : %d new trans.. %d old \n",Nme, M,begin_to_w[Nnodes]);
#endif
 if ((aux1=(int*)calloc(M+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in dlts_shuffle");
 ///// MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(osrc + begin_to_w[i], count_to_w[i], MPI_INT, 
						i, COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
	MPI_Irecv(aux1 + begin_from_w[i], count_from_w[i], MPI_INT,
						i, COLLAPSE_TAG, lts->comm,request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);
 free(osrc); osrc=aux1; aux1 = NULL;
 if ((aux4=(int*)calloc(M+1, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in dlts_shuffle");
 ///// MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(lab + begin_to_w[i], count_to_w[i], MPI_INT, i,	
						COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
	MPI_Irecv(aux4 + begin_from_w[i],count_from_w[i], MPI_INT, i, 
						COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);
 free(lab);lab=aux4;aux4=NULL;
 if((aux2=(int*)calloc(M, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in dlts_shuffle");
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(wdest + begin_to_w[i], count_to_w[i], MPI_INT, i, 
						COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
	MPI_Irecv(aux2 + begin_from_w[i],count_from_w[i], MPI_INT, i, 
						COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);
 free(wdest);wdest=aux2;aux2=NULL;

 if((aux3=(int*)calloc(M, sizeof(int)))==NULL)
	Fatal(1,error,"out of memory in dlts_shuffle");
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(odest + begin_to_w[i], count_to_w[i], MPI_INT, i, 
						COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
	MPI_Irecv(aux3 + begin_from_w[i], count_from_w[i], MPI_INT,
						i, COLLAPSE_TAG, lts->comm, request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);
 free(odest);odest=aux3;aux3=NULL;

 // sort on wdest
 for(i=0;i<=Nnodes;i++) count_to_w[i]=count_from_w[i]=begin_to_w[i] = 0;
 // sort4(wdest, lab, osrc, odest, M, Nnodes, begin_to_w);
 // sort on wdest
 ComputeCount(wdest, M, begin_to_w);
 Count2BeginIndex(begin_to_w, Nnodes);
 SortArray(&osrc,M,wdest,begin_to_w,Nnodes);
 SortArray(&lab,M,wdest,begin_to_w,Nnodes);
 SortArray(&odest,M,wdest,begin_to_w,Nnodes);

 // exchange mappings of destinations 
 //       exchange counts
 for(i=0;i<Nnodes;i++) count_to_w[i] = begin_to_w[i+1]-begin_to_w[i];
 MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,lts->comm);
 begin_from_w[0]=0; 
 for(i=1;i<=Nnodes;i++) begin_from_w[i]=begin_from_w[i-1]+count_from_w[i-1];
#ifdef DEBUG
 Warning(info,"%d: i'm the new destination of %d transitions, %d local",
 				 Nme, begin_from_w[Nnodes], begin_from_w[Nme]);
#endif
 aux1=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 aux2=(int*)calloc(begin_from_w[Nnodes], sizeof(int));
 if ((!aux1) || (!aux2))
	Fatal(1,error,"out of memory in dlts_shuffle");

 //       send odest / receive aux1  
 ///// MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(odest + begin_to_w[i],
						count_to_w[i], MPI_INT, i, COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(aux1 + begin_from_w[i],
						count_from_w[i], MPI_INT,i,COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);
#ifdef DEBUG
 Warning(info,"%d: sent odest, received aux1",Nme);
 for(i=0;i<begin_from_w[Nnodes];i++){
	assert(aux1[i] >= 0);
	assert(aux1[i] < N);
 }
#endif
 //      aux2 := wmap[aux1]
 for(i=0;i<begin_from_w[Nnodes];i++)
	aux2[i] = wmap[aux1[i]]; 
#ifdef DEBUG
 Warning(info,"%d: aux2 filled ",Nme);
#endif
 //      send aux2 / receive wdest
 MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(aux2 + begin_from_w[i],
						count_from_w[i], MPI_INT, i, COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(wdest + begin_to_w[i],
						count_to_w[i], MPI_INT,i,COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array); 
#ifdef DEBUG
 Warning(info,"%d: sent aux2, received wdest",Nme);
#endif
 //      aux2 := omap[aux1]
 for(i=0;i<begin_from_w[Nnodes];i++)
	aux2[i] = omap[aux1[i]];
 //      send aux2 / receive odest
 ///// MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++){
	MPI_Isend(aux2 + begin_from_w[i],
						count_from_w[i], MPI_INT, i, COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
	MPI_Irecv(odest + begin_to_w[i],
						count_to_w[i], MPI_INT,i,COLLAPSE_TAG, lts->comm,
						request_array + aux);
	aux++;
 }
 MPI_Waitall(aux, request_array, status_array);   
#ifdef DEBUG
 Warning(info,"%d: sent aux2, received odest",Nme);
#endif
 ///// MPI_Barrier(lts->comm);
 free(aux1); free(aux2);
#ifdef DEBUG
 Warning(info,"%3d: now re-creating lts..",Nme);
#endif
 // re-create lts
 if ((aux1=(int*)calloc(Nnodes+1, sizeof(int))) == NULL)
	Fatal(1,error,"out of memory in dlts_shuffle");
 ComputeCount(wdest,M,aux1);
 Count2BeginIndex(aux1, Nnodes);
 SortArray(&osrc,M,wdest,aux1,Nnodes);
 SortArray(&lab,M,wdest,aux1,Nnodes);
 SortArray(&odest,M,wdest,aux1,Nnodes);
 free(wdest);

 for(i=0;i<Nnodes;i++){
	lts->transition_count[Nme][i] = aux1[i+1]-aux1[i];
	lts->src[Nme][i] = osrc + aux1[i];
	lts->label[Nme][i] = lab + aux1[i];
	lts->dest[Nme][i] = odest + aux1[i];
 }
 free(aux1);
}


















/****************************************************


void dlts_clear_useless_transitions(dlts_t lts)


*****************************************************/

void dlts_clear_useless_transitions(dlts_t lts){

 int i, j,k, N, M, MU, MU_all, Nme, Nnodes;
 int *aux1, *aux2, *aux3;

 // init
 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme); 
 N = lts->state_count[Nme];
 if (Nme==0) Warning(info, "NOW COUNTING USELESS STATES and TRANSITIONS");
 // Warning(info,"%d: %d states",Nme,lts->state_count[Nme]);

 // sort
 /*
 aux1=(int*)calloc(N+1, sizeof(int));
 aux2=(int*)calloc(lts->label_count+1, sizeof(int));
 for (i = 0; i < Nnodes; i++){
	for(j = 0; j <= N; j++) aux1[j] = 0;
	M = lts->transition_count[Nme][i];
	
	ComputeCount(lts->src[Nme][i], M, aux1);
	Count2BeginIndex(aux1, N);
	
	//	Warning(info,"%d %d eee. %d transitions, %d states. aux1: %d %d %d",
	//			me, i, M, N, aux1[0], aux1[1], aux1[N]);

	//	SortArray(lts->dest[Nme]+i, M,
	//						lts->src[Nme][i], aux1, N);

	//	Warning(info,"%d %d iii. %d transitions, %d states. aux1: %d %d %d",
	//					me,i,M, N, aux1[0], aux1[1], aux1[N]);

	//	SortArray(lts->label[Nme]+i, M, 
	//					lts->src[Nme][i], aux1, N);

	// Warning(info,"%d %d ooo. %d transitions, %d states. aux1: %d %d %d",
	//			me,i,M, N,
	//			aux1[0], aux1[1], aux1[N]);
	
	//	for(j=0;j<N;j++)
	//	 for(k=aux1[j];k<aux1[j+1];k++)
	//		lts->src[Nme][i][k]=j;	
	
 */
	/*
		for(j=0;j<N;j++){
		for(k=0;k<lts->label_count+1;k++) aux2[k]=0;
		sort2(lts->label[Nme][i]+aux1[j], lts->dest[Nme][i]+aux1[j], 
		aux1[j+1] - aux1[j], lts->label_count, aux2);
		for(k=0;k<lts->label_count;k++)
		simplesort(lts->dest[Nme][i]+aux1[j]+aux2[k], aux2[k+1]-aux2[k]);
		}
	*/
 /*
 ///// MPI_Barrier(lts->comm);
 if (Nme==0) Warning(info,"%d DONE!!",i);
 }
 
 free(aux1);
 free(aux2);
 
 ///// MPI_Barrier(lts->comm);
 */
 // count
 MU = 0;
 for (i=0;i<Nnodes;i++){
	aux1=lts->src[Nme][i];
	aux2=lts->label[Nme][i];
	aux3=lts->dest[Nme][i];
	if ((aux1[0]==aux3[0])&&(aux2[0]==lts->tau)) MU++;
	for(j=1;j<lts->transition_count[Nme][i];j++)
	 if ((aux1[j]==aux3[j])&&(aux2[j]==lts->tau)) MU++;
	 else if ((aux1[j]==aux1[j-1])&&(aux2[j]==aux2[j-1])&&(aux3[j]==aux3[j-1])) MU++;
 }
 
 ///// MPI_Barrier(lts->comm);
 MPI_Allreduce(&MU, &MU_all, 1, MPI_INT, MPI_SUM, lts->comm);
 if (Nme==0) Warning(info,"at least %d transitions should be eliminated",MU_all);
 (void) N;
}









void dlts_shrink(dlts_t lts, int* wscc, int* oscc, int** weight){
// input: wscc, oscc
//        !!! must be stable, i.e. if scc(x)=y then scc(y)=y
// if weight==NULL then it doesn't get allocated and computed...
// output: lts only with nonempty states and RE-NUMBERED
//        weight, i.e. 0 for non-head states
//                      for head states: the number of states "owned"
//   !!!!!!!!!!   destroys oscc
// also the root->seg and root->ofs are actualized

 int i,j, newN, aux, maxscc, nscc;
 int *newindex, *tmp;
 int *count_from_w, *begin, *count_to_w;
 int* ooo;

 MPI_Request* request_array;
 MPI_Status* status_array;
 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme); 
 if (Nme==0) Warning(info,"\n\nSHRINK");

 request_array = (MPI_Request*)calloc(2*Nnodes, sizeof(MPI_Request));
 status_array = (MPI_Status*)calloc(2*Nnodes, sizeof(MPI_Status));

 if (Nme == lts->root_seg){
	lts->root_seg = wscc[lts->root_ofs];
	lts->root_ofs = oscc[lts->root_ofs];
 }
 // build newindex
 newindex=(int*)calloc(lts->state_count[Nme], sizeof(int));
 j=0;
 for(i=0;i<lts->state_count[Nme];i++)
	if ((wscc[i]==Nme)&&(oscc[i]==i))
	 newindex[i] = j++;
	else 
	 newindex[i] = -1;
 newN = j;

 // update src[Nme][i] , forall i and dest[Nme][Nme]
 for(i = 0; i < Nnodes;i++)
	for(j=0;j<lts->transition_count[Nme][i];j++){
#ifdef DEBUG
	 assert(lts->src[Nme][i][j] < lts->state_count[Nme]);
	 assert( newindex[lts->src[Nme][i][j]] >= 0);
#endif
	 lts->src[Nme][i][j] = newindex[lts->src[Nme][i][j]];
	}
 for(j=0;j<lts->transition_count[Nme][Nme];j++){
#ifdef DEBUG
	assert(lts->dest[Nme][Nme][j] < lts->state_count[Nme]);
	assert( newindex[lts->dest[Nme][Nme][j]] >= 0);
#endif
	lts->dest[Nme][Nme][j] = newindex[lts->dest[Nme][Nme][j]];
 }

 //  Warning(info,"%3d: exchanging counts",Nme);

 // exchange dest -> tmp
 // counts 
 begin = (int*)calloc(Nnodes+1, sizeof(int)); 
 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 MPI_Alltoall(lts->transition_count[Nme],1,MPI_INT,count_from_w,1,MPI_INT,lts->comm);
 count_from_w[Nme]=0;
 begin[0]=0; 
 for(i=1;i<=Nnodes;i++) begin[i]=begin[i-1]+count_from_w[i-1];

 //  Warning(info,"%3d: dest -> tmp",Nme);
 // dest->tmp
 tmp=(int*)calloc(begin[Nnodes], sizeof(int));
 ///// MPI_Barrier(lts->comm);
 aux=0;
 for(i=0;i<Nnodes;i++)
	if (i != Nme){
	 //	 assert(lts->dest[Nme][i] != NULL);
	 //	 Warning(info,"%3d: send to %3d %12d(%p)\n%3d: receive from %3d %12d", 
	 //					 Nme, i, lts->transition_count[Nme][i], lts->dest[Nme][i], Nme,i,count_from_w[i]);
	 MPI_Isend(lts->dest[Nme][i], lts->transition_count[Nme][i], MPI_INT, 
						 i, COLLAPSE_TAG, lts->comm, request_array + aux);
	 aux++;
	 MPI_Irecv(tmp + begin[i], count_from_w[i], MPI_INT,
						 i, COLLAPSE_TAG, lts->comm,request_array + aux);
	 aux++;
	}
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);

 //  Warning(info,"%3d: transforming tmp",Nme);
 // transform tmp
 for(i=0;i<begin[Nnodes];i++){
#ifdef DEBUG
	assert(tmp[i] < lts->state_count[Nme]);
	assert(newindex[tmp[i]] >= 0);
#endif
	tmp[i] = newindex[tmp[i]];
 }

 //  Warning(info,"%3d: tmp -> dest",Nme);


 // tmp->dest
 aux=0;
 for(i=0;i<Nnodes;i++)
	if (i != Nme){
	 MPI_Isend (tmp + begin[i], count_from_w[i], MPI_INT, 
						 i, COLLAPSE_TAG, lts->comm, request_array + aux);
	 aux++;
	 MPI_Irecv (lts->dest[Nme][i], lts->transition_count[Nme][i], MPI_INT,
						 i, COLLAPSE_TAG, lts->comm, request_array + aux);
	 aux++;
	}
 MPI_Waitall(aux, request_array, status_array);
 ///// MPI_Barrier(lts->comm);

 free(tmp);  free(begin); free(count_from_w);

 // compute weight (and !! destroy oscc)
 if (weight != NULL){
	(*weight) = (int*)calloc(newN, sizeof(int));
	begin = (int*)calloc(Nnodes+1, sizeof(int)); 
	count_to_w = (int*)calloc(Nnodes+1, sizeof(int)); 
	count_from_w = (int*)calloc(Nnodes+1, sizeof(int)); 
	ComputeCount(wscc, lts->state_count[Nme], begin);
	for(i=0;i<Nnodes;i++) count_to_w[i]=begin[i];
	Count2BeginIndex(begin, Nnodes);
	//			Warning(info,"newN=%d, *weight is %p weight is %p w[0]=%d",newN, (*weight),weight, (*weight)[0]);
	ooo = oscc; SortArray(&ooo, lts->state_count[Nme], wscc, begin, Nnodes);
	
	i=0;
	MPI_Allreduce(lts->state_count + Nme, &i, 1, MPI_INT, MPI_MAX, lts->comm);
	MPI_Alltoall(count_to_w,1,MPI_INT,count_from_w,1,MPI_INT,lts->comm);
	tmp = (int*)calloc(i, sizeof(int));
	for(i=0;i<Nnodes;i++){
	 MPI_Isend(ooo + begin[i], count_to_w[i], MPI_INT, 
						 i, WEIGHT_TAG, lts->comm, request_array);
	 MPI_Irecv(tmp, count_from_w[i], MPI_INT,
						 i, WEIGHT_TAG, lts->comm,request_array+1);
	 MPI_Wait(request_array, status_array);	 
	 MPI_Wait(request_array+1, status_array+1);
	 for(j=0;j<count_from_w[i];j++){
#ifdef DEBUG
		assert(tmp[j]<lts->state_count[Nme]);
#endif
		(*weight)[newindex[tmp[j]]]++;
	 }
	}  
	free(tmp); free(begin);free(count_to_w); free(count_from_w); 	
 }

 free(newindex);
 
 // update state_count
 
 //Warning(info,"*weight is %p weight is %p",(*weight),weight);
 if (weight != NULL) {	
	maxscc=0; nscc=0; j=0;
	for(i=0;i<newN;i++) {
	 j = j+ (*weight)[i]; 
	 if ((*weight)[i] > maxscc) maxscc = (*weight)[i];
	 if ((*weight)[i] == 1) nscc++;
	 //	 assert ((*weight)[i] >= 1);
	}
	Warning(info,"%3d: old N is %12d (weight: %12d), new N is %12d", 
				 Nme, lts->state_count[Nme], j, newN); 
	///// MPI_Barrier(lts->comm);
	i=lts->state_count[Nme]; 
	MPI_Reduce(&i, &aux, 1, MPI_INT, MPI_SUM, 0, lts->comm);
	i=maxscc; MPI_Reduce(&i, &maxscc, 1, MPI_INT, MPI_MAX, 0, lts->comm);
	i=nscc; MPI_Reduce(&i, &nscc, 1, MPI_INT, MPI_SUM, 0, lts->comm);
	if (Nme==0){
	 Warning(info,"largest SCC has %12d Nnodes",maxscc);
	 Warning(info,"there are %12d trivial SCCs", nscc);
	}
 }
 lts->state_count[Nme] = newN;
 ///// MPI_Barrier(lts->comm);
 // Warning(info,"\noooooooo\n");
}





void dlts_mark(dlts_t lts, int L, char* mark){
// input: a label 
// output: mark = 1 for all states with an outgoing transition "label"
 
 int i,j,N;

 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme); 

 // N=lts->state_count[Nme];
 for(i=0;i<Nnodes;i++)
	for(j=0; j < lts->transition_count[Nme][i]; j++)
	 if (lts->label[Nme][i][j] == L)
		mark[lts->src[Nme][i][j]] = 1;

 j=0;
 for(i=0;i<lts->state_count[Nme];i++)
	if (mark[i]) j++;
 Warning(info,"%3d: %12d states marked (out of %12d)",
				 Nme, j, lts->state_count[Nme]);
}


void dlts_reachback(dlts_t lts, char* mark){
// input: mark =.. 0 for most states, 1 for a small set
// output: mark: 1 for all states reachable (on any labels) 
//          from the initial marked set

 int *bufin, *bufout;
 int *count_from_w;
 int i,j, max, changed, allchanged;

 MPI_Request reqsend, reqrecv;
 MPI_Status stsend, strecv;

 ///// MPI_Barrier(lts->comm);
 MPI_Comm_size(lts->comm, &Nnodes);
 MPI_Comm_rank(lts->comm, &Nme); 
 if (Nme==0) Warning(info,"\nREACH BACK");

 count_from_w = (int*)calloc(Nnodes+1, sizeof(int));
 MPI_Alltoall(lts->transition_count[Nme],1,MPI_INT,count_from_w,1,MPI_INT,lts->comm);

 max=1; for(i=0;i<Nnodes;i++) if (count_from_w[i]>max) max=count_from_w[i];
 bufin=(int*)calloc(max, sizeof(int));
 max=1; for(i=0;i<Nnodes;i++) if (lts->transition_count[Nme][i]>max) max=lts->transition_count[Nme][i];
 bufout=(int*)calloc(max, sizeof(int));

 allchanged = 1;
 while (allchanged > 0){
	changed = 0;
	for(i=0;i<Nnodes;i++){
	 // exchange destinations
	 MPI_Isend(lts->dest[Nme][i], lts->transition_count[Nme][i], MPI_INT, 
						 i, REACH_TAG, lts->comm, &reqsend);
	 MPI_Irecv(bufin, count_from_w[i], MPI_INT,
						 i, REACH_TAG, lts->comm,&reqrecv);
	 MPI_Wait(&reqsend, &stsend);	 
	 MPI_Wait(&reqrecv, &strecv);
	 // write marks
	 for(j=0;j<count_from_w[i];j++)
		bufin[j] = mark[bufin[j]];
	 // exchange back
	 MPI_Isend(bufin, count_from_w[i], MPI_INT, 
						 i, REACH_TAG, lts->comm, &reqsend);
	 MPI_Irecv(bufout, lts->transition_count[Nme][i], MPI_INT,
						 i, REACH_TAG, lts->comm, &reqrecv);
	 MPI_Wait(&reqsend, &stsend);	 
	 MPI_Wait(&reqrecv, &strecv); 
	 // mark new states
	 for(j=0;j<lts->transition_count[Nme][i];j++)
		if ((!mark[lts->src[Nme][i][j]]) && (bufout[j])){
		 mark[lts->src[Nme][i][j]] = 1;
		 changed++;
		}
	} // end for i
	MPI_Allreduce(&changed, &allchanged, 1, MPI_INT, MPI_SUM, lts->comm );
	//	if (Nme==0) Warning(info,"new reachable: %12d",allchanged);
 } // end while changed
 free(bufin); free(bufout); free(count_from_w);
 if (Nme==0) Warning(info,"\n\n");
}







void print_status(taudlts_t t){
 int nh,mh,nf,mf,n,m,sn,sm,i,j;
 ///// MPI_Barrier(t->comm);
 MPI_Comm_rank(t->comm, &Nme);
 MPI_Reduce(&Nhooked, &nh, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&Mhooked, &mh, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&Nfinal, &nf, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&Mfinal, &mf, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&Nleft, &n, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&Mleft, &m, 1, MPI_INT, MPI_SUM, 0, t->comm);
 MPI_Reduce(&(t->M), &sm, 1, MPI_INT, MPI_SUM, 0, t->comm);
 j=0;
 for(i=0;i<t->N;i++)
	if (t->begin[i]<t->begin[i+1]) j++;
 MPI_Reduce(&j, &sn, 1, MPI_INT, MPI_SUM, 0, t->comm);
 if (Nme==0) 
	Warning(info,"\n_____________\nleft:      N=%12d    M=%12d     \nfinal:     N=%12d    M=%12d    \non cycles: N=%12d    M=%12d  \n\"real\"     N=%12d    M=%12d\n_____________\n",n,m,nf,mf,nh,mh,sn,sm);
 ///// MPI_Barrier(t->comm);
}








