/*
**  permlex.c
*/
/*
**  September 9, 1997
**  this program implements Algorithm 2.14, 2.15 and 2.16
**  successor, ranking and unranking algorithm for permutations
**  in lex order
*/
/*
**  Compile with:
**	gcc -O -c permlex.c
**	gcc -O permlex.o -opermlex
**
**  Run with:
**	permlex n
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define mod %
#define true 1
#define false 0
 typedef int * permutation;
 permutation rho;
 int * f; /*factorials*/

void PrintPerm( int n, permutation T)
/*
**  Prints the permutation T to stdout.
*/
{
  int j;
  printf("[%d",T[1]);
  for(j=2;j<=n;j++) printf(",%d",T[j]);
  printf("]");
}

void PermLexSuccessor(int n ,permutation pi, int * flag)
/*
**  Algorithm 2.14
**
**  replaces pi by its successor
**  flag is set to false if pi is the last permutation
*/
{
 int h,i,j,t;
 pi[0] = 0;
 i = n-1;
 while (pi[i+1] < pi[i])
   i = i-1;
 if (i == 0)
  (*flag) = false;
 else
  {
   (*flag) = true;
   j = n;
   while (pi[j] < pi[i])
     j = j-1;
   t = pi[j];
   pi[j] = pi[i];
   pi[i] = t;
   for(h=i+1;h<=n;h=h+1)
     rho[h] = pi[h];
   for(h=i+1;h<=n;h=h+1)
     pi[h] = rho[n+i+1-h];     
 }
}

int PermLexRank(int n,permutation pi)
/*
**  Algorithm 2.15
**
**  Returns r, the rank of permutation pi.
*/
{
 int i,j,r;
 r = 0;
 for(j=1;j<=n;j=j+1) rho[j]=pi[j];
 for(j=1;j<=n;j=j+1)
 {
   r = r + (rho[j] - 1) * f[n-j];
   for(i=j+1;i<=n;i=i+1) 
   {
     if (rho[i] > rho[j])
       rho[i] = rho[i] - 1;
   }
 }
 return(r);
}



void PermLexUnrank(int n, int r,permutation pi)
/*
**  Algorithm 2.16
**
**  Returns the permutation pi that has rank r.
*/
{
 int i,j,d;
 pi[n] = 1;
 for(j=1;j<n;j=j+1)
 {
   d = (r mod f[j+1])/f[j];
   r = r - d * f[j];
   pi[n-j] = d+1;
   for(i=n-j+1;i<=n;i=i+1)
   {
     if(pi[i] > d)
        pi[i] = pi[i] + 1;
   }
  }
}

void Init(int n)
/*
**  Initialize factorial array f and allocate 
**  storage for  the permutation rho.
*/
{
 int j;
 f=(int *) calloc(n+1,sizeof(int));
 f[0] = 1;
 printf ("middle %d = %f\n", n, ceil((double) n/2));
 for(j=0;j<=ceil((double)  n/2)-1;j=j+1) {
	f[j] = j*1;
	f[n-j-1] = f[j];
	printf ("f[%d] = %d, f[%d] = %d\n", j, f[j], n-j-1, f[n-j-1]);
 }
 rho=(permutation) calloc(n+1,sizeof(int));
}
 

int main(int ac,char *av[])
{
  int n,r,i;
  int flag;
  permutation pi;

  if(ac!=2)
  {
    fprintf(stderr,"Usage: %s n\n",av[0]);
    exit(1);
  }
  n=atoi(av[1]);
  Init(n);
  pi=(permutation) calloc(n+1,sizeof(int));
  printf("Testing rank and unrank\n");
  for(i=0;i<f[n];i=i+1) 
  {
    PermLexUnrank(n,i,pi);
    printf("%4d: ",i);
    PrintPerm(n,pi);
    r=PermLexRank(n,pi);
    printf("	rank = %d\n",r);
  }
 printf("Testing the successor function\n");
 for(i=1;i<=n;i++)
   pi[i] = i; 
 flag = true;
 while(flag)
 {    
   PrintPerm(n,pi);
   printf("	");
   r=PermLexRank(n,pi);
   printf("%d\n",r);
   PermLexSuccessor(n,pi,&flag);
 }
 return(0);
}
