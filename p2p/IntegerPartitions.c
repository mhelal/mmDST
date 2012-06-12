/*
**  partitions.c
*/
/*
**  October 1, 1997
**  This program implements Algorithm 3.1 - 3.9
*/
/*
**  Compile with:
**	gcc -O -c IntegerPartitions.c
**	gcc -O IntegerPartitions.o -oIntegerPartitions
**
**  Run with:
**      IntegerPartitions
*/
#include<stdio.h>
#include<stdlib.h>
#define false 0;
#define true 1;
#define Min(x,y) ((x<y)?x:y)
  typedef struct partitionnode 
  {
    int m;
    int numparts;
    int parts[101];
  } partition;
  int Pn[40][40];
  int P[50];
  
void Ferrers(partition a)
/*
**  print out the Ferrers diagram of a partition
*/
{
 int i,j,k,n;
 printf("\n");
 n= a.numparts;
 for(i=1;i<=n;i=i+1)
 {
    j = a.parts[i];
    for(k=1;k<=j;k=k+1)
     printf("*");
    printf("\n");
 }
 printf("\n");
}

void output(partition a)
/*
**  print out the partition a
*/
{
 int i,j;
 printf("%d = ",a.m);
 for(i=1;i<=a.numparts;i=i+1)
 {
    j = a.parts[i];
    printf("%d",j);
    if (i != a.numparts)
      printf("+");
 }
}

void RecPartition(partition a, int m, int B, int N)
/*
**  See Algorithm 3.1
**  generate all partitions with largest part of size at most B
**  recursive procedure used in other procedures
*/
{
 int i;
 if(m==0)
 {
   a.numparts = N;   
   output(a); printf("\n");
 }
 else
 {
   for (i=1;i<=Min(B,m);i=i+1)
   {
     a.parts[N+1] = i;
     RecPartition(a,m-i,i,N+1);
   }
 }
}

void GenPartitions(int m)
/*
**  Algorithm 3.1
**  generate all partitions of m
*/
{
 partition a;
 a.m =  m;
 RecPartition(a,m,m,0);
}

void ConjPartition(partition a, partition b)
/*
**  Algorithm 3.2
**  conjugate a given partition a, producing b
*/ 
{ 
 int i,j,n_prime; 
 b.m =  a.m;

 n_prime =  a.parts[1];
 for (i=1;i<=n_prime;i=i+1) b.parts[i] = 1; 
 b.numparts = n_prime;
 for (j=2;j<=a.numparts;j=j+1)
 {
   for(i=1;i<=a.parts[j];i=i+1)
     b.parts[i] =  b.parts[i] + 1;
 }
}


void GenPartitions2(int m, int n)
/*
** Algorithm 3.3
** generate all partitions of m with largest part of size n
*/
{
 partition a;
 a.m =  m;
 a.parts[1] = n;
 RecPartition(a,m-n,n,1);
}

void RecPartition2(partition a, int m, int B, int N)
/*
**  See  Algorithm 3.3
**  generate all partitions with largest part of size at most B
**  recursive procedure used in other procedures
*/
{
 int i;
 partition b;
 if(m==0)
 {
   a.numparts = N;   
   output(b); printf("\n");
 }
 else
 {
   for (i=1;i<=Min(B,m);i=i+1)
   {
     a.parts[N+1] = i;
     RecPartition(a,m-i,i,N+1);
   }
 }
}
void GenPartitions3(int m,int n)
/*
**  Algorithm 3.4
**  generate all partitions of m with n parts
*/
{
 partition a;
 a.m =  m;
 a.parts[1] = n;
 RecPartition2(a,m-n,n,1);
}





void EnumPartitions(int m, int n )
/*
**  Algorithm 3.5
**  compute the partition numbers P[i] and P[i,j] for i <= m, j <= n
*/
{
 int i,j;
 Pn[0][0] =  1;
 P[0] = 1;
 for(i=1;i<=m;i=i+1) Pn[n][0]=0;
 for(i=1;i<=m;i=i+1)
 {
   for (j=1;j<=Min(i,m);j=j+1)
   {
       if (i < (2*j)) 
         Pn[i][j] =  Pn[i-1][j-1];
       else 
         Pn[i][j] =  Pn[i-1][j-1] + Pn[i-j][j];
     P[i] = P[i] + Pn[i][j];
   }
   P[i] = 0;
   for(j=1;j<=m;j++) P[i] = P[i] + Pn[i][j];
 }
}

void EnumPartitions2(int m)
/*
**  Algorithm 3.6
**  compute the partition numbers P[i] for i <= m
*/
{
 int i,j,sum,omega,omega2,sign;
 P[0] = 1;
 P[1] = 1;
 for (i=2;i<=m;i=i+1)
 {
   sign = 1;
   sum = 0;
   omega = 1;
   j = 1;   
   omega2 =  omega+j;
   while(omega <= i)
   {
     if(sign==1)
      sum =  sum + P[i-omega];
     else 
      sum =  sum - P[i-omega];
     if(omega2 <= i)
     {
       if(sign==1)
         sum =  sum + P[i-omega2];
       else 
         sum =  sum - P[i-omega2];
     }
     omega =  omega + 3*j+1;
     j = j+1;
     omega2 =  omega+j;
     sign =  0 - sign;
   }
   P[i] =  sum;
 }
}

void PartitionLexUnrank(int m, int n, int r, partition *a)
/*
**  Algorithm 3.9
**  find the partition of m into n parts having rank r
*/
{
  int i;
  EnumPartitions(m,n);
  (*a).m = m;
  (*a).numparts = n;
  for(i=1;i<=n;i=i+1)
   (*a).parts[i] = 0;
  while (m > 0)
  {
      if (r < Pn[m-1][n-1])
      {
          (*a).parts[n] = (*a).parts[n]+1;
          m = m-1;
          n = n-1;
      }
      else 
      {
         for (i=1;i<=n;i=i+1)
           (*a).parts[i] = (*a).parts[i] + 1;
         r = r - Pn[m-1][n-1];
         m = m - n;
      }
  }
}


int PartitionLexRank(int m,int n,partition a)
/*
**  Algorithm 3.8
**  find the rank of partition a,
**  where a is given in standard form
*/
{
  int i,r;
  partition b;
  EnumPartitions(m,n);
  b =  a;
  m =  b.m;
  n =  b.numparts;
  r = 0;
  while (m > 0)
  {
     if (b.parts[n] == 1)
     {
         m =  m-1;
         n =  n-1;
     }
     else
     {         
        for(i=1;i<=n;i=i+1)
          b.parts[i]  =  b.parts[i] - 1;
        r =  r + Pn[m-1][n-1];
        m =  m-n;         
     }
   }
   return(r);
}

void PartitionLexSuccessor(int m,int n,partition *a,int *flag)
/*
**  Algorithm 3.7
**  replaces the partition a by its successor,
**  where a is given in standard form
*/
{
  int i,j,d;

  i =  2;
  while( (i<=n) && ((*a).parts[1]<=((*a).parts[i]+1)) )
    i=i+1;
  if(i == (n+1))
  {
    *flag = false;
  }
  else
  {
     (*a).parts[i] = (*a).parts[i] + 1;
     d = -1;
     for(j=(i-1);j>=2;j=j-1)
     {
         d = d + (*a).parts[j] - (*a).parts[i];
         (*a).parts[j]  =  (*a).parts[i];
     }
     (*a).parts[1] =  (*a).parts[1] + d;
  }
}
/*
int main()
{ 
  int m,n,i,j,r,s,num;
  int flag;
  partition a;
  char junk;
  m=7;
  printf("Test of Algorithm 3.1 with m=%d\n\n",m);
  GenPartitions(m);
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);


  m=10;
  n=4;
  printf("Test of Algorithm 3.3 with m=%d n=%d \n\n",m,n);
  GenPartitions2(m,n);
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);


  m=10;
  n=4;
  GenPartitions3(m,n);
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);


  m =  30;
  n =  10;
  printf("Test of Algorithm 3.5 with m=%d n=%d \n\n",m,n);
  EnumPartitions(m,n);
  printf("TABLE 3.1\nA table of partition numbers\n");
  printf("                 ");
  for(j=1;j<=(n/2);j++) printf("     ");
  printf("p(i,j)\n");
  printf("      |  p(i) |  n=1");
  for(j=2;j<=n;j++) printf("%5d",j);
  printf("\n");
  printf("---------------"); for(j=1;j<=n;j++) printf("-----");
  printf("\n");
  for (i=1;i<=m;i=i+1)
  { 
   printf("i=%3d",i);
   printf(" | ");
   printf("%5d |",P[i]);
   for (j=1;j<=Min(i,n);j=j+1)
   {   
     printf("%5d",Pn[i][j]);
   }
   printf("\n");
  }
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);


  m =  20;
  printf("Test of Algorithm 3.6 with m=%d \n\n",m);
  EnumPartitions2(m);
  printf("  i  |"); for(i=1;i<=m;i=i+1) printf("%5d",  i ); printf("\n");
  printf("-----|"); for(i=1;i<=m;i=i+1) printf("-----"); printf("\n");
  printf("P[i] |"); for(i=1;i<=m;i=i+1) printf("%5d",P[i]); printf("\n");
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);

[2,1,3,5,4]     1

  m = 10;
  n = 4;
  printf("Test of Algorithm 3.7 with m=%d n=%d \n\n",m,n);
  a.m = m;
  a.numparts = n;
  a.parts[1] =  m-n+1;
  for(i=2;i<=n;i=i+1)  
    a.parts[i] = 1;
  flag =  true;
  while(flag)
  {     
     output(a);
     r=PartitionLexRank(m,n,a);
     printf(": rank = %d\n",r);     
     PartitionLexSuccessor(m,n,&a,&flag);     
  }
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);


  m = 13;
  n = 5;
  printf("Test of Algorithm 3.8 with m=%d n=%d \n\n",m,n);
  EnumPartitions(m,n);
  num =  Pn[m][n];
  for(r=0;r<num;r=r+1)
  {
     printf("%4d: ",r);
     PartitionLexUnrank(m,n,r,&a);
     output(a);
     s=PartitionLexRank(m,n,a);
     printf(" rank=%4d",s);
     printf("\n");
  }
  printf("\nEnd of test.  Press enter for next test\n"); scanf("%c",&junk);

  m = 10;
  n = 4;
  printf("\nTest of Algorithm 3.9 with m=%d n=%d \n\n",m,n);
  for(r=0;r<=8;r=r+1)
  {
    PartitionLexUnrank(m,n,r,&a);
    printf(" rank = %3d: ",r);
    output(a); printf("\n");
  }
  printf("\nEnd of all tests.\n\n");

  return(0);
}
*/

int main()
{ 
  partition a;
  int m, n, i, j, num, s; 
  m = 13;
  n = 5;
  printf("Test of Algorithm 3.8 with m=%d n=%d \n\n",m,n);
  EnumPartitions(m,n);
  num =  Pn[m][n];
  for(i=0;i<num;i++)
  {
  for(j=0;j<n;j++)
  {
     printf(" %4d",Pn[i][j]);
  }
     printf("\n");
  }
  return(0);
}

