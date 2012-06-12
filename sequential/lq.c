#include <stdio.h>
#include <stdlib.h>
#include "lq.h"


void enqueue(struct qnode **first,int n)
{
  struct qnode *p=*first,*q;
  q=(struct qnode *)mmalloc(sizeof(struct qnode));
  if(*first==NULL)
    {
      *first=q;
      q->data=n;
      q->next=NULL;
    }
  else
    {
      while(p->next!=NULL)
	{
	  p=p->next;
	}
      p->next=q;
      q->data=n;
      q->next=NULL;
    }
}

int dequeue(struct qnode **first)
{
  int temp;
  struct qnode *p=*first;
  if(*first==NULL)
    {
      //printf("\nqueue is empty\n");
      return(-1);
    }
 else
   {
     p=*first;
     *first=(*first)->next;
     temp=p->data;
     free(p);
     return(temp);
   }
}

void displayq(struct qnode *p)
{
  if(p==NULL)
    {
      printf("\nqueue is empty\n");
    }
  else
    {
      while(p!=NULL)
	{
          printf("\n%d",p->data);
          p=p->next;
	}
    }
  printf("\n");
}

void freeq (struct qnode  **first) 
{
  struct qnode *p;
  
  while(*first!=NULL)
    {
      p=*first;
      *first=(*first)->next;
      free(p);
    }
  
}

