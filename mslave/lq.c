#include <stdio.h>
#include <stdlib.h>
#include "lq.h"
#include "globals.h"
#include "utils.h"


int inQueue (struct qnode *p, int element) {
  int found;
  found = 0;
  if(p!=NULL)
    {
      while(p!=NULL)
	{
          if (p->data == element)
	    found = 1;
          p=p->next;
	}
    }
  
  return found;
}

void enqueue(struct qnode **first, int n)
{
  struct qnode *p=*first,*q;
  if (inQueue((*first), n) == 0) {
    q=(struct qnode *)malloc(sizeof(struct qnode));
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
}


int dequeue(struct qnode **first)
{
  int temp;
  struct qnode *p=*first;
 
  if(*first==NULL)
    {
      mprintf(1, "\nqueue is empty\n", threadnum);
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
	char msg [100];
  if(p==NULL)
    {
      mprintf(1, "\nqueue is empty\n", threadnum);
    }
  else
    {
      while(p!=NULL)
	{
          sprintf(msg,"\n%d",p->data);
					mprintf (1, msg, threadnum);
          p=p->next;
	}
    }
  mprintf(1, "\n", threadnum);
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

