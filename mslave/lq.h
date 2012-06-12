#ifndef LQH_
#define LQH_

struct qnode
{
  int data;
  struct qnode *next;
};
struct qnode * processors;

void enqueue(struct qnode **first,int n);
int dequeue(struct qnode **first);
void displayq(struct qnode *p);
void freeq (struct qnode  **first); 


#endif
