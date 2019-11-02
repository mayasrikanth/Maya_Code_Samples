#ifndef QUEUE_H
#define QUEUE_H
#include "virtualmem.h"

/* Note to user: feel free to change the typedef size of age_t. */
typedef uint8_t age_t;
/*! A single node of a queue of resident "pages". */
typedef struct _queuenode {
    /* page member in the QueueNode. */
    page_t page;
    /* Age member in QueueNode. */
    age_t age;
    struct _queuenode *prev;
    struct _queuenode *next;
} QueueNode;



/*!
 * A queue of Queuenodes representing loaded pages.  This type is used to
 * keep track of pages in various states within the user-space
 * threading library.
 */
typedef struct _queue {
    QueueNode *head;  /*!< The first "page" in the queue. */
    QueueNode *tail;  /*!< The last "page" in the queue. */
} Queue;


int aging_queue_empty(Queue *queuep);
void aging_queue_append(Queue *queuep, page_t page);
page_t aging_queue_take(Queue *queuep);
int queue_remove_aging(Queue *queuep, QueueNode *nodep);
QueueNode *aging_lowest_queue(Queue *queuep);


#endif /* QUEUE_H */
