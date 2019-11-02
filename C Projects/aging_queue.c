#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "aging_queue.h"



/*!
 * Returns true (1) if the specified queue is empty.  Otherwise, returns
 * false (0).
 */
int aging_queue_empty(Queue *queuep) {
    assert(queuep != NULL);
    return (queuep->head == NULL);
}

/* Traverse the queue of all resident pages and  and return the QueueNode
 * representing the page with the lowest age value.
 */
QueueNode *aging_lowest_queue(Queue *queuep) {
    QueueNode *lowest = queuep->head;
    age_t lowest_age = queuep->head->age;
    QueueNode *curr = queuep->head->next;
    while (curr != NULL) {
        if(curr->age < lowest_age) {
            lowest_age = curr->age;
            lowest = curr;
        }
        curr = curr->next;
    }
    return lowest;

}


/*!
 * Add the QueueNode representing the page to front of the queue.
 * If the queue is empty, add the singleton element.  Otherwise, add the
 * element as the tail.
 */
void aging_queue_append(Queue *queuep, page_t page) {
    QueueNode *nodep = (QueueNode *) malloc(sizeof(QueueNode));
    if (nodep == NULL) {
        fprintf(stderr, "Couldn't allocate QueueNode\n");
        abort();
    }

    nodep->page = page;
    age_t age = 1;
    /* Setting topmost bit to 1. */
    age = (age << (8 * sizeof(age_t) -1));
    nodep->age = age;

    if(queuep->head == NULL) {
        nodep->prev = NULL;
        nodep->next = NULL;
        queuep->head = nodep;
        queuep->tail = nodep;
    }
    else {
        queuep->tail->next = nodep;
        nodep->prev = queuep->tail;
        nodep->next = NULL;
        queuep->tail = nodep;
    }
}


/*!
 * Get the first "page" from the queue.  Returns NULL if the queue is empty.
 */
page_t aging_queue_take(Queue *queuep) {
    QueueNode *nodep;
    page_t page;

    assert(queuep != NULL);

    /* Return NULL if the queue is empty */
    if(queuep->head == NULL)
        return -1;

    /* Go to the final element */
    nodep = queuep->head;
    if(nodep == queuep->tail) {
        queuep->head = NULL;
        queuep->tail = NULL;
    }
    else {
        nodep->next->prev = NULL;
        queuep->head = nodep->next;
    }

    page = nodep->page;
    free(nodep);
    return page;
}


/*!
 * Remove a specified node from a queue.
 *
 * Returns 1 if the node was found in the queue, or 0 if the node was not
 * found in the queue.
 *
 * NOTE:  DO NOT use this operation in an assertion, e.g.
 *            assert(queue_remove(somequeuep, somethreadp));
 *        Assertions may be compiled out of a program, and if they are, any
 *        side-effects in the assertion's test will also be compiled out.
 */
int queue_remove_aging(Queue *queuep, QueueNode *nodep) {

    assert(queuep != NULL);

    if (nodep->prev != NULL)
        nodep->prev->next = nodep->next;
    if (nodep->next != NULL)
        nodep->next->prev = nodep->prev;

    /* Reset head and tail pointers */
    if(queuep->head == nodep)
        queuep->head = nodep->next;
    if(queuep->tail == nodep)
        queuep->tail = nodep->prev;

    /* Delete the node. */
    free(nodep);

    /* We removed a node so return 1. */
    return 1;
}
