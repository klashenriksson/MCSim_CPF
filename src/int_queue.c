#include <stdlib.h>
#include <stdio.h>

#include "int_queue.h"

int_queue_t int_queue_create(size_t max_elements) {
    int* array = (int*)malloc(max_elements * sizeof(int));

    int_queue_t queue;
    queue.array = array;
    queue.back = array;
    queue.front = array;
    queue.max_size = max_elements;
    return queue;
}

int int_queue_size(int_queue_t* q)
{
    return (int)(q->back - q->front);
}

int int_queue_front(int_queue_t* q)
{
    return *q->front;
}

int int_queue_push_back(int_queue_t* q, int i)
{
    if ((q->back - q->array) >= q->max_size)
    {
        printf("Shit! Hit max size..");
    }
    *(q->back) = i;
    q->back += 1;
    //TODO maybe wrap here if needed..

    return 0;
}

int int_queue_pop_front(int_queue_t* q)
{
    if ((q->front - q->array) >= q->max_size)
    {
        printf("Shit! Hit max size..");
    }
    int val = *q->front;
    q->front += 1;
    //TODO: maybe wrap here if needed...
    return val;
}

void int_queue_destroy(int_queue_t* q)
{
    free(q->array);
    q->array = NULL;
    q->back = NULL;
    q->front = NULL;
    q->max_size = 0;
}

void int_queue_empty(int_queue_t* q)
{
    q->front = q->array;
    q->back = q->array;
}