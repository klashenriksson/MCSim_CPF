typedef struct int_queue {
    int* array;
    int max_size;
    int* front;
    int* back;
} int_queue_t;

int_queue_t int_queue_create(size_t max_elements);
int int_queue_size(int_queue_t* q);
int int_queue_front(int_queue_t* q);
int int_queue_push_back(int_queue_t* q, int i);
int int_queue_push_back_array(int_queue_t* q, int* array, int len);
int int_queue_pop_front(int_queue_t* q);
int int_queue_pop_front_array(int_queue_t* q, int* arr, int len);
void int_queue_destroy(int_queue_t* q);
void int_queue_empty(int_queue_t* q);