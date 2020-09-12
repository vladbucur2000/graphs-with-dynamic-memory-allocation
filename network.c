# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <ctype.h>
# include <limits.h>
# define NODES_MAX 50 // number of nodes
# define inf INT_MAX
# define SELECTED_NODE 1 // source node for dijkstra
# define SELECTED_NODE_BFS 1 // source node for BFS

struct pair {
    int x, pos;
};

struct input {
    bool edg;
    int val; 
};

struct Node {
    int vertex, cost; //vertex number and the cost of the node
    struct Node *next; 
};

struct Vector {
    struct Node *first; // an array of nodes (adjacency list)
};

struct Graph {
    int nrVer; //number of vertices
    struct Vector *adj; //an array of Vector (an array of adjacency lists) 
    bool *check;
};

struct Queue { // FIRST IN FIRST OUT 
    int front, back, size, capacity; //position for the first and last element in the array, size and allocated capacity
    int *v; //array
};

struct HeapNode {
    int vertex, cost;
};

struct Heap {
    struct HeapNode *array; // array with HeapNodes
    int size;
};

struct Stack {
    int *array;
    int size, capacity;
};

typedef struct Graph Graph;
typedef struct Node Node;
typedef struct Vector Vector;
typedef struct input input;
typedef struct pair pair;
typedef struct Queue Queue;
typedef struct Heap Heap;
typedef struct HeapNode HeapNode;
typedef struct Stack Stack;

void invalidInput() {
    printf("Invalid input.\n");
    exit(1);
}

int max(int x, int y) {
    if (x > y) return x;
    return y;
}

//converts chars from the input in integers 
pair charToInt(const char s[], int j) {
    bool OK = true;
    int n = strlen(s), x = 0;

    while (j < n && OK) {
        ++j;
        if (isdigit(s[j])) x = x * 10 + (s[j] - '0');
            else OK = false;
    }

    pair ans;
    ans.pos = j, ans.x = x;

    return ans;
}

//==========STACK FUNCTIONS==========
Stack *createStack(int capacity) {
    Stack *st = (Stack *) malloc(sizeof(Stack));
    st -> array = (int *) malloc (sizeof(int) * capacity);
    st -> size = 0;
    st -> capacity = capacity;

    for (int i = 0; i < st -> capacity; ++i) st -> array[i] = 0;
    
    return st;
}

void freeStack(Stack *st) {
    free (st -> array);
    free (st);
}

bool isStackFull(Stack *st) {
    if (st -> size == st -> capacity - 1) return true;
    return false;
}

bool isStackEmpty(Stack *st) {
    if (!(st -> size)) return true;
    return false;
}

void pushStack(Stack *st, int val) {
    if (isStackFull(st)) exit(1);

    st -> size++;
    st -> array[st -> size] = val;
}

void popStack(Stack *st) {
    if (isStackEmpty(st)) exit(1);

    st -> size--;
}

//==========HEAP MIN FUNCTIONS==========

//alocates memory to heap structure and array of HeapNodes
Heap *createHeap (int capacity) {
    Heap *pq = (Heap *)malloc(sizeof(Heap));
    pq -> size = 0;
    pq -> array = (HeapNode *) malloc(sizeof(HeapNode) * capacity);

    for (int i = 0; i < capacity; ++i) //initialisation of the heap array
         pq -> array[i].cost =  pq -> array[i].vertex = 0;
    
    return pq;
}

void freeHeap (Heap *pq) {
    free(pq -> array);
    free (pq);
}

//create HeapNode adding cost and vertex
HeapNode createHeapNode(int x, int y) {
    HeapNode node;
    node.cost = x, node.vertex = y;
    return node;
}

//swap HeapNodes
void swapHeap(Heap *pq, int x, int y) {
    HeapNode aux = (pq -> array[x]);
    pq -> array[x] = (pq -> array[y]);
    pq -> array[y] = aux;
}

//update minheap when inserting a new element
void heapUp(Heap *pq, int x) {
    if (x <= 1) return;
    if (pq -> array[x].cost < pq -> array[x/2].cost) {
        swapHeap(pq, x, x/2);
        heapUp(pq, x/2);
    }
}

//add HeapNode to minheap
void addHeap(Heap *pq, HeapNode new) {
    pq -> size++;
    int n = pq -> size;
    pq -> array[n] = new;
    heapUp(pq, n);
}

//update minheap after deleting the first/smallest element
void heapDown(Heap *pq, int n, int size) {
    int left = 2 * n, right = 2 * n + 1;
    int cost = pq -> array[n].cost;
    if (left <= size && cost > pq -> array[left].cost) {
                swapHeap(pq, n, left);
                heapDown(pq, left, size);
    } else if (right <= size && cost > pq -> array[right].cost) {
                swapHeap(pq, n, right);
                heapDown(pq, right, size);
    } 
}

//deletes the first/smallest element
void popHeap(Heap *pq) {
    swapHeap(pq, 1, pq -> size);
    pq -> size--;
    heapDown(pq, 1, pq -> size);
}

//returns the first/smallest element in the minheap
HeapNode getMin(Heap *pq) {
    return pq->array[1]; 
}

//==========QUEUE FUNCTIONS==========

//allocates memory to queue structure and array
Queue *createQueue(int capacity) {
    Queue *q = (Queue *) malloc(sizeof(Queue));
    q -> capacity = capacity; 
    q -> front = 0, q -> size = 0, q -> back = capacity - 1;
    q -> v = (int *) malloc(sizeof(int) * (q -> capacity));
    return q;
}

// checks if queue size > capacity
bool fullQueue(Queue *q) { 
    if (q -> capacity == q -> size) return true;
    return false;
}

// check if queue size == 0
bool emptyQueue(Queue *q) { 
    if (!(q -> size)) return true;
    return false;
}

//push element in queue
Queue *pushBack(Queue *q, int val) { 
    if (fullQueue(q)) {
        printf("ERROR queue full\n");
        exit(1);
    }

    int mod = q -> capacity;
    q -> back = (q -> back + 1) % mod;
    q -> v[q -> back] = val;
    q -> size++;

    return q;
}

//get the first element waiting in the queue
int getFront(Queue *q) {
    return q -> v[q -> front];
}

//removes the first element waiting in the queue
Queue *popFront(Queue *q) {
    if (emptyQueue(q)) {
        printf("ERROR queue empty\n");
        exit(1);
    }
    int mod = q -> capacity;
    q -> front = (q -> front + 1) %  mod;
    q -> size--;

    return q;
}

void freeQueue(Queue *q) {
    free (q -> v);
    free (q);
}

// ========== GRAPH FUNCTIONS ==========

void freeGraph(Graph *G) {
    int n = G -> nrVer;
    for (int i = 0; i <= n; ++i) {
        Node *x = G -> adj[i].first;
        while (x != NULL) {
            Node *cpy = x;
            x = x -> next;
            free(cpy);
        }
    }
    free(G -> adj);
    free (G -> check);
    free(G);
}

//alocates memory for all the nodes in the graph
Graph *createGraph(int nr) {
   
    Graph *G = (Graph *)malloc(sizeof(Graph));
    G -> nrVer = nr;
    G -> adj = malloc(sizeof(Vector) * (nr+1));

    for (int i = 0; i <= nr; ++i)
      G -> adj[i].first = NULL;

    G -> check = malloc(sizeof(int) * (nr + 1));

  return G;
}

//adds an edge to the graph
Graph *addEdge(int i, int j, int cost, Graph *G) {
    
    Node *node = (Node *)malloc(sizeof(Node));
    node -> vertex = j;
    node -> cost = cost;
    node -> next = G -> adj[i].first;
    G -> adj[i].first = node;
    
    return G;
}

//adds the edges given in the input
Graph *addEdges(input ok[NODES_MAX + 5][NODES_MAX + 5], int n, Graph *graph, int reverse) {

    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            if (ok[i][j].edg) {
                if (!reverse) graph = addEdge(i, j, ok[i][j].val, graph);
                    else graph = addEdge(j, i, ok[i][j].val, graph);
            }

    return graph;
}

//from input => graph (adjacency lists)
Graph *readAndMakeGraph(const char s[], int reverse) {
  
  int n = strlen(s), i, j, x, y, z, nrNodes = 0;
  input ok[NODES_MAX + 5][NODES_MAX + 5];

  //initialisation for the "ok" which holds the values of the input 
  // ok[x][y].edg keeps if there is it and edges between node x and node y (true if it is, false if it isn't)
  // ok[x][y].val keeps the cost of the x-y edge. it is initialised with 1.

  for (i = 1; i <= NODES_MAX; ++i) 
    for (j = 1; j <= NODES_MAX; ++j)
        ok[i][j].edg = false, ok[i][j].val = 1;

  for (i = 0; i < n; ++i) {
      x = 0, y = 0, z = 0;

      j = i - 1;
      pair per = charToInt(s, j);
      x = per.x, j = per.pos; 
 
      if (s[j] != '-') invalidInput();
      
      per = charToInt(s, j);
      y = per.x, j = per.pos;

      if (s[j] == '/') {
          per = charToInt(s, j);
          z = per.x, j = per.pos; 
      }

      if (j != n && s[j] != ',') invalidInput();
      //if (x == y) invalidInput();
      
      ok[x][y].edg = true;
      if (z) ok[x][y].val = z;

      nrNodes = max(nrNodes, x);
      nrNodes = max(nrNodes, y);
    
      i = j;
      if (s[i] != ',' && i != n) invalidInput();
  }
  
   Graph *graph = createGraph(nrNodes);
   graph = addEdges(ok, nrNodes, graph, reverse);
   
   return graph;
}

//classic bfs with queue
void bfs(Graph *G, int node) {
    bool select[G -> nrVer + 1];
    
    for (int i = 1; i <= G -> nrVer; ++i) select[i] = false; // initialisation of select
    select[node] = true;

    Queue *q = createQueue(G -> nrVer + 1); //create queue 
    q = pushBack(q, node); // push source node

    int lvl[G -> nrVer + 1];
    for (int i = 1; i <= G -> nrVer; ++i) lvl[i] = -1; // initialisation of lvl (level array) -> (-1) means no edges to reach from node to i
    lvl[node] = 0;

    while (!emptyQueue(q)) {
        int x = getFront(q);
        Node *current = G -> adj[x].first;

        while (current != NULL) {

           if (!select[current -> vertex]) {
                select[current -> vertex] = true; 
                q = pushBack(q, current -> vertex);
                lvl[current -> vertex] = lvl[x] + 1;
           }
           current = current -> next;
        }
        
        q = popFront(q);
    }

    for (int i = 1; i <= G -> nrVer; ++i)
        if (lvl[i] == -1) printf("from %d -> %d you don't have a path\n", SELECTED_NODE_BFS, i);
            else printf("from %d -> %d you have to go through %d edges\n", SELECTED_NODE_BFS, i, lvl[i]);
    
    freeQueue(q);
}

//finds the minimum number of edges to go through from SELECTED_NODE_BFS to the others -(without a cost per edge)
void shortestPath(Graph *G, int n) {
    printf("Shortest path from node %d to the other nodes:\n", SELECTED_NODE_BFS);
    bfs(G, SELECTED_NODE_BFS);
}

//classic dijkstra algorithm with minheap
void dijkstra(Graph *G, int dp[G -> nrVer + 1], bool OK) {
    Heap *pq = createHeap(G -> nrVer + 1);
    
    for (int i = 1; i <= G -> nrVer; ++i)
        dp[i] = inf;
    
    dp[SELECTED_NODE] = 0;
    
    HeapNode x = createHeapNode(0, SELECTED_NODE);
    addHeap(pq, x);

    while (pq -> size) {
        HeapNode x = getMin(pq);
        int cost = x.cost, vertex = x.vertex;
        if (cost == dp[vertex]) {
            Node *current = G -> adj[vertex].first;
            while (current != NULL) {
                if (dp[current -> vertex] > cost + current -> cost) {
                    dp[current -> vertex] = cost + current -> cost;
                    HeapNode x = createHeapNode(dp[current -> vertex], current -> vertex);
                    addHeap(pq, x);
                }
                current = current -> next;
            }
        }
        popHeap(pq);
    }
    
    if (OK) {
    printf("Minimum cost path from node %d to the other nodes:\n", SELECTED_NODE);
    for (int i = 1; i <= G -> nrVer; ++i)
        if (dp[i] == inf) printf("from %d -> %d you don't have a path\n", SELECTED_NODE, i);
            else printf("from %d -> %d minimum cost is %d\n",SELECTED_NODE, i, dp[i]);
    }

    freeHeap(pq);
}

void dfs(Graph *x, Stack *st, int node, int var) {
    x -> check[node] = true;

    Node *current = x ->adj[node].first;

    while (current != NULL) {
        if (! (x -> check[current -> vertex]) ) dfs(x, st, current -> vertex,  var);
        current = current -> next;
    }

    if (var == 1) pushStack(st, node);
}

void memsetcheck(Graph *x,  Graph *y) {
    for (int i = 1; i <= x -> nrVer; ++i)
        x -> check[i] = y -> check[i] = false;
}

int components(Graph *x, Graph *y) {
    int components = 0;
    Stack *st = createStack(x -> nrVer + 1);

    memsetcheck(x, y);

    for (int i = 1; i <= x -> nrVer; ++i)
        if (!(x -> check[i])) dfs(x, st, i, 1);

    for (int i = y -> nrVer; i >= 1; --i) 
        if (!(y -> check[st -> array[i]])) {
            dfs(y, st,  st -> array[i], 2); 
            components++;
        }
        
    freeStack(st);

    return components;
}

//==========TESTING===========
//============================

void assert(int line, bool b) {
  if (b) return;
  printf("The test on line %d fails.\n", line);
  exit(1);
}

void testMax() {
    assert(__LINE__, max(1,5) == 5);
    assert(__LINE__, max(5,5) == 5);
    assert(__LINE__, max(13,1924) == 1924);
    assert(__LINE__, max(-1, 3) == 3);
    assert(__LINE__, max(7,20) == 20);
}

void testCharToInt() {
    assert(__LINE__, charToInt("123",-1).x == 123 && charToInt("123", -1).pos == 3);
    assert(__LINE__, charToInt("2423432", -1).x == 2423432 && charToInt("2423432", -1).pos == 7);
    assert(__LINE__, charToInt("z321", 0).x == 321 && charToInt("321", 0).pos == 3);
    assert(__LINE__, charToInt("sad11",2).x == 11 && charToInt("sad11", 2).pos == 5);
}

void testCreateHeap() {
    Heap *x = createHeap(5);
   
    for (int i = 0; i < 5; ++i)
        assert(__LINE__, x -> array[i].cost == 0 && x -> array[i].vertex == 0);

    assert(__LINE__, x -> size == 0);
    
    freeHeap(x);
}

void testCreateHeapNode() {
    assert(__LINE__, createHeapNode(1,1).vertex == 1 && createHeapNode(1,1).cost == 1);
    assert(__LINE__, createHeapNode(1,5).vertex == 5 && createHeapNode(1,5).cost == 1);
    assert(__LINE__, createHeapNode(3,4).vertex == 4 && createHeapNode(3,4).cost == 3);
}

void testSwapHeap() {
    Heap *x = createHeap(5);
    HeapNode node = createHeapNode(2, 5), node2 = createHeapNode (3, 4);

    x -> array[1] = node, x -> array[3] = node2;
    swapHeap(x, 1, 3);

    assert(__LINE__, x -> array[1].cost == node2.cost && x -> array[3].cost == node.cost && x -> array[1].vertex == node2.vertex && x -> array[3].vertex == node.vertex);

    freeHeap(x);
}

void testAddHeap() {
    Heap *x = createHeap(4);
    HeapNode node = createHeapNode(1,3), node2 = createHeapNode(5,1), node3 = createHeapNode(10,2);
    
    addHeap(x, node2);
    assert(__LINE__, x -> array[1].cost == 5 && x -> array[1].vertex == 1);
    
    addHeap(x, node3);
    assert(__LINE__, x -> array[2].cost == 10 && x -> array[2].vertex == 2);

    addHeap(x, node);
    assert(__LINE__, x -> array[1].cost == 1 && x -> array[1].vertex == 3);

    freeHeap(x);
}

void testPopHeap() {
    Heap *x = createHeap(4);
    HeapNode node = createHeapNode(7,3), node2 = createHeapNode(1,5), node3 = createHeapNode(3,2);
    addHeap(x, node), addHeap(x, node2), addHeap(x, node3);

    popHeap(x);
    assert(__LINE__, x -> array[1].cost == 3 && x -> array[1].vertex == 2);
    
    popHeap(x);
    assert(__LINE__, x -> array[1].cost == 7 && x -> array[1].vertex == 3);

    freeHeap(x);
}

void testGetMin() {
    Heap *x = createHeap(5);
    HeapNode node = createHeapNode(1,3), node2 = createHeapNode(5,1), node3 = createHeapNode(0,7);
    
    addHeap(x, node2), addHeap(x, node);
    assert(__LINE__, getMin(x).cost == 1 && getMin(x).vertex == 3);
    
    addHeap(x, node3);
    assert(__LINE__, getMin(x).cost == 0 && getMin(x).vertex == 7);

    freeHeap(x);
}

void testCreateQueue() {
    Queue *q = createQueue(5);

    assert(__LINE__, q -> size == 0);
    assert(__LINE__, q -> capacity == 5);
    assert(__LINE__, q -> front == 0);
    assert(__LINE__, q -> back == 4);
    
    freeQueue(q);
}

void testPushBack() {
    Queue *q = createQueue(5);
    q = pushBack(q, 3);
    q = pushBack(q, 5);
    q = pushBack(q, 1);
    assert(__LINE__, q -> v[q -> back] == 1);
    freeQueue(q);
}

void testFullQueue() {
    Queue *q = createQueue(2);
    q = pushBack(q, 3);
    q = pushBack(q, 5);
    assert(__LINE__, fullQueue(q));
    freeQueue(q);
}

void testPopFront() {
    Queue *q = createQueue(3);
    q = pushBack(q, 3), q = pushBack(q, 5), q = pushBack(q, 6);
    q = popFront(q);
    assert(__LINE__, q -> v[q -> front] == 5);
    freeQueue(q);
}

void testEmptyQueue() {
    Queue *q = createQueue(3);
    q = pushBack(q, 3); 
    assert(__LINE__, !emptyQueue(q));
    q = popFront(q);
    assert(__LINE__, emptyQueue(q));
    freeQueue(q);
}

void testGetFront() {
    Queue *q = createQueue(3);
    q = pushBack(q, 1), q = pushBack(q, 2), q = pushBack(q, 3);
    assert(__LINE__, getFront(q) == 1);
    q = popFront(q);
    assert(__LINE__, getFront(q) == 2);
    q = popFront(q);
    assert(__LINE__, getFront(q) == 3);
    freeQueue(q);
}

void testCreateGraph() {
    Graph *x = createGraph(5);
    assert(__LINE__, x -> nrVer == 5);

    for (int i = 1; i <= x -> nrVer; ++i)
        assert(__LINE__, x -> adj[i].first == NULL);

    freeGraph(x);
}

void testReadAndMakeGraph() {
    Graph *x = readAndMakeGraph("1-2,1-3,2-4,2-5", 0);
    assert(__LINE__, x -> nrVer == 5);

    Node *current = x -> adj[1].first;
    assert(__LINE__, current -> vertex == 3);
    current = current -> next;
    assert(__LINE__, current -> vertex == 2); 

    current = x -> adj[2].first;
    assert(__LINE__, current -> vertex == 5);
    current = current -> next;
    assert(__LINE__, current -> vertex == 4);

    freeGraph(x);
}

bool check(int n, int x[n + 1], int y[]) {
    if (n != 5) return false;
    for (int i = 1; i <= n; ++i)
        if (x[i] != y[i]) return false;

    return true;
}

void testDijkstra() {
    Graph *x = readAndMakeGraph("1-2,1-4/2,4-3/4,2-3/2,4-5/3,3-5/6",0);
    int dp[x -> nrVer + 1];
    dijkstra(x, dp, false);

    assert(__LINE__, check(x -> nrVer, dp, (int[]){0,0,1,3,2,5}));
    
    freeGraph(x);
}

void testMemsetCheck() {
    Graph *x = readAndMakeGraph("1-2", 0);
    Graph *y = readAndMakeGraph("1-2", 1);
    memsetcheck(x, y);
    
    for (int i = 1; i <= x -> nrVer; ++i)
        assert(__LINE__, !(x -> check[i]) && !(y -> check[i]));

   freeGraph(x);
   freeGraph(y);
}

void testComponents() {
    Graph *x = readAndMakeGraph("2-1,1-3,3-2,1-4,4-5", 0);
    Graph *y = readAndMakeGraph("2-1,1-3,3-2,1-4,4-5", 1);
    assert(__LINE__, components(x, y) == 3);
    freeGraph(x), freeGraph(y);
   
    x = readAndMakeGraph("1-2,2-1", 0);
    y = readAndMakeGraph("1-2,2-1", 1);
    assert(__LINE__, components(x, y) == 1);
    freeGraph(x), freeGraph(y);

    x = readAndMakeGraph("1-2,2-6,6-7,7-6,3-1,3-4,2-3,4-5,5-4,6-5,5-8,8-7", 0);
    y = readAndMakeGraph("1-2,2-6,6-7,7-6,3-1,3-4,2-3,4-5,5-4,6-5,5-8,8-7", 1);
    assert(__LINE__, components(x, y) == 2);
    freeGraph(x), freeGraph(y);
}

void testCreateStack() {
    Stack *st = createStack(5);
    assert(__LINE__, st -> capacity == 5);
    assert(__LINE__, st -> size == 0);

    for (int i = 0; i < st -> capacity; ++i) 
        assert(__LINE__, st -> array[i] == 0);

    freeStack(st);
}

void testPushStack() {
    Stack *st = createStack(5);
    pushStack(st, 1), pushStack(st, 2);
    assert(__LINE__, st -> size == 2);
    assert(__LINE__, st -> array[1] == 1 && st -> array[2] == 2);
    freeStack(st);
}

void testPopStack() {
    Stack *st = createStack(5);
    pushStack(st, 1), pushStack(st, 2);
    popStack(st);
    assert(__LINE__, st -> size == 1);
    popStack(st);
    assert (__LINE__, st -> size == 0);
    freeStack(st);
}

void testIsStackEmpty() {
    Stack *st = createStack(5);
    pushStack(st, 1), popStack(st);
    assert(__LINE__, isStackEmpty(st) == true);
    pushStack(st, 1);
    assert(__LINE__, isStackEmpty(st) == false);
    freeStack(st);
}

void testIsStackFull() {
    Stack *st = createStack(3);
    pushStack(st, 1), pushStack(st, 1);
    assert(__LINE__, isStackFull(st) == true);
    popStack(st);
    assert(__LINE__, isStackFull(st) == false);
    freeStack(st);
}

void testStack() {
    testCreateStack();
    testPushStack();
    testPopStack();
    testIsStackEmpty();
    testIsStackFull();
}

void testHeap() {
   testCreateHeap();
   testCreateHeapNode();
   testSwapHeap();
   testAddHeap(); // tests both addHeap and heapUp
   testPopHeap(); // tests both popHeap and heapDown
   testGetMin();
}

void testQueue() {
   testCreateQueue();
   testPushBack();
   testFullQueue();
   testPopFront();
   testEmptyQueue();
   testGetFront();   
}

void testGraph() {
    testCreateGraph();
    testReadAndMakeGraph(); //also tests addEdge and addEdges
    testDijkstra();
    testMemsetCheck();
    testComponents(); //also tests dfs
}

void test() {
    testMax();
    testCharToInt();
    testHeap();
    testQueue();
    testStack();
    testGraph();
    printf("All tests passed.\n");
}

int main (int n, char *args[n]) {
   setbuf(stdout, NULL);
   if (n == 1) {
        test();
        return 0;
    }

    if (n == 2) {
       Graph *graph = readAndMakeGraph(args[1], 0);

       printf("Number of nodes: %d\n", graph -> nrVer);
       
       shortestPath(graph, graph -> nrVer);

       int dp[graph -> nrVer + 1]; // minimum cost from selected vertex to the other ones
       dijkstra(graph, dp, true);
       
       Graph *graph2 = readAndMakeGraph(args[1], 1); // transpose of the graph

       printf("The number of strongly connected components: %d\n", components(graph, graph2));

       freeGraph(graph);
       freeGraph(graph2);
       
       return 0;
    } 
    
     invalidInput();

    return 0;
}
