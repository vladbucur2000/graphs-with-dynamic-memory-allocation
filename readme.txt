The input must contains the edges and weights: a-b or a-b/x (in the first case the cost of a-b is 1 by default).

My program counts the number of nodes, calculates the shortest path from 
SELECTED_NODE_BFS (which is defined as 1, but can be changed) to the other nodes (using BFS algorithm). 

It also calculates the minimum cost path from SELECTED_NODE (which is defined as 
1, but can be changed) to the other nodes (DIJKSTRA algorithm). To get the best complexity O(number of edges * log(number of nodes)) I used a minHeap to keep the cost of the edges.

I also used DFS to find out how many strongly connected components are in the directed graph. 

