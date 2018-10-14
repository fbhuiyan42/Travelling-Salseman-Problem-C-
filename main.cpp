#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include<algorithm>
#include<vector>
#include<fstream>
#include <iomanip>
#include <limits.h>

using namespace std;
typedef vector<int> Tour;
ofstream fout("1005081_ratio.csv");
ofstream fout2("1005081_compare.csv");
//ofstream f("input.txt");
struct Point
{
    double x,y;
};

double euclideanDistance(struct Point a,struct Point b)
{
    return sqrt( ((a.x - b.x)*(a.x - b.x)) + ((a.y - b.y)*(a.y - b.y)));
}

double score(const Point *point,const vector<int>& tour,int V)
{
    double result = euclideanDistance(point[tour[0]],point[tour[V-1]]);
    for(int i = 1; i < V; ++i)
    {
       result = result + euclideanDistance(point[tour[i - 1]], point[tour[i]]);
    }
    //cout<<endl;
    return result;
}

vector<int> get_default_tour(const Point *point,int V)
{
    vector<int> tour(V);
    for(int i = 0; i < V; ++i)
    {
      tour[i] = i;
    }
    return tour;
}

int brute_force(const Point *point,int V)
{
    Tour tour(get_default_tour(point,V));
    Tour best_tour(tour);
    double best_score = score(point,tour,V);
    /*for (int i = 0; i < V; i++)
    {
        cout<<tour[i]<<" ";
    }
    cout <<"  score = " << best_score << endl;*/

    while (next_permutation(tour.begin()+1, tour.end()))
    {
        double s = score(point,tour,V);
        if (s < best_score)
        {
            best_score = s;
            best_tour = tour;
        }
        /*for (int i = 0; i < V; i++)
        {
            cout<<tour[i]<<" ";
        }
        cout<< "  score = " << s << endl;*/
    }
    cout << "TSP : ";
    for (int i = 0; i < V; i++)
    {
        cout<<best_tour[i]<<" ";
    }
    cout<<best_tour[0];
    cout<< "\nCost = " << best_score << endl;
    fout <<best_score<<",";
    return best_score;
}

void Task_1()
{
    ifstream fin ("1005081_input.txt");
    cout<<"Exact Exponential Method \n";
    int caseno =1;
    ofstream file_out("1005081_Task1.csv");
    file_out << "x, y" << endl;
    int V;
    //for(int V=2;caseno<11;V=V+1)
     while (fin >> V )
    {
        cout<<"\n----------CASE "<<caseno<<" ------------\n\n";
        caseno++;
        Point *point;
        point=(Point*)malloc(sizeof(Point)*V);
        clock_t start, stop;
        //srand(1);
        start = clock();
        for(int i=0;i<V;i++)
        {
            fin>>point[i].x;
            fin>>point[i].y;
            //point[i].x = rand()%10 ;
            //point[i].y = rand()%10 ;
        }
        int cost=brute_force(point,V);
        stop = clock();
        printf("\nfor input size %d execution time is %lf seconds\n", V,(((double)stop)-start)/CLOCKS_PER_SEC);
        file_out << V << ", " << (((double)stop)-start)/CLOCKS_PER_SEC << endl;
    }
    fout <<endl;
    fin.close();
    return ;
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int minKey(int key[], bool mstSet[],int V)
{
   // Initialize min value
   int min = INT_MAX, min_index;

   for (int v = 0; v < V; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;

   return min_index;
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation
int primMST(vector<vector<int> > &graph,int V)
{
     int parent[V]; // Array to store constructed MST
     int key[V];   // Key values used to pick minimum weight edge in cut
     bool mstSet[V];  // To represent set of vertices not yet included in MST

     // Initialize all keys as INFINITE
     for (int i = 0; i < V; i++)
        key[i] = INT_MAX, mstSet[i] = false;

     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST

     // The MST will have V vertices
     for (int count = 0; count < V-1; count++)
     {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet,V);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < V; v++)

           // graph[u][v] is non zero only for adjacent vertices of m
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if graph[u][v] is smaller than key[v]
          if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v])
             parent[v]  = u, key[v] = graph[u][v];
     }
     int cost=0;
     for (int i = 1; i < V; i++)
        cost=cost+graph[i][parent[i]];
     return cost;
}

int createMST(const Point *point,int v,int upto)
{
    int V=v-upto;
    int graph[V][V];
    int row=0;
    int i,j;
    for(i=upto;i<v;i++)
    {
        int column=0;
        for(j=upto;j<v;j++)
        {
            if(i!=j)
            {
                graph[row][column]=euclideanDistance(point[i],point[j]);
                column++;
            }
        }
        row++;
    }
    vector<vector<int> > vec;
    for (int i = 0; i < V; i++)
    {        //creating row
        vec.push_back(vector<int>());
    }
    for (int n = 0; n < V; n++)
    {        //creating columns for the rows
        for (int m = 0; m < V; m++)
        {
            vec[m].push_back(0);
        }
    }
    for (int m = 0; m < V; m++)
    {        //storing data
        for (int n = 0; n < V; n++)
        {
            vec[m][n] = graph[m][n];
        }
    }
    int cost=primMST(vec,V);
    return cost;
}

int lowerbound(const Point *point,int *p,int V,int i)
{
    double MSTcost,costA,costB,total;
    MSTcost=createMST(point,V,i);
    costA=99999;
    costB=99999;
    if(i>1)
    {
        for (int j = i; j < V; j++)
        {
            if(euclideanDistance(point[p[i-2]], point[p[j]])< costA)
            {
                costA= euclideanDistance(point[p[i-2]], point[p[j]]);
            }
        }
        for (int j = i; j < V; j++)
        {
            if(euclideanDistance(point[p[i-1]], point[p[j]])< costB)
            {
                costB= euclideanDistance(point[p[i-1]], point[p[j]]);
            }
        }
    }
    total= MSTcost+costA+costB+euclideanDistance(point[p[i-1]], point[p[i-2]]) ;
    return total;
}

void swaps(int *x,int *y)
{
    int temp;
    temp=*x;
    *x=*y;
    *y=temp;
}

int findTour(const Point *point,int *p,int V,int i,double w,int *TreeSoFar,double *bestsofar,int *pruned)
{
    int j;
    double t;
    double lower;
    if (i == V-1)
    {
        t = w + euclideanDistance(point[p[0]], point[p[V-1]]);
        if (t < *bestsofar)
        {
            memcpy(TreeSoFar, p, sizeof(int)*V);
            *bestsofar = t;
        }
        return *bestsofar;
    }
    else
    {
        if(i>1)
        {
            lower=lowerbound(point,p,V,i);
            if (lower >= *bestsofar)
            {
                pruned++;
                return *bestsofar;
            }
        }
    }
    for (j = i+1; j < V; j++)
    {
        swaps(&p[i+1], &p[j]);
        t = findTour(point, p, V, i+1,w + euclideanDistance(point[p[i]], point[p[i+1]]),TreeSoFar,bestsofar,pruned);
        swaps(&p[i+1], &p[j]);
    }
    return *bestsofar;
}

int branch_bound(const Point *point,int V)
{
  int *p, *TreeSoFar, i;
  double t, bestsofar = 999;
  int pruned=0;

  p = (int*)malloc(sizeof(int)*V);
  for (i = 0; i < V; i++)
  {
    p[i] = i;
  }
  TreeSoFar = (int*)malloc(sizeof(int)*V);
  findTour(point, p, V, 0,0.0, TreeSoFar, &bestsofar,&pruned);
  cout<<"TSP : ";
  for(int i=0;i<V;i++) cout<<TreeSoFar[i]<<" ";
  cout<<TreeSoFar[0];
  cout<<endl;
  cout<<"pruned "<<pruned<<" times";
  printf("\ncost = %lf \n",bestsofar);
  fout2 <<bestsofar<<",";
  return bestsofar;
}


void Task_2()
{
    ifstream fin ("1005081_input.txt");
    cout<<"Branch and Bounding Method \n";
    int caseno =1;
    ofstream file_out("1005081_Task2.csv");
    file_out << "x, y" << endl;
    int V;
    //for(int V=2;caseno<11;V=V+1)
    while (fin >> V )
    {
        cout<<"\n----------CASE "<<caseno<<" ------------\n\n";
        caseno++;
        Point *point;
        point=(Point*)malloc(sizeof(Point)*V);
        clock_t start, stop;
        srand(0);
        start = clock();
        for(int i=0;i<V;i++)
        {
            fin>>point[i].x;
            fin>>point[i].y;
            //point[i].x = rand()%10 ;
            //point[i].y = rand()%10 ;
        }
        int cost=branch_bound(point,V);
        stop = clock();
        printf("\nfor input size %d execution time is %lf seconds\n", V,(((double)stop)-start)/CLOCKS_PER_SEC);
        file_out << V << ", " << (((double)stop)-start)/CLOCKS_PER_SEC << endl;
    }
    fout2<<endl;
    fin.close();
    return;
}


struct Tree
{
    int V, E;
    struct Edge* edge;
};

struct Edge
{
    int src, dest, weight;
};

struct Tree* createTree(int V, int E)
{
    struct Tree* tree = (struct Tree*) malloc( sizeof(struct Tree) );
    tree->V = V;
    tree->E = E;

    tree->edge = (struct Edge*) malloc( tree->E * sizeof( struct Edge ) );

    return tree;
}

// a structure to represent a connected, undirected and weighted graph
struct Graph
{
    // V-> Number of vertices, E-> Number of edges
    int V, E;

    // graph is represented as an array of edges. Since the graph is undirected, the edge from src to dest
    // is also edge from dest to src. Both are counted as 1 edge here.
    struct Edge* edge;
};

// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
    struct Graph* graph =(struct Graph*) malloc( sizeof(struct Graph) );
    graph->V = V;
    graph->E = E;
    graph->edge = (struct Edge*) malloc( graph->E * sizeof( struct Edge ) );
    return graph;
}

// A structure to represent a subset for union-find
struct subset
{
    int parent;
    int rank;
};

// A utility function to find set of an element i (uses path compression technique)
int find(struct subset subsets[], int i)
{
    // find root and make root as parent of i (path compression)
    if (subsets[i].parent != i)
    {
        subsets[i].parent = find(subsets, subsets[i].parent);
    }
    return subsets[i].parent;
}

// A function that does union of two sets of x and y (uses union by rank)
void Union(struct subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Attach smaller rank tree under root of high rank tree
    // (Union by Rank)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;

    // If ranks are same, then make one as root and increment its rank by one
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Compare two edges according to their weights.
// Used in qsort() for sorting an array of edges
int myComp(const void* a, const void* b)
{
    struct Edge* a1 = (struct Edge*)a;
    struct Edge* b1 = (struct Edge*)b;
    return a1->weight >= b1->weight;
}

bool notVisited(int node,int index,int result[])
{
    for(int i=0;i<index;i++)
    {
        if(node==result[i]) return false;
    }
    return true;
}

double preOrder(struct Tree* tree, int root,int index,int *Result,const Point *point,Tour  tour,double &temp)
{
    tour[index]=root;
    Result[index++]= root;
    printf("%d ",Result[index-1]);
    for(int i=0;i< tree->E;i++)
    {
        if(tree->edge[i].src==root && notVisited(tree->edge[i].dest,index-1,Result)) preOrder(tree,tree->edge[i].dest,index,Result,point,tour,temp);
        else if(tree->edge[i].dest==root && notVisited(tree->edge[i].src,index-1,Result)) preOrder(tree,tree->edge[i].src,index,Result,point,tour,temp);
    }
    double cost=0;
    if(index>0)
    {
        cost= score(point,tour,tree->V);
        if(temp>cost)
        {
            cost=temp;
        }
        temp=cost;
    }
    return cost;
}

// The main function to construct MST using Kruskal's algorithm
double KruskalMST(struct Graph* graph,const Point *point)
{
    int V = graph->V;
    struct Edge result[V];  // Tnis will store the resultant MST
    int e = 0;  // An index variable, used for result[]
    int i = 0;  // An index variable, used for sorted edges

    // Step 1:  Sort all the edges in non-decreasing order of their weight
    // If we are not allowed to change the given graph, we can create a copy of
    // array of edges

    qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);

    // Allocate memory for creating V ssubsets
    struct subset *subsets =
        (struct subset*) malloc( V * sizeof(struct subset) );

    // Create V subsets with single elements
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }
    // Number of edges to be taken is equal to V-1
    while (e < V - 1)
    {
        // Step 2: Pick the smallest edge. And increment the index
        // for next iteration
        struct Edge next_edge = graph->edge[i++];
        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);
        // If including this edge does't cause cycle, include it
        // in result and increment the index of result for next edge
        if (x != y)
        {
            result[e++] = next_edge;
            Union(subsets, x, y);
        }
        // Else discard the next_edge
    }
    // print the contents of result[] to display the built MST
    //printf("Following are the edges in the constructed MST\n");
    int k=0;
    struct Tree* tree = createTree(graph->V, e);
    struct Tree* sortedtree = createTree(graph->V, e);
    for (i = 0; i < e; ++i)
    {
        tree->edge[k].src = result[i].src;
        tree->edge[k].dest = result[i].dest;
        tree->edge[k].weight = result[i].weight;
        k++;
    }
    cout<<"MST :\n";
    for (int i = 0; i < e; ++i)
    {
        printf("%d -- %d == %d\n", tree->edge[i].src, tree->edge[i].dest,tree->edge[i].weight);
    }
    int Result[tree->V];
    for(int i=0;i<tree->V;i++) Result[i]=0;
    printf("TSP : ");
    Tour tour(get_default_tour(point,tree->V));
    double temp=0;
    double cost=preOrder(tree,tree->edge[0].src,0,Result,point,tour,temp);
    printf("%d ",tree->edge[0].src);
    cout<< "\nCost = " << cost << endl;
    fout <<cost<<",";
    fout2 <<cost<<",";
    return cost;
}

void Task_3()
{
    ifstream fin ("1005081_input.txt");
    //ifstream fin ("input.txt");
    cout<<"Greedy 2-Approximation Algorithm \n";
    int E,k,V;
    int caseno =1;
    ofstream file_out("1005081_Task3.csv");
    file_out << "x, y" << endl;
    //for(int V=2;caseno<11;V=V+1)
    while (fin >> V )
    {
        cout<<V;
        cout<<"\n----------CASE "<<caseno<<" ------------\n\n";
        caseno++;
        E=V*(V-1);
        k=0;
        Point *point;
        point=(Point*)malloc(sizeof(Point)*V);
        clock_t start, stop;
        srand(1);
        start = clock();
        //f<< V <<endl;
        for(int i=0;i<V;i++)
        {
            fin>>point[i].x;
            fin>>point[i].y;
            //point[i].x = rand()%10 ;
            //point[i].y = rand()%10 ;
            //f<< point[i].x <<" "<<point[i].y<< endl;
        }
        //f<<endl;
        struct Graph* graph = createGraph(V, E);
        //cout<<"Constructed Graph :\n";
        for(int i=0;i<V;i++)
        {
            for(int j=0;j<V;j++)
            {
                if(i!=j)
                {
                    graph->edge[k].src = i;
                    graph->edge[k].dest = j;
                    graph->edge[k].weight = euclideanDistance(point[i], point[j]);
                    //cout<<graph->edge[k].src<<"--"<<graph->edge[k].dest<<" == "<<graph->edge[k].weight<<"\n";
                    k++;
                }
            }
        }
        double cost=KruskalMST(graph,point);
        stop = clock();
        printf("\nfor input size %d execution time is %lf seconds\n", V,(((double)stop)-start)/CLOCKS_PER_SEC);
        file_out << V << ", " << (((double)stop)-start)/CLOCKS_PER_SEC << endl;
    }
    fin.close();
    return;
}

int main()
{
    Task_1();
    Task_2();
    Task_3();
    return 0;
}
