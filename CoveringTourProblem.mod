/*********************************************
 * OPL 12.8.0.0 Model
 * Author: mikke
 * Creation Date: 11 Jun 2019 at 14:40:20
 *********************************************/

 
 
int NumV = ...;
int NumW = ...;
int n = NumV+1;


range nodes = 0..n;
range Yrange = 0..NumV;
range Vnodes = 1..NumV;
range Wnodes = 0..NumW;
//range customers = 1..NumCustomers;


// Edges between each pair (i,j), i != j, of "customer nodes" and all nodes
tuple edge {int i; int j;}
setof(edge) C_Edges = { <i,j> | ordered i,j in Vnodes};
setof(edge) G_Edges = { <0,i> | i in Vnodes} union C_Edges;
setof(edge) Edges   = G_Edges union { <i,n> | i in Vnodes};

//int Demands[customers] = ...;
int Distance[nodes][nodes] = ...;
//float Profits[customers] = ...;


int Clusters[Wnodes][Yrange] = ...;


// Decision variables: binary variables for each edge (i,j), i < j
dvar boolean x[Edges];


// Decision variables: Binary variables y(i) equal to one if customer i is visited
dvar boolean y[Yrange];

// Flow variables f(i,j) for each i,j, i != j. If edge (i,j) is in the
// tour, f(i,j) is the demand covered by the vehicle on the path from
// the depot up to and including node j. On the contrary, f(j,i)
// is the empty space left then for the vehicle, i.e. f(j,i)=Q-f(i,j)
dvar float f[nodes][nodes] in 0..n;

// Expression representing length of the tour
dexpr float routeCost= sum ( <i,j> in Edges ) Distance[i][j]*x[<i,j>];

/*****************************************************************************
 *
 * MODEL
 * 
 *****************************************************************************/

// Objective
minimize routeCost;

subject to {
   // Each customer node that is visited is adjacent to 2 edges in the solution
   forall (i in Vnodes)
     sum (<i,j> in Edges) x[<i,j>] + sum (<j,i> in Edges) x[<j,i>] - 2*y[i]== 0;

   // Outflow of source depot node 0
   sum( j in Vnodes ) ( f[0][j] - y[j]) == 0;

   // Inflow to source depot node 0
   sum( j in Vnodes ) ( f[j][0] + y[j] )== n-1;

   // Outflow sink depot node
   sum( j in Vnodes) f[n][j] == n-1;

   // "Flow conservation (inflow=outflow) at each customer node
   forall( i in Vnodes)
     sum( j in nodes: j != i ) (f[j][i]-f[i][j]) - 2*y[i] == 0;

   // If the tour goes from i to j, f(i,j) + f(j,i) equals vehicle's capacity
   forall ( <i,j> in Edges )
     f[i][j] + f[j][i] == (n-1)*x[<i,j>];
     
   // No loops
   forall ( i in nodes )
     f[i][i] == 0;
   
   forall ( j in Wnodes)
     sum(i in Yrange) Clusters[j][i]*y[i] >= 1;          
   
};

{edge} edges = { e | e in Edges : x[e]==1 }; 
{int} servedCusts = { i | i in Vnodes : y[i] == 1 };
