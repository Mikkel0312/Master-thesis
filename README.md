# Master-thesis
Source code for my Master's thesis: "Using route length approximations to solve the Covering Tour Problem".

In my Master's thesis I solved the Covering Tour Problem using 4 different heuristic solution methods:
1. A Cost-allocation method
2. An "Online"/"Offline" regression method
3. A mathematical program with different TSP route length approximations as the objective function
4. A minimum-distance-maximum-covering(MDMC) method.

The code for the first three heuristics is in "ConstructionCTP.py".
The code for the fourth heuristic is in "MDMC.py"

I also developed an improvement heuristic, which improves an initial solution. The code for this algorithm is in "ImproveCTP.py".

To compare the solutions obtained by the heuristics against the optimal solution, I solved the CTP to optimality using CPLEX and Baldacci's formulation of the CTP. The python code for this is in "CTP_TSPLIB.py" and the OPL model is in "CoveringTourProblem.mod". 

The code for the computational tests is in "Tables.py".

In the "CTP_main.py" is a script, that generates a random set of points and solves the CTP on those points. Here you specify which of the four algorithms to use, and also specify if the improvement heuristic should be used on not. It prints the objective value and plots the solution.

I used a Python interface for Concorde(http://www.math.uwaterloo.ca/tsp/concorde.html), which loads CPLEX. In order to run the code, you must have CPLEX installed on your machine.


