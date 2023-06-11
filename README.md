# Intersection Optimizer
## Summary
Monte Carlo simulation to determine the expected times and various sub properties of the intersection L crossing problem.  
  
The intersection L crossing problem is as follows. You are coming to a 4 way intersection with stop lights in both directions.
You must cross both directions (ex. moving from Nw to the SE corner of the intersection). You can either make both crossings
at the intersection, or you can choose to cross the street immediately (ex there exists a pedestrian crossing or a green light where
you currently are), and make only one crossing at the intersection. In what situation is either method optimal?  

optimizer.py can be run to perform a Monte Carlo simulation of the problem. math_plotter.py plots a theoretical derivation of the
expected value in both cases. Intersection Optimization.pdf is the report summarizing the problem and our results

## optimizer.py
To run optimizer.py, ensure you have python 3, scipy, numpy, and matplotlib installed on your machine. The code can be run by
calling `python optimizer.py` in the root directory. The code will perform the monte carlo simulation for the specified stoplight timings
and plot results obtained over the entire range of simulation. 

The code sweeps over a range of "on" times for both street lights and performs a full Monte Carlo simulation for each permutation.
to change this range, modify the `minT`, `maxT`, and `increment` variables. There are several other problem properties that can be
specified which are listed below:  
- `walkSpeed`: how fast you are assumed to walk in m/s. Affects the time it takes to cross the intersection
- `crossWidth`: width of the intersection in m. Affects the time it takes to cross the intersection
- `deadTime`: time between when one crossing light switches from green to red and the other light switches from red to green.
set by default to 6 seconds to simulate Toronto stop lights. Can also be used to simulate a time in which there is not enough
time available to cross near the end of a light's on period.
- `numSamples`: number of samples used for the Monte Carlo simulation

## math_plotter.py
math_plotter.py ran be run, and works similarly to optimizer.py. Instead of running a monte carlo simulation for each set of stop
light parameters, it determines the expected value using the theoretical derivation obtained in the paper.

## Results
The following plot is obtained for mean times of both strategies, assuming an initial travel direction of north to south on the west 
side of the road. strategy A refers to the method of making both crossings at the intersection, while strategy B is the method of
making the west to east crossing immediately, and then crossing north to south at the intersection. Results are obtained from 10,000
Monte Carlo Samples  
  
![Matplotlib output of a countour plot and heatmap showing the expected value times for both strategies mentioned above](/Results.png)




