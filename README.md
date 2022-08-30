# coupled_ice_hydrology
 
Summarizing different attempts to couple hydrology and ice flow

This document provides a rough outline to each script I have developed and what assumptions it makes. It will describe failure points in addition to initial outputs. 
1.	SSA_new_coeff.m
This is the SSA code adjusted with the new sliding law term and non-dimensionalization coefficients. I set N equal to100 kPa, 100 total time steps over 10 years, and a grid of 200 points in space. One thing I immediately noticed was that with such a “fine” time scale, the solver either “stalls” or returns a solution with “potential inaccuracies”. This is a problem inherent to the system I believe, for if I use A. Robel’s code as is and keep the number of time steps as 100 and the end time as 10 years, it returns the same errors. 
 <img width="256" alt="image" src="https://user-images.githubusercontent.com/25989736/187533367-17117dca-f50b-4da2-a8fb-74b72ea68ee2.png">
<img width="256" alt="image" src="https://user-images.githubusercontent.com/25989736/187533391-4cef0b89-95ed-4819-b4b8-b0082fae60e8.png">

One thing I did to evaluate the errors qualitatively was run it with larger time steps and see if it returns the same results. They do seem consistent, however it always seems like the first few steps are ineffective, but then the later steps are solved normally. Also, running A. Robel’s code on its own gives similar issues without changing the time step. 
Results (all with 100 steps):
    <img width="629" alt="image" src="https://user-images.githubusercontent.com/25989736/187533577-5b4fc0c0-038a-4d53-8a4c-d88ebc6bfa7f.png">

2.	hydro_realistic_h.m
This script takes in the initial steady state output for glacier thickness from SSA_new_coeff.m and uses that to determine a value for psi. It needs to interpolate between whatever grid that the ice model uses and its own grid (a uniform grid which accounts for stretching of the grounding line). Using a set timestep of 0.05 (nondimensionalized) and 1000 total steps, this simulation runs for 38 or so years. We also see that it reaches a steady state throughout the system quite quickly, except for the channel exit area which grows over time (this is expected). 
 <img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533630-ee45f0e5-d7c7-4bd7-90b7-a78ec8a3d9f6.png">

3.	hydro_realistic_hu.m
On top of the same outputs taken in by hydro_realistic_h, this model also takes in the steady state velocity. Using the exact same setup, these are the results:
 <img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533643-6b0d2538-d740-42e9-b5b3-7487039617c6.png">

Worryingly, it seems almost the exact same as the previous results. 

We know there is a velocity component that is advecting the surface area:
 <img width="190" alt="image" src="https://user-images.githubusercontent.com/25989736/187533660-6d4cd548-0710-4d09-9043-909886e7db7c.png">


Looking more closely at the differences between the two:
 <img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533687-7fcc7de0-3b9b-46f7-84c7-5b2a69133d91.png">

We see there is negligible difference throughout the channel with the exception of near/at the channel terminus, where the model accounting for velocity yields lower values. We see that the channel exit area does get smaller compared to that of the model without velocity, but it doesn’t reach a steady equilibrium. A steady state does exist – setting time derivatives to zero gives:
 <img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533701-5623464b-fe14-405a-8c4c-5c2576ceb2d1.png">

We see that the flow and effective pressure in the non-steady models very rapidly reach a steady state, but then the channel area takes a while. 
I think one issue with this is that xg = xg_old, so that term is always zero
4.	one_way_hydro_from_ice.m
This model is only one directionally coupled as well, except this time it also resolves the ice flow model independently on each time step. It first calculates a steady state ice flow, then on the first time step, recalculates ice flow (because we need two values for grounding line for the dxg/dt term in the hydro equations). Then using this new time step’s ice, we solve for hydrology. Interestingly enough, this seems to reach a steady state channel area quicker, and it is lower than what the steady state is determined to be in the previous example. The other values seem about the same. 
 
<img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533854-796ef146-2545-4253-a74e-a623d8a8654f.png">

5.	two_way_loose.m 
This model is the same process as above, except it feeds in the effective pressure from the hydrology equations back into the ice flow equations. I also have a setting that can toggle whether or not to use steady state hydrology. When I do not use steady state hydrology (I leave the time derivative in there, the solver fails to determine a solution on the next time step’s ice equations (exceeds max number of iterations). The same issue happens with steady state steps (with the same time scale). It could be that the initial solve without hydrology is just too far from the actual solution. These solutions were attempted on 1000 steps over 40 years, like previously.

HOWEVER, if we do time steps closer to that of the ice model (100 steps over 1000 years), we see that the model does indeed generate a solution (this is using steady state hydrology). It seems that with an active channel below, the grounding line recedes quite a bit, while thickness accumulates faster. Furthermore, it seems the glacier is overall slower. I’m a bit iffy with how the model spins up since it doesn’t use hydrology those first two steps. 
 
<img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533906-80bc76c1-61cb-4baa-ad28-0aa2657b14f6.png">

6.	combined.m
This groups all the equations together to solve at once. It takes a while, and for smaller time steps, the solver stalls like the ice sheet model. It solves for an initial steady state before calculating transient evolution like the ice sheet model. We can see that the initial condition is vastly different than before. This is 100 time steps over 1 year. 

 <img width="393" alt="image" src="https://user-images.githubusercontent.com/25989736/187533940-fd42db58-fa46-4659-b6e7-f00e9e79055a.png">

Seems pretty much the same as just 10 steps:
 
<img width="398" alt="image" src="https://user-images.githubusercontent.com/25989736/187533962-8ec6d375-efd5-40b2-8e30-fab8b9cf6344.png">

For some reason, the hydrology equations here can now handle larger time steps.  Doing 100 time steps over 1000 years yields this, which seems to be consistent with the earlier results. 
 
 <img width="468" alt="image" src="https://user-images.githubusercontent.com/25989736/187533984-1a08c6d1-4b1d-4cff-816e-dd3f7659cc88.png">


7.	combined_steady_hydro.m
Now I want to investigate if just calculating steady state hydrology every time step yields the same results. It does yield pretty much the same results (need to quantify) but this makes sense. Since th0 is a lot smaller than t0, the dt for hydro is a lot larger than the dt for ice. As a result, this makes the time derivative terms close to zero even when we’re calculating transience, making it the same as the steady state. 
