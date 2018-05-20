# N-body Gravity Simulator
This was a numerical methods project for grad school in Computational Physics. My python script models the orbits of celestial bodies based on gravitational forces. The initial conditions of the planets are from NASA. 

## Methods
To calculate the trajectories of the planets, I used the Leapfrog integration method. Leapfrog, a 2nd order method, is useful for mechanics because of its time reversability. If you integrate n steps, and then integrate backward by n steps, you arrive at your initial position. Leapfrog is also symplectic in that it conserves energy. Runge-Kutta, which is a 4th order method, is also popular but it does not conserve energy and the system will drift over time. 

I probed the solar system to see if it is maximally packed. In simulation, you can add an Earth-sized planet between Mars and Jupiter. Considering the mass of asteroid belt, this makes sense. 

## Results
We expect the total energy and total angular momentum to be conserved. The roundoff error of the change in energy and momentum cancels out and nicely oscillates over time. 

![Energy Figure Not Found!](/img/NBodyOrbit40Energy.png)

![Momentum Figure Not Found!](/img/NBodyOrbit40Momentum.png)

Over 40 years, note the magnitude of the change in energy and momentum is quite low. 

![Orbit Figure Not Found!](/img/NBodyOrbit40.png)

I included Pluto in the simulation, even though it's not a planet. 

## Additional Comments
It's interesting exploring 3-body systems with the gravity simulator. Perturbing unstable equilibria produces chaotic orbits. 


