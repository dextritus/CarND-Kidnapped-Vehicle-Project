# CarND-Kidnapped-Vehicle-Project
In this project, a particle filter was implemented for a car driving around a map with 42 landmarks. The motion model of the car was given by the bicycle model. Measurement observations were received at each time step from the Udacity simulator. 

The general steps of the algorithm were as follows:

* based on initial GPS position and heading data, the filter is initialized with ```num_particles``` particles according to a normal distribution with a certain standard deviation for x, y and heading.
* the position of each particle is updated in the ```prediction``` step according to the bicycle model; additional "process" noise is inserted;
* in the ```updateWeights``` method:
	* the observations from the car coordinate systems are transformed in the map coordinate system
	* the closest landmark is found for each observation, based on the minimum distance between the observation and all landmarks
	* the closest landmark for each observation is used to update the weight of the particle (product of probabilities from all observations)
* ```num_particles``` are resampled with replacement according to the weights obtained previously 
* the best particle is used as the best possible car location

A number of 30 particles is used in this project. A higher number of particles can decrease the RMS error of the heading and x and y coordinates, however the computational time increases. Furthermore, a larger number of particles could be required if the map and/or the motion of the vehicle becomes more complex. 