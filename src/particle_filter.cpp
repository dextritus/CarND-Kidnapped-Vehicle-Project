/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// number of particles can be low for this map and vehicle Multivariate_normal_distribution
	num_particles = 30;

	//engine to generate numbers from distribution
	std::default_random_engine rnd_eng;

	std::normal_distribution<double> gauss_x(x, std[0]);
	std::normal_distribution<double> gauss_y(y, std[1]);
	std::normal_distribution<double> gauss_theta(theta, std[2]);

	for (int i = 0; i<num_particles; i++) {
		Particle part;

		part.x = gauss_x(rnd_eng);
		part.y = gauss_y(rnd_eng);
		part.theta = gauss_theta(rnd_eng);

		part.weight = 1.0/num_particles;
		weights.push_back(1.0/num_particles);
	    //don't forget the id 
	    part.id = i;

	    //add the particle to the the vector of particles
	    particles.push_back(part);
	}	

	//done, don't come back here
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine rnd_eng;

	std::normal_distribution<double> gauss_x(0.0, std_pos[0]);
	std::normal_distribution<double> gauss_y(0.0, std_pos[1]);
	std::normal_distribution<double> gauss_theta(0.0, std_pos[2]);
	
	//! if yaw rate is 0, different equation
	for (auto &part : particles) {
		double this_theta = part.theta;

		if (abs(yaw_rate) > 0.0001) {
			part.x += velocity / yaw_rate * (sin(this_theta + yaw_rate*delta_t) - sin(this_theta));
			part.y += velocity / yaw_rate * (cos(this_theta) - cos(this_theta + yaw_rate*delta_t));
			part.theta += yaw_rate*delta_t;
		}

		else {
			part.x +=  velocity * cos(this_theta) * delta_t;
			part.y +=  velocity * sin(this_theta) * delta_t;
			//theta does not change
		}
		part.x += gauss_x(rnd_eng); 
		part.y += gauss_y(rnd_eng); 
		part.theta += gauss_theta(rnd_eng); 
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // did not use this method
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	//   denominator of the multivariate gaussian distribution
	const double multi_den = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);

	//   sum of all particle weights
	double tot_weight = 0.0;
	
    for (int p = 0; p<num_particles; p++) {

    	vector<int> associations;

    	vector<double> sense_x;
    	vector<double> sense_y;

    	vector<LandmarkObs> observations_map;	

		for (const auto obs:observations) {
			LandmarkObs observation_map;

			observation_map.x = particles[p].x + obs.x * cos(particles[p].theta) - obs.y*sin(particles[p].theta);
			observation_map.y = particles[p].y + obs.x * sin(particles[p].theta) + obs.y*cos(particles[p].theta);
			observations_map.push_back(observation_map);
		}

		//reinitialize particle weight
		particles[p].weight = 1.0;

		//find closest landmark to each observation
		for (int i = 0; i<observations_map.size(); i++) {

			double min_dist = sensor_range;
			int association = -1;
			//here, what happens if none of the landmarks are in the sensor range?
			//associate it to 0 
			for (int j = 0; j<map_landmarks.landmark_list.size(); j++) {
				double land_x = map_landmarks.landmark_list[j].x_f;
				double land_y = map_landmarks.landmark_list[j].y_f;

				double dist = sqrt(pow(observations_map[i].x - land_x, 2) + pow(observations_map[i].y - land_y, 2));
				//only look at landmarks at least as close as the sensor range; save the closest one
				if (dist < min_dist) {
					min_dist = dist; 
					association = j;
				}
			}

			if (association != -1) {
				double obs_x = observations_map[i].x;
				double obs_y = observations_map[i].y;
				double mu_x = map_landmarks.landmark_list[association].x_f;
				double mu_y = map_landmarks.landmark_list[association].y_f;
				long double prob = multi_den*exp( -(pow(obs_x - mu_x,2)/pow(std_landmark[0], 2) + pow(obs_y - mu_y,2)/pow(std_landmark[1], 2)));

				particles[p].weight *= prob;
				tot_weight += particles[p].weight;
			}
			associations.push_back(association + 1); //id of the map starts of 1
			sense_x.push_back(observations_map[i].x);
			sense_y.push_back(observations_map[i].y);
		}
		//if weight was never updated, disregard this particle
		if (particles[p].weight == 1) particles[p].weight = 0.0;
		SetAssociations(particles[p], associations, sense_x, sense_y);
		weights[p] = particles[p].weight;
	}
	//normalize weights of all  particles
	for (auto& part:particles) part.weight/=tot_weight;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::default_random_engine rnd_eng;
	//generate a number from 0 to weights.size()-1, with probability given by the weights themselves
	discrete_distribution<int> discrete_dist(weights.begin(), weights.end());

	vector<Particle> new_particles;

	for (int i = 0; i<num_particles; i++) {
		new_particles.push_back(particles[discrete_dist(rnd_eng)]);
	}

	particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
