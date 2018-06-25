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

#define EPS 0.00001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    if (is_initialized) return;

    num_particles = 100; // initialize arbitrary number of particles

    // normal distributions
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // generate particles
    for (int i=0; i<num_particles; i++) {
        Particle particle; // new instance
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0; // initial weight
        particles.push_back(particle); // add to collection
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // normal distributions
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    // calculate new state for all particles
    for (int i=0; i<num_particles; i++) {
        const double theta = particles[i].theta;
        // predict state based on our motion model
        if (fabs(yaw_rate) < EPS) {
            particles[i].x += velocity * delta_t * cos(theta);
            particles[i].y += velocity * delta_t * sin(theta);
        } else {
            particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate*delta_t) - sin(theta));
            particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + yaw_rate*delta_t));
            particles[i].theta += yaw_rate * delta_t;
        }
        // include noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    const unsigned int n_observations = observations.size();
    const unsigned int n_predictions = predicted.size();

    for (unsigned int i = 0; i < n_observations; i++) {
        double min_dist = numeric_limits<double>::max();
        int map_id = -1;
        // find predictions based on closest distance
        for (unsigned int j = 0; j < n_predictions; j++) {
            double x_dist = observations[i].x - predicted[j].x;
            double y_dist = observations[i].y - predicted[j].y;
            double dist = x_dist*x_dist + y_dist*y_dist;

            if (dist < min_dist) {
                min_dist = dist;
                map_id = predicted[j].id;
            }
        }
        observations[i].id = map_id;
    }
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

    for (int i = 0; i < num_particles; i++) {
        const double x = particles[i].x;
        const double y = particles[i].y;
        const double theta = particles[i].theta;
        // find the list of landmarks within the sensors range
        const double sensor_range_2 = sensor_range * sensor_range; // using squares for comparision below, just to avoid sqrt
        vector<LandmarkObs> landmarks_in_range;
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            float landmark_x = map_landmarks.landmark_list[j].x_f;
            float landmark_y = map_landmarks.landmark_list[j].y_f;
            int id = map_landmarks.landmark_list[j].id_i;
            const double dx = x - landmark_x;
            const double dy = y - landmark_y;
            if (dx*dx + dy*dy < sensor_range_2) {
                landmarks_in_range.push_back(LandmarkObs{id,  landmark_x,  landmark_y});
            }
        }

        // transform observation coordinates
        vector<LandmarkObs> mapped_observations;
        for (unsigned int j = 0; j < observations.size(); j++) {
            // transform based on transformation matrix logic
            const double map_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
            const double map_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
            mapped_observations.push_back(LandmarkObs{observations[j].id, map_x, map_y});
        }

        // associate mapped observations to landmarks in range
        dataAssociation(landmarks_in_range, mapped_observations);

        // reset and calculate weight
        particles[i].weight = 1.0;
        for (unsigned int j = 0; j < mapped_observations.size(); j++) {
            const double x_obs = mapped_observations[j].x;
            const double y_obs = mapped_observations[j].y;
            int landmark_id = mapped_observations[j].id;

            double landmark_x, landmark_y;
            for (unsigned int k = 0; k < landmarks_in_range.size(); k++) {
                if (landmarks_in_range[k].id == landmark_id) {
                    landmark_x = landmarks_in_range[k].x;
                    landmark_y = landmarks_in_range[k].y;
                    break;
                }
            }
            const double dx = x_obs - landmark_x;
            const double dy = y_obs - landmark_y;
            // weight based on Gaussian equation
            const double weight = ((1/(2*M_PI*std_landmark[0]*std_landmark[1]))
                                    * exp(-(dx*dx/(2*std_landmark[0]*std_landmark[1])
                                            + (dy*dy/(2*std_landmark[0]*std_landmark[1])))
                                        )
                                  );
            if (weight == 0) { // avoid zero
                particles[i].weight *= EPS;
            } else {
                particles[i].weight *= weight;
            }
        }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // get weights
    vector<double> weights;
    double max_weight = numeric_limits<double>::min();
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
        if (particles[i].weight > max_weight) {
            max_weight = particles[i].weight;
        }
    }

    // distributions
    uniform_real_distribution<double> dist_weights(0.0, max_weight); // for weight
    uniform_int_distribution<int> dist_index(0, num_particles-1); // for selecting particle

    // resample based on Sebastian's logic of distributing weights in a circle
    // get random index
    int index = dist_index(gen);
    double beta = 0.0; // set to min
    // the circle/wheel
    vector<Particle> resampled_particles;
    for (int i = 0; i < num_particles; i++) {
        beta += dist_weights(gen) * 2.0; // consider a range of current + 2*w
        while (beta > weights[index]) {
            beta -= weights[index]; // keep moving forward
            index = (index + 1) % num_particles; // circle around
        }
        resampled_particles.push_back(particles[index]);
    }
    particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
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
