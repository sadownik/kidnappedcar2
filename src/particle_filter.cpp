/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles
  
    std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // TODO: Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; ++i) {
    
    Particle newParticle;
    
    newParticle.id = i;
    newParticle.x = dist_x(gen);
    newParticle.y = dist_y(gen);
    newParticle.theta = dist_theta(gen);
    newParticle.weight = 1.0;
    // Print your samples to the terminal.
    //std::cout << "Sample " << i + 1 << " " << newParticle.x << " " << newParticle.y << " " << newParticle.theta << std::endl;
    
    particles.push_back(newParticle);
    weights.push_back(newParticle.weight);
  }
  
    is_initialized = true; 
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;

  
  for (int i = 0; i < particles.size(); ++i) {
    
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;
    
    //double theta_f = delta_t * yaw_rate + particle_theta;
    //particle_y = particle_y + velocity /(yaw_rate)* (cos(particle_theta) - cos(theta_f));
    //particle_x = particle_x + velocity /(yaw_rate)* (sin(theta_f) - sin(particle_theta));
    //particle_theta = theta_f;
    
	  double pred_x;
	  double pred_y;
	  double pred_theta;
	  //Instead of a hard check of 0, adding a check for very low value of yaw_rate
	  if (fabs(yaw_rate) < 0.0001) {
	    pred_x = particle_x + velocity * cos(particle_theta) * delta_t;
	    pred_y = particle_y + velocity * sin(particle_theta) * delta_t;
	    pred_theta = particle_theta;
	  } else {
	    pred_x = particle_x + (velocity/yaw_rate) * (sin(particle_theta + (yaw_rate * delta_t)) - sin(particle_theta));
	    pred_y = particle_y + (velocity/yaw_rate) * (cos(particle_theta) - cos(particle_theta + (yaw_rate * delta_t)));
	    pred_theta = particle_theta + (yaw_rate * delta_t);
	  }
	  
	  normal_distribution<double> dist_x(pred_x, std_pos[0]);
	  normal_distribution<double> dist_y(pred_y, std_pos[1]);
	  normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
    
    /*
    normal_distribution<double> dist_x(particle_x, std_pos[0]);
    normal_distribution<double> dist_y(particle_y, std_pos[1]);
    normal_distribution<double> dist_theta(particle_theta, std_pos[2]);
    */
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);    
    particles[i].theta = dist_theta(gen);
    
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  
  //std::cout << observations.size()<< std::endl;

  
  for (int i = 0; i < particles.size(); ++i) {
    particles[i].sense_x.erase (particles[i].sense_x.begin(),particles[i].sense_x.end());
    particles[i].sense_y.erase (particles[i].sense_y.begin(),particles[i].sense_y.end());
    for (int p = 0; p < observations.size(); ++p) {
    	float x_obs_map = observations[p].x * cos(particles[i].theta) - observations[p].y * sin(particles[i].theta) + particles[i].x; 
    	float y_obs_map = observations[p].x * sin(particles[i].theta) + observations[p].y  * cos(particles[i].theta) + particles[i].y; 
        particles[i].sense_x.push_back(x_obs_map);
        particles[i].sense_y.push_back(y_obs_map);
      
    }
  }  
  

  /*  
  for (int i = 0; i < map_landmarks.landmark_list.size(); ++i) {
  std::cout << map_landmarks.landmark_list[i].x_f << "   " << map_landmarks.landmark_list[i].y_f << std::endl;
  }
  std::cout << "-----end landmarks------" <<  std::endl;   
  */
  //observations.size()

    for (int t = 0; t < particles.size(); ++t) {
    particles[t].associations.erase (particles[t].associations.begin(),particles[t].associations.end());
    for (int i = 0; i < observations.size(); ++i) {
      float min = sensor_range*sqrt(2);
      int index = 0; // get the index of the landmark list to correlate best with the observation

      for (int p = 0; p < map_landmarks.landmark_list.size(); ++p) {
        float distance = sqrt(pow(particles[t].sense_x[i] - map_landmarks.landmark_list[p].x_f, 2.0) + pow(particles[t].sense_y[i] - map_landmarks.landmark_list[p].y_f, 2.0));
          if (distance<min) {
              min=distance;
              index = p+1;
          }
      }
      particles[t].associations.push_back(index);
    }

  }
  /*
    for (int i = 0; i < observations.size(); ++i) {
  std::cout << particles[0].sense_x[i] << "   " << particles[0].sense_y[i] <<  std::endl;
    
  std::cout <<  map_landmarks.landmark_list[particles[0].associations[i]-1].x_f << "   " << map_landmarks.landmark_list[particles[0].associations[i]-1].y_f <<  std::endl;
   std::cout << "-----------" <<  std::endl; 
    
  }

  
  
  
  
  
  
  for (int i = 0; i <  particles[0].associations.size(); ++i) {
  std::cout << particles[0].associations[i]<< std::endl;
  }
  std::cout << "-----end assos------" <<  std::endl;  
  */
  /**
     NOW THE UPDATE STEP
  */
  /*
  float weight_normalizer = 0.0;
    const double x_denom = 2 * std_landmark[0] * std_landmark[0];
  const double y_denom = 2 * std_landmark[1] * std_landmark[1];
  const double a = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  for (int i = 0; i < particles.size(); ++i) { // loop through particles to calculate each w
    
    particles[i].weight = 1.0;

    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];

    double denominator = (1.0/(sqrt(2.0 * M_PI * sigma_x * sigma_y)));
    
    	for (int t = 0; t < observations.size(); ++t) {
    		float trans_obs_x = particles[i].sense_x[t];
    		float trans_obs_y = particles[i].sense_y[t];
          	float pred_landmark_x = map_landmarks.landmark_list[particles[i].associations[t]-1].x_f;
          	float pred_landmark_y = map_landmarks.landmark_list[particles[i].associations[t]-1].y_f;
          
          
          
          
            float nn_x = map_landmarks.landmark_list[particles[i].associations[t]-1].x_f;
      		float nn_y = map_landmarks.landmark_list[particles[i].associations[t]-1].y_f;
      
      // Calculate multi-variate Gaussian distribution
      		double x_diff = trans_obs_x - nn_x;
      		double y_diff = trans_obs_y - nn_y;
      		double b = ((x_diff * x_diff) / x_denom) + ((y_diff * y_diff) / y_denom);
      		particles[i].weight *= a * exp(-b);

          

          
    		//particles[i].weight *= denominator * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y))));
          //std::cout << number  <<  std::endl; 
          
          
        }
    weight_normalizer += particles[i].weight;
    //std::cout << number  <<  std::endl; 
  
  }
  
*/
  
 float weight_normalizer = 0.0;
    const double x_denom = 2 * std_landmark[0] * std_landmark[0];
  const double y_denom = 2 * std_landmark[1] * std_landmark[1];
  const double a = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  for (int i = 0; i < particles.size(); ++i) { // loop through particles to calculate each w
    
    particles[i].weight = 1.0;

    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

    //double denominator = (1.0/(sqrt(2.0 * M_PI * sigma_x * sigma_y)));
    
    	for (int t = 0; t < observations.size(); ++t) {
    		float obs_x = particles[i].sense_x[t];
    		float obs_y = particles[i].sense_y[t];
          	float mu_x = map_landmarks.landmark_list[particles[i].associations[t]-1].x_f;
          	float mu_y = map_landmarks.landmark_list[particles[i].associations[t]-1].y_f;
      
      // Calculate multi-variate Gaussian distribution
            double gauss_norm;
            gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

            // calculate exponent
            double exponent;
            exponent = (pow(obs_x - mu_x, 2) / (2 * pow(sig_x, 2)))
                         + (pow(obs_y - mu_y, 2) / (2 * pow(sig_y, 2)));

            // calculate weight using normalization terms and exponent
            double weight;
            particles[i].weight *= gauss_norm * exp(-exponent);
          

          
    		//particles[i].weight *= denominator * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y))));
          //std::cout << number  <<  std::endl; 
          
          
        }
    weight_normalizer += particles[i].weight;
    //std::cout << number  <<  std::endl; 
  
  }
    
  
 /*  
for (int i = 0; i < particles.size(); ++i) {
  std::cout << particles[i].weight  <<  std::endl;
  }
  std::cout << "-----end particles weights------" <<  std::endl;  
 
*/

  

  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= weight_normalizer;
    weights[i] = particles[i].weight;
    if(std::isnan(particles[i].weight)){
              particles[i].weight = 0.0;
    }
    if(std::isinf(particles[i].weight)){
              particles[i].weight = 1;
    }                         
  }

  for (int t = 0; t < particles.size(); ++t) {
    //std::cout<<particles[t].weight << std::endl;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
 
 std::vector<Particle> particles_new;
 std::default_random_engine gen;

 std::vector<int> weights;
 
    for(int i=0; i<particles.size(); ++i) {
        weights.push_back(particles[i].weight);
      
    }
  //std::cout<<weights.size() << std::endl;

 std::discrete_distribution<> d(weights.begin(), weights.end());
  
 for (int i = 0; i < particles.size(); ++i) {
  int number = d(gen); 
  std::cout << number << std::endl;
  particles_new.push_back(particles[number]);
  }
  //std::cout << "------end resampled weights-----" <<  std::endl;  
 particles = particles_new;
 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}