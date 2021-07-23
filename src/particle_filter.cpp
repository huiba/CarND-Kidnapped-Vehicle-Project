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
#include <string>
#include <vector>
#include <map>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  normal_distribution<double> x_dist{x, std[0]};
  normal_distribution<double> y_dist{y, std[1]};
  normal_distribution<double> theta_dist{theta, std[2]};

  for (unsigned int i =0; i < num_particles; ++i) {
    Particle par;
    par.id = i;
    par.x  = x_dist(rnd_gen);
    par.y  = y_dist(rnd_gen);
    par.theta  = theta_dist(rnd_gen);
    par.weight = 1;
    particles.push_back(par);
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
  //std::default_random_engine rnd_gen;
  normal_distribution<double> x_dist{0, std_pos[0]};
  normal_distribution<double> y_dist{0, std_pos[1]};
  normal_distribution<double> theta_dist{0, std_pos[2]};
  for (auto& p : particles) {
    if (std::fabs(yaw_rate) >= 0.00001) {
      p.x += (velocity / yaw_rate) * (std::sin(p.theta + yaw_rate*delta_t) - std::sin(p.theta)) ;
      p.y += (velocity / yaw_rate) * (std::cos(p.theta) - std::cos(p.theta + yaw_rate*delta_t)) ;
      p.theta += yaw_rate * delta_t ;
    }
    else {
      p.x += velocity * delta_t * cos(p.theta);
      p.y += velocity * delta_t * sin(p.theta);
    }

    p.x += x_dist(rnd_gen);
    p.y += y_dist(rnd_gen);
    p.theta += theta_dist(rnd_gen);
  }
}

//void ParticleFilter::prediction(double delta_t, double std_pos[],
//                                double velocity, double yaw_rate) {
//    /**
//     * TODO: Add measurements to each particle and add random Gaussian noise.
//     * NOTE: When adding noise you may find std::normal_distribution
//     *   and std::default_random_engine useful.
//     *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
//     *  http://www.cplusplus.com/reference/random/default_random_engine/
//     */
//    //Normal distributions for sensor noise
//    normal_distribution<double> disX(0, std_pos[0]);
//    normal_distribution<double> disY(0, std_pos[1]);
//    normal_distribution<double> angle_theta(0, std_pos[2]);
//
//    for (int i = 0; i < num_particles; i++) {
//        if (fabs(yaw_rate) >= 0.00001) {
//            particles[i].x +=
//                    (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
//            particles[i].y +=
//                    (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
//            particles[i].theta += yaw_rate * delta_t;
//        } else {
//            particles[i].x += velocity * delta_t * cos(particles[i].theta);
//            particles[i].y += velocity * delta_t * sin(particles[i].theta);
//        }
//        // Add noise
//        particles[i].x += disX(rnd_gen);
//        particles[i].y += disY(rnd_gen);
//        particles[i].theta += angle_theta(rnd_gen);
//
//    }
//
//}

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
  

  double sum_w = 0;
  weights.clear();
  for (auto& p : particles) {
    // 1. transform to map coordinate system
    double x = p.x;
    double y = p.y;
    double theta = p.theta;
    // calculate basis vector in map coordinate system
    double u1 = std::cos(theta);
    double u2 = std::sin(theta);
    double v1 = -u2;
    double v2 = u1;
    vector<LandmarkObs> transformed_observations;
    for (auto& ob : observations) {
      auto t_x = ob.x * u1 + ob.y * v1 + x; 
      auto t_y = ob.x * u2 + ob.y * v2 + y; 
      transformed_observations.push_back(LandmarkObs{ob.id, t_x, t_y});
    }
    // 2. Associate observations with landmarks
    // key: obs id, value: landmark
    std::map<int, LandmarkObs> associated_landmarks;
    // filter out invalid landmarks
    std::vector<LandmarkObs> valid_landmarks; 
    for (auto& lm : map_landmarks.landmark_list) {
      if (std::fabs(x  - lm.x_f) <= sensor_range && std::fabs(y - lm.y_f) <= sensor_range) {
        valid_landmarks.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
      }
    }
    // find associations
    for (auto& ob : transformed_observations) {
      double min_dist = std::numeric_limits<double>::max();
      int mid = -1;
      double mx;
      double my;
      for (auto& vlm : valid_landmarks) {
        auto d = dist(ob.x, ob.y, vlm.x, vlm.y);
        if (d < min_dist) {
          min_dist = d;
          mid = vlm.id;
          ob.id = vlm.id;
          mx = vlm.x;
          my = vlm.y;
        }
      }
      associated_landmarks[mid] = LandmarkObs{mid, mx, my};
    }
    // 3. calculate weights
    double w = 1.0;
    double gaussian_normalizer = 2 * M_PI * std_landmark[0] * std_landmark[1];
    for (auto& ob : transformed_observations) {
      auto x_normalized = (ob.x - associated_landmarks[ob.id].x) / std_landmark[0];
      auto y_normalized = (ob.y - associated_landmarks[ob.id].y) / std_landmark[1];
      double w_t = std::exp(-(x_normalized * x_normalized + y_normalized * y_normalized) / 2.0) / gaussian_normalizer;
      if (w_t == 0) {
        w *= 0.00001; 
      } else {
        w *= w_t;
      }
    }
    sum_w += w;
    p.weight = w;
    weights.push_back(p.weight);
  }
}

void ParticleFilter::resample() {
  //std::default_random_engine rnd_gen;
  std::discrete_distribution<> d{weights.begin(), weights.end()};
  std::vector<Particle> new_particles;
  for (unsigned int i=0; i < particles.size(); ++i) {
    int re_id = d(rnd_gen);
    auto p = particles[re_id];
    new_particles.push_back(p);
  }
  particles = new_particles;
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