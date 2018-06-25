# CarND-Kidnapped-Vehicle: Particle Filters

Term2-Project3: Kidnapped Vehicle with Particle Filters

The goal of this project is to make the kidnapped vehicle recognize itself in
the given 2D world. This project impelemntes a Particle Filter for localization
of car in this 2D world. The particle filters is applied to the Markov
Localization algorithm.

The inputs to the program constitutes of a map with 2D coordiantes & properties
of the fixed objects and noisy data from multiple sensors with reference to the 
car's coordinate space. A simple cyclic model is assumed to model the motion of 
the vehicle. Using Markov Localization with Particle filters, the car is able 
to identify its location on this map.

This project works with the Udacity's simulator, which reads data from our
implementation and calculates the error in our predictions. If the errors are
under the considerable limit, the implementation is termed as successful.

![error](https://raw.githubusercontent.com/nitheeshkl/CarND-Kidnapped-Vehicle-Project/master/error.png)
