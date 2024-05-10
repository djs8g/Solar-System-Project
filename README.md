## Solar-System-Project
This is my final project for ASTR 5470 at the University of Virginia. The plan is to create a Solar System Simulator that will emulate the Earth's Solar System with the ability to modify the model and study the results.

### ASTR 5470 – Computational Astrophysics
Final Project Proposal by Darryl J. Schnellenberger

The basis of the project is to create a 2D model and eventually a 3D model. This all depends on if I am able to stabilize the 2D model and get at least one 3D animation tool to work on my laptop.

The solar system model is being developed with real gravitational physics in the model. Various functions are being added to the finalprojectfunctions.py file, which includes integrators (Runge-Kutta 2 and Leap Frog), N-Body models, various gravitational mechanics models (accelerations, forces, vectors, etc.), positional updates, randomization function (to randomize the starting positions), and various animation and movement functions.

Once all of this is working, the idea is to be able to modify the solar objects' configurations, in particular to play around with the mass of Jupiter and watch what happens to the orbits of the other planets and the Sun (requiring the Sun's position to be modifiable). While the model is still not stable and producing even a 2D model, primarily due to issues with animation attempts, including matplotlib, Turtle, and PyGames.

Some functions have also been written for the potential add-on portion of capturing an interstellar object that enters the solar system, but these remain unused and untested at this time. The capture problem will include determining whether an object's trajectory will intersect with the sphere of influence, also known as the Hill sphere, of the larger body, in this case the Sun. The Hill sphere is the region around a celestial body where its gravity dominates over that of other bodies.

An assumption in this project is that the solar system developed as-is into its current configuration. This does not try to form the solar system, although I have great interest in trying to develop something that will do such a thing.

This model defines classes for the star and planetary objects, with various methods to update changing orbital parameters. Constants are also defined in the finalprojectconstants.py file, many of which are included from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/.

There are also plans to create a user interface, rather than just making changes directly in the code, but that portion is delayed indefinitely.

This project currently utilizes Python code and scripts. The accompanying Jupyter Notebook is scrapped for now given that a stable model is not available. The code is currently being run independently from a Unix/Linux Kernel.  Outputs will include 2-D (and 3-D) animations of a working solar system model, plus graphs and plots to show calculations from planned upgrades to the initial model. Eventually there are plans to write some code in C and compare the calculations.

As of now, you will only need to have Anaconda or some other working python interface to run the available functions/methods. Access to libraries such as matplotlib and numpy are also needed. It would also be useful to have Turtle installed and working, since one of the model attempts is currently running with Turtle. I have been unable to get PyGames working, nor have I been able to get other graphics interfaces to work, but if I do I will update this appropriately for any that may be used. Jupyter Notebook or Jupyter Lab should also be available for when the model is working and I am able to lay out a useful set of steps with these environments.
the submission of this proposal but have not been successful yet.

