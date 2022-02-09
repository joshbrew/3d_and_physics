# 3D and Physics Math Tests - JS
 Some 3D and rendering math in js. 
 
 - Basic physics (with AABB and octrees, but not implemented yet in the collision checks), 
 - graph node trees for translating/rotating meshes 
 - 3D math e.g. basic linear algebra and calculus. Efficient nearest neighbor search.
 - pinhole camera model

Largely untested, more for exercise but eventually will get used by me (and thus cleaned up & demo'd). 

The camera matrix math in this version is based on using 2D matrices which is not quite as efficient but is conceptually easier. Meshes are still 1D arrays e.g. `[p0x,p0y,p0z,p1x,p1y,p1z,p2x,p2y,p2z]`

Bonus jank Double Pendulum simulation with 3 different derivations rendering side by side. They all have... issues... 
This is being made to explore more of the math of complex systems and learn lagrangian mechanics.

The web does not have good options for efficient pluggable 2d and 3d physics models to go on top of, say, a ThreeJS renderer. There are some things to do to leverage the gpu for simulating large physics body collections that I cannot find any JS implementations for that are developer oriented (if not behind a paywall or deeply integated in a larger system). I prefer few to no dependencies and otherwise just keeping my object formats fairly standard and general across libraries so it's easier to mix and match based on my needs. Maybe there will be a standard some day for all of the expected systems and data structures we are constantly redeveloping in each language etc, since right now that is mainly managed by what is actually efficient - which gets really esoteric really fast.
