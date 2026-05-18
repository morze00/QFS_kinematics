# QFS_kinematics
Event generator for (p,2p) and (p,pn) reactions
Original code was written in FORTRAN by Leonid Chulkov.
Kinematics generator for QFS reactions in direct and inverse kinematics:

 A(p,2p)B, A(p,pn)B, A(n,np)B

Use --direct-kinematics, otherwise inverse kinematics will be used by default.

In inverse kinematics the program simulates quasi-elastic scattering
of a target proton "i" (being at rest) on a cluster/nucleon "a" bound in the projectile nucleus A

 A --> i  ==> B + (i -> a) =>  B + 1 + 2

 The kinematics of the residual heavy spectator "B", the knocked-out outgoing cluster (2)
 and the scattered target nucleon (1) is stored in the output tree.

 Modify info.hh to change beam parameters, masses, internal momentum spread, etc.

  Pass --direct-kinematics flag to simulate in direct kinematics the scattering of
 of a beam-like proton "i" on a cluster/nucleon "a" bound in target nucleus A which is at rest.

 See output tree variables for the kinematics of the particles
 1, 2 and B in the lab and CM frame and other kinematical variables
