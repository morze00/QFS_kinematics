# QFS kinematics generator
Event generator for (p,2p) and (p,pn) quasi-free scattering reactions in direct and inverse kinematics.
Original code was written in FORTRAN by Leonid Chulkov.
The code can be used to generate kinematics for the following QFS reactions:

A(p,2p)B, A(p,pn)B, A(n,np)B

where A and B=(A-1) are the initial and final-state nucleus.

In inverse kinematics the program simulates quasi-elastic scattering
of a target proton "i" (being at rest) on a cluster/nucleon "a" bound in the projectile nucleus A

A --> i  ==> B + (i -> a) =>  B + 1 + 2

The kinematics of the residual heavy spectator "B", the knocked-out outgoing cluster (2)
and the scattered target nucleon (1) is stored in the output tree.

Modify info.hh to change beam parameters, masses, internal momentum spread, etc.

Pass '--direct-kinematics' flag to simulate in direct kinematics the scattering of
of a beam-like proton "i" on a cluster/nucleon "a" bound inside the target nucleus A, which is at rest.
By default, inverse kinematics and (p,2p) reactions are assumed (if no flag is set).

See output tree variables for the kinematics of the particles
1, 2 and B in the lab and CM frames. 

Compilation requires and existing ROOT installation.
Compile the code using 'make' command in the terminal. 

Example usage:

<pre>
./qfs --root=output.root --max-events=50000 --cugnon --ppn")
</pre>
./qfs --root=output.root --direct-kinematics
