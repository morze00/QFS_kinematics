# QFS kinematics generator
Event generator for (p,2p) and (p,pn) quasi-free scattering reactions in direct and inverse kinematics.
Original code was written in FORTRAN by Leonid Chulkov.
The code can generate kinematics of the following QFS reactions:
<pre>
A(p,2p)B
A(p,pn)B
A(n,np)B
</pre>
where A and B=(A-1) are the initial and final-state nucleus.

In inverse kinematics the program simulates quasi-elastic scattering
of a target proton "i" (being at rest) on a cluster/nucleon "a" bound in the projectile nucleus A

A --> i  ==> B + (i -> a) =>  B + 1 + 2

Pass '--direct-kinematics' flag to simulate in direct kinematics the scattering of
of a beam-like nucleon "i" on a cluster/nucleon "a" bound inside the target nucleus A, which is at rest.
By default, inverse kinematics and (p,2p) reactions are assumed (if no flags are set).

The kinematics of the residual heavy spectator "B", the knocked-out outgoing cluster (2)
and the scattered target nucleon (1) are stored in the output tree (in the lab and CM frames).

***--> Modify info.hh to change beam parameters, masses, internal momentum spread, etc.***

The code implements parameterized NN differential cross section dsigma/dt in the NN center-of-mass frame
using prescription by Cugnon et al.  NIM in Physics Research B 111(1996) 215-220:

https://doi.org/10.1016/0168-583X(95)01384-9

Pass '--cugnon' flag to use the parameterization, otherwise isotropic NN scattering will be used by default.

Use '--ppn' or '--nnp' flag to simulate QFS scattering of a bound neutron (or proton) "a" from proton (or neutron) probe "i".

 

Compilation requires an existing ROOT installation.
Compile the code using 'make' command in the terminal. 

Example usage:

<pre>
make     --> compile the code (requires ROOT installation)
./qfs -h    -->see available options and usage
./qfs --root=output.root --max-events=50000 --cugnon --ppn     -->Simulate 50k (p,pn) events in inverse-kinematics with Cugnon NN parameterization 
./qfs --root=output.root --direct-kinematics      -->Simulate 100k (p,2p) events in direct kinematics asumming isotropic NN cross sections
</pre>
