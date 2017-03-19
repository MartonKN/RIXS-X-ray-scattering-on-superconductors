# RIXS-X-ray-scattering-on-superconductors

This MATLAB code ran on the Harvard cluster, and for the system sizes we considered, often needed a few hundred to thousand nodes. The routine distribute_SRIXS_tasks.m distributes the input parameters between jobs evenly. Then, once the first part of the calculation is done, one needs to add up all the RIXS intensities using 'sum_XAS_contributions.m' and 'sum_IRIXS_contributions.m'.

Paper corresponding to this work:
Phys. Rev. B 94, 165127 (2016)
