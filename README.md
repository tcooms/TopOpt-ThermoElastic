# Thermoelastic topology optimization code for 2D problems in Matlab
This repository contains Matlab files for thermoelastic topology optimization considering steady-state and transient heat conduction.


The file top_tml_shc.m solves a topology optimization problem for a simply-supported beam subjected to a mechanical point load on the top and a thermal boundary condition at the bottom edge. More information about the theory and implementation can be found in:

Ooms, T., Vantyghem, G., Thienpont, T., Van Coile, R., De Corte W. Compliance-based topology optimization of structural components subjected to thermo-mechanical loading. Struct Multidisc Optim (2023). https://doi.org/10.1007/s00158-023-03563-3

The included results can be replicated with the developed MATLAB code in the supplementary material and executed by a command of the following form:

top_tml_shc(L,h,t,z,Vf,rmin,pE,pk,pb)

where the variables refer to the corresponding parameters discussed in the paper. 

For example, the optimization with the default parameters can be solved by entering the following command:

top_tml_shc(1200,400,10,10,0.4,3,3,3,3)

A modification of the original top_tml_shc.m code is made for design-independent thermal loads such as a uniform temperature difference. The code top_tml_uniform.m is added to this repository and can be used in a similar manner as described above.

The MATLAB files mmasub.m and subsolv.m for using the MMA algorithm are freely available on http://www.smoptit.se/ under the GNU General Public License (GPLv3). One should reference them correctly in line 3 of the MATLAB code to carry out the optimization procedure.
