Documentation for HyperTRIBE
============================
For source code please visit `GitHub <https://github.com/rosbashlab/HyperTRIBE>`_


About HyperTRIBE
----------------
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. This is an improved version of a previously developed technique called TRIBE (Targets of RNA-binding proteins Identified By Editing). HyperTRIBE couples an RBP to the catalytic domain of the Drosophila RNA editing enzyme ADAR and expresses the fusion protein in vivo. As the RBP-ADARcd (catalytic domain) fusion protein lacks the RNA recognition features of ADAR, the specificity of the RBP should determine the editing specificity of the fusion protein.  RBP targets are marked with novel RNA editing events and identified by sequencing RNA. This repository provides the necessary computational pipeline needed to analyze these RNA sequencing libraries.

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.

**For previous version of this software (TRIBE), please visit** `GitHub <https://github.com/rosbashlab/TRIBE>`_

Contents:

.. toctree::
   :maxdepth: 2

   installation
   run
   outputs
   example
   mariadb

Support
-------
Please use Github issues to bring up any errors that occur with software.


