Documentation for HyperTRIBE
============================
For source code please visit `GitHub <https://github.com/rosbashlab/HyperTRIBE>`_


About HyperTRIBE
----------------
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. This is an improved version of a previously developed technique called TRIBE (Targets of RNA-binding proteins Identified By Editing). TRIBE expresses a fusion protein consisting of a queried RBP and the catalytic domain from RNA editing enzyme ADAR (ADARcd), which marks target RNA transcripts by converting adenosine to inosine near the RBP binding sites. In spite of its usefulness, TRIBE is constrained by a low editing efficiency and editing-sequence bias from the ADARcd. So, we developed HyperTRIBE by incorporating a previously characterized hyperactive mutation, E488Q, into the ADARcd. This strategy increases the editing efficiency and reduce sequence bias, which dramatically increased sensitivity of this technique without sacrificing specificity. HyperTRIBE provides a more powerful strategy to identify RNA targets of RBPs with an easy experimental and computational protocol at low cost. 

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. Rna (2017)

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016)

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


