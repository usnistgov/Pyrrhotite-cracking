NIST Data Publication:
2D Finite Element Model for Quasistatic Crack Formation  in Multiphase
Materials Driven by Internal Phases Expansion
Version 1.0.0
DOI: https://doi.org/10.18434/mds2-4002

If you used this code for your work, please cite the dataset and the associated NIST 
Internal Report (NIST IR 8591):

* O.H.A Dehwah, E.J. Garboczi, N.S. Martys and S.S. Watson, (2025). 2D Finite Element
Based Model of Quasistatic Crack Formation Driven by Expansion of Internal Phases
in Multi-phase Materials, Version 1.0.0. National Institute of Standards and Technology,
https://doi.org/10.18434/mds2-4002 (Accessed: 12/17/2025)

*  O.H.A Dehwah, E.J. Garboczi, N.S. Martys and S.S. Watson, (2025). Quasistatic Crack Formation
in Multiphase Materials Driven by Internal Phases Expansion Mechanisms: Application to 
Cement-Based-Materials. (National Institute of Standards and Technology, Gaithersburg, MD), 
NIST Internal Report (IR) NIST IR 8591. https://doi.org/10.6028/NIST.IR.8591

-------
Authors
-------

  Osamah H.A Dehwah
    National Institute of Standards and Technology
    Materials and Structural Systems Division
    
  Edward J. Garboczi
    National Institute of Standards and Technology
    Applied Chemicals and Materials Division
    
  Nicos S. Martys
    National Institute of Standards and Technology
    Materials and Structural Systems Division
    
  Stephanie S. Watson
    National Institute of Standards and Technology
    Materials and Structural Systems Division


Contact:
  Osamah Dehwah
    osamah.dehwah@nist.gov

-----------
Description
-----------

This dataset contains the source code for modeling crack formation caused by the
expansion internal phases, driven by deterioration mechanisms such as
pyrrhotite oxidation or Alkali Silica Reaction (ASR) in concrete. It also
includes simulation results and the  necessary input files for the various 
Fortran components. A comprehensive  description is available in the NIST IR: 
Quasistatic Crack Formation in Multiphase  Materials Driven by Internal
Phases Expansion Mechanisms: Application to  Cement-Based-Materials, DOI:
https://doi.org/10.6028/NIST.IR.8591.


--------------
Data Use Notes
--------------

This data is publicly available according to the NIST statements of
copyright, fair use and licensing; see
https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software

You may cite the use of this data as follows:
Dehwah, Osamah H.A, Garboczi, Edward J., Martys, Nicos S., Watson, Stephanie S.
(2025), 2D Finite Element Based Model of Quasistatic Crack Formation Driven by
Driven by Expansion of Internal Phases  in Multi-phase Materials, Version 1.0.0,
National Institute of Standards and Technology,
https://doi.org/10.18434/mds2-4002 (Accessed: [12-18-2026])

-------------
Data Overview
-------------

Files included in this publication:

 Crack Formation_main_codes.zip

    Main Codes
    This is the main codes:  To run simulation: thermal2d-critical-stress-
      cracking.f To continue/restart simulation: thermal2d-stress-
      cracking_restart.f To obtain elastic properties: thermal2d-moduli.f
    Format: Fortran

 10_cracks_restart.zip

    Restart results example
    This file contains results for the restart simulation for 10 cracks per
      iteration.
      
  Simulation_results_10&100cracks.zip

    Results_Example
    This file contains simulation results for 10 and 100 cracks per iteration.
      Images can be combined to create a GIF image.

  10_cracks_Thermal.zip

    Elastic Properties
    This file shows some results of the elastic properties obtained using
      "thermal2d-moduli.f". Also, it contains script.py and output_Thermal.py,
      which are used to run the thermal2d-moduli.f and to perform analysis
      (post-processing) of the results.
    Format: Fortran

* Detailed file descriptions can be found in the provided reference, NIST Internal Report:
Quasistatic Crack Formation in Multiphase Materials Driven by Internal Phases Expansion Mechanisms: 
Application to Cement-Based-Materials.

---------------
Version History
---------------

1.0.0 (this version)
  initial release

