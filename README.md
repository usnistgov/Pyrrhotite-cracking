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

  Osamah H.A Dehwah,
    National Institute of Standards and Technology,
    Materials and Structural Systems Division
    
  Edward J. Garboczi,
    National Institute of Standards and Technology,
    Applied Chemicals and Materials Division
    
  Nicos S. Martys,
    National Institute of Standards and Technology,
    Materials and Structural Systems Division
    
  Stephanie S. Watson,
    National Institute of Standards and Technology,
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

 Crack Formation_main_codes

    Crack formation main codes: This is the main codes:  
    To run simulation: thermal2d-critical-stress-cracking.f 
    To continue/restart simulation: thermal2d-stress-cracking_restart.f 
    To obtain elastic properties: thermal2d-moduli.f
    Format: Fortran

 Simulation Results Examples: This folder contains the results of the codes. Three zip folders:

 Simulation_results_10&100cracks.zip

    Simulation Results Example: This folder contains simulation results for 10 and 100 cracks per iteration.
      Images can be combined to create a GIF image.
      
 10_cracks_restart.zip

    Restart results example: This folder contains results for the restart simulation for 10 cracks per iteration. '##.gif' files: Visual representations of crack patterns. 'fort.10##' files: Data files containing the microstructure associated with each gif. 'micro.in': The primary input file. 'phasemod_values': Output file recording the modulus for each phase type. Further technical details regarding these formats are available in Sections 4.2 (Input Data Files) and 4.3 (FEM Program manual: B.	To restart the code) of the Internal Report.
      
  10_cracks_Thermal.zip

    Elastic Properties: This folder contains results for elastic properties obtained using thermal2d-moduli.f. While this example uses five distinct microstructure inputs (fort files), the process is customizable to any number of inputs. The included script.py and output_Thermal.py automate the execution of thermal2d-moduli.f and perform post-processing. 'script.py': Organizes fort.10## files into individual directories and executes the Fortran code. output_Thermal.py: Performs post-processing and generates two identical output files (output_csv and output_txt) for user convenience. Results are reported in GPa, though the output units will consistently match the units used for the input phases in the main code. For further details, refer to the NIST Internal Report, Section 4.3 (FEM Program Manual: C. Physical/Mechanical Properties).
   

* Detailed file descriptions can be found in the provided reference, NIST Internal Report:
Quasistatic Crack Formation in Multiphase Materials Driven by Internal Phases Expansion Mechanisms: Application to Cement-Based-Materials.

---------------
Sources of Variability 
---------------
-A very small tolerance on the norm of the residual is specified to ensure that the final solution is effectively independent of the tolerance itself. Since the appropriate tolerance value is sensitive to the system size, careful selection is required. In this code, a typical tolerance on the order of nx × ny × 10⁻¹² =5.47×10^(-7)  has been used, where nx and ny denote the number of elements in the x and y directions, respectively.

-The code represents the displacement field using a piecewise linear approximation, which would make the normal strain and stress continues at interfaces, but they are not forced to be. 

* More details on sources of variability can be found in the associated NIST Internal Report, Section 3.5.

---------------
Version History
---------------

1.0.0 (this version)
  initial release
  
---------------
Disclaimer
---------------
This data/work was created by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States.  This data/work may be subject to foreign copyright.

The data/work is provided by NIST as a public service and is expressly provided “AS IS.” NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST does not warrant or make any representations regarding the use of the data or the results thereof, including but not limited to the correctness, accuracy, reliability or usefulness of the data. NIST SHALL NOT BE LIABLE AND YOU HEREBY RELEASE NIST FROM LIABILITY FOR ANY INDIRECT, CONSEQUENTIAL, SPECIAL, OR INCIDENTAL DAMAGES (INCLUDING DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF BUSINESS INFORMATION, AND THE LIKE), WHETHER ARISING IN TORT, CONTRACT, OR OTHERWISE, ARISING FROM OR RELATING TO THE DATA (OR THE USE OF OR INABILITY TO USE THIS DATA), EVEN IF NIST HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

To the extent that NIST may hold copyright in countries other than the United States, you are hereby granted the non-exclusive irrevocable and unconditional right to print, publish, prepare derivative works and distribute the NIST data, in any medium, or authorize others to do so on your behalf, on a royalty-free basis throughout the world.

You may improve, modify, and create derivative works of the data or any portion of the data, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the data and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the data:  Data citation recommendations are provided at https://www.nist.gov/open/license.

Permission to use this data is contingent upon your acceptance of the terms of this agreement and upon your providing appropriate acknowledgments of NIST’s creation of the data/work.

---------------
Commercial Product Disclaimer
---------------
Certain equipment, instruments, software, or materials are identified in this data/code in order to specify the experimental procedure adequately.  Such identification is not intended to imply recommendation or endorsement of any product or service by NIST, nor is it intended to imply that the materials or equipment identified are necessarily the best available for the purpose.




