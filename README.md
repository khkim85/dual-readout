# dual-readout
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter.

## How-to
### Compile
After fetching the repository, do

    source setenv-cc7-gcc8.sh
    mkdir build
    cd build
    cmake3 ..
    make -j4

Copy excutable file in bin to each directory.

    cp bin/DRsim DRsim/

### Running GEANT4
#### 1. GEANT4 standalone particle gun
In build/DRsim,

    ./DRsim <run_macro> <filenumber> <filename>
    ./DRsim run_ele.mac 0 eltest

generates, `<filename>_<filenumber>.root`

#### 2. Using HepMC input
This requires the ROOT file generated from `Gen`. Assuming the name of the file `<filename>_<filenumber>.root`,

    ./DRsim run_hepmc.mac <filenumber> <filename>

### Reconstruction
This requires the ROOT file generated from `DRsim`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/Reco,

    ./Reco <filenumber> <filename>

### generation macros
generate box/ele_0.root with GEANT4 simulation

    source rungun.sh run_ele $(Process) box/ele

generate with condor jobs

    condor_submit runel.co

### Analysis
This requires the ROOT file generated from `Reco`. Assuming the name of the file `<filename>_<filenumber>.root`, in build/analysis,

    ./<your_analysis_program> <filenumber> <filename>

Merge files from condor jobs with name `box/ele_#.root`,# from 0 to 49.
to eltestdrsim.root, eltestreco.root

    ./merge ../box/ele_%.root 50 ../eltest

Process image with `eltestdrsim.root` and `eltestreco.root` to `eltest.root`

    ./process ../eltest ../eltest.root 0

eltest.root have "event" tree contains processed data by each event.

for python contents can be inspected by accessing tree

    import ROOT
    import numpy as np
    
    infile=ROOT.TFile("eltest.root","read")
    event=infile.Get("event")
    event.Print()
    
    event.GetEntry(0)
    print(event.E_DRcorr,np.array(image_ecore_s).shape,np.array(fiber_ecor_s).shape))
    
images have 168 * 168 shape, but imaging algorithm is still in progress.

### Precaution
Since GEANT4 takes very large amount of time per an event, P8ptcgun, DRsim and Reco are assumed to run a few events only per ROOT file. The executables can be run on parallel using `torque` or `condor`, and can be merged before analysis step using `hadd` from ROOT.

