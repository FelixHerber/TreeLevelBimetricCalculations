(* ::Package:: *)

<< "/opt/ohpc/pub/Mathematica/Applications/Slurm.m" (* Load the Slurm-awareness package. OLD: << Slurm`   *)
SlurmDebugLevel = 0;            (* Set the verbosity level (0 = none) *)
SlurmShowConfiguration[]        (* Show the used nodes and cores      *)
SlurmLaunchKernels[]  
$LoadAddOns={"FeynArts"};
<<FeynCalc`
$ParallelizeFeynCalc=True;
$VeryVerbose=0;

(*References used here:
F. Herber (2025) The Scattering Amplitudes of Bimetric Theory. 
[Master's thesis, Stockholm University]. Diva Portal. Link:
https://urn.kb.se/resolve?urn=urn:nbn:se:su:diva-248139*)

(*The paths to the bimetric theory modefile, where the amplitudes are stored and the code implementing the calculations.*)
hpcBimTheoryFAPath  = "/cfs/home/fehe3416/.Mathematica/Applications/FeynCalc/FeynArts/Models/xActExpansionOfBimetricTheory_FA/xActExpansionOfBimetricTheory_FA";
savePath = "/cfs/home/fehe3416/SavedAmplitudes";
codeFile = "/cfs/home/fehe3416/BimetricAmplitudes.wl";

(* ------------------------------------------------------------ *)
(* Main Code*)
(* ------------------------------------------------------------ *)

ParallelEvaluate[
  Needs["FeynCalc`"];
  ScalarProduct[A,B] = 7;
  Print["Scalar Product should be 7 = ", Contract[FV[A,m]FV[B,m]]];
]

(*Import the functions in BimetricAmplitudes.wl*)
Get[codeFile];
ParallelEvaluate[Get[codeFile]];

(*Check if BimetricAmplitudes.wl has been loaded properly*)
ParallelEvaluate[Names["Global`SplitSingleRoot"]]
ParallelEvaluate[Names["Global`ExpandDummyPolsAndSetMandelstam"]]

(*We only need to consider the processes GG->GG, GG->MM, GM->MM, MM->MM. Can be seen looking at scatteringProcesses.*)
processes = {{1,1,1,1}, {1,1,2,2}, {1,2,2,2}, {2,2,2,2}};
amplitudeStructures = Generate2To2AmpStructures[processes, hpcBimTheoryFAPath];

(*This is the list of distinct bimetric helicity scattering amplitudes we can have since the expansion of 
  bimetric theory has time reversal symmetry and parity symmetry (see Section 7.2 in F. Herber (2025)). 
  See: DerivationOfDistinctProcesses.nb for where these come.*)
scatteringProcesses = {{{1,2},{1,2},{1,2},{1,2}},{{1,2},{1,2},{1,2},{1,-2}},{{1,2},{1,2},{1,-2},{1,-2}},{{1,2},{1,2},{2,2},{2,2}},{{1,2},{1,2},{2,2},{2,1}},{{1,2},{1,2},{2,2},{2,0}},{{1,2},{1,2},{2,2},{2,-1}},{{1,2},{1,2},{2,2},{2,-2}},{{1,2},{1,2},{2,1},{2,1}},{{1,2},{1,2},{2,1},{2,0}},{{1,2},{1,2},{2,1},{2,-1}},{{1,2},{1,2},{2,1},{2,-2}},{{1,2},{1,2},{2,0},{2,0}},{{1,2},{1,2},{2,0},{2,-1}},{{1,2},{1,2},{2,0},{2,-2}},{{1,2},{1,2},{2,-1},{2,-1}},{{1,2},{1,2},{2,-1},{2,-2}},{{1,2},{1,2},{2,-2},{2,-2}},{{1,2},{1,-2},{2,2},{2,2}},{{1,2},{1,-2},{2,2},{2,1}},{{1,2},{1,-2},{2,2},{2,0}},{{1,2},{1,-2},{2,2},{2,-1}},{{1,2},{1,-2},{2,2},{2,-2}},{{1,2},{1,-2},{2,1},{2,1}},{{1,2},{1,-2},{2,1},{2,0}},{{1,2},{1,-2},{2,1},{2,-1}},{{1,2},{1,-2},{2,0},{2,0}},{{1,2},{2,2},{2,2},{2,2}},{{1,2},{2,2},{2,2},{2,1}},{{1,2},{2,2},{2,2},{2,0}},{{1,2},{2,2},{2,2},{2,-1}},{{1,2},{2,2},{2,2},{2,-2}},{{1,2},{2,2},{2,1},{2,1}},{{1,2},{2,2},{2,1},{2,0}},{{1,2},{2,2},{2,1},{2,-1}},{{1,2},{2,2},{2,1},{2,-2}},{{1,2},{2,2},{2,0},{2,0}},{{1,2},{2,2},{2,0},{2,-1}},{{1,2},{2,2},{2,0},{2,-2}},{{1,2},{2,2},{2,-1},{2,-1}},{{1,2},{2,2},{2,-1},{2,-2}},{{1,2},{2,2},{2,-2},{2,-2}},{{1,2},{2,1},{2,2},{2,2}},{{1,2},{2,1},{2,2},{2,1}},{{1,2},{2,1},{2,2},{2,0}},{{1,2},{2,1},{2,2},{2,-1}},{{1,2},{2,1},{2,1},{2,1}},{{1,2},{2,1},{2,1},{2,0}},{{1,2},{2,1},{2,1},{2,-1}},{{1,2},{2,1},{2,0},{2,0}},{{1,2},{2,1},{2,0},{2,-1}},{{1,2},{2,1},{2,-1},{2,-1}},{{1,2},{2,0},{2,2},{2,2}},{{1,2},{2,0},{2,2},{2,1}},{{1,2},{2,0},{2,2},{2,0}},{{1,2},{2,0},{2,1},{2,1}},{{1,2},{2,0},{2,1},{2,0}},{{1,2},{2,0},{2,0},{2,0}},{{1,2},{2,-1},{2,2},{2,2}},{{1,2},{2,-1},{2,2},{2,1}},{{1,2},{2,-1},{2,1},{2,1}},{{1,2},{2,-2},{2,2},{2,2}},{{2,2},{2,2},{2,2},{2,2}},{{2,2},{2,2},{2,2},{2,1}},{{2,2},{2,2},{2,2},{2,0}},{{2,2},{2,2},{2,2},{2,-1}},{{2,2},{2,2},{2,2},{2,-2}},{{2,2},{2,2},{2,1},{2,1}},{{2,2},{2,2},{2,1},{2,0}},{{2,2},{2,2},{2,1},{2,-1}},{{2,2},{2,2},{2,1},{2,-2}},{{2,2},{2,2},{2,0},{2,0}},{{2,2},{2,2},{2,0},{2,-1}},{{2,2},{2,2},{2,0},{2,-2}},{{2,2},{2,2},{2,-1},{2,-1}},{{2,2},{2,2},{2,-1},{2,-2}},{{2,2},{2,2},{2,-2},{2,-2}},{{2,2},{2,1},{2,2},{2,1}},{{2,2},{2,1},{2,2},{2,0}},{{2,2},{2,1},{2,2},{2,-1}},{{2,2},{2,1},{2,1},{2,1}},{{2,2},{2,1},{2,1},{2,0}},{{2,2},{2,1},{2,1},{2,-1}},{{2,2},{2,1},{2,0},{2,0}},{{2,2},{2,1},{2,0},{2,-1}},{{2,2},{2,1},{2,-1},{2,-1}},{{2,2},{2,0},{2,2},{2,0}},{{2,2},{2,0},{2,1},{2,1}},{{2,2},{2,0},{2,1},{2,0}},{{2,2},{2,0},{2,0},{2,0}},{{2,2},{2,-1},{2,1},{2,1}},{{2,1},{2,1},{2,1},{2,1}},{{2,1},{2,1},{2,1},{2,0}},{{2,1},{2,1},{2,1},{2,-1}},{{2,1},{2,1},{2,0},{2,0}},{{2,1},{2,1},{2,0},{2,-1}},{{2,1},{2,1},{2,-1},{2,-1}},{{2,1},{2,0},{2,1},{2,0}},{{2,1},{2,0},{2,0},{2,0}},{{2,0},{2,0},{2,0},{2,0}}};
scatteringProcessesFlatten = Map[Function[x, Flatten[x,1]], scatteringProcesses];

Print["Starting computations of all amplitudes."]
bimetricAmpsList = ParallelMap[Function[x, Calc2To2BimAmpFromStructuresSerial[x, amplitudeStructures, savePath]], scatteringProcessesFlatten];
Print["Computation finished."]
bimetricAmps = {scatteringProcessesFlatten, bimetricAmpsList} // Transpose;
DumpSave[
 FileNameJoin[{savePath, "allBimetricAmplitudes.mx"}],
  bimetricAmps
 ];
