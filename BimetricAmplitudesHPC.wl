(* ::Package:: *)

<< "/opt/ohpc/pub/Mathematica/Applications/Slurm.m" (* Load the Slurm-awareness package.  *)
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

(* ------------------------------------------------------------ *)
(* Paths                                                        *)
(* ------------------------------------------------------------ *)

(* Directory containing this HPC script *)
baseDirectory = DirectoryName[$InputFileName];

(* Project-local paths *)
savePath =
 FileNameJoin[
  {baseDirectory, "SavedAmplitudes"}
 ];

codeFile =
 FileNameJoin[
  {baseDirectory, "BimetricAmplitudes.wl"}
 ];

hpcBimTheoryFAPath =
 FileNameJoin[
  {
   baseDirectory,
   "xActExpansionOfBimetricTheory_FA",
   "xActExpansionOfBimetricTheory_FA"
   }
  ];

(* Check that all required files/directories exist *)
If[!DirectoryQ[savePath],
 Print[
  "ERROR: SavedAmplitudes directory not found at: ",
  savePath
  ];
 Abort[];
 ];

If[!FileExistsQ[codeFile],
 Print[
  "ERROR: BimetricAmplitudes.wl not found at: ",
  codeFile
  ];
 Abort[];
 ];

If[!FileExistsQ[hpcBimTheoryFAPath <> ".mod"] &&
   !FileExistsQ[hpcBimTheoryFAPath],
 Print[
  "ERROR: Bimetric FeynArts model file not found at: ",
  hpcBimTheoryFAPath
  ];
 Abort[];
 ];

(*Post-processes amplitudes with LeafCount >= 1000. The numerator is
rewritten to expose a recurring parameter combination, after which
FullSimplify is attempted with a time limit of 600 seconds per amplitude.
The simplified result is retained only if its LeafCount is smaller.
This post-processing is kept separate from BimetricAmplitudes.wl since
symbolic simplification times can vary considerably between amplitudes.*)
PostSimplify[amp_]:=Module[{ampTogether,ampNumerator,ampDenominator,intermediate,result,betaPlaceHolder},

If[LeafCount[amp]<1000,Return[amp]];

ampTogether=Together[amp];
ampNumerator=Numerator[ampTogether];
ampDenominator=Denominator[ampTogether];

(*Rewrite BetaBim1 in terms of the recurring combination c BetaBim1-c^3 BetaBim3. This combination occurs often enough in many of the amplitudes to yield considerable simplifications.*)
ampNumerator=ampNumerator/. {BetaBim1->(c^3*BetaBim3+betaPlaceHolder)/c};
intermediate=TimeConstrained[FullSimplify[ampNumerator, Assumptions->s>0],600,ampNumerator];
result=intermediate/. {betaPlaceHolder->c*BetaBim1-c^3*BetaBim3};
result=result/ampDenominator;
If[LeafCount[result]<LeafCount[amp],result,amp]]

(* ------------------------------------------------------------ *)
(* Main Code*)
(* ------------------------------------------------------------ *)

ParallelEvaluate[
  Needs["FeynCalc`"];
  ScalarProduct[A,B] = 7;
  Print["Scalar Product should be 7 = ", Contract[FV[A,m]FV[B,m]]];
]

(*Distribute PostSimplify across the kernels*)
DistributeDefinitions[PostSimplify];

(*Import the functions in BimetricAmplitudes.wl*)
Get[codeFile];
ParallelEvaluate[Get[codeFile]];

(*Check if BimetricAmplitudes.wl has been loaded properly*)
ParallelEvaluate[Names["Global`SplitSingleRoot"]]
ParallelEvaluate[Names["Global`ExpandDummyPolsAndSetMandelstam"]]

(*This is the list of distinct bimetric helicity scattering amplitudes we can have since the expansion of 
  bimetric theory has time reversal symmetry and parity symmetry together with swap symmetry. 
  See: GetOrbitsBinAmps.nb for where these come.*)
scatteringProcesses = {{{1,2},{1,2},{1,2},{1,2}},{{1,2},{1,2},{1,2},{1,-2}},{{1,2},{1,2},{1,-2},{1,-2}},{{1,2},{1,2},{2,2},{2,2}},{{1,2},{1,2},{2,2},{2,1}},{{1,2},{1,2},{2,2},{2,0}},{{1,2},{1,2},{2,2},{2,-1}},{{1,2},{1,2},{2,2},{2,-2}},{{1,2},{1,2},{2,1},{2,1}},{{1,2},{1,2},{2,1},{2,0}},{{1,2},{1,2},{2,1},{2,-1}},{{1,2},{1,2},{2,1},{2,-2}},{{1,2},{1,2},{2,0},{2,0}},{{1,2},{1,2},{2,0},{2,-1}},{{1,2},{1,2},{2,0},{2,-2}},{{1,2},{1,2},{2,-1},{2,-1}},{{1,2},{1,2},{2,-1},{2,-2}},{{1,2},{1,2},{2,-2},{2,-2}},{{1,2},{1,-2},{1,2},{1,-2}},{{1,2},{1,-2},{2,2},{2,2}},{{1,2},{1,-2},{2,2},{2,1}},{{1,2},{1,-2},{2,2},{2,0}},{{1,2},{1,-2},{2,2},{2,-1}},{{1,2},{1,-2},{2,2},{2,-2}},{{1,2},{1,-2},{2,1},{2,1}},{{1,2},{1,-2},{2,1},{2,0}},{{1,2},{1,-2},{2,1},{2,-1}},{{1,2},{1,-2},{2,0},{2,0}},{{1,2},{2,2},{1,2},{2,2}},{{1,2},{2,2},{1,2},{2,1}},{{1,2},{2,2},{1,2},{2,0}},{{1,2},{2,2},{1,2},{2,-1}},{{1,2},{2,2},{1,2},{2,-2}},{{1,2},{2,2},{1,-2},{2,2}},{{1,2},{2,2},{1,-2},{2,1}},{{1,2},{2,2},{1,-2},{2,0}},{{1,2},{2,2},{1,-2},{2,-1}},{{1,2},{2,2},{1,-2},{2,-2}},{{1,2},{2,2},{2,2},{2,2}},{{1,2},{2,2},{2,2},{2,1}},{{1,2},{2,2},{2,2},{2,0}},{{1,2},{2,2},{2,2},{2,-1}},{{1,2},{2,2},{2,2},{2,-2}},{{1,2},{2,2},{2,1},{2,1}},{{1,2},{2,2},{2,1},{2,0}},{{1,2},{2,2},{2,1},{2,-1}},{{1,2},{2,2},{2,1},{2,-2}},{{1,2},{2,2},{2,0},{2,0}},{{1,2},{2,2},{2,0},{2,-1}},{{1,2},{2,2},{2,0},{2,-2}},{{1,2},{2,2},{2,-1},{2,-1}},{{1,2},{2,2},{2,-1},{2,-2}},{{1,2},{2,2},{2,-2},{2,-2}},{{1,2},{2,1},{1,2},{2,1}},{{1,2},{2,1},{1,2},{2,0}},{{1,2},{2,1},{1,2},{2,-1}},{{1,2},{2,1},{1,2},{2,-2}},{{1,2},{2,1},{1,-2},{2,2}},{{1,2},{2,1},{1,-2},{2,1}},{{1,2},{2,1},{1,-2},{2,0}},{{1,2},{2,1},{1,-2},{2,-1}},{{1,2},{2,1},{2,2},{2,2}},{{1,2},{2,1},{2,2},{2,1}},{{1,2},{2,1},{2,2},{2,0}},{{1,2},{2,1},{2,2},{2,-1}},{{1,2},{2,1},{2,2},{2,-2}},{{1,2},{2,1},{2,1},{2,1}},{{1,2},{2,1},{2,1},{2,0}},{{1,2},{2,1},{2,1},{2,-1}},{{1,2},{2,1},{2,1},{2,-2}},{{1,2},{2,1},{2,0},{2,0}},{{1,2},{2,1},{2,0},{2,-1}},{{1,2},{2,1},{2,0},{2,-2}},{{1,2},{2,1},{2,-1},{2,-1}},{{1,2},{2,1},{2,-1},{2,-2}},{{1,2},{2,1},{2,-2},{2,-2}},{{1,2},{2,0},{1,2},{2,0}},{{1,2},{2,0},{1,2},{2,-1}},{{1,2},{2,0},{1,2},{2,-2}},{{1,2},{2,0},{1,-2},{2,2}},{{1,2},{2,0},{1,-2},{2,1}},{{1,2},{2,0},{1,-2},{2,0}},{{1,2},{2,0},{2,2},{2,2}},{{1,2},{2,0},{2,2},{2,1}},{{1,2},{2,0},{2,2},{2,0}},{{1,2},{2,0},{2,2},{2,-1}},{{1,2},{2,0},{2,2},{2,-2}},{{1,2},{2,0},{2,1},{2,1}},{{1,2},{2,0},{2,1},{2,0}},{{1,2},{2,0},{2,1},{2,-1}},{{1,2},{2,0},{2,1},{2,-2}},{{1,2},{2,0},{2,0},{2,0}},{{1,2},{2,0},{2,0},{2,-1}},{{1,2},{2,0},{2,0},{2,-2}},{{1,2},{2,0},{2,-1},{2,-1}},{{1,2},{2,0},{2,-1},{2,-2}},{{1,2},{2,0},{2,-2},{2,-2}},{{1,2},{2,-1},{1,2},{2,-1}},{{1,2},{2,-1},{1,2},{2,-2}},{{1,2},{2,-1},{1,-2},{2,2}},{{1,2},{2,-1},{1,-2},{2,1}},{{1,2},{2,-1},{2,2},{2,2}},{{1,2},{2,-1},{2,2},{2,1}},{{1,2},{2,-1},{2,2},{2,0}},{{1,2},{2,-1},{2,2},{2,-1}},{{1,2},{2,-1},{2,2},{2,-2}},{{1,2},{2,-1},{2,1},{2,1}},{{1,2},{2,-1},{2,1},{2,0}},{{1,2},{2,-1},{2,1},{2,-1}},{{1,2},{2,-1},{2,1},{2,-2}},{{1,2},{2,-1},{2,0},{2,0}},{{1,2},{2,-1},{2,0},{2,-1}},{{1,2},{2,-1},{2,0},{2,-2}},{{1,2},{2,-1},{2,-1},{2,-1}},{{1,2},{2,-1},{2,-1},{2,-2}},{{1,2},{2,-1},{2,-2},{2,-2}},{{1,2},{2,-2},{1,2},{2,-2}},{{1,2},{2,-2},{1,-2},{2,2}},{{1,2},{2,-2},{2,2},{2,2}},{{1,2},{2,-2},{2,2},{2,1}},{{1,2},{2,-2},{2,2},{2,0}},{{1,2},{2,-2},{2,2},{2,-1}},{{1,2},{2,-2},{2,2},{2,-2}},{{1,2},{2,-2},{2,1},{2,1}},{{1,2},{2,-2},{2,1},{2,0}},{{1,2},{2,-2},{2,1},{2,-1}},{{1,2},{2,-2},{2,1},{2,-2}},{{1,2},{2,-2},{2,0},{2,0}},{{1,2},{2,-2},{2,0},{2,-1}},{{1,2},{2,-2},{2,0},{2,-2}},{{1,2},{2,-2},{2,-1},{2,-1}},{{1,2},{2,-2},{2,-1},{2,-2}},{{1,2},{2,-2},{2,-2},{2,-2}},{{2,2},{2,2},{2,2},{2,2}},{{2,2},{2,2},{2,2},{2,1}},{{2,2},{2,2},{2,2},{2,0}},{{2,2},{2,2},{2,2},{2,-1}},{{2,2},{2,2},{2,2},{2,-2}},{{2,2},{2,2},{2,1},{2,1}},{{2,2},{2,2},{2,1},{2,0}},{{2,2},{2,2},{2,1},{2,-1}},{{2,2},{2,2},{2,1},{2,-2}},{{2,2},{2,2},{2,0},{2,0}},{{2,2},{2,2},{2,0},{2,-1}},{{2,2},{2,2},{2,0},{2,-2}},{{2,2},{2,2},{2,-1},{2,-1}},{{2,2},{2,2},{2,-1},{2,-2}},{{2,2},{2,2},{2,-2},{2,-2}},{{2,2},{2,1},{2,2},{2,1}},{{2,2},{2,1},{2,2},{2,0}},{{2,2},{2,1},{2,2},{2,-1}},{{2,2},{2,1},{2,2},{2,-2}},{{2,2},{2,1},{2,1},{2,1}},{{2,2},{2,1},{2,1},{2,0}},{{2,2},{2,1},{2,1},{2,-1}},{{2,2},{2,1},{2,1},{2,-2}},{{2,2},{2,1},{2,0},{2,0}},{{2,2},{2,1},{2,0},{2,-1}},{{2,2},{2,1},{2,0},{2,-2}},{{2,2},{2,1},{2,-1},{2,-1}},{{2,2},{2,1},{2,-1},{2,-2}},{{2,2},{2,0},{2,2},{2,0}},{{2,2},{2,0},{2,2},{2,-1}},{{2,2},{2,0},{2,2},{2,-2}},{{2,2},{2,0},{2,1},{2,1}},{{2,2},{2,0},{2,1},{2,0}},{{2,2},{2,0},{2,1},{2,-1}},{{2,2},{2,0},{2,1},{2,-2}},{{2,2},{2,0},{2,0},{2,0}},{{2,2},{2,0},{2,0},{2,-1}},{{2,2},{2,0},{2,0},{2,-2}},{{2,2},{2,0},{2,-1},{2,-1}},{{2,2},{2,-1},{2,2},{2,-1}},{{2,2},{2,-1},{2,2},{2,-2}},{{2,2},{2,-1},{2,1},{2,1}},{{2,2},{2,-1},{2,1},{2,0}},{{2,2},{2,-1},{2,1},{2,-1}},{{2,2},{2,-1},{2,1},{2,-2}},{{2,2},{2,-1},{2,0},{2,0}},{{2,2},{2,-1},{2,0},{2,-1}},{{2,2},{2,-1},{2,-1},{2,-1}},{{2,2},{2,-2},{2,2},{2,-2}},{{2,2},{2,-2},{2,1},{2,1}},{{2,2},{2,-2},{2,1},{2,0}},{{2,2},{2,-2},{2,1},{2,-1}},{{2,2},{2,-2},{2,0},{2,0}},{{2,1},{2,1},{2,1},{2,1}},{{2,1},{2,1},{2,1},{2,0}},{{2,1},{2,1},{2,1},{2,-1}},{{2,1},{2,1},{2,0},{2,0}},{{2,1},{2,1},{2,0},{2,-1}},{{2,1},{2,1},{2,-1},{2,-1}},{{2,1},{2,0},{2,1},{2,0}},{{2,1},{2,0},{2,1},{2,-1}},{{2,1},{2,0},{2,0},{2,0}},{{2,1},{2,0},{2,0},{2,-1}},{{2,1},{2,-1},{2,1},{2,-1}},{{2,1},{2,-1},{2,0},{2,0}},{{2,0},{2,0},{2,0},{2,0}}};
scatteringProcessesFlatten = Map[Function[x, Flatten[x,1]], scatteringProcesses];

(*We only need to consider certain sectors like GG->GG, GG->MM, GM->MM. Can be seen looking at scatteringProcesses.*)
processes = DeleteDuplicates[Map[Function[x,{x[[1]][[1]],x[[2]][[1]],x[[3]][[1]],x[[4]][[1]]}],scatteringProcesses]];
amplitudeStructures = Generate2To2AmpStructures[processes, hpcBimTheoryFAPath];

Print["Starting computations of all amplitudes."];
bimetricAmpsList = ParallelMap[Function[x, Calc2To2BimAmpFromStructuresSerial[x, amplitudeStructures, savePath]], scatteringProcessesFlatten];
Print["Computation finished. Polishing the results"];
bimetricAmpsSimplifiedList = ParallelMap[PostSimplify, bimetricAmpsList];
Print["Polishing finished"];
bimetricAmps = {scatteringProcessesFlatten, bimetricAmpsSimplifiedList} // Transpose;
Put[
 bimetricAmps,
 FileNameJoin[{savePath, "allBimetricAmplitudes.wl"}]
];
