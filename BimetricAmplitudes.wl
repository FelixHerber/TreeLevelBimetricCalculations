(* ::Package:: *)

(*References used here:
F. Herber (2025) The Scattering Amplitudes of Bimetric Theory. 
[Master's thesis, Stockholm University]. Diva Portal. Link:
https://urn.kb.se/resolve?urn=urn:nbn:se:su:diva-248139*)

(*The path to the bimetric theory modefile, provided this notebook sits in the same folder as it.
  Else, replace this with the path to the xActExpansionOfBimetricTheory_FA folder.
*)
bimTheoryFAPath = FileNameJoin[{NotebookDirectory[],"xActExpansionOfBimetricTheory_FA","xActExpansionOfBimetricTheory_FA"}];
(* ------------------------------------------------------------ *)
(* Error Handling *)
(* ------------------------------------------------------------ *)

(*Basic input validators for the public interface.
Particle labels: 1 = massless spin-2 field 
                 2=massive spin-2 field;
Allowed helicities:
                 particle 1: -2, +2 and "gauge" 
                 particle 2: -2,-1,0,1,2 

In the present bimetric model the only allowed masses are 0 and mFP.*)
GetParticleMass[particle_] := Switch[particle, 1, 0, 2, mFP]
ValidParticleQ[p_] := MemberQ[{1,2},p];
ValidHelicityQ[h_] := MemberQ[{-2,-1,0,1,2,"gauge"},h];
ValidMassQ[m_] := MemberQ[{0,mFP},m]

(*Generic length-four list validator, used for momenta,particles,helicities,and masses in 2->2 processes.*)
ValidFourListQ[x_]:=ListQ[x] && Length[x]== 4

(*These check if the particles have valid types, helcities and masses respectively.*)
ValidParticlesQ[particles_] := ValidFourListQ[particles] && And @@ Map[ValidParticleQ,particles]
ValidHelicitiesQ[helicities_] := ValidFourListQ[helicities] && And @@Map[ValidHelicityQ,helicities]
ValidMassesQ[masses_] := ValidFourListQ[masses] && And @@ Map[ValidMassQ,masses]

(*Checks if helicity is a valid helicity for the particle.*)
ValidHelicityForParticle[particle_,helicity_] := Module[{},
 Switch[particle, 
  1, 
   If[MemberQ[{2,-2,"gauge"}, helicity],
    Return[True],
                   
    Return[Failure["InvalidParticleHelicity",<|"Particle"-> particle,"Helicity"-> helicity,"AllowedHelicities"->{2,-2,"gauge"}|>]]
    ],
  2,
   If[MemberQ[{2,1,0,-1,-2}, helicity],              
    Return[True],
                   
    Return[Failure["InvalidParticleHelicity",<|"Particle"-> particle,"Helicity"-> helicity,"AllowedHelicities"->{2,1,0,-1,-2}|>]]
    ],
  _, 
   Return[Failure["InvalidParticle",<|"Particle"->particle,"AllowedParticles"-> {1,2}|>]]]
 ]

(*Validates the public process input {particle1,helicity1,particle2,helicity2,particle3,helicity3,particle4,helicity4}.*)
ValidProcessQ[process_] := Module[{particleData,validityData,failure},
 If[Not[ListQ[process] && Length[process] == 8],
  Return[Failure["InvalidProcess",<|"Process" -> process,"AllowedProcess" -> "List of length 8 of the form {particle1, helicity1, ..., particle4, helicity4}"|>]]
  ];

 particleData = Partition[process,2];
 validityData = ValidHelicityForParticle @@@ particleData;
 failure=SelectFirst[validityData,FailureQ,None];
 If[failure=!= None,
  Return[failure]
  ];
 True
]

(*Checks if the momenta is a four list with distinct symbols like momenta = {p1,p2,p3,p4}*)
ValidMomentaQ[momenta_] := ValidFourListQ[momenta] && DuplicateFreeQ[momenta] && VectorQ[momenta, MatchQ[#, _Symbol]&]

(*Checks if expr contains nested square roots.*)
NestedHalfIntegerPowerQ[expr_] :=! FreeQ[expr, Power[base_, Rational[_,2]] /; !FreeQ[base, Power[_, Rational[_,2]]]];

(*Looks through an expression and outputs the first time a failure is found if there is any failure.*)
FirstFailure[expr_] := FirstCase[expr, _Failure, None, {0,Infinity}]

(* ------------------------------------------------------------ *)
(* Raw Feynman Diagram Generation *)
(* ------------------------------------------------------------ *)

(* ::The Spin-2 Propagators::*)
MasslessSpin2Numerator[li1p1_,li1p2_,li2p1_,li2p2_,mom_] := I/2*(MetricTensor[li1p1,li2p2]*MetricTensor[li1p2,li2p1]+MetricTensor[li1p1,li2p1]*MetricTensor[li1p2,li2p2]-MetricTensor[li1p1,li1p2]*MetricTensor[li2p1,li2p2]);

(*Transverse-traceless projector, used as a building block in the massive spin-2 propagator.*)
G[inds1_,inds2_, mom_,mass_] := MetricTensor[inds1,inds2]-FV[mom,inds1]FV[mom,inds2]/mass^2;

(*Massive spin-2 propagator*)
MassiveSpin2Numerator[li1p1_,li1p2_,li2p1_,li2p2_,mom_,mass_] := I/2 * (G[li1p1,li2p2, mom,mass]G[li1p2,li2p1, mom,mass]+G[li1p2,li2p2, mom,mass]G[li1p1,li2p1, mom,mass] - 2/3 * G[li1p1,li1p2, mom,mass]G[li2p1,li2p2, mom,mass])// Expand;

(*Replaces all occurrences of FAPropagatorNumerator in amplitude with the spin-2 propagators*)
SetFAPropagatorNumerator[amplitude_, useSimplifiedMassive_ : False] := Module[{ampSetPropagatorNumerator},
 ampSetPropagatorNumerator = amplitude /. {
 FAPropagatorNumerator[li1p1_,li1p2_,li2p1_,li2p2_, mom_,mass_] :> If[PossibleZeroQ[mass],
  MasslessSpin2Numerator[li1p1,li1p2,li2p1,li2p2,mom], 
  If[useSimplifiedMassive,
   SimplifiedMassiveSpin2Numerator[li1p1,li1p2,li2p1,li2p2, mom,mass],
   MassiveSpin2Numerator[li1p1,li1p2,li2p1,li2p2, mom,mass]
   ]
  ]
 };
 ampSetPropagatorNumerator
 ]
(*Takes an expr which can be factored as polarization tensors times some expression and
  separates the polarization tensors in question from that expression.

  Symbolically: Takes an expression of the form: \[CurlyEpsilon]^{h_1}_{\[Mu]_1 \[Nu]_1}(p_1)... \[CurlyEpsilon]^{h_n}_{\[Mu]_n \[Nu]_n}(p_n)\[Times]\[CurlyEpsilon]^{*,r_1}_{\[Alpha]_1 \[Beta]_1}(q_1)...\[Times]X 
  and outputs {X,{\[CurlyEpsilon]_1,...,\[CurlyEpsilon]^*_n}}.*)
SeparatePolTensorsFromExpr[expr_] := Module[{polTensorsInExpr,polTensorsConjInExpr,polTensorsInExprTot,exprNoPolTensors},
 polTensorsInExpr = Cases[expr, PolarizationTensor[A__]];
 polTensorsConjInExpr = Cases[expr, Conjugate[PolarizationTensor][A__]];
 polTensorsInExprTot = Flatten[Join[polTensorsInExpr, polTensorsConjInExpr]];
 exprNoPolTensors = expr/. {PolarizationTensor[A__] -> 1, Conjugate[PolarizationTensor][A__] -> 1};
 Return[{exprNoPolTensors, polTensorsInExprTot}]
]
(*Convert Amplitude from FeynArts to FeynCalc, sets couplings and propagators.
  incomingMomenta, outgoingMomenta are lists of momenta and could have forms such as {p1},{p1,p2},{p1,p2,p3} 
  depending on the process under consideration.

  Note that M$FACouplings is a global variable that will be set as soon as a model is scanned, that is after Get2To2Amp is used.
*)
FeynArtsAmpToFeynCalcAmp[ampFA_, incomingMomenta_, outgoingMomenta_, dimension_ : 4, useSimplifiedMassive_ : False] := 
Module[{ampFCFAConvert,ampSetCouplings},
 ampFCFAConvert = FCFAConvert[ampFA, IncomingMomenta -> incomingMomenta, OutgoingMomenta -> outgoingMomenta, ChangeDimension -> dimension];
 ampSetCouplings = SetFAPropagatorNumerator[Map[FeynAmpDenominatorExplicit,ampFCFAConvert/.{M$FACouplings}[[1]]], useSimplifiedMassive];
 Return[ampSetCouplings]
]
(*Generates the tree-level FeynArts amplitude for a 2->2 process and converts it to FeynCalc notation.
  particles = {T1,T2,T3,T4}, where Ti is 1 or 2. 
  momenta = {p1,p2,p3,p4}, with p1,p2 incoming and p3,p4 outgoing.*)
Get2To2Amp[particles_, momenta_, modelFileFADir_ : bimTheoryFAPath]:= 
Module[{particle1Type, particle2Type, particle3Type, particle4Type, p1, p2, p3, p4, treeLevel2To2Process, ampFA, ampFC},
 If[Not[ValidMomentaQ[momenta]],
  
  Return[Failure["InvalidMomenta",<|"Momenta" -> momenta, "AllowedMomenta" -> "List of 4 distinct undefined symbols."|>]];
  ];
 If[Not[ValidParticlesQ[particles]],
   
   Return[Failure["InvalidParticles",<|"Particles"->particles,"AllowedParticles"-> {1,2}|>]];
  ];

 {particle1Type, particle2Type, particle3Type, particle4Type} = {particles[[1]], particles[[2]], particles[[3]], particles[[4]]};
 {p1, p2, p3, p4}={momenta[[1]], momenta[[2]], momenta[[3]], momenta[[4]]};

 treeLevel2To2Process = InsertFields[CreateTopologies[0, 2-> 2], 
                                     {T[particle1Type],T[particle2Type]}-> {T[particle3Type], T[particle4Type]}, 
                                     GenericModel -> {modelFileFADir,"Close"}, 
                                     Model -> {modelFileFADir,"Close"}, 
                                     InsertionLevel -> {Classes}];
 ampFA = CreateFeynAmp[treeLevel2To2Process];

 ampFC = FeynArtsAmpToFeynCalcAmp[ampFA, {p1,p2},{p3,p4}];
 Return[ampFC]
]

(* ------------------------------------------------------------ *)
(* Set Kinematics *)
(* ------------------------------------------------------------ *)

(*For a given particle with specified momentum, helicity and incoming/outgoing state, {p,h,+(-)I}, give the 
  relevant polarization vectors that will appear in an expression with one such particle.*)
GetRelevantVectors[particle_] := Block[{mom = particle[[1]], helicity = particle[[2]], incoming = particle[[3]]},
 Switch[helicity,
  2, Return[{Pol[mom, 1, incoming], mom}],
  1, Return[{Pol[mom, 1, incoming], Pol[mom, 0, incoming], mom}],
  0, Return[{Pol[mom, 1, incoming], Pol[mom, 0, incoming], Pol[mom, -1, incoming], mom}],
  -1, Return[{Pol[mom, -1, incoming], Pol[mom, 0, incoming],mom}],
  -2, Return[{Pol[mom, -1, incoming],mom}],
  "gauge", Return[{Pol[mom,1,incoming],mom}],
  _, Return[Failure["InvalidHelicity",<|"Helicity"->helicity,"AllowedHelicities"->{-2,-1,0,1,2,"gauge"}|>]]
 ]
]
(*This computes all possible scalar products that can occur between polarization vectors and (or) momentum vectors of the particles 
  with momenta as momenta = {p1,p2,p3,p4}, helicities as helicities = {h1,h2,h3,h4} and masses = {m1,m2,m3,m4} and puts the result in an association. 
  Allows fast conversion of scalar products to Mandelstam expressions*)
GetKinematicsAssociations[momenta_,helicities_,masses_]:=
Module[{particleData,incomingData,relevantPolVectors,relevantScalarProducts,ComputeSP,evaluatedSPs, relevantScalarProductsFinal},

 If[Not[ValidMomentaQ[momenta]],

  Return[Failure["InvalidMomenta",<|"Momenta"->momenta,"AllowedMomenta"-> "List of 4 distinct undefined symbols."|>]]
  ];
 If[Not[ValidHelicitiesQ[helicities]],

  Return[Failure["InvalidHelicities",<|"Helicities"->helicities,"AllowedHelicities"-> {-2,-1,0,1,2,"gauge"}|>]]
  ];
 If[Not[ValidMassesQ[masses]],

  Return[Failure["InvalidMasses",<|"Masses"->masses,"AllowedMasses"-> {0,mFP}|>]]
  ];

 incomingData = {I,I,-I,-I};
 particleData = Transpose[{momenta, helicities, incomingData}];
 relevantPolVectors = Map[GetRelevantVectors, particleData];
 relevantScalarProducts = With[{L=relevantPolVectors}, Tuples[Flatten[L], 2]];
 ComputeSP[vectorPairs_] := Block[{vector1 = vectorPairs[[1]], vector2 = vectorPairs[[2]]}, 
  ScalarProductInCM[vector1, vector2, masses, momenta]
  ];

 evaluatedSPs = Map[ComputeSP, relevantScalarProducts];
 With[{failure = FirstFailure[evaluatedSPs]}, If[failure=!=None,Return[failure]]];

 relevantScalarProductsFinal=Map[Function[x,Pair[Momentum[x[[1]]],Momentum[x[[2]]]]], relevantScalarProducts];
 AssociationThread[relevantScalarProductsFinal -> evaluatedSPs]
]
(*Constructs assumptions on the Mandelstam variables and related expressions for 2 -> 2 scattering. 
 The condition oneMinusCosThetaSquared > 0 is obtained from requiring 1 - Cos[theta]^2 >0 in the CM frame,
 after putting the expression over a common positive denominator and substituting for u. In this bimetric implementation,
 each external mass is either 0 or mFP.*)
KinematicAssumptions[masses_]:= 
Module[{m1=masses[[1]],m2=masses[[2]],m3=masses[[3]],m4=masses[[4]],kinematicAssumptions, positiveMasses,massPositivityConditions,
cos\[Theta],oneMinusCosThetaSquared},

 positiveMasses = Select[masses, Function[x, Not[PossibleZeroQ[x]]]];
 massPositivityConditions = Apply[And, Map[# > 0&, positiveMasses]];
(*This expression corresponds to 1-cos\[Theta]^2 which we want to be positive.
  More precisely, this is obtained from taking 1-cos\.08\[Theta]^2 and then using together to put everything on common denominator.
  That denominator is always positive so it can be disgarded. Then, one demands that the numerator is positive, 
  after having substituted for u = - s - t + m1^2 + m2^2 + m3^2 + m4^2.
*)
 oneMinusCosThetaSquared = -4 s (m2^4 m3^2 + m1^4 m4^2-m1^2 (-m4^4 + m2^2 (m3^2 + m4^2 - s) + m4^2 s + m3^2 (m4^2 - t) + m4^2 t + s t) + s (m3^2 (m4^2 - t) + t ( - m4^2 + s + t)) + m2^2 (m3^4 + (m4^2 - s) t -m3^2 (m4^2 + s + t)));
 kinematicAssumptions = Simplify[(s > m1^2 + m2^2) && (s > m3^2 + m4^2) && t \[Element] Reals && u \[Element] Reals && massPositivityConditions && oneMinusCosThetaSquared > 0];
 Return[kinematicAssumptions]
]
(*Compute the scalar product of vectors v1 and v2 which can be momenta (p) (but must then agree with the allowed momenta in 
  momenta = {p1,p2,p3,p4}) or polarization vectors Pol[p,h,+-I] using the CM frame. 
  The polarization vector Pol[p,h,+-I] is more general, 
  it describes the polarization vector for a particle with momentum p, helicity h = -1,0,1 and which can be ingoing I or outgoing -I.

  For the formulas in this function, see section 6.2 of my F. Herber (2025).
*)
ScalarProductInCM[v1_, v2_, masses_, momenta_, Pol_ : Pol] :=
Module[{m1 = masses[[1]], m2 = masses[[2]],m3 = masses[[3]],m4 = masses[[4]], matrixVectors, p1 = momenta[[1]],p2 = momenta[[2]],p3 = momenta[[3]],
        p4 = momenta[[4]], \[Eta], cos\[Theta], sin\[Theta], pi, pf, p1Mom,p2Mom, p3Mom, p4Mom, pol1Plus, pol1Minus, pol1Zero, pol2Plus, pol2Minus, 
        pol2Zero, pol3PlusC, pol3MinusC, pol3ZeroC, pol4PlusC, pol4MinusC, pol4ZeroC, simplifyAssumptions, vector1, vector2, result},

 (*The Minkowski Metric*)
 \[Eta]={{1, 0, 0, 0},{0, -1, 0, 0},{0, 0, -1, 0},{0, 0, 0, -1}};

 (*Cosine and sine of the scattering angle \[Theta] in the CM frame. *)
 cos\[Theta] =((m1 - m2) * (m1 + m2)*(m3 - m4)*(m3 + m4)+ s*(t-u))/ (Sqrt[m1^4 +(m2^2 - s)^2 - 2 * m1^2 * (m2^2 + s)] * Sqrt[m3^4 + (m4^2 - s)^2 - 2m3^2 * (m4^2 + s)]) // Simplify;
 sin\[Theta] = Sqrt[s * (-2 m1^2 (m2^2 (s-2 (m3^2 + m4^2)) - m3^2 (2 m4^2 + 2 s - t + u) + m3^4 - 2 m4^2 s - m4^2 t + m4^2 u + m4^4 + s^2) + m1^4 (s - 2 (m3^2 + m4^2)) - 2 m2^2 (- m3^2 (2 m4^2 + 2 s + t - u) + m3^4 + m4^2 (-2 s + t - u) + m4^4 + s^2) + m2^4 (s - 2 (m3^2 + m4^2)) + s (-2 m3^2 (m4^2 + s) + m3^4 + (m4^2 - s)^2 - (t - u)^2))]/(Sqrt[-2 * m1^2 * m2^2 - 2m1^2 * s + m1^4 - 2m2^2 * s + m2^4 + s^2]*Sqrt[-2m3^2 * m4^2 - 2m3^2 * s + m3^4 - 2m4^2 * s + m4^4 + s^2])// Simplify;

(*Magnitude of incoming and outgoing three-momenta in the CM frame*)
 pi = Sqrt[((m1^2 - m2^2)^2 - 2s*(m1^2 + m2^2) + s^2) / (4s)];
 pf = Sqrt[((m3^2 - m4^2)^2 - 2s*(m3^2 + m4^2) + s^2)/(4s)];

(*Incoming Four-Momenta*)
 p1Mom={Sqrt[m1^2 + pi^2], 0, 0, pi};
 p2Mom={Sqrt[m2^2 + pi^2], 0, 0, -pi};

(*Outgoing Four-Momenta*)
 p3Mom={Sqrt[m3^2 + pf^2], pf * sin\[Theta], 0, pf * cos\[Theta]} // Simplify;
 p4Mom={Sqrt[m4^2 + pf^2], - pf * sin\[Theta], 0,- pf * cos\[Theta]} // Simplify;

(*Minkowski Metric in mostly minus signature*)
 \[Eta] = {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0},{0, 0, 0, -1}};

(*Incoming Helicity Polarization Vectors*)
 pol1Plus = 1 / Sqrt[2]{0, 1, I, 0};
 pol1Minus = 1 / Sqrt[2]{0, 1, -I, 0};
 pol1Zero = If[PossibleZeroQ[m1], 
             0, 
             1/m1 * {pi, 0, 0, Sqrt[m1^2+pi^2]} // Simplify
             ];

 pol2Plus = -1/Sqrt[2]{0, 1, -I, 0};
 pol2Minus = -1/Sqrt[2]{0, 1, I, 0};
 pol2Zero = If[PossibleZeroQ[m2], 
             0,
             1 / m2 * {pi, 0, 0, - Sqrt[m2^2 + pi^2]}//Simplify
             ];

(*Outgoing Helicity Polarization Vectors*)
 pol3PlusC = 1 / Sqrt[2]{0, cos\[Theta], -I, -sin\[Theta]} // Simplify;
 pol3MinusC = 1/Sqrt[2]{0, cos\[Theta], I, -sin\[Theta]} // Simplify;
 pol3ZeroC = If[PossibleZeroQ[m3], 
              0, 
              1 / m3 * {pf,Sqrt[m3^2 + pf^2] * sin\[Theta], 0, Sqrt[m3^2 + pf^2] * cos\[Theta]} // Simplify];

 pol4PlusC = 1 / Sqrt[2]{0, -cos\[Theta], -I, sin\[Theta]} // Simplify;
 pol4MinusC = 1 / Sqrt[2]{0, -cos\[Theta], I, sin\[Theta]} // Simplify;
 pol4ZeroC = If[PossibleZeroQ[m4], 
              0, 
              1 / m4 * {pf, -Sqrt[m4^2 + pf^2] * sin\[Theta], 0, -Sqrt[m4^2 + pf^2] * cos\[Theta]} // Simplify];

 simplifyAssumptions = KinematicAssumptions[masses];


(*Conversion between abstract vectors and polarization vectors to explicit ones*)
 matrixVectors = <|p1 -> p1Mom, p2 -> p2Mom, p3 -> p3Mom, p4 -> p4Mom, Pol[p1,1,I] -> pol1Plus, Pol[p1, -1, I] -> pol1Minus, 
                 Pol[p1, 0, I] -> pol1Zero, Pol[p2,1,I] -> pol2Plus, Pol[p2, -1, I] -> pol2Minus, Pol[p2, 0, I] -> pol2Zero, 
                 Pol[p3, 1, -I] -> pol3PlusC, Pol[p3, -1, -I] -> pol3MinusC, Pol[p3, 0, -I] -> pol3ZeroC, Pol[p4, 1, -I] -> pol4PlusC, 
                 Pol[p4, -1, -I] -> pol4MinusC, Pol[p4, 0, -I] -> pol4ZeroC|>;

 vector1 = matrixVectors[v1];
 vector2 = matrixVectors[v2];

 If[MissingQ[vector1],
  
  ScalarProductInCM::unknownvec="Unknown vector `1` in ScalarProductInCM.";
  Return[Failure["UnknownVector",<|"MessageTemplate" -> ScalarProductInCM::unknownvec,"MessageParameters" -> {v1},"Vectors" -> {v1,v2}|>]];
  ];
 If[MissingQ[vector2],
  
  ScalarProductInCM::unknownvec="Unknown vector `2` in ScalarProductInCM.";
  Return[Failure["UnknownVector",<|"MessageTemplate" -> ScalarProductInCM::unknownvec,"MessageParameters" -> {v2},"Vectors" -> {v1,v2}|>]];
  ];

 result = FullSimplify[Transpose[vector1] . \[Eta] . vector2, Assumptions -> simplifyAssumptions];
 result
]

(*Takes expr which represents a piece of a 2->2 process involving particles with masses described by masses = {m1,m2,m3,m4}, 
  momenta as in momenta = {p1,p2,p3,p4} and outputs the result computed in the CM frame in terms of Mandelstam variables.
  This function can be very inefficient if expr contains many FeynCalc scalar products. 
*)
SetKinematics[expr_, masses_, momenta_, Pol_ : Pol] := FCE[FCI[expr]/.{
  Pair[Momentum[A_], Momentum[B_]]:>ScalarProductInCM[A,B, masses,momenta, Pol]
  }
 ];

(* ------------------------------------------------------------ *)
(* Go from Dummy Pols to Mandelstam *)
(* ------------------------------------------------------------ *)

(*Strategy:
1. Replace each spin-2 polarization tensor by a product of dummy polarization vectors PolA and PolB.
2. Perform all Lorentz contractions.
3. For each external particle,generate the possible spin-1 polarization-vector decompositions corresponding to its helicity.
4. Construct all combinations of these decompositions.
5. Apply each combination to the contracted expression. This reproduces the result that would have been obtained by expanding
   the polarization tensors before contraction,while avoiding the large intermediate expressions that occur in bimetric-theory vertices.*)

(*Takes an expression. Replaces all polarization tensors in it by the outer product of two dummy polarization vectors, 
  contracts metrics in the expression such that the expression consists of dummy polarization vectors in various scalar products
  and use on shell properties on the polarization tensors to remove excess terms. *)

SetDummyPols[expr_] := Module[{exprInDummyPols,exprContracted,exprUseOnShellProperties},
 exprInDummyPols = expr //. {PolarizationTensor[x_, mom_, A_, B_] -> FV[PolA[mom], A]FV[PolB[mom], B]};
 exprContracted = Contract[exprInDummyPols];
 exprUseOnShellProperties = FCI[exprContracted] /. {
  Pair[Momentum[mom_], Momentum[PolA[mom_]]] -> 0, Pair[Momentum[mom_], Momentum[PolB[mom_]]] -> 0,
  Pair[Momentum[mom_], Momentum[PolA[-mom_]]] -> 0, Pair[Momentum[mom_], Momentum[PolB[-mom_]]] -> 0,
  Pair[Momentum[PolA[mom_]], Momentum[PolB[mom_]]] -> 0
 };
 Return[exprUseOnShellProperties]
]

(*Returns a list of replacement-rule sets for dummy polarization vectors with mom as momentum. Each element corresponds to one term in
 the spin-2 polarization tensor decomposition. Example: helicity 0 returns three rule sets, corresponding
 to \[CurlyEpsilon](+1)\[CurlyEpsilon](-1),\[CurlyEpsilon](-1)\[CurlyEpsilon](+1),\[CurlyEpsilon](0)\[CurlyEpsilon](0) with their appropriate coefficients.*)
PolTensorExpansionRules[mom_, helicity_] := Module[{dummyPolReplaceRules, incoming, decompositionTerm, finalMom},
(*If the momentum is equal to -p, it is outgoing while if it is equal to p it is incoming;*)
 If[FreeQ[mom, -1], 
  incoming = I; 
  finalMom = mom, 
  incoming = -I; 
  finalMom = -mom];

 decompositionTerm[h1_, h2_, coeff_] := {Momentum[PolA[mom]]-> coeff * Momentum[Pol[finalMom, h1, incoming]], Momentum[PolB[mom]] -> coeff * Momentum[Pol[finalMom, h2, incoming]]};

 Switch[helicity,
  2, dummyPolReplaceRules = {decompositionTerm[1, 1, 1]},
  1, dummyPolReplaceRules = {decompositionTerm[1, 0, 1 / Sqrt[Sqrt[2]]], decompositionTerm[0, 1, 1 / Sqrt[Sqrt[2]]]},
  0, dummyPolReplaceRules = {decompositionTerm[1, -1, 1 / Sqrt[Sqrt[6]]], decompositionTerm[-1, 1, 1 / Sqrt[Sqrt[6]]], decompositionTerm[0, 0, I * Sqrt[2] / Sqrt[Sqrt[6]]]},
  -1, dummyPolReplaceRules = {decompositionTerm[-1, 0, 1 / Sqrt[Sqrt[2]]], decompositionTerm[0, -1, 1 / Sqrt[Sqrt[2]]]},
  -2, dummyPolReplaceRules = {decompositionTerm[-1, -1, 1]},
  "gauge", dummyPolReplaceRules = {{Momentum[PolA[mom]] -> 1 / Sqrt[Sqrt[2]] * Momentum[Pol[finalMom, 1, incoming]], Momentum[PolB[mom]] -> 1 /Sqrt[Sqrt[2]] * Momentum[mom]},{Momentum[PolA[mom]] -> 1/Sqrt[Sqrt[2]] * Momentum[mom], Momentum[PolB[mom]] -> 1 / Sqrt[Sqrt[2]] * Momentum[Pol[finalMom, 1, incoming]]}},
  _, Return[Failure["InvalidHelicity",<|"Helicity" -> helicity,"AllowedHelicities"->{-2, -1, 0, 1, 2, "gauge"}|>]]
  ];

 Return[dummyPolReplaceRules]
];

(*Takes an expression with dummy polarization vectors that is fully contracted and where those dummy polarization vectors 
  have momenta as specified in momenta = {p1,p2,p3,p4}. Depending on the helicities, expands these dummy polarization vectors
  in such a way as to reproduce the corresponding helicity and puts the expansion in a list.*)
ExpandPolTensors[exprWithDummyPols_, momenta_, helicities_] := Block[{particleData, ruleFamilies, allRuleSets, exprExpandedPolTensors},

 particleData=Transpose[{momenta, helicities}];

 (*The list of all replacement rules of the dummy polarizations*)
 ruleFamilies = PolTensorExpansionRules@@@particleData;
 With[{failure = FirstFailure[ruleFamilies]}, If[failure =!= None, Return[failure]]];

 allRuleSets = Map[Flatten, Tuples[ruleFamilies]];

 exprExpandedPolTensors = Map[exprWithDummyPols/. #&, allRuleSets];
 Return[exprExpandedPolTensors]
]
(*Takes diagramWithDummyPols, an expression which is written as contractions of the dummy polarizations, 
  where the dummy polarizations have momenta as described in momenta = {p1,p2,p3,p4}. 
  Then given the helicities and masses of the particles that appear in the process and the detailed kinematics of the process 
  [as described by kinematicsAssoc, giving all possible scalar products between the complete polarization vectors] and
   some assumptions, evaluates the dummy polarization vectors such that their helicities are described by helicities, 
   evaluates the scalar products to obtain expressions consisting only of Mandelstam variables. *)
ExpandDummyPolsAndSetMandelstam[diagramWithDummyPols_, momenta_, helicities_, masses_, kinematicsAssoc_, assumptions_] :=
Block[{signedMomenta, diagramInsertPolTensors, diagramSetKin, totMassSquared = Sum[j^2, {j, masses}], diagramSplitRoots, 
       diagramSimplified},
(*Must be used since the last two particles should be outgoing*)
 signedMomenta = {momenta[[1]], momenta[[2]], -momenta[[3]], -momenta[[4]]};

 Print["[ExpandDummyPolsAndSetMandelstam]: Fixing Polarization Tensors with ExpandPolTensors"];
 diagramInsertPolTensors = ExpandPolTensors[diagramWithDummyPols, signedMomenta, helicities];
 With[{failure = FirstFailure[diagramInsertPolTensors]}, If[failure =!= None, Return[failure]]];

 Print["[ExpandDummyPolsAndSetMandelstam]: Setting Kinematics"];
 diagramSetKin = diagramInsertPolTensors /. Normal[kinematicsAssoc];
 Print["[ExpandDummyPolsAndSetMandelstam]: Simplifying. "];

 diagramSplitRoots = Map[Function[x, SplitSingleRoot[x, assumptions]], diagramSetKin];
 With[{failure = FirstFailure[diagramSplitRoots]}, If[failure =!= None, Return[failure]]];

 diagramSimplified = Map[Function[x, Simplify[x /.{u -> - s - t + totMassSquared}, Assumptions -> assumptions]], diagramSplitRoots];
 Print["[ExpandDummyPolsAndSetMandelstam]: Finished Simplify"];
 diagramSimplified
]

(* ------------------------------------------------------------ *)
(* Algebraic Simplification *)
(* ------------------------------------------------------------ *)

(*Takes an expression that contains contractions of dummy polarization vectors. 
  Evaluates this expression and simplifies it, yielding an expression consisting only of Mandelstam variables.*)
SimplifyDummyPolExpr[dummyPolExpr_, momenta_, helicities_, masses_, kinematicsAssoc_]:=
Module[{listOfContributionsToExpr, rootExpandedExpr, simplifiedResult, assumptions, diagramSimplified},

 assumptions = KinematicAssumptions[masses];
 listOfContributionsToExpr = ExpandDummyPolsAndSetMandelstam[dummyPolExpr, momenta, helicities, masses, kinematicsAssoc, assumptions];
 With[{failure = FirstFailure[listOfContributionsToExpr]}, If[failure =!= None, Return[failure]]];

 rootExpandedExpr = SplitSingleRoot[Plus @@ listOfContributionsToExpr, assumptions];
 With[{failure = FirstFailure[rootExpandedExpr]}, If[failure =!= None, Return[failure]]];

 Print["[SimplifyDummyPolExpr]:  Starting simplification of all contributions to a diagram"];
 simplifiedResult = Simplify[rootExpandedExpr, Assumptions -> assumptions];
 Print["[SimplifyDummyPolExpr]: Completed Simplify. "];
 Return[simplifiedResult]
]
(*Takes a sum of diagrams. Tries to use the Mandelstam trick s+t+u = m1^2+m2^2+m3^2+m4^2 to write the result on 
  simplest possible form.*)
MandelstamDiagramsIntoAmp[diagrams_, masses_] := 
Module[{assumptions, ampSimplifyRoots, ampSimplify, ampTrickMandelstam},
 assumptions = KinematicAssumptions[masses];

 ampSimplifyRoots = SplitSingleRoot[Plus @@diagrams, assumptions];
 With[{failure = FirstFailure[ampSimplifyRoots]}, If[failure =!= None, Return[failure]]];
 Print["[MandelstamDiagramsIntoAmp]: Summing all diagrams together for the process, simplifying and using TrickMandelstam."];
 ampSimplify = Simplify[ampSimplifyRoots, Assumptions -> assumptions];

 ampTrickMandelstam = TrickMandelstam[PowerExpand[ampSimplify, s], {s, t, u, Plus @@ (masses^2)}];
 Print["[MandelstamDiagramsIntoAmp]: Completed"];
 Return[ampTrickMandelstam]
]

(* ------------------------------------------------------------ *)
(* Root Management Program *)
(* ------------------------------------------------------------ *)

(* Returns True if expr can be proven positive under the assumptions. Returns False if positivity cannot be proven.*)
IsPositiveDef[expr_, assumptions_] := TrueQ[Simplify[Element[expr, Reals] && expr > 0, Assumptions -> assumptions]]
(*Splits a list of factors into factors proven positive,
  factors proven negative, and factors whose sign cannot be determined under the assumptions.*)
SplitBySignFromList[inputList_, assumptions_] := Block[{positiveElements, negativeElements,realElements, otherElements},

 positiveElements = Select[inputList, Function[x, IsPositiveDef[x, assumptions]]];
 negativeElements = Select[inputList, Function[x, IsPositiveDef[-x, assumptions]]];
 realElements = Join[positiveElements, negativeElements];
 otherElements = DeleteCases[inputList, Alternatives @@ realElements];
 {positiveElements, negativeElements, otherElements}
]
(*For rootArg = N / D and p = rootPower/2, tries to split off any part of N / D that has a definite sign 
  under the assumptions;
  If D > 0, use (N / D) ^ p = N ^ p / D ^ p;
  If D < 0, rewrite N / D=(-N) / (-D), then use -D > 0;
  If N > 0, split only the numerator: (N / D) ^ p = N ^ p (1 / D) ^ p;
  If N < 0, rewrite N / D =(-N)(-1 / D), then use -N > 0;

  This avoids using the unsafe identity (1 / D) ^ p = 1 / D ^ p when D is complex or has undetermined sign.*)
RootTogether[rootArg_, rootPower_, assumptions_] := 
Block[{p, argInRootTogether, argInRootNumerator, argInRootDenominator, splitFractionRoot},

 p = Rational[rootPower, 2];
 argInRootTogether = Together[rootArg];
 argInRootNumerator = Numerator[argInRootTogether];
 argInRootDenominator = Denominator[argInRootTogether];

 Which[
  IsPositiveDef[argInRootDenominator, assumptions],
  splitFractionRoot = Power[argInRootNumerator, p] * 1 / Power[argInRootDenominator, p];,

  IsPositiveDef[-argInRootDenominator, assumptions],
  splitFractionRoot = Power[-argInRootNumerator, p] * 1 / Power[-argInRootDenominator, p];,

  IsPositiveDef[argInRootNumerator, assumptions],
  splitFractionRoot = Power[argInRootNumerator, p] * Power[1 / argInRootDenominator, p];,

  IsPositiveDef[-argInRootNumerator, assumptions],
  splitFractionRoot = Power[-argInRootNumerator, p] * Power[- 1 / argInRootDenominator, p];,

  True,
  splitFractionRoot = Power[argInRootTogether, p];
  ];
 splitFractionRoot
];
(*Factors rootArg and pulls out factors with definite sign as described by assumptions; Positive factors are pulled out directly. 
  Negative factors x are rewritten using -x > 0;, while the accumulated sign is left inside the remaining root.*)
RootFactor[rootArg_, rootPower_, assumptions_] := 
Block[{argInRootFactor, argInRootFactorInList, positiveFactors, negativeFactors, otherFactors,splitFactorRoot, splitPositiveFactors,
       splitNegativeFactors, splitOtherFactors, multiplyPosFactors, multiplyNegFactors, multiplyOtherFactors},

 argInRootFactor = Factor[rootArg];
 If[Head[argInRootFactor] === Times, 

(*Put Each Factor in a List*)
  argInRootFactorInList = List @@ argInRootFactor;
  {positiveFactors, negativeFactors, otherFactors} = SplitBySignFromList[argInRootFactorInList, assumptions];

  splitPositiveFactors = Map[Function[x, Power[x, Rational[rootPower, 2]]], positiveFactors];
  splitNegativeFactors = Map[Function[x, Power[-x, Rational[rootPower, 2]]], negativeFactors];
  splitOtherFactors = Join[{(-1)^(Length[negativeFactors])}, otherFactors];
  multiplyPosFactors = Apply[Times, splitPositiveFactors];
  multiplyNegFactors = Apply[Times, splitNegativeFactors];
  multiplyOtherFactors = Power[Apply[Times, splitOtherFactors], Rational[rootPower, 2]];

  splitFactorRoot = Times[multiplyPosFactors, multiplyNegFactors, multiplyOtherFactors];,

  splitFactorRoot = Power[rootArg, Rational[rootPower, 2]]
  ];
 splitFactorRoot
];
(*Applies RootTogether and RootFactor to all half-integer powers in expr.

Limitations: -expr must be a single expression, not a list.
             -nested half-integer powers such as Sqrt[Sqrt[x]] are not supported.
             -the transformation is conservative: factors are only pulled out when their sign can be proven under the assumptions.*)
SplitSingleRoot[expr_, assumptions_] := Block[{exprSplitFraction, exprSplitFactor},

 If[ListQ[expr],
  Return[Failure["ExpressionIsList", <|"MessageTemplate" -> "SplitSingleRoot cannot handle lists.", "Expression" -> expr|>]]
  ];
 If[NestedHalfIntegerPowerQ[expr],
  Return[Failure["NestedRoot", <|"MessageTemplate" -> "SplitSingleRoot cannot handle nested half-integer powers.", "Expression" -> expr|>]]
  ];

 exprSplitFraction = expr /. {Power[x_, Rational[y_,2]] :> RootTogether[x, y, assumptions]};
 exprSplitFactor=exprSplitFraction /. {Power[x_, Rational[y_, 2]] :> RootFactor[x, y, assumptions]};
 exprSplitFactor
]

(* ------------------------------------------------------------ *)
(* Automatic Bimetric Amplitude Calculation *)
(* ------------------------------------------------------------ *)

(*Computes the 2->2 tree-level amplitude for a bimetric scattering process. 
Input format:{particle1,helicity1,particle2,helicity2,particle3,helicity3,particle4,helicity4}.
The process is particle1+particle2->particle3+particle4.
Particle labels are 1 for the massless spin-2 field and 2 for the massive spin-2 field. 
The output is expressed in Mandelstam variables.

Uses ParallelMap, is not parallelizable
*)
Calc2To2BimAmp[process_] :=
Block[{validity = ValidProcessQ[process], particles, helicities, masses, ampWithDiagramsInList, simplifiedDiagramsInList, 
       kinematicAssoc, ampWithDummyPols, totAmp},

 If[FailureQ[validity], 
  Return[validity]
  ];

 particles = {process[[1]], process[[3]], process[[5]], process[[7]]};
 helicities = {process[[2]], process[[4]], process[[6]], process[[8]]};
 masses = Map[GetParticleMass, particles];

(*This below creates the structures behind the FeynCalc Amplitude.*)
 ampWithDiagramsInList=Get2To2Amp[particles, {p1, p2, p3, p4}];
 If[FailureQ[ampWithDiagramsInList],
  Return[ampWithDiagramsInList]
  ];

(*This code below sets dummy polarizations on each of the diagrams.*)
 ampWithDummyPols = ParallelMap[SetDummyPols, ampWithDiagramsInList];
 With[{failure = FirstFailure[ampWithDummyPols]}, If[failure =!= None, Return[failure]]];

 (*This evaluates the kinematics*)
 kinematicAssoc = GetKinematicsAssociations[{p1, p2, p3, p4}, helicities, masses];
 If[FailureQ[kinematicAssoc], Return[kinematicAssoc]];

 simplifiedDiagramsInList = ParallelMap[Function[x, SimplifyDummyPolExpr[x, {p1, p2, p3, p4}, helicities, masses, kinematicAssoc]], ampWithDummyPols];
 With[{failure = FirstFailure[simplifiedDiagramsInList]}, If[failure =!= None, Return[failure]]];

 totAmp = MandelstamDiagramsIntoAmp[simplifiedDiagramsInList, masses];
 totAmp
]

(* Generates the FeynArts/FeynCalc amplitude structures for a list of
   particle processes.

   Each particle process has the form {T1,T2,T3,T4}, where Ti is a particle
   label, so in {1, 2}. For example, {1,1,2,2} corresponds to
   particle 1 + particle 1 -> particle 2 + particle 2.

   The output is an Association mapping each particle process to the list of
   unsimplified diagram amplitudes. This is useful when many helicity
   amplitudes share the same external particle content, since the expensive
   FeynArts/FeynCalc generation only has to be done once. 
   *)
Generate2To2AmpStructures[particleProcesses_, filename_] := Block[{validParticles, ampStructures, 
                                                            ampStructuresDummyPols, dictOfAmps},
 
 validParticles = And @@ Map[ValidParticlesQ, particleProcesses];
 If[Not[validParticles],
  Return[Failure["InvalidParticles",<|"Particles" -> validParticles,"AllowedParticles"-> {1,2}|>]];
  ];
 Print["[Generate2To2AmpStructures]: Constructing FeynCalc amplitudes"];
 ampStructures = Map[Function[x, Get2To2Amp[x, {p1, p2, p3, p4}, filename]], particleProcesses];
 With[{failure = FirstFailure[ampStructures]}, If[failure =!= None, Return[failure]]];
 
 Print["[Generate2To2AmpStructures]: Setting dummy polarization vectors and contracting."];
 ampStructuresDummyPols = Map[Function[x, SetDummyPols[x]], ampStructures];
 dictOfAmps = AssociationThread[particleProcesses, ampStructuresDummyPols];
 dictOfAmps
]

(* Computes a 2->2 helicity amplitude using precomputed amplitude structures.

   process has the form
     {particle1, helicity1, particle2, helicity2,
      particle3, helicity3, particle4, helicity4}.

   ampStructureDict should be produced by Generate2To2AmpStructures.
   This function avoids regenerating the FeynArts/FeynCalc diagrams and is
   therefore preferable when computing many helicity amplitudes. 
   
   saveAmpFolder describes the folder in which the result should be saved as a .mx file or loaded from.
   If saveAmpFolder = None, the function will not save the result in a .mx file and neither load the result, just compute
   the amplitude for Mathematica.
   
   Uses ParallelMap to simplify multiple diagrams at once. Is not parallelizable.
   *)
Calc2To2BimAmpFromStructures[process_, ampStructureDict_, saveAmpFolder_:None]:= 
Module[{validity = ValidProcessQ[process], useCache, processIdentifier, saveFileName, particles, helicities, masses, 
        currentAmpStructure, kinematicsAssoc, amplitudeSplitRoots, simplifiedDiagrams, summedDiagrams, assumptions, 
        amplitudeSimplify, totMassSquared, amplitudeTrickMandelstam},
 
 If[FailureQ[validity], 
  Return[validity]
  ];                  
 
 particles = {process[[1]], process[[3]], process[[5]], process[[7]]};
 helicities = {process[[2]], process[[4]], process[[6]], process[[8]]};
 masses = {GetParticleMass[process[[1]]], GetParticleMass[process[[3]]], 
           GetParticleMass[process[[5]]], GetParticleMass[process[[7]]]};
 
 If[!KeyExistsQ[ampStructureDict, particles],
 Return[Failure["MissingAmplitudeStructure", <|"Particles" -> particles, 
                "AvailableParticleStructures" -> Keys[ampStructureDict]|>]]
  ]; 
 
 processIdentifier = StringReplace[ToString[process], {"{"-> "", "}"-> "", ","-> "_", " "-> "", "-"-> "Minus"}];
 useCache = saveAmpFolder =!= None;

 If[useCache,
  If[! DirectoryQ[saveAmpFolder],
    CreateDirectory[saveAmpFolder, CreateIntermediateDirectories -> True]
    ];

  saveFileName = FileNameJoin[{saveAmpFolder, "bimetricAmp_" <> processIdentifier <> ".mx"}];

  If[FileExistsQ[saveFileName],
   Print["[Calc2To2BimAmpFromStructures]: Loading saved amplitude for process ", process];
   Get[saveFileName];
   Return[amplitudeTrickMandelstam];
   ];
  ];

 currentAmpStructure = ampStructureDict[particles];
 kinematicsAssoc = GetKinematicsAssociations[{p1, p2, p3, p4}, helicities, masses];
 
 If[FailureQ[kinematicsAssoc], Return[kinematicsAssoc]];
  
 Print["[Calc2To2BimAmpFromStructures]: Starting calculation of process:", process];
 simplifiedDiagrams = ParallelMap[Function[diagram, SimplifyDummyPolExpr[diagram, {p1, p2, p3, p4}, helicities, masses, kinematicsAssoc]], currentAmpStructure];
 With[{failure = FirstFailure[simplifiedDiagrams]}, If[failure =!= None, Return[failure]]];
  
 Print["[Calc2To2BimAmpFromStructures]: Simplifying the sum of diagrams contributing to the process:", process];
  
 amplitudeTrickMandelstam = MandelstamDiagramsIntoAmp[simplifiedDiagrams, masses];
 With[{failure = FirstFailure[amplitudeTrickMandelstam]}, If[failure =!= None, Return[failure]]];
  
 If[useCache,
  Print["[Calc2To2BimAmpFromStructures]: Saving calculation of process ", process];
  DumpSave[saveFileName, amplitudeTrickMandelstam];
  ];
 
amplitudeTrickMandelstam
]

(* Same as Calc2To2BimAmpFromStructures but this is parallelizable, uses Map when simplifying each diagram.
   Good when outer loop is parallelized.
   *)
Calc2To2BimAmpFromStructuresSerial[process_, ampStructureDict_, saveAmpFolder_ : None]:= 
Module[{validity = ValidProcessQ[process], useCache, processIdentifier, saveFileName, particles, helicities, masses, 
        currentAmpStructure, kinematicsAssoc, amplitudeSplitRoots, simplifiedDiagrams, summedDiagrams, assumptions, 
        amplitudeSimplify, totMassSquared, amplitudeTrickMandelstam},
 
 If[FailureQ[validity], 
  Return[validity]
  ];                  
 
 particles = {process[[1]], process[[3]], process[[5]], process[[7]]};
 helicities = {process[[2]], process[[4]], process[[6]], process[[8]]};
 masses = {GetParticleMass[process[[1]]], GetParticleMass[process[[3]]], 
           GetParticleMass[process[[5]]], GetParticleMass[process[[7]]]};
 
 If[!KeyExistsQ[ampStructureDict, particles],
 Return[Failure["MissingAmplitudeStructure", <|"Particles" -> particles, 
                "AvailableParticleStructures" -> Keys[ampStructureDict]|>]]
  ]; 
 
 processIdentifier = StringReplace[ToString[process], {"{"-> "", "}"-> "", ","-> "_", " "-> "", "-"-> "Minus"}];
 useCache = saveAmpFolder =!= None;

 If[useCache,
  If[! DirectoryQ[saveAmpFolder],
    CreateDirectory[saveAmpFolder, CreateIntermediateDirectories -> True]
    ];

  saveFileName = FileNameJoin[{saveAmpFolder, "bimetricAmp_" <> processIdentifier <> ".wl"}];

  If[FileExistsQ[saveFileName],
   Print["[Calc2To2BimAmpFromStructures]: Loading saved amplitude for process ", process];
   Return[Get[saveFileName]];
   ];
  ];

 currentAmpStructure = ampStructureDict[particles];
 kinematicsAssoc = GetKinematicsAssociations[{p1, p2, p3, p4}, helicities, masses];
 
 If[FailureQ[kinematicsAssoc], Return[kinematicsAssoc]];
  
 Print["[Calc2To2BimAmpFromStructures]: Starting calculation of process:", process];
 simplifiedDiagrams = Map[Function[diagram, SimplifyDummyPolExpr[diagram, {p1, p2, p3, p4}, helicities, masses, kinematicsAssoc]], currentAmpStructure];
 With[{failure = FirstFailure[simplifiedDiagrams]}, If[failure =!= None, Return[failure]]];
  
 Print["[Calc2To2BimAmpFromStructures]: Simplifying the sum of diagrams contributing to the process:", process];
  
 amplitudeTrickMandelstam = MandelstamDiagramsIntoAmp[simplifiedDiagrams, masses];
 With[{failure = FirstFailure[amplitudeTrickMandelstam]}, If[failure =!= None, Return[failure]]];
  
 If[useCache,
  Print["[Calc2To2BimAmpFromStructures]: Saving calculation of process ", process];
 Put[amplitudeTrickMandelstam, saveFileName];
];
 
amplitudeTrickMandelstam
]
