(* ::Package:: *)

(* Lightweight browser for the precomputed bimetric helicity amplitudes.
   Place this file and BimetricAmplitudeBrowser.nb in the same directory as
   allBimetricAmplitudes.wl. No FeynArts or FeynCalc installation is needed. *)

ClearAll[
  A, AForm, AInfo, AVTheta, ToVTheta, SimplifyResult,
  Pretty, AvailableHelicities,
  LoadBimetricAmplitudes, AmplitudeDataFile,
  amplitudeBrowserDirectory, resolveAmplitudeFile,
  validAmplitudeDataQ, particleCode, validHelicityQ,
  makeAmplitudeKey, amplitudeKeyString, ensureAmplitudeData,
  swapMandelstamTU, symmetryOrbitCandidates,
  findStoredRepresentative, decodeAmplitudeKey,
  amplitudeLookup, legsFromRules, representativeInfoFromLegs,
  amplitudeLabel, amplitudeFormFromLegs, particleTypeset,
  kallenPolynomial, vThetaKinematics, toVThetaFromLegs
];

$BimetricAmplitudeData = {};
$BimetricAmplitudeIndex = <||>;
$BimetricAmplitudeFile = Missing["NotLoaded"];

LoadBimetricAmplitudes::nofile =
  "Could not find allBimetricAmplitudes.wl in `1`. Pass an explicit file name to LoadBimetricAmplitudes[file].";
LoadBimetricAmplitudes::ambiguous =
  "More than one file matching allBimetricAmplitudes*.wl was found in `1`: `2`. Rename the desired file to allBimetricAmplitudes.wl or pass its path explicitly.";
LoadBimetricAmplitudes::invalid =
  "The file `1` does not have the expected form {{key, amplitude}, ...}, with keys of length eight.";
LoadBimetricAmplitudes::duplicate =
  "The file `1` contains duplicate amplitude keys.";

amplitudeBrowserDirectory[] := Module[{directory},
  directory = Quiet @ Check[NotebookDirectory[], $Failed];
  If[StringQ[directory], directory, Directory[]]
];

resolveAmplitudeFile[Automatic] := Module[{directory, preferredFile, candidates},
  directory = amplitudeBrowserDirectory[];
  preferredFile = FileNameJoin[{directory, "allBimetricAmplitudes.wl"}];

  If[FileExistsQ[preferredFile], Return[preferredFile]];

  candidates = FileNames["allBimetricAmplitudes*.wl", directory];
  Which[
    Length[candidates] == 1,
      First[candidates],
    Length[candidates] == 0,
      Message[LoadBimetricAmplitudes::nofile, directory];
      Failure["AmplitudeFileNotFound", <|"Directory" -> directory|>],
    True,
      Message[LoadBimetricAmplitudes::ambiguous, directory, candidates];
      Failure["AmbiguousAmplitudeFiles", <|"Files" -> candidates|>]
  ]
];

resolveAmplitudeFile[file_String] := Module[{directory, candidate},
  directory = amplitudeBrowserDirectory[];
  candidate = If[FileExistsQ[file], file, FileNameJoin[{directory, file}]];
  If[
    FileExistsQ[candidate],
    ExpandFileName[candidate],
    Message[LoadBimetricAmplitudes::nofile, DirectoryName[candidate]];
    Failure["AmplitudeFileNotFound", <|"File" -> candidate|>]
  ]
];

validAmplitudeDataQ[data_] :=
  ListQ[data] && AllTrue[
    data,
    Function[row,
      ListQ[row] && Length[row] == 2 &&
      ListQ[row[[1]]] && Length[row[[1]]] == 8
    ]
  ];

(* Lookup treats a list as a list of several keys.  The stored eight-entry
   process keys are therefore encoded as InputForm strings in the index. *)
amplitudeKeyString[key_List] := ToString[key, InputForm];

LoadBimetricAmplitudes[file_: Automatic] := Module[{path, data, keys},
  path = resolveAmplitudeFile[file];
  If[FailureQ[path], Return[path]];

  data = Quiet @ Check[Get[path], $Failed];
  If[data === $Failed || !validAmplitudeDataQ[data],
    Message[LoadBimetricAmplitudes::invalid, path];
    Return[Failure["InvalidAmplitudeData", <|"File" -> path|>]]
  ];

  keys = First /@ data;
  If[!DuplicateFreeQ[keys],
    Message[LoadBimetricAmplitudes::duplicate, path];
    Return[Failure["DuplicateAmplitudeKeys", <|"File" -> path|>]]
  ];

  $BimetricAmplitudeData = data;
  $BimetricAmplitudeIndex = AssociationThread[
    amplitudeKeyString /@ keys,
    Last /@ data
  ];
  $BimetricAmplitudeFile = path;

  <|
    "File" -> FileNameTake[path],
    "Directory" -> DirectoryName[path],
    "AmplitudesLoaded" -> Length[data]
  |>
];

AmplitudeDataFile[] := $BimetricAmplitudeFile;

ensureAmplitudeData[] := Module[{loadResult},
  If[AssociationQ[$BimetricAmplitudeIndex] && Length[$BimetricAmplitudeIndex] > 0,
    Return[True]
  ];
  loadResult = LoadBimetricAmplitudes[];
  If[FailureQ[loadResult], loadResult, True]
];

particleCode[particle_] := Which[
  SameQ[particle, G] || SameQ[particle, "G"] || SameQ[particle, 1], 1,
  SameQ[particle, M] || SameQ[particle, "M"] || SameQ[particle, 2], 2,
  True, Missing["InvalidParticle", particle]
];

validHelicityQ[1, helicity_] := MemberQ[{-2, 2}, helicity];
validHelicityQ[2, helicity_] := MemberQ[Range[-2, 2], helicity];

makeAmplitudeKey[legs_List] := Module[{codes, helicities, invalidPairs},
  If[
    Length[legs] != 4 || !AllTrue[legs, MatchQ[#, {_, _Integer}] &],
    Return @ Failure[
      "InvalidLegs",
      <|
        "MessageTemplate" ->
          "Expected four legs of the form {particle, helicity}.",
        "Legs" -> legs
      |>
    ]
  ];

  codes = particleCode /@ legs[[All, 1]];
  If[AnyTrue[codes, MissingQ],
    Return @ Failure[
      "InvalidParticle",
      <|
        "MessageTemplate" -> "Particles must be G or M (equivalently 1 or 2).",
        "Legs" -> legs
      |>
    ]
  ];

  helicities = legs[[All, 2]];
  invalidPairs = Select[
    Transpose[{codes, helicities}],
    !validHelicityQ[#[[1]], #[[2]]] &
  ];
  If[invalidPairs =!= {},
    Return @ Failure[
      "InvalidHelicity",
      <|
        "MessageTemplate" ->
          "Allowed helicities are +/-2 for G and -2,-1,0,1,2 for M.",
        "InvalidParticleHelicityPairs" -> invalidPairs
      |>
    ]
  ];

  Flatten @ Transpose[{codes, helicities}]
];

(* A simultaneous t <-> u replacement, protected against sequential rule application. *)
swapMandelstamTU[expression_] := Module[{temporaryT, temporaryU},
  expression /. {t -> temporaryT, u -> temporaryU} /.
    {temporaryT -> u, temporaryU -> t}
];

(* The symmetry generators act on the ordered legs {A,B,C,D} as
     s_i: {B,A,C,D},       t <-> u, phase (-1)^(hA+hB+hC+hD),
     s_f: {A,B,D,C},       t <-> u, phase (-1)^(hA+hB+hC+hD),
     tau: {C,D,A,B},       phase (-1)^(hA+hB+hC+hD),
     pi:  hX -> -hX.
   The canonical products s_i^a s_f^b tau^c pi^d give all 16 elements of
   D8 x Z2.  Each candidate stores the relation from the requested amplitude
   to the transformed amplitude. *)
symmetryOrbitCandidates[legs_List] := Module[
  {flags, helicitySum},
  flags = SortBy[Tuples[{0, 1}, 4], Total];
  helicitySum = Total[legs[[All, 2]]];

  Map[
    Function[flag,
      Module[
        {initialExchange, finalExchange, timeReversal, parity,
         transformedLegs, transformations},

        {initialExchange, finalExchange, timeReversal, parity} = flag;
        transformedLegs = legs;
        transformations = {};

        If[initialExchange == 1,
          transformedLegs = transformedLegs[[{2, 1, 3, 4}]];
          AppendTo[transformations, "InitialExchange"]
        ];
        If[finalExchange == 1,
          transformedLegs = transformedLegs[[{1, 2, 4, 3}]];
          AppendTo[transformations, "FinalExchange"]
        ];
        If[timeReversal == 1,
          transformedLegs = transformedLegs[[{3, 4, 1, 2}]];
          AppendTo[transformations, "TimeReversal"]
        ];
        If[parity == 1,
          transformedLegs = ({#[[1]], -#[[2]]} &) /@ transformedLegs;
          AppendTo[transformations, "Parity"]
        ];

        <|
          "RepresentativeLegs" -> transformedLegs,
          "RepresentativeKey" -> makeAmplitudeKey[transformedLegs],
          "Transformations" -> transformations,
          "SwapTU" -> OddQ[initialExchange + finalExchange],
          "Phase" -> If[
            OddQ[initialExchange + finalExchange + timeReversal],
            (-1)^helicitySum,
            1
          ]
        |>
      ]
    ],
    flags
  ]
];

decodeAmplitudeKey[key_List] :=
  ({If[#[[1]] == 1, G, M], #[[2]]} &) /@ Partition[key, 2];

findStoredRepresentative[legs_List] := Module[
  {loadResult, requestedKey, normalizedLegs, candidates, representative,
   representativeAmplitude, transformedAmplitude},

  loadResult = ensureAmplitudeData[];
  If[FailureQ[loadResult], Return[loadResult]];

  requestedKey = makeAmplitudeKey[legs];
  If[FailureQ[requestedKey], Return[requestedKey]];
  normalizedLegs = Partition[requestedKey, 2];

  (* Amplitudes with exactly one M vanish because no G^n M vertex exists. *)
  If[Count[requestedKey[[{1, 3, 5, 7}]], 2] == 1,
    Return @ <|
      "RequestedLegs" -> normalizedLegs,
      "RequestedKey" -> requestedKey,
      "RepresentativeLegs" -> Missing["NotRequired"],
      "RepresentativeKey" -> Missing["NotRequired"],
      "Transformations" -> {"SingleMassiveFieldSelectionRule"},
      "SwapTU" -> False,
      "Phase" -> 1,
      "RepresentativeAmplitude" -> 0,
      "Amplitude" -> 0
    |>
  ];

  candidates = symmetryOrbitCandidates[normalizedLegs];
  representative = SelectFirst[
    candidates,
    KeyExistsQ[
      $BimetricAmplitudeIndex,
      amplitudeKeyString[#["RepresentativeKey"]]
    ] &,
    Missing["NotFound"]
  ];

  If[MissingQ[representative],
    Return @ Failure[
      "RepresentativeNotFound",
      <|
        "MessageTemplate" ->
          "No stored representative was found in the symmetry orbit of the requested process. The data file may be incomplete or use different conventions.",
        "Legs" -> legs,
        "Key" -> requestedKey
      |>
    ]
  ];

  representativeAmplitude = Lookup[
    $BimetricAmplitudeIndex,
    amplitudeKeyString[representative["RepresentativeKey"]]
  ];
  transformedAmplitude = If[
    TrueQ[representative["SwapTU"]],
    swapMandelstamTU[representativeAmplitude],
    representativeAmplitude
  ];
  transformedAmplitude = representative["Phase"] transformedAmplitude;

  Join[
    <|
      "RequestedLegs" -> normalizedLegs,
      "RequestedKey" -> requestedKey
    |>,
    representative,
    <|
      "RepresentativeAmplitude" -> representativeAmplitude,
      "Amplitude" -> transformedAmplitude
    |>
  ]
];

amplitudeLookup[legs_List] := Module[{result},
  result = findStoredRepresentative[legs];
  If[FailureQ[result], result, result["Amplitude"]]
];

legsFromRules[
  Rule[incomingParticles_List, outgoingParticles_List],
  Rule[incomingHelicities_List, outgoingHelicities_List]
] := Module[{particles, helicities},
  particles = Join[incomingParticles, outgoingParticles];
  helicities = Join[incomingHelicities, outgoingHelicities];
  If[
    Length[incomingParticles] != 2 || Length[outgoingParticles] != 2 ||
    Length[incomingHelicities] != 2 || Length[outgoingHelicities] != 2,
    Return @ Failure[
      "InvalidProcessRules",
      <|"MessageTemplate" -> "Both sides of each rule must contain two entries."|>
    ]
  ];
  MapThread[List, {particles, helicities}]
];

(* Raw amplitude lookup. The first two legs are incoming and the last two outgoing. *)
A[leg1_List, leg2_List, leg3_List, leg4_List] :=
  amplitudeLookup[{leg1, leg2, leg3, leg4}];

A[Rule[incoming_List, outgoing_List]] := Module[{legs},
  If[Length[incoming] != 2 || Length[outgoing] != 2,
    Return @ Failure[
      "InvalidProcessRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two legs."|>
    ]
  ];
  legs = Join[incoming, outgoing];
  amplitudeLookup[legs]
];

A[
  particleRule : Rule[_List, _List],
  helicityRule : Rule[_List, _List]
] := Module[{legs},
  legs = legsFromRules[particleRule, helicityRule];
  If[FailureQ[legs], Return[legs]];
  amplitudeLookup[legs]
];

representativeInfoFromLegs[legs_List] := Module[{result, representativeProcess},
  result = findStoredRepresentative[legs];
  If[FailureQ[result], Return[result]];

  representativeProcess = If[
    MissingQ[result["RepresentativeKey"]],
    Missing["NotRequired"],
    decodeAmplitudeKey[result["RepresentativeKey"]]
  ];

  <|
    "RequestedProcess" -> decodeAmplitudeKey[result["RequestedKey"]],
    "RepresentativeProcess" -> representativeProcess,
    "Transformations" -> result["Transformations"],
    "MandelstamMap" -> If[TrueQ[result["SwapTU"]], "t <-> u", "identity"],
    "Phase" -> result["Phase"],
    "RepresentativeAmplitude" -> result["RepresentativeAmplitude"],
    "Amplitude" -> result["Amplitude"]
  |>
];

(* Diagnostic form: report the representative and the symmetry relation used. *)
AInfo[leg1_List, leg2_List, leg3_List, leg4_List] :=
  representativeInfoFromLegs[{leg1, leg2, leg3, leg4}];

AInfo[Rule[incoming_List, outgoing_List]] := Module[{legs},
  If[Length[incoming] != 2 || Length[outgoing] != 2,
    Return @ Failure[
      "InvalidProcessRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two legs."|>
    ]
  ];
  legs = Join[incoming, outgoing];
  representativeInfoFromLegs[legs]
];

AInfo[
  particleRule : Rule[_List, _List],
  helicityRule : Rule[_List, _List]
] := Module[{legs},
  legs = legsFromRules[particleRule, helicityRule];
  If[FailureQ[legs], Return[legs]];
  representativeInfoFromLegs[legs]
];

(* ------------------------------------------------------------ *)
(* Centre-of-mass variables v and theta                         *)
(* ------------------------------------------------------------ *)

kallenPolynomial[x_, y_, z_] := x^2 + y^2 + z^2 - 2 x y - 2 x z - 2 y z;

(* General 2 -> 2 CM kinematics in the conventions of BimetricAmplitudes.wl.
   The input is the ordered particle list {A,B,C,D}. *)
vThetaKinematics[particles_List] := Module[
  {codes, masses, massSquares, lambdaIncoming, lambdaOutgoing,
   massDifferenceProduct, tMinusU, tExpression, uExpression},

  If[Length[particles] != 4,
    Return @ Failure[
      "InvalidParticles",
      <|"MessageTemplate" -> "Expected four ordered particles {A,B,C,D}."|>
    ]
  ];

  codes = particleCode /@ particles;
  If[AnyTrue[codes, MissingQ],
    Return @ Failure[
      "InvalidParticle",
      <|"MessageTemplate" -> "Particles must be G or M (equivalently 1 or 2)."|>
    ]
  ];

  masses = codes /. {1 -> 0, 2 -> mFP};
  massSquares = masses^2;
  lambdaIncoming = kallenPolynomial[s, massSquares[[1]], massSquares[[2]]];
  lambdaOutgoing = kallenPolynomial[s, massSquares[[3]], massSquares[[4]]];
  massDifferenceProduct =
    (massSquares[[1]] - massSquares[[2]])
    (massSquares[[3]] - massSquares[[4]]);

  tMinusU =
    (Sqrt[lambdaIncoming lambdaOutgoing] Cos[\[Theta]] -
      massDifferenceProduct)/s;
  tExpression = (Total[massSquares] - s + tMinusU)/2;
  uExpression = (Total[massSquares] - s - tMinusU)/2;

  {t -> tExpression, u -> uExpression}
];

Options[ToVTheta] = {SimplifyResult -> True};

(* Express an amplitude in {s,v,theta}, using
       v = Sqrt[1 - 4 mFP^2/s],
   and the CM scattering-angle convention of the calculation file. *)
ToVTheta[expression_, particles_List, OptionsPattern[]] := Module[
  {kinematicRules, massRule, assumptions, result},

  kinematicRules = vThetaKinematics[particles];
  If[FailureQ[kinematicRules], Return[kinematicRules]];

  massRule = mFP -> Sqrt[s (1 - v^2)]/2;
  assumptions =
    Element[{s, v, \[Theta]}, Reals] && s > 0 && 0 < v < 1 &&
    0 < \[Theta] < Pi;

  result = (expression /. kinematicRules) /. massRule;
  If[
    TrueQ[OptionValue[SimplifyResult]],
    FullSimplify[result, Assumptions -> assumptions],
    result
  ]
];

ToVTheta[
  expression_,
  Rule[incomingParticles_List, outgoingParticles_List],
  OptionsPattern[]
] := Module[{particles},
  If[Length[incomingParticles] != 2 || Length[outgoingParticles] != 2,
    Return @ Failure[
      "InvalidParticleRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two particles."|>
    ]
  ];
  particles = Join[incomingParticles, outgoingParticles];
  ToVTheta[
    expression,
    particles,
    SimplifyResult -> OptionValue[SimplifyResult]
  ]
];

Options[AVTheta] = Options[ToVTheta];

toVThetaFromLegs[legs_List, simplifyResult_] := Module[{amplitude, particles},
  amplitude = amplitudeLookup[legs];
  If[FailureQ[amplitude], Return[amplitude]];
  particles = legs[[All, 1]];
  ToVTheta[
    amplitude,
    particles,
    SimplifyResult -> simplifyResult
  ]
];

(* Convenience form: look up the amplitude and immediately convert it. *)
AVTheta[
  leg1_List, leg2_List, leg3_List, leg4_List,
  OptionsPattern[]
] := toVThetaFromLegs[
  {leg1, leg2, leg3, leg4},
  OptionValue[SimplifyResult]
];

AVTheta[
  Rule[incoming_List, outgoing_List],
  OptionsPattern[]
] := Module[{legs},
  If[Length[incoming] != 2 || Length[outgoing] != 2,
    Return @ Failure[
      "InvalidProcessRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two legs."|>
    ]
  ];
  legs = Join[incoming, outgoing];
  toVThetaFromLegs[legs, OptionValue[SimplifyResult]]
];

AVTheta[
  particleRule : Rule[_List, _List],
  helicityRule : Rule[_List, _List],
  OptionsPattern[]
] := Module[{legs},
  legs = legsFromRules[particleRule, helicityRule];
  If[FailureQ[legs], Return[legs]];
  toVThetaFromLegs[legs, OptionValue[SimplifyResult]]
];

(* Display-only replacements. Pretty does not alter the algebraic expression returned by A. *)
$AmplitudeDisplayRules = {
  HoldPattern[mFP] :> Subscript[m, Style["FP", FontSlant -> "Plain"]],
  HoldPattern[mPl] :> Subscript[m, Style["Pl", FontSlant -> "Plain"]],
  HoldPattern[BetaBim1] :> Subscript[\[Beta], 1],
  HoldPattern[BetaBim3] :> Subscript[\[Beta], 3]
};

Pretty[expression_] := TraditionalForm[expression /. $AmplitudeDisplayRules];

particleTypeset[1] := Style["G", FontSlant -> "Italic"];
particleTypeset[2] := Style["M", FontSlant -> "Italic"];

amplitudeLabel[legs_List] := Module[{key, particles, helicities},
  key = makeAmplitudeKey[legs];
  If[FailureQ[key], Return[key]];
  particles = key[[{1, 3, 5, 7}]];
  helicities = key[[{2, 4, 6, 8}]];

  Row[{
    Subsuperscript[
      Style["\[ScriptCapitalA]", FontSlant -> "Italic"],
      Row[{
        helicities[[1]], ",", helicities[[2]], ";",
        helicities[[3]], ",", helicities[[4]]
      }],
      Row[{
        particleTypeset[particles[[1]]], particleTypeset[particles[[2]]],
        " \[Rule] ",
        particleTypeset[particles[[3]]], particleTypeset[particles[[4]]]
      }]
    ],
    Row[{"(", Style["s", FontSlant -> "Italic"], ",",
      Style["t", FontSlant -> "Italic"], ",",
      Style["u", FontSlant -> "Italic"], ")"}]
  }]
];

amplitudeFormFromLegs[legs_List] := Module[{amplitude, label},
  amplitude = amplitudeLookup[legs];
  If[FailureQ[amplitude], Return[amplitude]];
  label = amplitudeLabel[legs];
  Row[{label, " = ", Pretty[amplitude]}]
];

AForm[leg1_List, leg2_List, leg3_List, leg4_List] :=
  amplitudeFormFromLegs[{leg1, leg2, leg3, leg4}];

AForm[Rule[incoming_List, outgoing_List]] := Module[{legs},
  If[Length[incoming] != 2 || Length[outgoing] != 2,
    Return @ Failure[
      "InvalidProcessRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two legs."|>
    ]
  ];
  legs = Join[incoming, outgoing];
  amplitudeFormFromLegs[legs]
];

AForm[
  particleRule : Rule[_List, _List],
  helicityRule : Rule[_List, _List]
] := Module[{legs},
  legs = legsFromRules[particleRule, helicityRule];
  If[FailureQ[legs], Return[legs]];
  amplitudeFormFromLegs[legs]
];

AvailableHelicities[particle1_, particle2_, particle3_, particle4_] :=
  Module[{loadResult, codes},
    loadResult = ensureAmplitudeData[];
    If[FailureQ[loadResult], Return[loadResult]];

    codes = particleCode /@ {particle1, particle2, particle3, particle4};
    If[AnyTrue[codes, MissingQ],
      Return @ Failure[
        "InvalidParticle",
        <|"MessageTemplate" -> "Particles must be G or M (equivalently 1 or 2)."|>
      ]
    ];

    Cases[
      $BimetricAmplitudeData,
      {key_, _} /; key[[{1, 3, 5, 7}]] === codes :>
        key[[{2, 4, 6, 8}]]
    ]
  ];

AvailableHelicities[Rule[incoming_List, outgoing_List]] := Module[{particles},
  particles = Join[incoming, outgoing];
  If[Length[incoming] != 2 || Length[outgoing] != 2,
    Return @ Failure[
      "InvalidParticleRule",
      <|"MessageTemplate" -> "The incoming and outgoing sides must each contain two particles."|>
    ]
  ];
  Apply[AvailableHelicities, particles]
];
