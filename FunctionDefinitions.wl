(* ::Package:: *)

(* ::Text:: *)
(*Install the necessary functionality.*)


If[PacletFind["EinsteinSummation"]=={},
ResourceFunction["PacletizeResourceFunction"]["EinsteinSummation"]
]
If[PacletFind["MatrixPartialTrace"]=={},
ResourceFunction["PacletizeResourceFunction"]["MatrixPartialTrace"]
]


(* ::Text:: *)
(*Density operator represented by Bloch vectors.*)


rhoQubit[n_]:=1/2(IdentityMatrix[2]+n . Array[PauliMatrix,3])


(* ::Subsection:: *)
(*Parallel Plotting*)


Options[ParallelDensityPlotData] = {
  "Resolution" -> {10^-2, 10^-2},        (* \[CapitalDelta]x, \[CapitalDelta]y snapping *)
  "InitialGrid" -> {25, 25},             (* seed points *)
  "MaxIters" -> 5,
  "Refinement" -> 3
}
ParallelDensityPlotData[f_, {x_, xmin_, xmax_}, {y_, ymin_, ymax_}, opts : OptionsPattern[]]:=Module[
  {
    res = OptionValue["Resolution"],
    init = OptionValue["InitialGrid"],
    maxIters = OptionValue["MaxIters"],
    refine = OptionValue["Refinement"],
    
    snap, cache = <||>, data, surrogate, proposePoints, evalBatch,
    seedPts, newPts, iter = 0
  },
  (*Function to enforce certain resolution to help with caching*)
  snap[{xx_,yy_}]:={Round[xx,res[[1]]],Round[yy,res[[2]]]};
  (*Takes a batch of points and evaluates the function on them in parallel*)
  evalBatch = 
    Function[pts,
      Module[{todo, vals},
        todo = Select[DeleteDuplicates[pts], ! KeyExistsQ[cache, #] &];
        If[todo === {}, Return[{}]];
          vals = ParallelMap[f @@ # &, todo];
          AssociateTo[cache,Thread[todo->vals]];
          todo
        ]
    ];
    
  seedPts = snap /@ Flatten[
    Table[{xx,yy},
      {xx,Subdivide[xmin,xmax,init[[1]]-1]},
      {yy,Subdivide[ymin,ymax,init[[2]]-1]}
    ],
    1
  ];
  (*First evaluate on the seed grid and cache the results*)
  evalBatch[seedPts];
  
  (*Supplies a dummy function to DensityPlot to leverage its adaptive sampling algorithm to suggest points*)
  proposePoints[s_] :=
    Module[{plt, pts},
      {plt, pts} = Reap[
        DensityPlot[
          s[x, y],
          {x, xmin, xmax}, {y, ymin, ymax},
          PlotPoints -> init, MaxRecursion -> refine,
          EvaluationMonitor :> Sow[snap[{x, y}], "pts"],
          PerformanceGoal -> "Quality"
          ],
          "pts"
        ];
      DeleteDuplicates @ First @ pts
    ];
  
  While[iter < maxIters,
	(*Get plot points from cache*)  
    data = KeyValueMap[{#1,#2}&,cache];
    (*Create a dummy function*)
    surrogate = Interpolation[N[data], InterpolationOrder->1];
    (*Get suggestions and evaluate function*)
    newPts = evalBatch[proposePoints[surrogate]];
    iter++;refine++;
  ];
  (*Return the cache*)
  cache
]


(* ::Subsection:: *)
(*Hamiltonian Evolution*)


(* ::Text:: *)
(*Definition of the 2-qubit Hamiltonian*)


Hamiltonian[\[Omega]_,\:2206\[Omega]_,coupling_:DiagonalMatrix[{JX,JY,JZ}]]:=ArrayFlatten[
	\[Omega]*TensorProduct[PauliMatrix[3],IdentityMatrix[2]]
	+(\[Omega]+\:2206\[Omega])*TensorProduct[IdentityMatrix[2],PauliMatrix[3]]
	+Sum[coupling[[i,j]]*TensorProduct[PauliMatrix[i],PauliMatrix[j]],{i,3},{j,3}]
]


(* ::Text:: *)
(*The system-actuator unitaries are obtained by exponentiating the total Hamiltonian.*)


SetAttributes[Propagator,HoldAll];
Propagator[t_,H_]:=MatrixExp[-I t H]
Propagator[t,HoldPattern[Hamiltonian[\[Omega]_,\[CapitalDelta]\[Omega]_,DiagonalMatrix[{JX_,JY_,JZ_}]]]]=Propagator[t,Hamiltonian[\[Omega],\[CapitalDelta]\[Omega],DiagonalMatrix[{JX,JY,JZ}]]]//FullSimplify;


(* ::Text:: *)
(*We use the simplified formula for the Choi matrix representation of unitaries. The transpose is necessary because Mathematica follows row-major convention whereas out mathematical formalism is formulated in a column-major convention.*)


SetAttributes[PropagatorChoiMatrix,HoldAll];
PropagatorChoiMatrix[t_:t,H_:Hamiltonian[]]:=TensorProduct[
	Flatten[Propagator[t,H]//Transpose],
	Flatten[Propagator[-t,Transpose[H]]//Transpose]
]
PropagatorChoiMatrix[t,HoldPattern[Hamiltonian[\[Omega]_,\[CapitalDelta]\[Omega]_,DiagonalMatrix[{JX_,JY_,JZ_}]]]]=PropagatorChoiMatrix[t,Hamiltonian[\[Omega],\[CapitalDelta]\[Omega],DiagonalMatrix[{JX,JY,JZ}]]]//FullSimplify;


(* ::Text:: *)
(*To link together the system-actuator unitaries we first need to reshape them into the appropriate tensor whose indices correspond to the tensor product of system-actuator input and output Hilbert spaces.*)


SetAttributes[PropagatorChoiTensor,HoldAll];
PropagatorChoiTensor[t_:t,H_:Hamiltonian[]]:=Fold[
	Partition[##]&,
	PropagatorChoiMatrix[t,H],
	{{2,2},{2,2},{2,2}}
]
PropagatorChoiTensor[t,HoldPattern[Hamiltonian[\[Omega]_,\[CapitalDelta]\[Omega]_,DiagonalMatrix[{JX_,JY_,JZ_}]]]]=PropagatorChoiTensor[t,Hamiltonian[\[Omega],\[CapitalDelta]\[Omega],DiagonalMatrix[{JX,JY,JZ}]]];


(* ::Subsection:: *)
(*Process tensor and SDP*)


(* ::Text:: *)
(*With that the link product can be written as a tensor contraction. To contract tensors in a readable way, we use the EinsteinSummation resource function. Now we can start gradually building up to the SDP. We start by combining together the initial state and the process tensor into ProcessMap. It takes as input the initial system-actuator state, which is assumed to be provided as a tensor product, a list of timings for individual unitaries, and the Hamiltonian.*)


SetAttributes[ProcessMap,HoldAll];
ProcessMap[rhoIn_,t_List,H_]:=EinsteinSummation[
	{
	(*The initial state is contracted with the tensor along the first four indices*)
	#[[1,1;;4]],
	Sequence@@#
	}&@Array[
	(*Constructs a list of all indices and then rule replacements are used to indicate contraction*)
		{
		(*Each tensor has eight indices*)
		Subscript[a,i,2#-1],
		Subscript[a,j,2#-1],
		Subscript[b,i,2#-1],
		Subscript[b,j,2#-1],
		Subscript[a,i,2#],
		Subscript[a,j,2#],
		Subscript[b,i,2#],
		Subscript[b,j,2#]
		}&,Length[t]
	]/.{
		(*Indices a_{_,2l+1} are replaced with a_{_,2l} to indicate the contraction*)
		Subscript[a,i_,k_]:>Subscript[a, i,k-1]/;k>1&&OddQ[k],
		(*After the last system-actuator unitary the actuator is traced out*)
		Subscript[b,j,2Length[t]]->Subscript[b, i,2Length[t]]},
	{
	rhoIn,
	Sequence@@(PropagatorChoiTensor[#,H]&/@t)
	}
]


(* ::Text:: *)
(*ProcessComb takes the ProcessCombMap and contracts it with the target state to compute the  overlap upon contraction with a tester tensor.*)


SetAttributes[ProcessComb,HoldAll];
ProcessComb[rhoIn_,rhoTgt_,t_List,H_]:=EinsteinSummation[
	{
	(Join@@Table[{Subscript[b, i,l],Subscript[b, j,l]},{l,2,2Length[t]-1}])~Join~{ai,aj},
	{aj,ai}
	},
	{
	ProcessMap[rhoIn,t,H],
	rhoTgt
	}
]


(* ::Text:: *)
(*ProcessMatrix takes ProcessComb and turns it into a matrix.*)


SetAttributes[ProcessCombMatrix,HoldAll];
ProcessCombMatrix[rhoIn_,rhoFin_,t_List,H_]:=
	Nest[ArrayFlatten,#,(Dimensions[#]//Length)/2]&[ProcessComb[rhoIn,rhoFin,t,H]];


(* ::Text:: *)
(*Next we need formulate the constraints. We will use Indexed to work with directly with the components of the deterministic superinstrument . For partial trace we use the resource function MatrixPartialTrace. The constraints are constructed and simplified once and cached as they only depend on the number of time steps.*)


CausalityConstraints[timeSteps_Integer]:=CausalityConstraints[timeSteps]=Module[
	{Tvar=Array[Indexed[T,{##}]&,{4,4}^(timeSteps-1)]},
	Table[
		ParallelMap[
			FullSimplify,
			KroneckerProduct[
				MatrixPartialTrace[Tvar,-2i;;-1,2],
				1/(2^(2i)) IdentityMatrix[2^(2i)]
			]
			-KroneckerProduct[
				MatrixPartialTrace[Tvar,-(2i-1);;-1,2],
				1/(2^(2i-1)) IdentityMatrix[2^(2i-1)]
			],
			{2},Method->"FinestGrained"
		]==0,
		{i,timeSteps-1}
	]
]


(* ::Text:: *)
(*Now we are ready to formulate the SDP. Mathematica has a built-in for semidefinite programs called SemidefiniteOptimization which by default uses the CSDP solver. The problem can be specified using either, scalar, vector, or matrix variables but it has to specified as a minimisation problem. We can easily turn our maximisation problem into a minimisation problem by minimising the negative of the objective function. By default the method only returns the solution but we can prompt it to return both the minimum value of the objective function and the solution. The option inheritance is added for additional convenience when experimenting with different settings. The output of SDPResult is an association that stores the maximal population in the target state under the key "Population", the optimal deterministic superinstrument under "Tester", and the final state of the system under "Final State".*)


Options[SDPResult]=Options[SemidefiniteOptimization];
SDPResult[rhoIn_,rhoFin_,t_?(VectorQ[#,NumericQ]&),H_?(MatrixQ[#,NumericQ]&),opt:OptionsPattern[]]:=Module[
	{
		timeSteps=Length[t],
		dim,
		Tvar,
		result,
		TesterTensor
	},
	dim=4^(timeSteps-1);
	Tvar=Array[Indexed[T,{##}]&,{dim,dim}];
	result=AssociationThread[
		{"Population","Tester"},
		SemidefiniteOptimization[
			-Re@Tr[Transpose[ProcessCombMatrix[rhoIn,rhoFin,t,H]] . T],
			{
				VectorGreaterEqual[{T,0},{"SemidefiniteCone",dim}],
				Tr[Tvar]==2^(timeSteps-1)
			}~Join~CausalityConstraints[timeSteps],
			T\[Element]Matrices[{dim,dim},Complexes,Hermitian[{1,2}]],
			{"PrimalMinimumValue","PrimalMinimizer"},
			Evaluate@FilterRules[{opt}~Join~Options[SDPResult],Options[SemidefiniteOptimization]]
		]
	];
	result["Population"]=-result["Population"];
	result["Tester"]=result["Tester"][[1]];
	TesterTensor=Nest[Partition[#,{2,2}]&,result["Tester"],2(timeSteps-1)-1];
	AssociateTo[result,
	"Final State"->Chop@EinsteinSummation[
		{
		#~Join~{ai,aj},
		#
		}&[Join@@Table[{Subscript[b, i,l],Subscript[b, j,l]},{l,2,2Length[t]-1}]],
		{
		ProcessCombMap[rhoIn,t,H],
		TesterTensor
		}]
	]
];


(* ::Subsection:: *)
(*Unitary Strategy*)


(* ::Text:: *)
(*We use the following polynomial parameterization of unitaries. The normalisation of the unitary is slightly redundant since we will explicitly enforce the optimisation to be over unit vectors but it makes sure the matrix is always unitary even if the numerical algorithm slightly violates the constraints.*)


Unitary[w_,x_,y_,z_,normalised_:True]:=If[normalised,1/Norm[{w,x,y,z}],1](
	w IdentityMatrix[2]
	-I x PauliMatrix[1]
	-I y PauliMatrix[2]
	-I z PauliMatrix[3]
);
UnitaryChoi[w_,x_,y_,z_,normalised_:True]:=TensorProduct[
	Flatten[Unitary[w,x,y,z,normalised]//Transpose],
	Flatten[Unitary[w,-x,y,-z,normalised]//Transpose]
]//Simplify;


(* ::Text:: *)
(*We define the objective function, which is not necessary but cleans up the code. To optimise over uncorrelated unitaries we can reuse the ProcessMatrix to construct the objective function by contracting it with the tensor product of the appropriate number of Choi matrices of the actuator unitaries.*)


SetAttributes[UnitaryObjectiveFunction,HoldAll];
UnitaryObjectiveFunction[rhoIn_,rhoTgt_,t_List,H_]:=With[
	{vars=ToExpression[Map[StringJoin["U",ToString[#]]&,Range[Length[t]-1]]]},
	{varsIndexed=Function[s,Array[Indexed[s,#]&,4]]/@vars},
	Tr[Dot[
		Transpose[ProcessCombMatrix[rhoIn,rhoTgt,t,H]],
		KroneckerProduct@@UnitaryChoi@@@varsIndexed
		]
	]//Re//ComplexExpand//Chop//Simplify[#,Assumptions->Element[Alternatives@@Flatten@vars,Reals]]&
]


(* ::Text:: *)
(*Similar to SDPResult the result of UnitaryStrategy is an association where the maximal population in the target state is stored under the key "Population" and the individual unitaries are stored under keys "U1", "U2", etc.*)


SetAttributes[UnitaryStrategy,HoldAll];
Options[UnitaryStrategy]=Options[NMaximize];
UnitaryStrategy[rhoIn_,rhoTgt_,t_List,H_,opt:OptionsPattern[]]:=Module[{max,sol},
	{max,sol}=With[{vars=ToExpression[Map[StringJoin["U",ToString[#]]&,Range[Length[t]-1]]]},
		NMaximize[
			UnitaryObjectiveFunction[rhoIn,rhoTgt,t,H],
			Element[Alternatives@@vars,Sphere[4]],
			Evaluate@FilterRules[{opt}~Join~Options[UnitaryStrategy],Options[NMaximize]]
		]
	];
	Association[
		"Population"->max,
		KeyValueMap[ToString[#1]->Unitary@@#2&,Association[sol]]
	]
]
