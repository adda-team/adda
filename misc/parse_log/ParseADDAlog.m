(* ::Package:: *)

(* :Title: Parsing ADDA logs *)

(* :Author: Maxim Yurkin (yurkin@gmail.com) *)

(* :Summary: This package performs parsing of ADDA log files *)

(* :Context: ParseADDAlog` *)

(* :Package Version: 1.1 *)

(* :Mathematica Version: 7.0 *)

BeginPackage["ParseADDAlog`"];

ParseLog::usage="ParseLog[\"fname\"] parses ADDA log file, extracting the values of all recognized entries. A full list of names of supported entries is available via EntryNames. The main part of parsing is based on table LogFormat inside the private package section.";
EntryNames::usage="EntryNames[] returns a full list of entry names, scanned by ParseLog.";
EntryIndex::usage="EntryIndex[\"name(s)\"] returns index of the entry, which can be furhter used to index data block returned by ParseLog. If a list of names is given, the function returns a list of indices.";
EntryValue::usage="EntryValue[data,\"name(s)\"] provides a direct method to extract entry values from data returned by ParseLog. If a list of names is given, the function returns a list of values.";

Begin["`Private`"];

(* matches number in usual or engineering format *)
NumberPattern={NumberString~~{"e","E"}~~NumberString,NumberString};
(* converts string in usual or engineering format into a string *)
ToNumber[str_]:=Module[{spl},
	spl=StringSplit[str,{"e","E"}];
	If[Length[spl]==1,
		ToExpression[spl[[1]]],
		ToExpression[spl[[1]]]10.^ToExpression[spl[[2]]]
	]
]
(* builds a complex number from strings for real and imaginary parts and sign *)
ToComplex[x_,y_,s_]:=If[s=="+",ToNumber[x]+ToNumber[y]I,ToNumber[x]-ToNumber[y]I];
SetAttributes[ToNumber,Listable];
SetAttributes[ToComplex,Listable];
(* matches number in usual or engineering format, and produces number instead of a string *)
GenNumberString=x:NumberPattern:>ToNumber[x];
(* matchs complex string and convert it to numeric format *)
Off[RuleDelayed::rhs];
ComplexPattern[x_,y_,s_]:=x:NumberPattern~~s:{"+","-"}~~y:NumberPattern~~"i";
On[RuleDelayed::rhs];
ComplexString=ComplexPattern[x,y,s]:>ToComplex[x,y,s];
(* matchs string "..x..x.." and transforms to numeric array of size 3 *)
GridString=x:NumberPattern~~"x"~~y:NumberPattern~~"x"~~z:NumberPattern:>ToNumber[{x,y,z}];
(* matchs vector writtern as "(..,..,..)" and transforms to numeric array of size 3 *)
VectorString="("~~x:NumberPattern~~","~~y:NumberPattern~~","~~z:NumberPattern~~")":>ToNumber[{x,y,z}];
(* matchs complex vector writtern as "(..,..,..)" and transforms to numeric array of size 3 *)
VectorComplexString=("("~~ComplexPattern[x1,y1,s1]~~","~~ComplexPattern[x2,y2,s2]~~","~~ComplexPattern[x3,y3,s3]~~")":>ToComplex[{x1,x2,x3},{y1,y2,y3},{s1,s2,s3}]);

(* Each row contains 1) name,  2) starting (True) or contained (False) 3) pre-string, and 4) pattern (or rule) to match. Due to performance
   considerations they should be in the same order as present in the file, but elements which are not present in the file will be skipped automatically *)
LogFormat={ 
	{"version",True,"Generated by ADDA v.",___},
	(* computer name is skipped in seq runs *)
	{"Ncores",True,"The program was run on: ",{x:NumberString~~" processors":>ToExpression[x],___->1}},
	{"command",True,"command: ","'"~~s___~~"'"->s},
	{"lambda",True,"lambda: ",GenNumberString},
	(* more specific parameters can be extracted from shape, possibly multi-line *)
	{"shape",True,"shape: ",___}, 
	(* granule information is skipped *)
	{"dipGrid",True,"box dimensions: ",GridString},
	(* multi-line when more than 1 value, also will not work for anisotropic ones *)
	{"refIndex",True,"refractive index: ",ComplexString},
	(* perfectly reflecting substrate is skipped *)
	{"subRefIndex",True,"Particle is placed near the substrate with refractive index ",ComplexString},
	{"subHeight",True,"  height of the particle center: ",GenNumberString},
	(* the ratios (given for non-cubical dipoles) are currently ignored *)
	{"dipSize",True,"Dipole size: ",{GridString,GenNumberString}},
	{"dpl",True,"Dipoles/lambda: ",GenNumberString},
	{"reqEps",True,"Required relative residual norm: ",GenNumberString},
	{"Ndipoles",True,"Total number of occupied dipoles: ",GenNumberString},
	(* values of occupied dipoles per domain are skipped *)
	{"sizeParameter",True,"Volume-equivalent size parameter: ",GenNumberString},
	{"incBeamType",True,"Incident beam: ",___},
	{"incCenter",True,"Incident beam center position: ",VectorString},
	{"incProp",True,"Incident propagation vector: ",VectorString},
	{"incPolY",True,"Incident polarization Y(par): ",VectorString},
	{"incPolX",True,"Incident polarization X(per): ",VectorString},
	{"reflProp",True,"Reflected propagation vector: ",VectorString},
	{"transProp",True,"Transmitted propagation vector: ",VectorString},
	{"orient",True,"Particle orientation: ",___},
	(* orient-modified propagation and polarization vectors are skipped *)
	(* which scattering matrices are calculated is skipped *)
	{"polRel",True,"Polarizability relation: ",___},
	{"scatQuan",True,"Scattering quantities formulae: ",___},
	{"intTerm",True,"Interaction term prescription: ",___},
	{"reflTerm",True,"Reflected Green's tensor formulae: ",___},
	{"FFT",True,"FFT algorithm: ",___},
	{"oclFFT",True,"OpenCL FFT algorithm: ",___},
	{"iterMethod",True,"Iterative Method: ",___},
	{"symmetry",True,"Symmetries: ",___},
	{"optimization",True,"Optimization is done for ",___},
	(* checkpoint information is skipped *)
	
	(* OpenCL, FFT and memory section *)
	{"OCLdevice",True,"Using OpenCL device ",___},
	{"OCLdeviceMemory",True,"Device memory: total - ",x:NumberString~~" MB, maximum object - "~~y:NumberString:>ToExpression[{x,y}]},
	{"fftGrid",True,"The FFT grid is: ",GridString},
	{"memoryFFTeach",True,"Memory usage for MatVec matrices",GenNumberString},
	{"memoryTotal",True,"Total memory usage: ",GenNumberString},
	{"memoryMaxEach",True,"Maximum memory usage of single processor: ",GenNumberString},
	{"memoryOCL",True,"OpenCL memory usage: peak total - ",x:NumberString~~" MB, maximum object - "~~y:NumberString:>ToExpression[{x,y}]},
	(* maximum-object OpenCL memory is skipped *)

	(* Only the two starting parameters are extracted from the iterative solution, the iteration progress itself is not yet supported. During 
	   orientation averaging, this section is not present at all - instead various orientations are listed. Parsing of the latter is also not yet supported *)
	{"polarizability_Y",True,"CoupleConstant: ",{VectorComplexString,ComplexString}},
	{"x0_Y",True,"x_0 = ",___},
	{"polarizability_X",True,"CoupleConstant: ",{VectorComplexString,ComplexString}},
	{"x0_X",True,"x_0 = ",___},
	
	(* timing section; all optional communication time are skipped, since they cannot be reliably parsed without much extra effort *)
	{"NsingleParticle",True,"Total number of single particle evaluations: ",GenNumberString},
	{"Niter",True,"Total number of iterations: ",GenNumberString},
	{"NMVP",True,"Total number of matrix-vector products: ",GenNumberString},
	{"NcalcEfield",True,"Total planes of E field calculation (each ",x:NumberString~~" points): "~~y:NumberString:>ToExpression[{x,y}]},
	{"timeWall",True,"Total wall time: ",GenNumberString},
	{"timeSinceInit",True,"Time since MPI_Init: ",GenNumberString},
	{"timeTotal",True,"Total time: ",GenNumberString},
	{"timeInit",False,"Initialization time: ",GenNumberString},
	{"timeInitOCL",False,"init OpenCL: ",GenNumberString},
	{"timeInitInteraction",False,"init interaction: ",GenNumberString},
	{"timeInitDmatrix",False,"init Dmatrix: ",GenNumberString},
	{"timeInitDmatrixComm",False,"communication: ",GenNumberString},
	{"timeInitFFT",False,"FFT setup: ",GenNumberString},
	{"timeInitParticle",False,"make particle: ",GenNumberString},
	{"timeGranule",False,"granule generator: ",GenNumberString},
	(* communication time for granule generator is skipped *)
	{"timeIF",False,"Internal fields: ",GenNumberString},
	{"timeIFsingleRun",False,"one solution: ",GenNumberString},
	{"timeIFsingleRunComm",False,"communication: ",GenNumberString},
	{"timeIFMVP",False,"matvec products: ",GenNumberString},
	{"timeIFMVPComm",False,"communication: ",GenNumberString},
	{"timeIFincBeam",False,"incident beam: ",GenNumberString},
	{"timeIFinitSolver",False,"init solver: ",GenNumberString},
	{"timeIFinitSolverComm",False,"communication: ",GenNumberString},
	{"timeIFoneIter",False,"one iteration: ",GenNumberString},
	{"timeIFoneIterComm",False,"communication: ",GenNumberString},
	{"timeIFoneIterMVP",False,"matvec products: ",GenNumberString},
	{"timeIFoneIterMVPComm",False,"communication: ",GenNumberString},
	{"timeSF",False,"Scattered fields: ",GenNumberString},
	{"timeSFonePlane",False,"one plane: ",GenNumberString},
	(* communication time for one plane is skipped *)
	{"timeSFallDir",False,"one alldir: ",GenNumberString},
	(* communication time for one alldir is skipped *)
	{"timeSFscatGrid",False,"one scat_grid: ",GenNumberString},
	(* communication time for one scat_grid is skipped *)
	{"timeOtherSQ",False,"Other sc.quantities: ",GenNumberString},
	{"timeOtherSQComm",False,"communication: ",GenNumberString},
	{"timeFileIO",True,"File I/O: ",GenNumberString},
	{"timeIntegration",True,"Integration: ",GenNumberString}
};
Nparams=Length[LogFormat];

ParseLog[fname_]:=Module[{FileLines,Nlines,vals,np,nl,nlLast,str,spl,case,found},
	FileLines=ReadList[fname,String];
	Nlines=Length[FileLines];
	vals=Table["None",{Nparams}];
	For[np=1;nl=1,np<=Nparams,np++,
		nlLast=nl;
		found=False;
		While[nl<=Nlines,
			(* following two cases can be combined in one, but this way it is easier to understand *)
			str=FileLines[[nl]];
			If[LogFormat[[np,2]], 
				(* check starting string *)
				spl=StringSplit[str,StartOfString ~~LogFormat[[np,3]]][[1]];
				nl++;
				If[StringLength[spl]<StringLength [str],
					found=True;
					case=StringCases[spl,LogFormat[[np,4]]];
					If[Length[case]>=1,vals[[np]]=case[[1]]];
					Break[];
				],
				(* check containing string *)
				spl=StringSplit[str,LogFormat[[np,3]]];
				nl++;
				If[Length[spl]>= 2,
					found=True;
					case=StringCases[spl[[2]],LogFormat[[np,4]]];
					If[Length[case]>=1,vals[[np]]=case[[1]]];
					Break[];
				]
			]
		];
		If[!found,nl=nlLast];
	];
	vals
];
EntryNames[]:=LogFormat[[All,1]];
EntryIndex[str_]:=Module[{tmp=Position[EntryNames[],str,1,1]},
	If[Length[tmp]==0,-1,tmp[[1,1]]]
];
SetAttributes[EntryIndex,Listable];
EntryValue[vals_,str_]:=vals[[EntryIndex[str]]];

End[];

EndPackage[];
