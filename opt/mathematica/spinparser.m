(* ::Package:: *)

(* ::Title:: *)
(*SpinParser Mathematica package*)


BeginPackage["SpinParser`"];


(* ::Chapter:: *)
(*Interface*)


(* ::Subsection:: *)
(*Symbols*)


Protect[Cutoff,Site,Reference,Verbose];
Protect[Component,XX,XY,XZ,YX,YZ,YZ,ZX,ZY,ZZ];


(* ::Subsection:: *)
(*Functions*)


GetLatticeBasis::usage="GetLatticeBasis[obsfile] 
Return the lattice basis for the lattice graph used in the calculation of the specified obsfile.

Args:
	obsfile (String): Path to the obsfile. 

Returns:
	List: List of coordinates of lattice basis sites. For an N-site basis, the return List has shape (N,3), with the second dimension storing the x, y, and z components of the position vector.";

GetLatticePrimitives::usage="GetLatticePrimitives[obsfile]
	Return the primitive lattice vectors for the lattice graph used in the calculation of the specified obsfile. 

	Args:
		obsfile (String): Path to the obsfile. 

	Returns:
		List: List of primitive lattice vectors. The return List has shape (3,3), with the first dimension enumerating the primitives, and the second dimension storing the x, y, and z components of the primitive.";

GetLatticeSites::usage="GetLatticeSites[obsFile]
	Get a list of all lattice sites within truncation range of a specified reference site. 

	Args:
		obsfile (String): Path to the obsfile.

	Options:
		Reference (Integer|List): Reference site. Can be either an integer `n` which refers to the n-th basis site, or a three-component list which specifies the real-space position of a basis site. Defaults to 0.
		Verbose (True|False): Defines whether the output should be verbose. Defaults to True.

	Returns:
		List|Association: If `verbose==False`, return a List of shape (N,3) where N is the number of sites within truncation range, and the second dimension stores the x, y, and z components of the lattice site position. 
			If the output is verbose, return an Association with keys `data`, containing the nonverbose data, and `reference`, which contains the real-space position of the specified reference site. ";

GetCorrelation::usage="GetCorrelation[obsfile]
	Obtain correlation data from an observables file.

	Args:
		obsfile (String): Path to the obsfile.

	Options:
		Cutoff (Real|List|All): Cutoff values at which to retrieve data. Can be a single cutoff value, a list of cutoff values, or All. Defaults to All.
		Site (List|All): Lattice sites at which to retrieve data. Can be a one-dimensional List of x,y,z real-space coordinates of a site, a two-dimensional N*3 List of multiple sites, or All. Defaults to All.
		Reference (Integer|List): Reference site for the correlation measurement. Can be either an integer `n` which refers to the n-th basis site, or a three-component list which specifies the real-space position of a basis site. Defaults to 0.
		Component (Symbol|All): Spin component of the correlation function to retrieve. Can be XX, XY, XZ, YX, ..., ZZ, or All. The latter is equivalent to XX+YY+ZZ. Defaults to All. 
		Verbose (True|False): Defines whether the output should be verbose. Defaults to True.

	Returns:
		List|Association: If `verbose==False`, return a List of shape (N,M) where N is the number of cutoff values selected and N is the number of sites selected.
			If the output is verbose, return an Association with keys `data` (contains the nonverbose data), `cutoff` (contains the selected cutoff values), `site` (contains the selected lattice sites), and `reference` (contains the real-space position of the specified reference site).";

GetStructureFactor::usage="GetStructureFactor[obsfile,momentum]
	Calculate the spin structure factor from the spin correlations stored in an observable file. 

	Args:
		obsfile (String): Path to the obsfile.
		momentum (List): Momentum points at which to compute the structure factor. Can be a one-dimensional List of three kx,ky,kz coordinates, or a two-dimensional List of multiple momentum points.

	Options:
		Cutoff (Real|List|All): Cutoff values at which to retrieve data. Can be a single cutoff value, a list of cutoff values, or All. Defaults to All.
		Component (Symbol|All): Spin component of the structure factor to retrieve. Can be XX, XY, XZ, YX, ..., ZZ, or All. The latter is equivalent to XX+YY+ZZ. Defaults to All. 
		Verbose (True|False): Defines whether the output should be verbose. Defaults to True.

	Returns:
		List|Association: If `verbose==False`, return a List of shape (N,M) where N is the number of cutoff values selected and N is the number of momentum points specified.
			If the output is verbose, return an Association with keys `data` (contains the nonverbose data), `cutoff` (contains the selected cutoff values), and `momentum` (momentum points).";
			
PlotCorrelationFlow::usage="PlotCorrelationFlow[obsFile]
	Plot the cutoff dependence of spin-spin correlations.

	Args:
		obsfile (String): Path to the obsfile.

	Options:
		Supports all plotting directives of ListLinePlot. Further options are listed below. 
		Site (List|All): Lattice sites for which to plot correlations. Can be a one-dimensional List of x,y,z real-space coordinates of a site, a two-dimensional N*3 List of multiple sites, or All. Defaults to All.
		Reference (Integer|List): Reference site for the spin-spincorrelation. Can be either an integer `n` which refers to the n-th basis site, or a three-component list which specifies the real-space position of a basis site. Defaults to 0.
		Component (Symbol|All): Spin component of the correlation function to retrieve. Can be XX, XY, XZ, YX, ..., ZZ, or All. The latter is equivalent to XX+YY+ZZ. Defaults to All.";
		
PlotStructureFactorFlow::usage="PlotStructureFactorFlow[obsFile,momentum]
	Plot the cutoff dependence of the momentum space spin structure factor.

	Args:
		obsfile (String): Path to the obsfile.
		momentum (List): Momentum points at which to plot the structure factor. Can be a one-dimensional List of three kx,ky,kz coordinates, or a two-dimensional List of multiple momentum points.

	Options:
		Supports all plotting directives of ListLinePlot. Further options are listed below. 
		Component (Symbol|All): Spin component of the structure factor to retrieve. Can be XX, XY, XZ, YX, ..., ZZ, or All. The latter is equivalent to XX+YY+ZZ. Defaults to All.";


(* ::Chapter:: *)
(*Implementation*)


(* ::Subchapter:: *)
(*Header*)


Begin["`Private`"];


(* ::Subchapter:: *)
(*Internal functions*)


IsNaN[x_]:=Round[x]==-9223372036854775808 (*workaround to identify numerical NaN values in Mathematica*)


(* ::Subchapter:: *)
(*External functions*)


(* ::Section:: *)
(*GetLatticeBasis*)


(* ::Input::Initialization:: *)
GetLatticeBasis[obsfile_]:=Module[{cmd,raw,processed},
cmd=ToString[StringForm["import spinparser.obs as o; o.getLatticeBasis(r'``')",obsfile]];
raw=ExternalEvaluate["Python",cmd];
processed=Normal[raw];
Return[processed];
]


(* ::Section:: *)
(*GetLatticePrimitives*)


(* ::Input::Initialization:: *)
GetLatticePrimitives[obsfile_]:=Module[{cmd,raw,processed},
cmd=ToString[StringForm["import spinparser.obs as o; o.getLatticePrimitives(r'``')",obsfile]];
raw=ExternalEvaluate["Python",cmd];
processed=Normal[raw];
Return[processed];
]


(* ::Section:: *)
(*GetLatticeSites*)


(* ::Input::Initialization:: *)
Options[GetLatticeSites]={Reference->0,Verbose->True};
GetLatticeSites[obsfile_,opts:OptionsPattern[]]:=Module[{argReference,argVerbose,cmd,raw,processed},
argReference=ExportString[OptionValue[Reference],"PythonExpression"];
argVerbose=If[OptionValue[Verbose]==True,"True","False"];
cmd=ToString[StringForm["import spinparser.obs as o; o.getLatticeSites(r'``', reference=``, verbose=``)",obsfile,argReference,argVerbose]];
raw=ExternalEvaluate["Python",cmd];
processed=If[raw[[0]]===Association,Normal/@raw,Normal@raw];
Return[processed];
]


(* ::Section:: *)
(*GetCorrelation*)


(* ::Input::Initialization:: *)
Options[GetCorrelation]={Cutoff->All, Site->All, Reference->0,Component->All,Verbose->True};
GetCorrelation[obsfile_,opts:OptionsPattern[]]:=Module[{argCutoff,argSite,argReference,argComponent,argVerbose,cmd,raw,processed},
argCutoff=If[OptionValue[Cutoff]===All,"\"all\"",ExportString[OptionValue[Cutoff],"PythonExpression"]];
argSite=If[OptionValue[Site]===All,"\"all\"",ExportString[OptionValue[Site],"PythonExpression"]];
argReference=ExportString[OptionValue[Reference],"PythonExpression"];
argComponent=Switch[OptionValue[Component],All,"\"all\"",XX,"\"XX\"",XY,"\"XY\"",XZ,"\"XZ\"",YX,"\"YX\"",YY,"\"YY\"",YZ,"\"YZ\"",ZX,"\"ZX\"",ZY,"\"ZY\"",ZZ,"\"ZZ\""];
argVerbose=If[OptionValue[Verbose]===True,"True","False"];
cmd=ToString[StringForm["import spinparser.obs as o; o.getCorrelation(r'``', cutoff=``, site=``, reference=``, component=``, verbose=``)",obsfile,argCutoff,argSite,argReference,argComponent,argVerbose]];
raw=ExternalEvaluate["Python",cmd];
processed=If[raw[[0]]===Association,Normal/@raw,Normal@raw];
Return[processed];
]


(* ::Section:: *)
(*GetStructureFactor*)


(* ::Input::Initialization:: *)
Options[GetStructureFactor]={Cutoff->All,Component->All,Verbose->True};
GetStructureFactor[obsfile_,momentum_,opts:OptionsPattern[]]:=Module[{argMomentum,argCutoff,argComponent,argVerbose,cmd,raw,processed},
argMomentum=ExportString[momentum,"PythonExpression"];
argCutoff=If[OptionValue[Cutoff]===All,"\"all\"",ExportString[OptionValue[Cutoff],"PythonExpression"]];
argComponent=Switch[OptionValue[Component],All,"\"all\"",XX,"\"XX\"",XY,"\"XY\"",XZ,"\"XZ\"",YX,"\"YX\"",YY,"\"YY\"",YZ,"\"YZ\"",ZX,"\"ZX\"",ZY,"\"ZY\"",ZZ,"\"ZZ\""];
argVerbose=If[OptionValue[Verbose]===True,"True","False"];
cmd=ToString[StringForm["import spinparser.obs as o; o.getStructureFactor(r'``', ``, cutoff=``, component=``, verbose=``)",obsfile,argMomentum,argCutoff,argComponent,argVerbose]];
raw=ExternalEvaluate["Python",cmd];
processed=If[raw[[0]]===Association,Normal/@raw,Normal@raw];
Return[processed];
]


(* ::Section:: *)
(*PlotCorrelationFlow*)


Options[PlotCorrelationFlow]=Join[{Site->All,Reference->0,Component->All},Options[ListLinePlot]];
PlotCorrelationFlow[obsfile_,opts:OptionsPattern[]]:=Module[{raw,data,labels,plot},
(*fetch data*)
raw=GetCorrelation[obsfile,Cutoff->All,Verbose->True,Evaluate[FilterRules[{opts},{Site,Reference,Component}]]];
data=Table[{raw["cutoff"][[c]],raw["data"][[c,s]]},{s,1,Length[raw["site"]]},{c,1,Length[raw["cutoff"]]}];
data=Select[data,(Not[Or@@IsNaN/@Flatten[#]])&];
(*generate labels*)
labels=Table[ToString[s],{s,raw["site"]}];
(*generate plot*)
plot=ListLinePlot[data,Evaluate[FilterRules[{opts},Options[ListLinePlot]]],PlotLegends->labels];
Return[plot];
]


(* ::Section:: *)
(*PlotStructureFactorFlow*)


Options[PlotStructureFactorFlow]=Join[{Component->All},Options[ListLinePlot]];
PlotStructureFactorFlow[obsfile_,momentum_,opts:OptionsPattern[]]:=Module[{raw,data,labels,plot},
(*fetch data*)
raw=GetStructureFactor[obsfile,momentum,Cutoff->All,Verbose->True,Evaluate[FilterRules[{opts},{Component}]]];
data=Table[{raw["cutoff"][[c]],raw["data"][[c,m]]},{m,1,Length[raw["momentum"]]},{c,1,Length[raw["cutoff"]]}];
data=Select[data,(Not[Or@@IsNaN/@Flatten[#]])&];
(*generate labels*)
labels=Table[ToString[m],{m,raw["momentum"]}];
(*generate plot*)
plot=ListLinePlot[data,Evaluate[FilterRules[{opts},Options[ListLinePlot]]],PlotLegends->labels];
Return[plot];
]


(* ::Chapter:: *)
(*Scratch*)


(* ::Subchapter:: *)
(*Footer*)


End[];
EndPackage[];
