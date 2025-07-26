```mathematica
(* ::Package:: *)

(* ::Title:: *)
(*BRUE Path Cost Region Plot (Mathematica Version)*)


(* ::Text:: *)
(*This script is a Mathematica implementation of the Python script 'plot_region_path_cost.py'.*)
(*It performs the following functions:*)
(*1. Defines link/path cost functions.*)
(*2. Constructs three types of feasible regions (reg, reg2, reg3) and an equilibrium region (regEqm).*)
(*3. Calculates the min/max and midpoints of each path cost over the feasible region.*)
(*4. Samples and plots the regions.*)
(*5. Exports the plot to a PDF file.*)


(* ::Section:: *)
(*1. Basic Constants*)


rho = 15;
sigma = 0.02;
e = 15.0; (* Error upper bound *)
eps = 10^-8; (* Buffer for strict inequalities *)

(* Boundary conditions *)
f1Min = 1000.0;
f1Max = 8000.0;
f2Min = 0.0;
f2Max = 7000.0;
totalMax = 10000.0;


(* ::Section:: *)
(*2. Link/Path Cost Functions*)


link1[f1_] := 18.0 * (1 + 0.15 * (f1 / 3600.0)^4);
link2[f2_] := 22.5 * (1 + 0.15 * (f2 / 3600.0)^4);
sharedDenom[f1_, f2_] := (10000.0 - f1 - f2) / 1800.0;
link3[f1_, f2_] := 12.0 * (1 + 0.15 * sharedDenom[f1, f2]^4);
link5[f1_, f2_] := 2.4 * (1 + 0.15 * sharedDenom[f1, f2]^4);
link8[f1_, f2_] := 12.0 * (1 + 0.15 * sharedDenom[f1, f2]^4);

(* Perceived cost function *)
rgCost[rawTime_] := rawTime + rho * (1 - Exp[-sigma * rawTime]);

(* Path cost functions *)
path1[f1_, f2_] := rgCost[link1[f1]];
path2[f1_, f2_] := rgCost[link2[f2]];
path5[f1_, f2_] := rgCost[link3[f1, f2] + link5[f1, f2] + link8[f1, f2]];

maxPath[f1_, f2_] := Max[path1[f1, f2], path2[f1, f2], path5[f1, f2]];
minPath[f1_, f2_] := Min[path1[f1, f2], path2[f1, f2], path5[f1, f2]];


(* ::Section:: *)
(*3. Constraints and Region Definitions*)


(* Box constraints *)
boxConstraints[f1_, f2_] := And[
    f1Min <= f1 <= f1Max,
    f2Min <= f2 <= f2Max,
    f1 + f2 <= totalMax
];

(* Order constraints for BS_0^zeta *)
orderConstraints[f1_, f2_] := And[
    path1[f1, f2] - path2[f1, f2] + eps <= 0,
    path1[f1, f2] - path5[f1, f2] + eps <= 0,
    path2[f1, f2] - path5[f1, f2] + eps <= 0
];

(* Pairwise constraints for BS_0^zeta *)
pairConstraints[f1_, f2_] := And[
    Abs[path1[f1, f2] - path2[f1, f2]] <= e,
    Abs[path1[f1, f2] - path5[f1, f2]] <= e,
    Abs[path2[f1, f2] - path5[f1, f2]] <= e
];

(* Feasible region for BS_0^zeta (corresponds to 'mask_reg' in Python) *)
reg[f1_, f2_] := And[
    boxConstraints[f1, f2],
    orderConstraints[f1, f2],
    pairConstraints[f1, f2]
];

(* Feasible region for S_0^zeta (corresponds to 'mask_reg2' in Python) *)
reg2[f1_, f2_] := And[
    boxConstraints[f1, f2],
    (maxPath[f1, f2] - minPath[f1, f2]) <= e
];


(* ::Section:: *)
(*4. Calculate Min/Max/Midpoints*)


(* Find min/max values of path costs within the feasible region 'reg' *)
{r1Min, r1Max} = {NMinValue[{path1[f1, f2], reg[f1, f2]}, {f1, f2}], NMaxValue[{path1[f1, f2], reg[f1, f2]}, {f1, f2}]};
{r2Min, r2Max} = {NMinValue[{path2[f1, f2], reg[f1, f2]}, {f1, f2}], NMaxValue[{path2[f1, f2], reg[f1, f2]}, {f1, f2}]};
{r5Min, r5Max} = {NMinValue[{path5[f1, f2], reg[f1, f2]}, {f1, f2}], NMaxValue[{path5[f1, f2], reg[f1, f2]}, {f1, f2}]};

(* Calculate midpoints *)
midPath1 = (r1Min + r1Max) / 2;
midPath2 = (r2Min + r2Max) / 2;
midPath5 = (r5Min + r5Max) / 2;

(* Feasible region for RS_0^zeta (corresponds to 'mask_reg3' in Python) *)
reg3[f1_, f2_] := And[
    reg[f1, f2],
    path1[f1, f2] <= midPath1,
    path2[f1, f2] <= midPath2,
    path5[f1, f2] <= midPath5
];


(* ::Section:: *)
(*5. Calculate Equilibrium Line*)


(* This function replicates the logic from 'calculate_equilibrium_line' in Python *)
calculateEquilibriumLine[leftBoundaryX_, upperLimitX_, upperLimitY_] := Module[
    {moneyMax, moneyMin, moneyRange, weights},
    moneyMax = Max[upperLimitY];
    moneyMin = Min[upperLimitY];
    moneyRange = moneyMax - moneyMin;
    
    weights = If[moneyRange > 0,
        0.7 + 0.1 * ((upperLimitY - moneyMin) / moneyRange),
        ConstantArray[0.5, Length[upperLimitY]]
    ];
    
    leftBoundaryX + (upperLimitX - leftBoundaryX) * weights
];

upperLimitX = {midPath1, midPath2, midPath5};
upperLimitY = {20.0, 15.0, 2.0}; (* Monetary costs *)
leftBoundaryX = {r1Min, r2Min, r5Min};

eqmLimitX = calculateEquilibriumLine[leftBoundaryX, upperLimitX, upperLimitY];

Print["Left Boundary X: ", leftBoundaryX];
Print["Upper Limit X: ", upperLimitX];
Print["Upper Limit Y (Monetary Cost): ", upperLimitY];
Print["Equilibrium Line X: ", eqmLimitX];

(* Feasible region for T_eqm (corresponds to 'mask_eqm' in Python) *)
regEqm[f1_, f2_] := And[
    reg[f1, f2],
    path1[f1, f2] <= eqmLimitX[[1]],
    path2[f1, f2] <= eqmLimitX[[2]],
    path5[f1, f2] <= eqmLimitX[[3]]
];


(* ::Section:: *)
(*6. Plotting*)


plot = Show[
    RegionPlot[
        reg2[f1, f2], {f1, f1Min, f1Max}, {f2, f2Min, f2Max},
        PlotStyle -> Directive[Blue, Opacity[0.3]], BoundaryStyle -> None,
        PlotPoints -> 100
    ],
    RegionPlot[
        reg[f1, f2], {f1, f1Min, f1Max}, {f2, f2Min, f2Max},
        PlotStyle -> Directive[Red, Opacity[0.3]], BoundaryStyle -> None,
        PlotPoints -> 100
    ],
    RegionPlot[
        reg3[f1, f2], {f1, f1Min, f1Max}, {f2, f2Min, f2Max},
        PlotStyle -> Directive[Green, Opacity[0.3]], BoundaryStyle -> None,
        PlotPoints -> 100
    ],
    RegionPlot[
        regEqm[f1, f2], {f1, f1Min, f1Max}, {f2, f2Min, f2Max},
        PlotStyle -> Directive[Purple, Opacity[0.5]], BoundaryStyle -> None,
        PlotPoints -> 100
    ],
    PlotRange -> {{f1Min, f1Max}, {f2Min, f2Max}},
    AxesLabel -> {"f_1", "f_2"},
    BaseStyle -> {FontFamily -> "Computer Modern"},
    GridLines -> Automatic,
    PlotLegends -> Placed[
        LineLegend[
            {Directive[Blue, Opacity[0.3]], Directive[Red, Opacity[0.3]], Directive[Green, Opacity[0.3]], Directive[Purple, Opacity[0.5]]},
            {"Region of S_0^\[Zeta]", "Region of BS_0^\[Zeta]", "Region of RS_0^\[Zeta]", "Region of T_eqm"}
        ],
        {Right, Top}
    ]
];

(* Display and Export the plot *)
Print[plot];
Export[FileNameJoin[{DirectoryName[$InputFileName], "ee_time_point_mathematica.pdf"}], plot];
Print["\[check] Plotting complete, saved to ee_time_point_mathematica.pdf"];

```
