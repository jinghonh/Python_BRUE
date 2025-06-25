# Boundedly Rational User Equilibria (BRUE): Mathematical Formulation and Solution Sets

Xuan Di $^{a,\ast}$ , Henry X. Liu $^{a}$ , Jong- Shi Pang $^{b}$ , Xuegang (Jeff) Ban $^{c}$

$^{a}$ Department of Civil Engineering, University of Minnesota, Twin Cities   $^{b}$ Department of Industrial and Enterprise Systems Engineering, University of Illinois at Urbana- Champaign   $^{c}$ Department of Civil and Environmental Engineering, Rensselaer Polytechnic Institute

# Abstract

Boundedly rational user equilibria (BRUE) represent traffic flow distribution patterns where travellers can take any route whose travel cost is within an 'indifference band' of the shortest path cost. Those traffic flow patterns satisfying the above condition constitute a set, named the BRUE solution set. It is important to obtain all the BRUE flow patterns, because it can help predict the variation of the link flow pattern in a traffic network under the boundedly rational behavior assumption. However, the methodology of constructing the BRUE set has been lacking in the established literature. This paper fills the gap by constructing the BRUE solution set on traffic networks with fixed demands connecting multiple OD pairs. According to the definition of the  $\epsilon$ - BRUE, where  $\epsilon$  is the indifference band for the perceived travel cost, we formulate the  $\epsilon$ - BRUE problem as a nonlinear complementarity problem (NCP), so that a BRUE solution can be obtained by solving a BRUE- NCP formulation. To obtain the whole BRUE solution set encompassing all BRUE flow patterns, we firstly propose a methodology of generating various path combinations which may be utilized under the boundedly rational behavior assumption. We find out that with the increase of the indifference band, the path set that contains boundedly rational equilibrium flows will be augmented, and the critical values of indifference bands to augment these path sets can be identified by solving a family of mathematical programs with equilibrium constraints (MPEC) sequentially. After these utilized path sets are attained, the BRUE solution set can be obtained when we assign all traffic demands to these utilized paths. Various numerical examples are given to illustrate our findings.

$\mathfrak{S}$  2012 The Authors. Published by Elsevier Ltd. Selection and/or peer review under responsibility of Delft University of Technology

# Keywords:

Bounded Rationality User Equilibria (BRUE), Indifference band,  $\epsilon$ - BRUE

# 1. Introduction

In the static traffic assignment problem, the traditional 'perfect rationality' (PR) route choice paradigm (Wardrop 1952) assumes that given the travel cost of each road, travellers only take the shortest paths (i.e., the least utility paths). However, this assumption is very restrictive in reality. Many empirical studies, including simulation experiments (Nakayama et al. 2001), stated preference surveys (Avineri and Prashker 2004), and GPS vehicle trajectory data (Morikawa et al. 2005, Zhu 2011) showed that in practice drivers do

not always choose the shortest paths, and the classical Wardrop user equilibrium (UE) assignment model cannot give accurate prediction of traffic flow patterns. Thus, many other behavioral models have been developed to relax the PR assumption.

As opposed to 'rationality as optimization,' Herbert Simon, in 1957, proposed that people are boundedly rational in their decision- making process. This is either because people lack accurate information, or they are incapable of obtaining the optimized decision due to the complexity of the situations. They tend to seek a satisfactory choice solution instead. Since then, the concept of 'bounded rationality' (BR) was studied extensively in Economics and Psychology literature (See Conlisk 1996 for a detailed review).

Mahmassani and Chang (1987) introduced the BR assumption for the first time in modeling pre- trip departure time selection for a single bottleneck. Since then, a large body of literature incorporated bounded rationality in various transportation models, such as route choice behavior (Han and Timmermans 2006; Nakayama et al. 2001), traffic safety (Sivak 2002), hyperpath assignment (Fonzone and Bell 2010), transportation planning (Gifford and Checherita 2007; Khisty and Arslan 2005), traffic policy making (Marsden et al. 2012) and so on. Some researchers also incorporated some thresholds to the discrete choice model to capture the impact of inertia on travellers' choices (Cantillo et al. 2006, 2007). All these studies indicated that the BR assumption plays a very important role in transportation modeling.

When the BR assumption is used to model drivers' route choice behavior, there are two aspects regarding the boundedly rational route choice process. Some studies suggested that travellers do not take the shortest paths because they are not capable of perceiving actual travel costs due to limited cognitive capacity, or it is too costly to search information about all alternative paths (Gabaix et al. 2006; Gao et al. 2011).

On the other hand, some studies assumed that all path cost information are available to travellers through some information system, but they will not switch to shorter paths due to the existence of inertia, which was quantified by a term named 'indifference band' (Mahmassani and Chang 1987). A series of experiments were conducted by Mahmassani and his colleges to validate this BR behavioral assumption and calibrate the values of indifference bands (Hu and Mahmassani 1997; Jayakrishnan et al. 1994a; Mahmassani and Chang 1987; Mahmassani and Jayakrishnan 1991; Mahmassani and Liu 1999; Srinivasan and Mahmassani 1999). These experiments were conducted on an interactive simulator- DYNASMART, incorporating pre- trip departure time, route choices and en- route path switching decisions. Subjects, as travellers, could change paths en- route at each node and also adjust their departure- time choices the next day based on the previous days' travel experiences. Travellers were assumed to follow the BR behavioral rule in decision- making processes, i.e., they would only switch routes when the improved trip time exceeded some indifference bands. The values of these indifference band depended on individual characteristics and network performances. Lu and Mahmassani (2008) further studied the impact of the congestion pricing on drivers' behavior within the boundedly rational behavioral framework.

In this study, we will adopt the second aspect of bounded rationality. Within this BR framework, travel costs are assumed to be deterministically flow- dependent, and travellers can perceive travel costs accurately but some indifference bands exist due to inertia to switch routes. When traffic flow patterns stabilize to some equilibrium, called 'boundedly rational user equilibria' (BRUE), travellers can take any route whose travel time is within an indifference band of the shortest path cost (Guo and Liu 2011; Lou et al. 2010). Indifference bands vary among origin- destination (OD) pairs. By introducing one parameter (indifference band) for each OD pair, the BR framework relaxes the restrictive PR assumption that travellers only take the shortest paths at equilibrium.

According to Ben- Akiva et al. (1984), travellers' route choice behavior is regarded as a two- stage process: path set generation (i.e., a path choice set is generated between origin and destination according to route characteristics) and traffic assignment (i.e., the traffic demands are mapped to these generated paths based on certain traffic assignment criteria). Accordingly, in this paper, we will study how to generate boundedly rational path sets first, and then assign traffic demands to these paths based on the BRUE condition. Many path generation algorithms were proposed in the existing literature, such as K- shortest path algorithm, labeling approach (Ben- Akiva et al. 1984), parametric least- generalized cost path algorithm (Mahmassani et al. 2005), doubly stochastic choice set generation (Bovy and Fiorenzo- Catalano 2007), pareto paths generation (Wang et al. 2009), and so on. We will propose a different path generation algorithm by solving a sequence of mathematical programs with equilibrium constraints.

After the BR path sets are generated, traffic demands will be assigned to these paths. One of the difficulty about the BR traffic assignment is that BRUE are generally not unique. The BRUE solution set includes all possible traffic flow equilibrium patterns over a traffic network. Obtaining all BRUE flow patterns offers urban planners an opportunity to investigate variability of traffic flows in a network. However, most existing literature on BRUE tried to avoid obtaining complete BRUE solutions due to its complex properties. Lou et al. (2010) for the first time studied some topological properties of the BRUE solution set, i.e., nonconvexity, and implemented congestion pricing based on the BR behavioral assumption. However, it did not provide a systematic approach of solving the BRUE set. This paper constructs the BRUE solution set by focusing on networks with fixed demand. Proposing the methodology of obtaining the BRUE set and exploring the fundamental mathematical properties of BRUE will serve as a building block for BRUE related applications, such as BR- related congestion pricing and other network design problems.

Following the two- stage route choice process, the rest of the paper is organized as follows: In Section 2, the  $\epsilon$  - BRUE is defined and formulated as a nonlinear complementarity problem (NCP). In Section 3, the BRUE- related acceptable path set is defined, and its structure is studied. In Section 4, how to obtain the acceptable path set is presented. In Section 5, we will construct the BRUE path flow solution set based on the acceptable path set. Some examples are given to illustrate the structure of the BRUE path flow solution set. Conclusions and future work are discussed in Section 6.

# 2. Definition of  $\epsilon$  -BRUE and nonlinear complementarity formulation

The traffic network is represented by a directed graph that includes a set of consecutively numbered nodes,  $\mathcal{N}$  , and a set of consecutively numbered links,  $\mathcal{L}$  Let  $w$  denote the OD pair set connected by a set of simple paths (composed of a sequence of distinct nodes),  $\mathcal{P}^w$  , through the network. The traffic demand for OD pair  $W$  is  $d^w$  .Let  $f_{i}^{w}$  denote the flow on path  $i\in \mathcal{P}$  for OD pair  $W$  , then the path flow vector is  $\mathbf{f} = \{f_i^w\}_{i\in \mathcal{P}^w}$  . The feasible path flow set is to assign the traffic demand on the feasible paths:  $\mathcal{F}\triangleq \{\mathbf{f}:\mathbf{f}\geqslant \mathbf{0},\sum_{i\in \mathcal{P}}f_i^w = d^w,\forall w\in \mathcal{W}\}$  . Denote  $x_{a}$  as the link flow on link  $a$  , then the link flow vector is  $\mathbf{x} = \{x_{a}\}_{a\in \mathcal{L}}$  . Each link  $a\in \mathcal{L}$  is assigned a cost function of the link flow, written as  $t(\mathbf{x})$  . Let  $\delta_{a,i}^{w} = 1$  if link  $a$  is on path  $i$  connecting OD pair  $W$  , and O if not; then  $\Delta \triangleq \{\delta_{a,i}^{w}\}_{a\in \mathcal{L},i\in \mathcal{P}}^{w}$  denotes the link- path incidence matrix. Therefore  $f_{i}^{w} = \sum_{a}c_{a,i}^{w}x_{a}$  , and it can be rewritten in a vector form as  $\mathbf{x} = \Delta \mathbf{f}$  . Denote  $C_i^w (\mathbf{f})$  as the path cost on path  $i$  for OD pair  $W$  , then the path cost vector  $C(\mathbf{f})\triangleq \{C_i^w (\mathbf{f})\}_{i\in \mathcal{P}}^{wW}$  . So  $C(\mathbf{f}) = \Delta^T t(\mathbf{x})$  under the additive path cost assumption.

In this paper, we assume the link cost is separable, continuous and linear with respect to its own link flow, i.e.,  $t(\mathbf{x}) = H\mathbf{x}$  , where  $H$  is the Jacobian matrix of the link cost. Then the path cost can be computed as:  $C(\mathbf{f}) = \Delta^T t(\mathbf{x}) = \Delta^T H\Delta \mathbf{f}\triangleq A\mathbf{f}$

# 2.1. Definition

Mahmassani and Chang (1987) defined a BRUE flow vector as the one 'whenever all users' perceived travel costs on their selected routes are constrained within their respective indifference bands'. Following this line, Lou et al. (2010) and Guo and Liu (2011) specialize this description in a mathematical way as follows.

Definition 2.1. For a given nonnegative vector  $\epsilon = (\epsilon^{w})_{i\in W},\epsilon^{w}\geqslant 0$  , a feasible path flow vector  $\mathbf{f}\in \mathcal{F}$  is said to be a  $\epsilon$  - boundedly rational user equilibrium (BRUE) path flow pattern, denoted by  $\mathbf{f}_{BRUE}^{\epsilon}$  , if

$$
f_{i}^{w} > 0\Rightarrow C_{i}^{w}(\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w},\forall i\in \mathcal{P},\forall w\in \mathcal{W} \tag{1}
$$

This definition says that, for a one path flow pattern which is a boundedly rational user equilibrium, travellers only pick any route that is within a given indifference band  $\epsilon$  of the shortest path.

Remark. - Equation (1) gives a necessary condition judging whether a flow pattern is BRUE, and is equivalent to the following condition:

$$
C_i^w (\mathbf{f}) > \min_{j\in \mathcal{P}^w}C_j^w (\mathbf{f}) + \epsilon^w\Rightarrow f_i^w = 0,
$$

But the inverse is not always true:

$$
f_{i}^{w} = 0\neq C_{i}^{w}(\mathbf{f}) > \min_{j\in \mathcal{P}^{w}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w}.
$$

In other words, an unused path may have lower cost than a used one, which will never happen in the UE setting. Therefore, if  $C_i^w (\mathbf{f})\leqslant \min_{j\in \mathcal{P}^w}C_j^w (\mathbf{f}) + \epsilon^w$  , then  $f_{i}^{w}\geqslant 0$

When  $\epsilon = 0$  the BRUE definition is reduced to:

$$
\mathcal{F}_{BRUE}^{0}\triangleq \mathcal{F}_{UE} = \{\mathbf{f}\in \mathcal{F}:f_{i}^{w} > 0\Rightarrow C_{i}^{w}(\mathbf{f}) = \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}),\forall i\in \mathcal{P},\forall w\in \mathcal{W}\} ; \tag{2}
$$

The path flow set satisfying the above definition is called the UE path flow set. Based on the UE path flow set  $\mathcal{F}_{UE}$  , the UE shortest path set  $\mathcal{P}_{UE}$  can be defined as:

$$
\mathcal{P}_{UE} = \{i\in \mathcal{P}:C_i^w (\mathbf{f}) = \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}),\forall \mathbf{f}\in \mathcal{F}_{UE}\} . \tag{3}
$$

Note. The UE shortest path from  $\mathcal{P}_{UE}$  may carry flow or may have no flow on it.

Usually the  $\epsilon$  - BRUE is non- unique. Denote a set containing all path flow patterns satisfying Definition (1) as the  $\epsilon$  - BRUE path flow solution set:

$$
\mathcal{F}_{BRUE}^{\epsilon}\triangleq \{\mathbf{f}\in \mathcal{F}:f_i^w >0\Rightarrow C_i^w (\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w,\forall i\in \mathcal{P},\forall w\in \mathcal{W}\} ; \tag{4}
$$

Proposition 2.1. If the link cost function is continuous, the  $\epsilon$  - BRUE solution  $(\epsilon \geqslant 0)$  is non- empty.

Proof. First, Patriksson (1994) showed that, when the link cost function is continuous, UE exists. Let  $\mathbf{f}\in \mathcal{F}_{UE}$  be the UE path flow pattern, when  $\epsilon \geq 0$

$$
f_{i}^{w} > 0\Rightarrow C_{i}^{w}(\mathbf{f}) = \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w},\forall i\in \mathcal{P},\forall w\in \mathcal{W}
$$

So f is also a  $\epsilon$  - BRUE (  $\epsilon \geq 0)$  , i.e.,  $\mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon}$  . Then  $\mathcal{F}_{UE}\subseteq \mathcal{F}_{BRUE}^{\epsilon}$  . In other words, UE must be contained in the BRUE set. Given the continuous link cost function, at least one BRUE flow pattern exists, and therefore  $\mathcal{F}_{BRUE}^{\epsilon}\neq \emptyset$

Note.  $\mathcal{F}_{UE}$  may be non- unique if the link cost function is not strictly monotone. In spite of it non- uniqueness, it is still contained in the  $\epsilon$  - BRUE set.

# 2.2. BRUE-NCP formulation

Given the continuous link cost function, the  $\epsilon$  - BRUE must exist. The next question is, how we can compute the equilibrium solutions. Since UE is a special case of the BRUE, we will start with the UE. Based on the Wardrop's first principle, the UE can be solved from a nonlinear complementarity problem (NCP). For all  $i\in \mathcal{P}^w$  and all  $w\in \mathcal{W}$  ..

$$
\begin{array}{rl} & 0\leqslant f_i^w\perp C_i^w (\mathbf{f}) - \pi^w\geqslant 0,\\ & 0\leqslant \pi^w\perp d^w -\sum_{i\in \mathcal{P}}f_i^w\geqslant 0. \end{array} \tag{5b}
$$

where  $\pi^w$  is the shortest path cost for the OD pair  $\mathcal{W}$

Similarly, BRUE can be formulated as a NCP as well, but some changes should be made.

Proposition 2.2 (BRUE NCP). Given  $\epsilon^w (w\in \mathcal{W})$  , and some  $\pmb {\rho} = (\rho_{i}^{w})_{i}^{w}$  , where  $0\leq \rho_{i}^{w}\leq \epsilon^{w}$  . A feasible path flow vector  $\mathbf{f}\in \mathcal{F}$  is a  $\epsilon$  - BRUE path flow pattern if and only if it solves the following NCP,  $\forall i\in$ $\mathcal{P}^{w}$ $\forall w\in \mathcal{W}$  ..

$$
\begin{array}{rl} & {(\mathcal{NCP}(\pmb {\rho}))\qquad 0\leqslant f_i^w\perp C_i^w (\mathbf{f}) + \rho_i^w -\pi^w\geqslant 0,}\\ & {\qquad 0\leqslant \pi^w\perp d^w -\sum_{i\in \mathcal{P}}f_i^w\geqslant 0.} \end{array} \tag{6b}
$$

where the physical meaning of  $\pi^w$  is the maximal path cost within the band  $\epsilon^w$  for OD pair  $w$

Proof. In the fixed- demand case, (6b) is reduced to the nonnegative constraint on  $\pi^w\geqslant 0$

(1) To prove the necessary part. Let  $\mathbf{f}$  be a feasible flow pattern and let  $(\rho ,\pi)$  be a pair such that Equation 
(6) holds, then  $0\leqslant \rho_{i}\leqslant \epsilon ,\pi^{w}\geqslant 0$  for all  $w\in W$  and  $i\in P^w$

Moreover,  $C_i^w (\mathbf{f}) + \rho_i^w - \pi^w\geqslant 0$  and  $0\leqslant \rho_{i}^{w}\leqslant \epsilon^{w}$  indicate that

$$
\pi^w\leqslant C_i^w (\mathbf{f}) + \rho_i^w\leqslant C_i^w (\mathbf{f}) + \epsilon^w,
$$

thus,

$$
\pi^w\leqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w.
$$

For  $f_{i}^{w} > 0$  (6a) holds iff  $C_i^w (\mathbf{f}) + \rho_i^w - \pi^w = 0$  ,i.e.,

$$
f_{i}^{w} > 0\Rightarrow C_{i}^{w}(\mathbf{f}) = \pi^{w} - \rho_{i}^{w}\leqslant \pi^{w}\leqslant \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w},
$$

which satisfies Definition (2.1), so f is a  $\epsilon$  - BRUE path flow pattern.

(2) To prove the sufficient part, suppose  $\mathbf{f}$  is a BRUE flow pattern. For  $i\in \mathcal{P}^w$  and  $w\in W$  , define

$$
\pi^w\triangleq \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w\geqslant 0. \tag{7}
$$

and

$$
\begin{array}{r}\rho_i^w\triangleq \left\{ \begin{array}{ll}\pi^w -C_i^w (\mathbf{f}),\mathrm{if}C_i^w (\mathbf{f})\leqslant \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon^w;\\ 0,\mathrm{if}C_i^w (\mathbf{f}) > \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon^w; \end{array} \right. \end{array} \tag{8}
$$

Equation (6) holds automatically for any BRUE flow pattern. It suffices to show that  $0\leqslant \rho_{i}^{w}\leqslant \epsilon^{w},\pi^{w}\geqslant$  0. Since f is a BRUE flow pattern, it follows that, by Definition (2.1), if  $f_{i}^{w}\geqslant 0,C_{i}^{w}(\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w}$  if  $C_i^w (\mathbf{f}) > \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w$  , then  $f_{i}^{w} = 0$

When  $f_{i}^{w}\geqslant 0$  , then  $C_i^w (\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w$

$$
C_i^w (\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w = \pi^w,
$$

which yields  $\rho_{i}^{w}\geqslant 0$  from (8). Again, by (8),

$$
\rho_{i}^{w} = \pi^{w} - C_{i}^{w}(\mathbf{f}) = \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w} - C_{i}^{w}(\mathbf{f})\leqslant \epsilon^{w},
$$

showing that  $\rho_{i}^{w}\leqslant \epsilon^{w}$  . So  $0\leqslant \rho_{i}^{w}\leqslant \epsilon^{w}$  if  $f_{i}^{w} > 0$

When  $C_i^w (\mathbf{f}) > \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w$  , i.e.,  $f_{i}^{w} = 0$  , then  $\rho_{i}^{w} = 0$  by (8).

In summary,  $0\leqslant \rho_{i}^{w}\leqslant \epsilon^{w}$  and  $\pi^w\geqslant 0$  for  $i\in P^w$  and  $w\in W$

Remark. - Compare the UE- NCP with the BRUE- NCP formulation, there is one additional term  $\rho_{i}^{w}$  in BRUE- NCP, we call it 'indifference function.' If  $C_i^w (\mathbf{f})\geqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon^w$  , then  $\rho_{i}^{w} = 0$  ; if  $C_i^w (\mathbf{f}) =$ $\min_{j\in \mathcal{P}}C_j^w (\mathbf{f})$  , then  $\rho_{i}^{w} = \epsilon^{w}$

The meaning of  $\pi^w$  is different in the two settings. Regarding UE,  $\pi^w$  is the shortest travel time; while in BRUE, its value is equal to the shortest travel cost plus  $\epsilon^w$  (see Equation (7)).  $\pi^w$  is a function of the flow pattern, so we call it the 'maximum path cost' for a specific BRUE flow pattern. For some flow pattern, there may not exist a path with the exact cost of  $\pi^w$ . If some flow pattern happens to have a path with the cost of  $\pi^w$ , this path may or may not carry flows.

From BRUE NCP (6) and its proof, we have another conclusion:

Corollary 1. If  $\mathbf{f} \in \mathcal{F}_{BRUE}^{\epsilon}$ , there must exist at least one vector pair  $(\rho , \pi)$  satisfying NCP (6). Moreover, one value of  $(\rho , \pi)$  is determined by (7) and (8).

By substituting one indifference function  $\rho \in \mathbb{R}_+^q$  into the BRUE- NCP, one BRUE path flow pattern can be obtained. The BRUE- NCP formulation provides an approach of solving one BRUE solution. In the following we will show how to construct the complete BRUE solution sets out of a specific solution.

# 3. Monotonically non-decreasing acceptable path set

In the last section, we show that, when the indifference band is zero, the BRUE is equivalent to the UE flow pattern, and travellers will only take shortest paths. When the indifference band gradually increases, some paths which are too costly to take under the UE may be utilized under the BRUE. In this section, we will discuss the relationship between the indifference band and the number of utilized paths.

# 3.1. Monotonically non-decreasing property

All feasible paths for one particular BRUE flow pattern can be classified into three categories:

Definition 3.1. Given a- BRUE flow pattern  $\mathbf{f} \in \mathcal{F}_{BRUE}^{\epsilon}$ , the total feasible paths could have three statuses: acceptable, zero- acceptable and unacceptable. The acceptable path carries flow, while its cost is within the shortest cost plus the band; the zero- acceptable path is acceptable in terms of the cost, but carries no flow; and the unacceptable path is longer than the shortest cost plus the band.

$$
\begin{array}{r}a^{\epsilon}(\mathbf{f}) = \{i\in \mathcal{P}:f_i > 0,C_i^w (\mathbf{f})\leqslant \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon^w,\forall w\} ;\\ 0^{\epsilon}(\mathbf{f}) = \{i\in \mathcal{P}:f_i = 0,C_i^w (\mathbf{f})\leqslant \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon^w,\forall w\} ;\\ u^{\epsilon}(\mathbf{f}) = \{i\in \mathcal{P}:f_i = 0,C_i^w (\mathbf{f}) > \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon^w,\forall w\} . \end{array} \tag{9b}
$$

There are three properties of the above three path sets:

Proposition 3.1. (1) Given one BRUE flow pattern  $\mathbf{f}$ , the union of paths at these three status is the feasible path set:  $a^{\epsilon}(\mathbf{f}) \cup 0^{\epsilon}(\mathbf{f}) \cup u^{\epsilon}(\mathbf{f}) = \mathcal{P}$

(2) Given one BRUE flow pattern  $\mathbf{f}$  we can always find at least one path which is acceptable or zeroacceptable:  $a^\epsilon (\mathbf{f})\cup 0^\epsilon (\mathbf{f})\neq 0$  
(3) Given one BRUE flow pattern f, if  $\begin{array}{r}\epsilon \geqslant \max_{j\in \mathcal{P}}C_j^w (\mathbf{f}) - \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}), \end{array}$  all feasible paths are either acceptable or zero-acceptable:  $u^{\epsilon}(\mathbf{f}) = \emptyset$

Proof. (1) It is obvious from the definition.

(2) Since the shortest path  $i = \arg \min_{j\in \mathcal{P}}C_j^w (\mathbf{f})$  always exists, so  $i\in a^{\epsilon}(\mathbf{f})\cup 0^{\epsilon}(\mathbf{f})\subseteq \mathcal{P}$

(3)  $\begin{array}{r}\max_{j\in \mathcal{P}}C_j^w (\mathbf{f}) - \min_{j\in \mathcal{P}}C_j^w (\mathbf{f})\leqslant \epsilon \end{array}$  implies  $C_i^w (\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_j^w (\mathbf{f}) + \epsilon ,\forall i\in \mathcal{P}$  , thus no path is unaceceptable, i.e.,  $u^{\epsilon}(\mathbf{f}) = \emptyset$

Definition (3.1) divides all the feasible paths for one BRUE flow pattern into three classes. Each status notation indicates the dependency of the path status on  $\epsilon$  and the specific BRUE flow pattern. The following proposition will discuss the relationship between the path status and the value of  $\epsilon$ .

Proposition 3.2. Given  $\mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon}$ $if0\leq \epsilon < \epsilon^{\prime}$  , then  $a^{\epsilon}(\mathbf{f})\subseteq a^{\epsilon^{\prime}}(\mathbf{f}),0^{\epsilon}(\mathbf{f})\subseteq 0^{\epsilon^{\prime}}(\mathbf{f}).$

Proof. It suffices to show that, if an arbitrary path  $i$  is acceptable or zero- acceptable under the indifference band  $\epsilon$  , it is also acceptable or zero- acceptable when the indifference band is  $\epsilon^{\prime} > \epsilon$ $\forall i\in a^{\epsilon}(\mathbf{f})$  , we know  $f_{i} > 0$  , and

$$
C_{i}^{w}(\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{w}< \min_{j\in \mathcal{P}}C_{j}^{w}(\mathbf{f}) + \epsilon^{'}\Rightarrow i\in a^{\prime}(\mathbf{f}).
$$

Therefore  $a^{\epsilon}(\mathbf{f})\subseteq a^{\epsilon^{\prime}}(\mathbf{f})$  . Similarly  $0^{\epsilon}(\mathbf{f})\subseteq 0^{\epsilon^{\prime}}(\mathbf{f})$

Proposition (3.2) says, for a  $\epsilon$  - BRUE flow pattern, the status of a path depends on the value of  $\epsilon$  . Since its total feasible paths are fixed, those unused paths under smaller  $\epsilon$  can be utilized with bigger  $\epsilon$  . Therefore, the bigger  $\epsilon$  is, the more paths are acceptable or zero- acceptable, and the less paths are unacceptable.

In the following, we will discuss the impact of the value of  $\epsilon$  on the size of the  $\epsilon$  - BRUE flow set.

Proposition 3.3 (Monotonically non- decreasing flow set). If  $0\leq \epsilon < \epsilon^{\prime}$  , then  $\mathcal{F}_{BRUE}^{\epsilon}\subseteq \mathcal{F}_{BRUE}^{\epsilon^{\prime}}$

Proof. It suffices to show that  $\forall \mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon}\Rightarrow \mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon^{\prime}}$

$$
\begin{array}{rl}\forall \mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon} & \Rightarrow \{f_i^w >0\Rightarrow C_i^w (\mathbf{f})\leqslant \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon \\ & < \underset {j\in \mathcal{P}}{\min}C_j^w (\mathbf{f}) + \epsilon ',\forall i\in \mathcal{P},\forall w\in \mathcal{W}\} \\ & \Rightarrow \mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon '} & \Rightarrow \mathcal{F}_{BRUE}^{\epsilon}\subseteq \mathcal{F}_{BRUE}^{\epsilon '} \end{array}
$$

Proposition (3.3) indicates that when  $\epsilon$  increases, more flow patterns will be included in the  $\epsilon$  - BRUE flow set.

When the  $\epsilon$  - BRUE flow set exists and is non- unique, each flow pattern has different combination of acceptable, zero- acceptable and unacceptable paths. A largest  $\epsilon$  - acceptable path set' contains all acceptable paths for every flow pattern in the  $\epsilon$  - BRUE flow set, mathematically:

$$
\mathcal{P}_l^\epsilon = \bigcup_{\mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon}}a^\epsilon (\mathbf{f}). \tag{10}
$$

The largest  $\epsilon$  - acceptable path set shares the similar property as the  $\epsilon$  - BRUE flow set has:

Proposition 3.4 (Monotonically non- decreasing path set). If  $0\leq \epsilon < \epsilon^{\prime}$ $\mathcal{P}_l^{\epsilon}\subseteq \mathcal{P}_l^{\epsilon^{\prime}}$  , where  $\mathcal{P}_l^{\epsilon}$  is defined in (10).

Proof.

$$
\begin{array}{rl}\mathcal{P}_l^\epsilon & = \bigcup_{\mathbf{f}\in \mathcal{F}_{BRUE}^\epsilon}a^\epsilon (\mathbf{f})\subseteq \bigcup_{\mathbf{f}\in \mathcal{F}_{BRUE}^\epsilon}a^\epsilon (\mathbf{f})\mathrm{(Proposition~(3.3))}\\ & \subseteq \bigcup_{\mathbf{f}\in \mathcal{F}_{BRUE}^{\epsilon'}}a^{\epsilon '}(\mathbf{f})\mathrm{(Proposition~(3.2))} = \mathcal{P}_l^{\epsilon '} \end{array}
$$

When  $\epsilon$  varies from zero to infinity, the minimum number of paths the largest  $\epsilon$  - acceptable path set contains is the UE shortest paths when  $\epsilon = 0$  , i.e.,  $\mathcal{P}_{BRUE}^{0}\triangleq \mathcal{P}_{UE}$  . The maximum number of paths the largest  $\epsilon$  - acceptable path set contains is all feasible paths, meaning all feasible paths will be utilized if the indifference band is too large. Then we have the following corollary:

Corollary 2.  $\mathcal{P}_{UE}\subseteq \mathcal{P}_l^{\epsilon}\subseteq \mathcal{P}$

Given  $\epsilon$  , the largest  $\epsilon$  - acceptable path set (defined in Equation (10)) is a set of all acceptable paths under the  $\epsilon$  - BRUE set. It is possible that some acceptable paths for one  $\epsilon$  - BRUE flow pattern are not acceptable for other flow patterns; vice versa. This necessaries the exploration of the interior structure of the  $\epsilon$  - acceptable path set. Proposition (3.4) provides us with one approach of analyzing its structure by varying values of  $\epsilon$

The path set is a set of finite paths, while  $\epsilon$  is treated as a continuous parameter for the time being. Start with the UE path set when  $\epsilon = 0$  : provided the topology of a network and the link cost functions,  $\mathcal{P}_{UE}$  can be determined by some established algorithms, e.g., column generation algorithm (Patriksson, 1994), gradient projection algorithm (Jayakrishnan et al., 1994b), or maximum entropy algorithm (Bell and Iida, 1997). According to Proposition (3.4), when  $\epsilon$  is gradually increased, more paths will be included, and we should be able to identify those acceptable paths one by one, until all alternative paths are included. This offers the theoretical foundation for deriving different combinations of acceptable paths by varying  $\epsilon$  subsequently.

# 3.2. Definition of  $\epsilon$  -acceptable path set

It is assumed that there are  $n$  alternative paths for OD pair  $W$  ,i.e.,  $\mathcal{P} = \{1,\dots ,n\}$  and  $|\mathcal{P}| = n$  , where  $|\mathcal{P}|$  is the cardinality of set  $\mathcal{P}$  . Among these  $n$  paths, there are  $p$  shortest paths at the UE, i.e.,  $\mathcal{P}_{UE} = \{1,\dots ,p\}$  and  $|\mathcal{P}_{UE}| = p\leqslant n$

Definition 3.2. Assuming there exists a unique sequence of finite critical points  $\epsilon_{j}^{*w}$ $(j = 1,\dots ,J)$  , with  $\epsilon_0^* = 0,\epsilon_{j + 1}^* = \infty$  dividing the nonnegative real line into  $(J + 1)$  intervals:  $[0,\infty) = [0,\epsilon_1^*)\cup [\epsilon_1^*,\epsilon_2^*)\dots \cup$ $[\epsilon_j^*,\infty) = \bigcup_{j = 0}^J [\epsilon_j^*,\epsilon_{j + 1}^*)$  . A sequence of critical points are defined as:

$$
\begin{array}{rcl}\epsilon_1^* & \triangleq & \inf \{\mathcal{P}_{UE}\subset \mathcal{P}_l^\epsilon \} ;\\ & & \vdots \\ \epsilon_j^* & \triangleq & \inf \{\mathcal{P}_l^{\epsilon_{j - 1}^\epsilon}\subset \mathcal{P}_l^\epsilon \} ;\\ & & \vdots \\ \epsilon_j^* & \triangleq & \inf \{\mathcal{P}_l^\epsilon = \mathcal{P}\} . \end{array} \tag{11}
$$

The largest  $\epsilon$  - acceptable path set will remain the same until  $\epsilon$  reaches these values, i.e., for  $\epsilon_{j}^{*}\leq \epsilon_{1}< \epsilon_{2}<$ $\epsilon_{j + 1}^{*},\mathcal{P}_{l}^{\epsilon_{j}^{*}} = \mathcal{P}_{l}^{\epsilon_{1}} = \mathcal{P}_{l}^{\epsilon_{2}}\subset \mathcal{P}_{l}^{\epsilon_{j + 1}^{*}}$

![](images/981483f9bac2967cc51ec81d399d454766bb1957eaf34aae74af91e1bb3d12f7.jpg)  
Fig. 1: Monotonically non-decreasing property illustration

A newly added path' is defined as the path which is unacceptable under  $\epsilon_{j - 1}^{*}$  but acceptable when  $\epsilon = \epsilon_{j}^{*}$

$$
r_j^*\triangleq \{i\in \mathcal{P}:i\in \mathcal{P}_l^{\epsilon_j^*},i\notin \mathcal{P}_l^{\epsilon_{j - 1}^*}\} , \tag{12}
$$

Given the current acceptable path set  $\mathcal{P}_l^{\epsilon_{j - 1}^*}$ , after newly acceptable paths are identified, the updated acceptable path set is:

$$
\mathcal{P}_l^{\epsilon_j^*}\triangleq \mathcal{P}_l^{\epsilon_{j - 1}^*}\bigcup r_j^*. \tag{13}
$$

Remark. There may exist two or more paths added at the same time, so  $r_j^*$  should be treated as a path set.

Provided a fixed indifference band  $\epsilon$ , and  $\epsilon_{I}^{*}\leqslant \epsilon < \epsilon_{I + 1}^{*}$ . Let  $K = \min \{I,J\}$ . The  $\epsilon$ - acceptable path set is defined as:

$$
\mathcal{P}^{\epsilon} = \{\mathcal{P}_{UE},\mathcal{P}_I^{\epsilon_1^*},\dots ,\mathcal{P}_I^{\epsilon_K^*}\} . \tag{14}
$$

In other words, the  $\epsilon$ - acceptable path set is composed of  $K$  acceptable path sets with  $\mathcal{P}_{UE}\subset \mathcal{P}_I^{\epsilon_1^*}\subset \dots \subset \mathcal{P}_I^{\epsilon_K^*}\subseteq \mathcal{P}$ .

# 4. Generation of the  $\epsilon$ -acceptable path set

In the last section, we introduced a BR- related path set: the largest  $\epsilon$ - acceptable path set. In this section, we will explore how to generate this path set. Definition (11) says that, the largest  $\epsilon$ - acceptable path set includes more paths when  $\epsilon$  increases to some critical values. Thus, a mathematical program with equilibrium constraint (MPEC) can be developed to solve these critical values:

$$
\begin{array}{rl} & 0\leqslant f_i\perp C_i(\mathbf{f}) + \rho_i - \pi \geqslant 0,\forall i\in \mathcal{P},\\ & \sum_{r\in \mathcal{P}}f_r = d,\\ & d - \sum_{j\in \mathcal{P}_l^{\epsilon_j^* -1}}f_j > 0,\\ & 0\leqslant \rho_i\leqslant \epsilon_j,\forall i\in \mathcal{P},\\ & f_i + C_i(\mathbf{f}) + \rho_i - \pi >0,\forall i\in \mathcal{P}. \end{array} \tag{15b}
$$

(15a- 15b) is to guarantee the path flow pattern is a feasible BRUE; (15c) tries to push a small amount of flow from the acceptable path set  $\mathcal{P}_I^{\epsilon_{j - 1}^*}$  to some newly acceptable path if  $\epsilon$  is increased a little bit; (15d) sets the bounds for the indifference function; (15e) ensures the strict complementarity condition in (15a) (Cottle et al. 2009).

If MPEC (15) is solvable, optimal solutions  $(\mathbf{f}^{*},\pmb{\rho}^{*},\pi^{*},\epsilon_{j}^{*})$  will be obtained. The newly added path  $r_j^*$  can be derived from  $\mathbf{f}^{*}$ . It is the path that is excluded from  $\mathcal{P}_I^{\epsilon_{j - 1}^*}$ , but begins to carry a very small amount of flow in  $\mathbf{f}^{*}$ .

(15c) and (15e) are inequalities without equal sign, which defines an open set causing non- attainability of the optimal solution. So we introduce a parameter  $0< \delta \leqslant 1$  such that  $d - \sum_{r\in \mathcal{P}_I^{\epsilon_{j - 1}^*}}f_r\geqslant \delta d$  and  $f_{i} + C_{i}(\mathbf{f})+$ $\rho_{i} - \pi \geqslant \delta$ . We call this modified version as  $\delta$ - MPEC, and we will solve this version in practice by giving

$\delta$  a very small value  $\delta = 0.01$  works well). Rewrite the  $\delta$ - MPEC in a compact form:

$$
\begin{array}{rl}(\delta \text{-MPEC}) & \min \qquad \epsilon_j\\ s.t.\qquad \mathbf{1}^T\mathbf{f} = d,\\ & \quad -\mathbf{v}^T\mathbf{f}\geqslant (\delta -1)d,\\ & \quad \epsilon_j\mathbf{1} - \rho \geqslant \mathbf{0},\\ & \quad \mathbf{f} + C(\mathbf{f}) + \rho -\pi \geqslant \delta ,\\ & \quad \rho ,\pi ,\epsilon_j\geqslant \mathbf{0},\\ & \quad \mathbf{f}\geqslant \mathbf{0},\\ & \quad C(\mathbf{f}) + \rho -\pi \mathbf{1}\geqslant \mathbf{0},\\ & \quad \mathbf{f}^T (C(\mathbf{f}) + \mathbf{f}^T\rho -\pi \mathbf{1}) = \mathbf{0}. \end{array} \tag{16h}
$$

where  $\mathbf{v}$  is an vector of the same dimension with the path flow vector, with the  $i^{th}$  component equal to one if  $f_{i} \in \mathcal{P}^{\epsilon_{j - 1}}$ , and zero otherwise.  $\mathbf{1}$  is a vector of 1.

In practice, this MPEC problem can be solved by GAMS software (General Algebraic Modeling System, see Rosenthal and Brooke 2007).

# 4.1. Solving critical points sequentially

When  $\epsilon_1^*$  is achieved, we include the corresponding path  $r_1^*$  into the UE shortest path set and get  $\mathcal{P}_l^{\epsilon_1^*}$ . Next we are interested in finding the critical point  $\epsilon_2^*$  based on current  $\mathcal{P}_l^{\epsilon_1^*}$ . We can solve the above MPEC again by replacing  $\mathcal{P}_l^{\epsilon_1^*}$  in (15b) with  $\mathcal{P}_l^{\epsilon_2^*}$ , and  $\epsilon_2^*$  will be obtained. Similarly,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ ,  $\epsilon_3^*$ , are able to be computed sequentially. This procedure will not stop until all feasible paths are included into the BRUE acceptable set or the critical value reaches the given  $\epsilon$ .

The above procedure provides a method of obtaining the critical points and the order of adding new paths to the acceptable path sets. We will illustrate how to implement this procedure on a small network in the following.

Example 4.1. A four- link network connecting one OD pair in parallel with demand 2, the link cost for each link is 1,  $x_{2} + 1.5$ ,  $x_{3} + 3$ ,  $x_{4} + 3$ . The UE is  $x_{1} = 2$ ,  $x_{2} = x_{3} = x_{4} = 0$ . Four paths are numbered as path 1, 2, 3, 4.

![](images/5bce2c11ae8d5cac34a7cd73ae8f8a17e410b3d8b4cd48536678b52061c3f6a1.jpg)  
Fig. 2: Single-OD pair network illustration

Solving MPEC, we have  $\epsilon_0^* = 0$ ,  $\epsilon_1^* = 0.5$ ,  $\epsilon_2^* = 2$ ,  $\epsilon_3^* = \infty$ . There are three cases for the largest  $\epsilon$ - acceptable path sets:

(1)  $0 \leqslant \epsilon < 0.5$ :  $\mathcal{P}_l^{\epsilon} = \{1\} , \mathbf{f} = [2, 0, 0, 0]$ ; 
(2)  $0.5 \leqslant \epsilon < 2$ :  $\mathcal{P}_l^{\epsilon} = \{1, 2\} , \mathbf{f} = [2, 0, 0, 0]$ ; 
(4)  $\epsilon \geqslant 2$ :  $\mathcal{P}_l^{\epsilon} = \{1, 2, 3, 4\} , \mathbf{f} = [2, 0, 0, 0]$ .

If  $\epsilon$  is calibrated from empirical data as 1.5, then  $\mathcal{P}_l^{1.5} = \{1, 2\}$ . Therefore,  $\mathcal{P}^{1.5} = \{\{1\} \{1, 2\} \}$ .

# 4.2.  $\epsilon$  -acceptable path set for multiple OD pairs

For a network with total  $W$  OD pairs, let  $\epsilon_{j}^{*w}$  be the critical point for OD pair  $w\in \mathcal{W},j = 0,1,\dots ,J^{w}$  then  $\epsilon_{j}^{*}\triangleq \{\epsilon_{j}^{*w}\}$  is a set of critical points for all OD pairs, assuming  $\epsilon_0^{*w} = 0,w\in \mathcal{W}$

Given  $\epsilon = (\epsilon^{1},\dots ,\epsilon^{W})$  , if  $\epsilon_{j}^{W}\leqslant \epsilon^{1}< \epsilon_{n_{K_{1} + 1}}^{1},\dots ,\epsilon_{n_{K_{W}}}^{W}\leqslant \epsilon^{W}< \epsilon_{n_{K_{W} + 1}}^{W}$  , then the largest  $\epsilon$  - acceptable set for multiple OD pairs is the union of the largest  $\epsilon$  - acceptable set for each OD path pair:

$$
\mathcal{P}_l^{\epsilon} = \mathcal{P}_l^{(\epsilon^1,\dots ,\epsilon^w,\dots ,\epsilon^W)} = \mathcal{P}_l^{\epsilon_{K_1}^{*1}}\bigcup \dots \bigcup \mathcal{P}_l^{\epsilon_{n_{K_1}}^{*W}}. \tag{17}
$$

where,  $\mathcal{P}_l^{\epsilon_{n_{K_w}}^{*w}}$  is the largest  $\epsilon_{n_{K_w}}^{*w}$  - acceptable set for OD pair  $W$  and can be solved by MPEC (15).

The  $\epsilon$  - acceptable path set is various combinations of largest  $\epsilon_{n_{K_w}}^{*w}$  - acceptable sets when varying  $\epsilon_{j}^{*w}$  among every critical point across all OD pairs:

$$
\mathcal{P}^{\epsilon} = \mathcal{P}^{(\epsilon^{1},\dots ,\epsilon^{w},\dots ,\epsilon^{W})} = \left\{\mathcal{P}_{l}^{(\epsilon^{1},\dots ,\epsilon_{j}^{*w},\dots ,\epsilon^{W})}\right\}_{\epsilon_{j}^{*w}},w\in \mathcal{W},j = 0,1,\dots ,J^{w}. \tag{18}
$$

Now we will discuss how to solve each critical point for every OD pair. Regarding one OD pair  $\nu \in \mathcal{W}$  the same approach of computing  $\epsilon_{j}^{*y}$  as mentioned in Equation (15) can be adopted. The only difference is, path costs need the information of path flows across all OD pairs. So path flows  $\mathbf{f}^w,w\in \mathcal{W},w\neq \nu$  are parameters when Equation (15) for  $\nu$  is calculated. In other words,  $\epsilon_{j}^{*y}$  is a function of  $\mathbf{f}^w,w\neq \nu$  . But we only want those  $\mathbf{f}^w,w\neq \nu$  such that  $\epsilon_{j}^{*y}$  can be achieved. So accordingly, we modify Equation (15) to accommodate this distinction. For any OD pair  $\nu$  and the  $j^{th},j = 1,\dots ,J^{w}$  critical point  $\epsilon_{j}^{y}$  can be computed as:

$$
\begin{array}{rl} & {\quad \mathrm{min}\qquad \epsilon_{j}^{y}}\\ & {\quad \mathrm{s.t.}}\\ & {\quad 0\leqslant f_{i}^{\nu}\perp C_{i}^{\nu}(\mathbf{f}) + \rho_{i}^{\nu} - \pi^{\nu}\geqslant 0,\forall i\in \mathcal{P}^{\nu},}\\ & {\quad \sum_{j\in \mathcal{P}^{\nu}}f_{j}^{w} = d^{w},\forall w\in \mathcal{W},}\\ & {\quad d^{\nu} - \sum_{j\in \mathcal{P}_{i}^{j - 1}}f_{j}^{\nu} > 0,}\\ & {\quad 0\leqslant \rho_{i}^{\nu}\leqslant \epsilon_{j}^{\nu},\forall i\in \mathcal{P}^{\nu},}\\ & {\quad 0\leqslant f_{i}^{\nu} + C_{i}^{\nu}(\mathbf{f}) + \rho_{i}^{\nu} - \pi^{\nu} > 0,\forall i\in \mathcal{P}^{\nu}.} \end{array} \tag{19e}
$$

The algorithm of calculating  $\epsilon$  - acceptable path set for multiple OD pairs can thus be summarized as follows:

1. Calculate  $\epsilon_{j}^{*y}(j = 1,\dots ,J^{w},\nu \in \mathcal{W})$  from Equation (19), and obtain the  $\epsilon^{\parallel}$  -acceptable path set for OD pair  $\nu$  
2. As a by-product in solving Equation (19), one feasible path flow pattern  $\mathbf{f}^{*w},w\neq \nu$  is attained simultaneously.  $\forall w\neq \nu$  , denote the longest used path for OD pair  $W$  as  $p^w = \left\{\max_{i\in \mathcal{P}^w}C_j^w (\mathbf{f}),f_i^w >0\right\}$  and compute  $C_p^w (\mathbf{f})$  . Then the  $i^{th}(i = 1,\dots ,J^{w})$  critical point for OD pair  $w\epsilon_{i}^{*w} = C_{p}^{w}(\mathbf{f}) - \min_{j\in \mathcal{P}^{w}}C_{j}^{w}(\mathbf{f})$  
3. After obtaining the critical points  $\epsilon_{i}^{*w}$ $(i = 1,\dots ,J^{w},w\neq \nu)$  , if  $\epsilon_{i - 1}^{*w}< C_p^w (\mathbf{f}) - \min_{j\in \mathcal{P}^w}C_j^w (\mathbf{f})< \epsilon_i^{*w}$  then  $\mathbf{f}^{*w}(w\neq \nu)$  is also a BRUE path flow pattern, and the path  $p^w$  is acceptable when  $\epsilon_{j}^{*y}< \epsilon_{j}^{*y}$  and  $\epsilon_{i}^{w}< \epsilon_{i}^{*w}$  ; or else path  $p^w$  is unacceptable; 
4. Combine the acceptable paths under various combinations of critical points among all OD pairs, the  $\epsilon$  -acceptable path set  $\mathcal{P}^{\epsilon}$  is obtained.

The following example will illustrate how to construct the  $\epsilon$  - acceptable path set for multiple OD pairs.

Example 4.2. The topology of the test network, two OD demands and link cost functions are illustrated above each link in Figure (3). Red curves on the right indicate six paths, denoted by the number of links it passes along: 1- 3, 1- 4, 2- 3, 2- 4, 1, 2. The first four paths belong to OD pair 1, the rest two belong to OD pair 2. The equilibrium path flow pattern is  $f_{1} = 0, f_{2} = 1, f_{3} = 0, f_{4} = 0, f_{5} = 1, f_{6} = 0$ .

![](images/bb195c62e672a3adc69077e5c1580df00790034c5678cc62d9d7f05f4af3a5b1.jpg)  
Fig. 3: Two-OD pair network and paths illustration

Based on Equation (19),  $\epsilon_1^{*1}, \epsilon_2^{*1}$  and  $\epsilon_1^{*2}$  can be computed separately. The critical points for each OD pair are:

For OD pair 1,

(1)  $0 \leqslant \epsilon^{1} < 4$ :  $\mathcal{P}^{\epsilon^{1}} = \{2\}$   
(2)  $4 \leqslant \epsilon^{1} < 8$ :  $\mathcal{P}^{\epsilon^{1}} = \{\{2\} , \{1,2\} \}$   
(3)  $8 \leqslant \epsilon^{1} < 12$ :  $\mathcal{P}^{\epsilon^{1}} = \{\{2\} , \{1,2\} , \{1,2,4\} \}$

For OD pair 2,

(1)  $0 \leqslant \epsilon^{2} < 8$ :  $\mathcal{P}^{\epsilon^{2}} = \{5\}$   
(2)  $\epsilon^{2} \geqslant 8$ :  $\mathcal{P}^{\epsilon^{2}} = \{\{5\} , \{5,6\} \}$

Combing two OD pairs, the overall  $\epsilon$ - acceptable path set under different combination of critical points is:

$$
\begin{array}{r l}&{\mathcal{P}^{\epsilon}=\left\{\begin{array}{l l}{\{2,5\},0\leqslant\epsilon^{1}\leqslant4,0\leqslant\epsilon^{2}\leqslant8;}\\ {\{\{2,5\},\{2,5,6\}\},0\leqslant\epsilon^{1}\leqslant4,\epsilon^{2}\geqslant8;}\\ {\{\{2,5\},\{1,2,5\}\},4\leqslant\epsilon^{1}\leqslant8,0\leqslant\epsilon^{2}\leqslant8;}\\ {\{\{2,5\},\{1,2,5\},\{2,5,6\},\{1,2,5,6\}\},4\leqslant\epsilon^{1}\leqslant8,\epsilon^{2}\geqslant8;}\\ {\{\{2,5\},\{1,2,5\},\{1,2,4,5\}\},8\leqslant\epsilon^{1}\leqslant12,0\leqslant\epsilon^{2}\leqslant8;}\\ {\{\{2,5\},\{1,2,5\},\{1,2,4,5\},\{2,5,6\},\{1,2,5,6\},\{1,2,4,5,6\}\},8\leqslant\epsilon^{1}\leqslant12,\epsilon^{2}\geqslant8;}\\ {\{\{2,5\},\{1,2,5\},\{1,2,4,5\},\{1,2,3,4,5\}\},\epsilon^{1}\geqslant12,0\leqslant\epsilon^{2}\leqslant8;}\\ {\{\{2,5\},\{1,2,5\},\{1,2,4,5\},\{2,5,6\},\{1,2,5,6\},\{1,2,4,5,6\},\{1,2,3,4,5\},\{1,2,3,4,5\},\{1,2,3,4,5\},\{\epsilon^{1}\geqslant12,\epsilon^{2}\geqslant8.}\end{array}\right.}\end{array} \tag{20}
$$

All acceptable path sets are also illustrated in the following figure: the numbers in each block displays acceptable path numbers for certain  $(\epsilon^{1}, \epsilon^{2})$  pair. The left bottom block for  $(\epsilon^{1} < 4, \epsilon^{2} < 8)$  is the UE shortest path set. How to get acceptable paths in the block to its right for  $(\epsilon^{1} \geqslant 4, \epsilon^{2} < 8)$  is explained in the following, and other blocks follow the same line of reason. As described above in step (2) and (3), one path flow pattern  $f_{1} = 0.001, f_{2} = 0.9999, f_{3} = 0, f_{4} = 0, f_{5} = 1, f_{6} = 1$  is attained when  $\epsilon^{1} = 4$  is solved from Equation (19). For OD pair 1, its path flow increases from 0 to a positive number 0.001, meaning that path 1 will start to carry flows if  $\epsilon^{1} > 4$ ; the utilized path 5 has the cost of 2 and the unused path 6 has the cost of 10, then their cost difference is 8. When  $\epsilon^{2} < 8$ , only path 5 is acceptable for OD pair 2. Therefore, when  $\epsilon^{1} \geqslant 4$  and  $\epsilon^{2} < 8$ , only paths 1,2 (connecting OD pair 1) and 5 (connecting OD pair 2) are acceptable.

![](images/414504f2f115e7804ea9b007519d756ad3fc69dd4ae9dc968011e4c9a971ea2d.jpg)  
Fig. 4: Acceptable paths for critical point pairs  $(\epsilon^{1},\epsilon^{2})$

By far we have proposed how to solve the  $\epsilon$ - acceptable path set for both single OD pair and multiple OD pairs. The following will discuss the methodology of constructing the  $\epsilon$ - BRUE path flow set.

# 5. Construction of  $\epsilon$ -BRUE path flow set

Generally the  $\epsilon$ - BRUE set is non- convex (Lou et al., 2010), so it is not easy to construct it directly. If we can decompose the whole BRUE set into small subsets, which are easier to study, then constructing the whole set can be reduced to constructing each subset. Based on this idea, the key step is to explore the interior structure of the  $\epsilon$ - BRUE set and identify these simpler subsets.

In Section (4), we analyze the interior structure of the  $\epsilon$ - acceptable path set: as the indifference band gradually increases, more paths will begin to carry flows. Correspondingly, the  $\epsilon$ - BRUE set can be decomposed into subsets as well, with only acceptable paths carrying flows for each subset. Denote the  $k^{th}$  subset as  $\mathcal{F}_k^\epsilon$ , and  $K$  is the total number of largest  $\epsilon$ - acceptable path sets. Mathematically,  $\mathcal{F}_{BRUE}^\epsilon$  is the union of  $K + 1$  disjoint subsets:

$$
\begin{array}{rl} & {\mathcal{F}_i^\epsilon \bigcap \mathcal{F}_j^\epsilon = \emptyset ,i,j = 0,\dots ,K,i\neq j;}\\ & {\mathcal{F}_{BRUE}^\epsilon = \bigcup_{k = 0}^K\mathcal{F}_k^\epsilon .} \end{array} \tag{21b}
$$

According to the largest acceptable path set defined in Equation (10), the  $k^{th}$ $\epsilon$ - BRUE path flow subset is defined as:

$$
\begin{array}{rcl}{\mathcal{F}_0^\epsilon} & = & {\{\mathbf{f}\in \mathcal{F}_{BRUE}^\epsilon :a^\epsilon (\mathbf{f})\subseteq \mathcal{P}_{UE}\} ,}\\ {\mathcal{F}_k^\epsilon} & = & {\{\mathbf{f}\in \mathcal{F}_{BRUE}^\epsilon :\mathcal{P}_l^{\epsilon_{k^{-1}}}\subset a^\epsilon (\mathbf{f})\subseteq \mathcal{P}_l^{\epsilon_k}\} ,k = 1,\dots ,K} \end{array} \tag{22}
$$

where  $\mathcal{P}_l^{\epsilon_k^*}$  is the largest  $\epsilon_k^*$ - acceptable path set defined in (13).

# 5.1.  $\epsilon$ -BRUE path flow set for one OD pair

Equation (22) defines the  $k^{th}$  subset of the  $\epsilon$ - BRUE path flow set. In this section, we will explore how to construct each subset. Define a sequence of sets  $S_k^\epsilon , k = 0,\dots ,K$ , and assign all travel demands to paths from the associated largest path set, then we get:

$$
\begin{array}{rl}S_0^\epsilon \triangleq \{\mathbf{f}\in \mathcal{F}: & \forall i\in \mathcal{P}_{UE}:f_i,f_j\geqslant 0,|C_i(\mathbf{f}) - C_j(\mathbf{f})\leqslant \epsilon ;\\ & \forall i\notin \mathcal{P}_{UE}:f_i = 0\} . \end{array}
$$

$$
\begin{array}{rl}S_k^{\epsilon}\triangleq \{\mathbf{f}\in \mathcal{F}:\forall i\in \mathcal{P}_l^{\epsilon_{k - 1}}:f_i\geqslant 0;\exists i\in \mathcal{P}_l^{\epsilon_k^*}\backslash \mathcal{P}_l^{\epsilon_{k - 1}^*}:j_i > 0;\\ \forall i,j\in \mathcal{P}_l^{\epsilon_k^*}:|C_i(\mathbf{f}) - C_j(\mathbf{f})|\leqslant \epsilon ;\\ \forall i\notin \mathcal{P}_l^{\epsilon_k^*}:f_i = 0\} ,k = 1,\dots ,K. \end{array} \tag{23}
$$

where  $\mathcal{P}_l^{\epsilon_k^*}$  is defined in (13).

The set  $S_0^{\epsilon}$  contains all feasible path flow patterns where only the UE shortest paths carry flows and the cost difference between any two shortest paths are within the band  $\epsilon$ . For the set  $S_0^{\epsilon}, k = 1, \dots , K$ , the newly acceptable paths will carry flows, while those paths belonging to the  $\epsilon_{k - 1}^{\epsilon}$ - largest acceptable paths can carry flow or not, and the cost difference between any two acceptable or zero- acceptable paths are within the band  $\epsilon$ . By definition, each subset is disjoint, i.e.,  $\bigcap_{k = 0}^{K} S_{k}^{\epsilon} = \emptyset$ ; and it is a subset of the  $\epsilon$ - BRUE set, i.e.,  $S_{k}^{\epsilon} \subseteq \mathcal{F}_{BRUE}^{\epsilon}$ . The following proposition will show that  $\mathcal{F}_k^{\epsilon}$  and  $S_{k}^{\epsilon}$  are equivalent sets.

# Proposition 5.1.

$$
\mathcal{F}_k^{\epsilon} = S_k^{\epsilon},k = 0,\dots ,K. \tag{24}
$$

where  $\mathcal{F}_k^{\epsilon}$  is defined in (22),  $S_{k}^{\epsilon}$  is defined in (23).

Proof. (1)  $\mathcal{F}_k^{\epsilon} \subseteq S_k^{\epsilon}, k = 0, \dots , K$ .

$\forall \mathbf{f}\in \mathcal{F}_k^{\epsilon}$  , by definition,  $a^{\epsilon}(\mathbf{f})\subseteq \mathcal{P}^{\epsilon_{i}}$  , i.e.,  $\forall i,j\in a^{\epsilon}(\mathbf{f}),C_{i}(\mathbf{f})\leqslant \min_{j\in \mathcal{P}}C_{j}(\mathbf{f}) + \epsilon \leqslant C_{j}(\mathbf{f}) + \epsilon \Rightarrow C_{i}(\mathbf{f}) - C_{j}(\mathbf{f})\leqslant \epsilon$ $\epsilon$  Similarly,  $C_j(\mathbf{f}) - C_i(\mathbf{f})\leqslant \epsilon ,$  In summary,  $|C_{i}(\mathbf{f}) - C_{j}(\mathbf{f})|\leqslant \epsilon ,\forall i,j\in \mathcal{P}_{l}^{\epsilon_{k}^{*}}\Rightarrow \mathbf{f}\in \mathcal{S}_{0}^{\epsilon}$ $(2)S_{k}^{\epsilon}\subseteq \mathcal{F}_{k}^{\epsilon},k = 0,\dots ,K$  When  $k = 0$  ..

$\forall i\in a^{\epsilon}(\mathbf{f}),\mathbf{f}\in S_{0}^{\epsilon}$  , we need to show  $i\in \mathcal{P}_{UE}$  . Assume  $i\notin \mathcal{P}_{UE}$  , since  $\mathbf{f}\in S_0^{\epsilon}$  , then  $f_{i} = 0$  .On the other hand, because  $i\in a^{\epsilon}(\mathbf{f}),f_{i} > 0$  , contradicted with  $f_{i} = 0$  .Thus  $\forall i\in a^{\epsilon}(\mathbf{f}),i\in \mathcal{P}_{UE}$  , i.e.,  $a^{\epsilon}(\mathbf{f})\subseteq \mathcal{P}_{UE},\forall \mathbf{f}\in S_0^{\epsilon}\subseteq \mathcal{F}_{BRUE}^{\epsilon}\Rightarrow \mathbf{f}\in \mathcal{F}_{0}^{\epsilon}$  . So  $\mathcal{F}_0^{\epsilon}\subseteq \mathcal{S}_0^{\epsilon}$  . Combing with the result from (1),  $\mathcal{F}_0^{\epsilon} = S_0^{\epsilon}$

When  $k = 1$  , similarly,  $a^{\epsilon}(\mathbf{f})\subseteq \mathcal{P}^{\epsilon_{1}}$ $\forall \mathbf{f}\in S_1^{\epsilon}$  , and  $\mathcal{P}_{UE}\subset a^{\epsilon}(\mathbf{f})$  as  $\exists i\in \mathcal{P}_l^{\epsilon_i}\backslash \mathcal{P}_{UE}:f_i > 0$  , i.e., at least one newly added path needs to carry flow. Because  $S_{1}^{\epsilon}\subseteq \mathcal{F}_{BRUE}^{\epsilon}\backslash S_{0}^{\epsilon}$  , and  $\mathcal{F}_0^{\epsilon} = S_0^{\epsilon}$  , therefore  $\mathcal{P}_{UE}\subset a^{\epsilon}(\mathbf{f})\subseteq \mathcal{P}_l^{\epsilon_i}$ $\forall \mathbf{f}\in \mathcal{P}_{BUE}^{\epsilon}\backslash \mathcal{F}_0^{\epsilon}$  , so  $\mathbf{f}\in \mathcal{F}_1^{\epsilon}$

We can repeat this proof similarly for  $k = 2, \dots , K$ . Therefore  $S_{k}^{\epsilon} \subseteq \mathcal{F}_{k}^{\epsilon}, k = 0, \dots , K$ .

In conclusion,  $\mathcal{F}_k^{\epsilon} = S_k^{\epsilon}, k = 0, \dots , K$ .

Proposition (5.1) shows that by constructing each flow subset as in Equation (23), then it is equivalent to the definition of the subset in Equation (22). The union of these subsets constitutes the  $\epsilon$ - BRUE set:

# Corollary 3.

$$
\mathcal{F}_{BRUE}^{\epsilon} = \bigcup_{k = 1}^{K}S_{k}^{\epsilon}. \tag{25}
$$

where  $S_{j}^{\epsilon}$  is defined in Equation (23). In summary, Proposition (24) and Corollary (3) provide the methodology of constructing the  $\epsilon$ - BRUE set. The following example will illustrate this methodology.

Example 5.1. The topology of the test network, the OD demand between nodes 1- 4 and link cost functions are illustrated in Figure (5), and  $\epsilon = 15$ . Red lines display four paths: 1- 3- 4 (path 1), 1- 3- 2- 4 (path 2), 1- 2- 3- 4 (path 3), 1- 2- 4 (path 4). The equilibrium path flow pattern is [2, 2, 0, 2], i.e., path 1, 2 and 4 are utilized under UE. Substitute  $\mathcal{P}_l^{\epsilon_0^*} = \{1, 2, 4\}$ , path costs and the demand into MPEC (16), we obtain  $\epsilon_1^* = 6.5$ ,  $\mathbf{f} = [1.5, 3, 0, 1.5]$ ,  $C(\mathbf{f}) = [96.5, 103, 103, 96.5]$ . In other words, if  $\epsilon_1^* > 6.5$ , path 3 is utilized as well. Since  $\epsilon = 15 > \epsilon_1^*$ , we know  $\mathcal{P}^{\epsilon = 15} = \{\{1, 2, 4\} , \{1, 2, 3, 4\} \}$ .

![](images/4503c4ae69283bc0bea844c4e9877ca4e0ec8fedf0ab3f8ae75d77d1986f9085.jpg)  
Fig. 5: Single OD pair network illustration

Due to the flow conservation of the fixed demand, its BRUE solution set can be characterized by the first three paths. The whole BRUE solution set is shown in Figure (6), composed of a 3- path yellow subset and a 4- path magenta subset. Each subset satisfies Equation (24):

$$
\begin{array}{rcl}\mathcal{F}_0^{\epsilon = 15} & = & \{\mathbf{f}\in \mathcal{F}:\forall i\in \{1,2,4\} :f_i,f_j\geqslant 0,|C_i(\mathbf{f}) - C_j(\mathbf{f})|\leqslant 15,f_1 + f_2 + f_4 = 6;f_3 = 0\} ;\\ \mathcal{F}_1^{\epsilon = 15} & = & \{\mathbf{f}\in \mathcal{F}:\forall i\in \{1,2,4\} :f_i\geqslant 0;f_3 > 0;f_1 + f_2 + f_3 + f_4 = 6;\forall i,j\in \{1,2,3,4\} :|C_i(\mathbf{f}) - C_j(\mathbf{f})|\leqslant 15\} . \end{array}
$$

In  $\mathcal{F}^{\epsilon = 15}$ , path 3 does not carry flow, so the 3- path subset is on the bottom of the  $(f_1,f_2,f_3)$  coordinates; In  $\mathcal{F}^{\epsilon = 15}$ , path 3 begins to carry flow and  $f_3 > 0$ . Either the 3- path subset or the 4- path subset is convex, but their union is not convex, which is consistent with results in Lou et al. (2010).

![](images/1d210ed6f720031cee8b61922a1e38827d69871e1a773877d0d49e11ad939b97.jpg)  
Fig. 6: BRUE solution set illustration composed of two pieces

# 5.2.  $\epsilon$ -BRUE path flow set for multiple OD pairs

After knowing the methodology of obtaining all acceptable paths sets for multiple OD pairs in Section (4.2), it is not difficulty to generalize the methodology of constructing the  $\epsilon$ - BRUE set for a single OD pair to multiple OD pairs. For a network with multiple OD pairs, the  $\epsilon$ - BRUE set is the union of all subsets where demands are assigned to all acceptable paths across OD pairs:

$$
\mathcal{F}_{BRUE}^{\epsilon} = \bigcup_{k = 1}^{K}S_{k}^{\epsilon^{*}}, \tag{26}
$$

where,

$$
\begin{array}{rl}S_1^{\epsilon}\triangleq \{\mathbf{f}\in \mathcal{F}: & \forall i,j\in \mathcal{P}_{UE}:f_i^w,f_j^w\geqslant 0,|C_i^w (\mathbf{f}) - C_j^w (\mathbf{f})|\leqslant \epsilon^w;\\ & \forall i\notin \mathcal{P}_{UE}:f_i^{w} = 0,w\in W\} . \end{array}
$$

$$
\begin{array}{rl}S_k^{\epsilon}\triangleq \{\mathbf{f}\in \mathcal{F}: & \forall i\in \mathcal{P}_{l^k -2}^{\epsilon_{k - 2}}:f_i^w\geqslant 0;\exists i\in \mathcal{P}_{l^k -1}^{\epsilon_{k - 1}}\backslash \mathcal{P}_{l^k -2}^{\epsilon_{k - 2}}:f_i^w >0;\\ & \forall i,j\in \mathcal{P}_{l^k}^{\epsilon_{k - 1}}:|C_i^w (\mathbf{f}) - C_j^w (\mathbf{f})|\leqslant \epsilon^w;\\ & \forall i\notin \mathcal{P}_{l^k}^{\epsilon_{k - 1}}:f_i^w = 0,w\in W\} ,k = 2,\dots ,K. \end{array}
$$

where,  $K$  is the total number of acceptable path sets;  $\mathcal{P}_l^{\epsilon_k^*}$  is defined in Equation (17).

# 6. Conclusion

The concept of the boundedly rational user equilibrium (BRUE) was proposed in the 1980s, and the boundedly rational behaviour assumption was validated extensively in the 1990s by various empirical studies and traffic experiments. However, no mathematical properties of the BRUE have been studied thoroughly, due to its non- uniqueness and non- convexity. This paper fills the theoretical gap by solving the BRUE solution set and studying its mathematical properties. A nonlinear complementarity problem (NCP) is formulated to solve one BRUE path flow pattern. The NCP formulation indicates that under the boundedly rational behaviour assumption, the path cost is perturbed by an indifference function related to its corresponding path flow. Since the BRUE set is generally non- unique, how to construct this set is the main task of this paper.

Before constructing the BRUE flow set, its associated BRUE path set is explored first. All feasible paths for one BRUE flow pattern can be divided into three categories: acceptable, zero- acceptable and unacceptable. According to the monotonically non- decreasing property of the BRUE set, as the value of the indifference band increases, some paths which are not utilized will be taken, and thus the path set that contains the equilibrium flow, named the acceptable path set, will be augmented. The critical values of the indifference band to augment the path set can be identified by sequentially solving a class of mathematical programs with equilibrium constraints. After the acceptable path sets are obtained, the whole BRUE flow set can be decomposed into several subsets containing different numbers of acceptable paths. Each BRUE flow subset can be obtained when traffic demands are assigned to the corresponding acceptable paths.

In summary, this paper proposed a systematic methodology of obtaining BRUE solutions for networks with fixed demands connecting multiple OD pairs. It can help predict BR link flow patterns in a network, which guides planners to make network design decisions accordingly when travellers behave boundedly rational.

In the future, we would like to utilize the proposed methodology of constructing the BRUE set to study network design problems within the bounded rationality framework. The classical network design problem is usually formulated as a bilevel program: the upper level is the decision made to either enhance capacities of the established links, execute congestion pricing, or add new links to an existing road network; the lower level is a equilibrium problem, describing how travellers distribute among a new network topology. Due

to the existence of the indifference band, travellers may respond differently to one traffic design proposal, causing difficulty of predicting the BRUE link flow patterns and evaluating the efficiency of the proposal. Therefore a new network design framework needs to be established to accommodate bounded rationality route choice behavior.

# Acknowledgement

This research was funded by the National Science Foundation under grants EFRI- 1024604, EFRI- 024647, EFRI- 1024984, and CMMI- 09969600. The views are those of the authors alone.

# References

Avineri, E., Prashker, J., 2004. Violations of expected utility theory in route- choice stated preferences: Certainty effect and inflation of small probabilities. Transportation Research Record, 1894 (1), 222- 229. Bell, M., Iida, Y., 1997. Transportation network analysis. Ben- Akiva, M., Bergman, M., Daly, A., Ramaswamy, R., 1984. Modeling inter- urban route choice behaviour. In: Proceedings from the ninth international symposium on transportation and traffic theory, VNU Science Press, Utrecht, Netherlands. pp. 299- 330. Bovy, P., Fiorenzo- Catalano, S., 2007. Stochastic route choice set generation: behavioral and probabilistic foundations. Transportmetrica 3 (3), 173- 189. Cantillo, V., de Dios Ortúzar, J., Williams, H. C., 2007. Modeling discrete choices in the presence of inertia and serial correlation. Transportation Science, 41 (2), 95- 205. Cantillo, V., Heydecker, B., de Dios Ortúzar, J., 2006. A discrete choice model incorporating thresholds for perception in attribute values. Transportation Research Part B, 40 (9), 807- 825. Conlisk, J., 1996. Why bounded rationality? Journal of Economic literature, 34 (2), 669- 700. Cottle, R., Pang, J., Stone, R., 2009. The linear complementarity problem. No. 60. Society for Industrial Mathematics. Fonzone, A., Bell, M. G., 2010. Bounded rationality in hyperpath assignment: the locally rational traveller model. Gabaix, X., Laibson, D., Moloch, G., Weinberg, S., 2006. Costly information acquisition: Experimental analysis of a boundedly rational model. The American Economic Review, 1043- 1068. Gao, S., Freijinger, E., Ben- Akiva, M., 2011. Cognitive cost in route choice with real- time information: An exploratory analysis. Procedia- Social and Behavioral Sciences, 17, 136- 149. Gifford, J. L., Checherita, C., 2007. Bounded rationality and transportation behavior: lessons for public policy. In: Transportation Research Board 86th Annual Meeting. No. 07- 2451. Guo, X., Liu, H., 2011. Bounded rationality and irreversible network change. Transportation Research Part B, 45 (10), 1606- 1618. Han, Q., Timmermans, H., 2006. Interactive learning in transportation networks with uncertainty, bounded rationality, and strategic choice behavior: Quantal response model. Transportation Research Record, 1964 (1), 27- 34. Hu, T., Mahmassani, H., 1997. Day- to- day evolution of network flows under real- time information and reactive signal control. Transportation Research Part C, 5 (1), 51- 69. Jayakrishnan, R., Mahmassani, H., Hu, T., 1994a. An evaluation tool for advanced traffic information and management systems in urban networks. Transportation Research Part C, 2 (3), 129- 147. Jayakrishnan, R., Tsai, W., Prashker, J., Rajadhyaksha, S., 1994b. A faster path- based algorithm for traffic assignment. Transportation Research Record, 75- 75. Khisty, C. J., Arslan, T., 2005. Possibilities of steering the transportation planning process in the face of bounded rationality and unbounded uncertainty. Transportation Research Part C, 13 (2), 77- 92. Lou, Y., Yin, Y., Lawphongpanich, S., 2010. Robust congestion pricing under boundedly rational user equilibrium. Transportation Research Part B, 44 (1), 15- 28. Lu, C., Mahmassani, H., 2008. Modeling user responses to pricing: Simultaneous route and departure time network equilibrium with heterogeneous users. Transportation Research Record, 2085 (1), 124- 135. Mahmassani, H., Chang, G., 1987. On boundedly rational user equilibrium in transportation systems. Transportation Science, 21 (2), 89- 99. Mahmassani, H., Jayakrishnan, R., 1991. System performance and user response under real- time information in a congested traffic corridor. Transportation Research Part A, 25 (5), 293- 307. Mahmassani, H., Liu, Y., 1999. Dynamics of commuting decision behaviour under advanced traveller information systems. Transportation Research Part C, 7 (2- 3), 91- 107. Mahmassani, H., Zhou, X., Lu, C., 2005. Toll pricing and heterogeneous users: approximation algorithms for finding bicriterion time- dependent efficient paths in large- scale traffic networks. Transportation Research Record, 1923 (1), 28- 36. Marsden, G., Frick, K. T., May, A. D., Deakin, E., 2012. Bounded rationality in policy learning amongst cities: lessons from the transport sector. Environment and Planning A, 44 (4), 905- 920. Morikawa, T., Miwa, T., Kurauchi, S., Yamamoto, T., Kobayashi, K., 2005. Driver's route choice behavior and its implications on network simulation and traffic assignment. Simulation Approaches in Transportation Analysis, 341- 369. Nakayama, S., Kitamura, R., Fujii, S., 2001. Drivers' route choice rules and network behavior: Do drivers become rational and homogeneous through learning? Transportation Research Record, 1752 (1), 62- 68. Patriksson, M., 1994. The traffic assignment problem: models and methods. VSP Intl Science.

Rosenthal, R., Brooke, A., 2007. GAMS: A User's Guide. GAMS Development Corporation.Sivak, M., 2002. How common sense fails us on the road: Contribution of bounded rationality to the annual worldwide toll of one million traffic fatalities. Transportation Research Part F, 5 (4), 259- 269. Srinivasan, K., Mahmassani, H., 1999. Role of congestion and information in trip- makers' dynamic decision processes: Experimental investigation. Transportation Research Record, 1676 (1), 44- 52. Wang, J., Raith, A., Ehrgott, M., 2009. Telling analysis with bi- objective traffic assignment. In: Multiple Criteria Decision Making for Sustainable Energy and Transportation Systems: Proceedings of the 19th International Conference on Multiple Criteria Decision Making, Auckland, New Zealand, 7th- 12th January 2008. Vol. 634. Springer, p. 117. Wardrop, J., 1952. Some theoretical aspects of road traffic research. Proceedings of the Institution of Civil Engineers, 1, 325- 378. Zhu, S., 2011. The roads taken: Theory and evidence on route choice in the wake of the I- 35w mississippi river bridge collapse and reconstruction. Ph.D. thesis, University of Minnesota.