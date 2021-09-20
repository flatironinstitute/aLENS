# Introduction

The protein property implemented in this package is the two-stage explicit KMC model.
The reference for this model is:  
**Blackwell, R. et al. Physical determinants of bipolar mitotic spindle assembly and stability in fission yeast. Science Advances 3, e1601603 (2017).**  
Details on the physical protein model can be found in **Section 1.3 of the Supplementary Material** of this paper.
Note that there are many typos, errors, and inconsistencies in that document.
In this markdown file, the notations are cleaned up, some symbols are changed, and the expressions to compute the four rates/probabilities are listed here.

This document explains the coding convention and program behaviors, corresponding to the physical model.

# Units

All calculations in the code are based on $\mu m, pN, s$ units.

- $k_BT=0.00411 pN\cdot \mu m$ at room temperature 300K.
- MT diameter $D=0.025 \mu m = 25 nm$.
- Water viscosity $\eta = 0.00089 pN\cdot s \cdot \mu m^{-2} = 0.00089 Pa \cdot s$.
- Stress unit $1 Pa = 1 pN\cdot \mu m^{-2}$.

# The struct `ProteinType`

This is the struct holding the intrinsic properties of a protein. Different `ProteinType`s are initialized from the configuration file `ProteinConfig.yaml`.
Currently, the `ProteinType` does not change over time during the simulation.
This can be made changeable in the future if necessary.
The code performs no check on the configuration.
It is the user's duty to provide meaningful settings for `ProteinType`.

If `fixedEnd0: true` is set, the `freeNumber: 0` must be set, otherwise the code may crash, because a free protein cannot have a fixed and bound end.
However, if `fixedEnd0: false` is set, `fixedLocationPerMT` can be arbitrarily set, which only means the initial position of each motor protein.

Variables in `ProteinType` are set by `ProteinConfig.hpp/cpp` routines by reading a `yaml` file named `ProteinConfig.yaml`.
The per head and per protein properties are directly set as POD types.
The pointer to the lookup table `LUTablePtr` is set and refreshed after the MPI particle exchange stage and before the KMC stage within each timestep, according to the `tag` of each `ProteinType`, to make sure valid LUTables are used.
The `LUTablePtr`s are merely users of those LookupTables.

## The list of protein properties.

- **Lower case letters are reserved for protein properties, unless there is no such property for filaments**
- **Upper case letters are reserved for filament properties**
- **In real code, all units are converted to $\mu m$, $pN$, and $s$**

### **Per head (end) property**

| Symbol     | Name                                                | Number       | Dimension                                                                               | Kinesin-5                                     | `ProteinType` member variable |
| ---------- | --------------------------------------------------- | ------------ | --------------------------------------------------------------------------------------- | --------------------------------------------- | ----------------------------- |
| $K_a$      | Association constant                                | $10\sim 100$ | $(\mu M)^{-1}$ per binding site                                                         | $90.9$ (To be checked)                        | `Ka[2]` for two heads         |
| $K_e$      | Association constant                                | $10\sim 100$ | $(\mu M)^{-1}$ per binding site if `useBindVol=True`, otherwise number per binding site | same as $K_a$ by default if `useBindVol=True` | `Ke[2]` for two heads         |
| $k_{o}^s$  | Singly bound turnover rate                          | $O(0.1)$     | $s^{-1}$                                                                                | $0.11$                                        | `ko_s[2]` for two heads       |
| $k_{o}^d$  | Doubly bound turnover rate                          | $O(0.1)$     | $s^{-1}$                                                                                | same as $k_o^s$                               | `ko_d[2]` for two heads       |
| $v_m$      | Max velocity                                        | $10\sim1000$ | $nm/s$                                                                                  | $-50$ on polar aligned MT pair                | `vmax[2]` for two heads       |
| $v_{m,AP}$ | Max velocity when doubly bound to anti-parallel MTs | $10\sim1000$ | $nm/s$                                                                                  | $-8$ on polar aligned MT pair                 | `vmaxAP[2]` for two heads     |
| $d^s$      | diffusivity when singly bound to MT                 | $10^{-3}$    | $um^2/s$                                                                                | 0                                             | `diffBoundS[2]` for two heads |
| $d^d$      | diffusivity when doubly bound to MT                 | $10^{-3}$    | $um^2/s$                                                                                | 0                                             | `diffBoundD[2]` for two heads |

### **Per protein property**

| Symbol     | Name                                    | Number          | Dimension               | Kinesin-5                  | `ProteinType` member variable |
| ---------- | --------------------------------------- | --------------- | ----------------------- | -------------------------- | ----------------------------- |
| $\epsilon$ | Effective binding site density along MT | $125\sim 1650$  | site $\cdot \mu m^{-1}$ | $1625$                     | `eps`                         |
| $\kappa$   | Spring constant                         | $O(100)$        | $pN/\mu m$              | $300$                      | `kappa`                       |
| $\lambda$  | Unbinding load sensitivity              | $O(0 \sim 0.4)$ | dimensionless           | $0.25822$ (To be verified) | `lambda`                      |
| $\ell_0$   | Rest length                             | $O(50)$         | $nm$                    | $53$                       | `freeLength`                  |
| $r_c$      | Capture radius                          | **Note**        | $\mu m$                 | **Note**                   | `rc`                          |
| $f_s$      | Stall force                             | $O(1 \sim 5)$   | $pN$                    | $5$                        | `fstall`                      |
| $d_0$      | Diffusivity unbound                     | $O(1)$          | $\mu m^2 /s$            | $4.5$                      | `diffUnbound`                 |

**Note**: The singly bound stage capture radius $r_c$ can be set in two different ways in ProteinConfig.yaml;

- $r_c = D/2 = 0.0125$. This makes the singly bound stage behavior the same as the old model.
- $r_c = (D + \ell_0)/2$. This makes the singly bound stage behavior follow the same convention as the new energy exponential factor including $-D$ in the doubly bound stage.

**Note 2**: The singly bound stage capture radius $r_c$ from the config yaml file is not used, if the protein diffusion length scale within one step is larger than the $r_c$ set by the user. In this case, $r_c=\sqrt{6d_0\Delta t}$, where $\Delta t$ is the timestep.

### **Protein behavior flags.**

- Flag for one end fixed, `bool fixedEnd0`.
- Flag for walking off the ends, `bool walkOff`.
- Flag for using bind volume in the $S\rightarrow D$ step, `bool useBindVol`.
- Flag for choosing how to calculate protein stretching length in the $S\rightarrow D$ step, `lookupType=` 1 or 2.

When `fixedEnd0==true`, the code is implemented as:

```c++
    if (type.fixedEnd0) {
        // end0 is fixed, set probability of changing binding status to zero
        type.ko_s[0] = 0;
        type.ko_d[0] = 0;
    }
```

# The struct `ProteinBindStatus`

This is the struct holding the position and binding (dynamic) information of a protein.
Each protein is tracked as a single point if at unbound or singly bound state (U or S), but tracked as a spring with two ends if at doubly bound state (D).

**Important** This is the only data that is updated by the KMC calculations.

Therefore, the datafields `double pos[3]` and `double posEndBind[2][3]` should be carefully handled.
More on this in the `ProteinData` section.

The definitions of default values are:

```cpp
constexpr int ID_UB = -1;
constexpr double NAND = std::numeric_limits<double>::quiet_NaN();
```

Here is a list of `ProteinBindStatus` data fields.  
| variable | default value | description |
| ---------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `double pos[3]` | `{NAND,NAND,NAND}` | The (lab frame) center position of this protein. Always valid |
| `bool changeBind[2]` | `{true, true}` | If `false`, update binding information only and bypass the KMC step. |
| `int idBind[2]` | `{ID_UB, ID_UB}` | The gid of each bound MT. Set to `ID_UB` if unbound |
| `int rankBind[2]` | `{ID_UB, ID_UB}` | The MPI rank of each bound MT. $\in[0,nProcs-1]$, otherwise set to `ID_UB`. |
| `double lenBind[2]` | `{NAND, NAND}` | The length of each bound MT. |
| `double distBind[2]` | `{NAND, NAND}` | The distance to the center of each bound MT. $\in [-lenBind/2,lenBind/2]$, positive is towards the plus end. |
| `double posEndBind[2][3]` | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (lab frame) location of each protein end (head) when bound. |
| `double centerBind[2][3]` | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (lab frame) location of the center of each bound MT |
| `double directionBind[2][3]` | `{{NAND,NAND,NAND},{NAND,NAND,NAND}}` | The (normalized) orientation vector of each bound MT |

# The Binding/Unbinding probabilities and rates in KMC

The U-S-D process follows Poisson process.
Our algorithm assumes that multiple binding/unbinding events do not occur in
the same time step $\delta t$. As $\delta t$ becomes large, this approximation
breaks down.

To find whether the current time step is appropriate we state that the
probability of two events occurring in a time $t$ must be

$$\max \{ P\left(C(t) \cup B(t')|A(0)\right) \} < \delta$$

for some small value $\delta$ where $A, B$, and $C$ are states at times $t$ and
$t'$. The $P(C(t) \cup B(t')|A(0)) = P(C(t)|B(t'))P(B(t')|A(0))$ and $P(B(t)|A(0))= 1 - \exp[-k_{A \rightarrow B} t]$ as it follows a Poisson process for a
single event occurring. The maximum probability for a double event occurring
happens when $t'=t'_{max}$ which is found by solving

$$\frac{dP(C(t) \cup B(t')|A(0))}{dt'}\Bigg|_{t'_{max}} = 0$$

$$k_{A\rightarrow B} \left( e^{k_{B\rightarrow C}(t-t')} - 1 \right) - k_{B\rightarrow C} \left( e^{k_{A\rightarrow B} t'} - 1 \right)  = 0.$$

When $t'_{max}$ is substituted into the total probability, you can determine
whether or not $t$ exceeds your tolerance for the $A\rightarrow B \rightarrow C$
process. While you can not analytically calculate $t'_{max}$, it is possible to
numerically compute this value. An iterative approach can also be used to calculate
the largest possible time step for a given $\delta$.

There are four unique processes that you must consider:
$U\rightarrow S \rightarrow U$, $U\rightarrow S \rightarrow D$,
$S\rightarrow D \rightarrow S$, and $D\rightarrow S \rightarrow U$.
The processes
$D\rightarrow S \rightarrow D =  S\rightarrow D \rightarrow S$ and
$S\rightarrow U \rightarrow S = U\rightarrow S \rightarrow U$ probabilistically
and, therefore, are not necessary to compute.

The KMC rates and probabilities depend on both properties and geometry.
In the following, the superscript $i=0,1$ denotes the two heads of a protein.
For the geometric dependent part, uppercase symbols are for MT, like $L,D$ for length and diameter.
Lowercase symbols are for protein, like $\ell,r_c$ for length and capture radius.

## Between Unbound and Singly Bound Stages

| Probability                              | Rate                                                                                                                                |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| $P_{on}^{s,i} = k_{on}^{s,i}\Delta t$    | $k_{on}^{s,i} = \left( \dfrac{L_c}{2r_c} \right) \left(\dfrac{3 K_a^i }{4\pi r_c^3} \right) \left( 2r_c\epsilon\right) k_{o}^{s,i}$ |
| $P_{off}^{s,i} = k_{off}^{s,i} \Delta t$ | $k_{off}^{s,i}=k_{o}^{s,i}$                                                                                                         |

$L_c$ is the length of MT within the capture-sphere with radius $r_c$ of protein.
$\Delta t$ is the discretized timestep size.

## Between Singly Bound and Doubly Bound Stages

During this step the key is to compute the protein stretch energy, i.e., the stretch length.
Since we track proteins by their bind locations along the center lines of bound filaments, one protein's length is from one bound center line location to another center line location.
We call this length `ProteinEndEndLength()`, this is simply the length $\ell$ between two `posEndBind[3]` along both filament **center line**s.
However, if we directly use this length subtracting the free length $\ell_0$ as the streching length, we over estimate the protein stretching because in reality binding happens between filament surfaces.
Therefore, we must compute another `ProteinForcceLength()`, based on `ProteinEndEndLength()` and some fixes, to compute the stretch length.
We implement two methods of computing this, controlled by the flag `lookupType`.
For a given singly bound geometry, we first find the projection of the bound head to the candidate filament center line, and call the distance from the bound head to this projection point as $d$.
Then, starting from this projection point, the candidate doubly bound site location along the center line is denoted by $s$.
The candidate `ProteinEndEndLength()`$=\sqrt{d^2+s^2}$.
We use $\ell_f(s)$ to denote `ProteinForceLength()`.

- `lookupType=0`, $\ell_f(s)=\sqrt{d^2+s^2}-D$. This is more general, but has larger error when $s$ is large, i.e., binding is oblique relative to the filament.
- `lookupType=1`, $\ell_f(s)=(1-D/d)\sqrt{d^2+s^2}$. This is accurate for parallel filaments, but is singular for some geometries such as $d\approx0$.

Note: $V_b$ is the statistical volume an unbound end explores when the other protein end is attached to a MT. This calculated by taking the integral of the Boltzmann factor over all possible lengths $r$ of the protein.

- `lookupType=0`, $V_b = \int \exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(r-\ell_0 - D\right)^2\right]dr^3$.
- `lookupType=1`, $V_b = \int \exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(r-\ell_0\right)^2\right]dr^3$.

| Probability                                            | Rate                                                                                                                                                                                                         |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| $\displaystyle P_{on}^{d,i} = k_{on}^{d,i} \Delta t$   | $\displaystyle k_{on}^{d,i} =  \left\{D\epsilon \left(\dfrac{K_e^i}{V_b} \right) \int\exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(\ell_f(s)-\ell_0\right)^2\right] \frac{1}{D} ds \right\}  k_{o}^{d,i}$ |
| $\displaystyle P_{off}^{d,i} = k_{off}^{d,i} \Delta t$ | $\displaystyle k_{off}^{d,i} = \exp\left[\frac{\lambda \kappa \beta}{2} \left(\ell_f(s)-\ell_0\right)^2\right] k_{o}^{d,i}$                                                                                     |

### Exceptions to binding volume interpretation

Experiments are capable of finding $K_a^i$ but finding $K_e^i$ can be difficult as geometric considerations become more important. However, _in vitro_ assays are able to find ratios between singly and doubly bound proteins between pairs of MTs. With this knowledge it is possible to calculate an effective $K_e^i$ if the numbers between singly and doubly bound proteins are comparable and both protein ends are the same e.g. Kinesin-5 or Ase1.

If `useBindVol==false`, the $\left(K_e^i/V_b\right)$ is replaced by simply $K_e^i$.
This also changes the units of $K_e^i$ to number per binding site.
The unit conversion is handled within the code.

You can approximate the new $K_e^i$ if you know the ratio of singly bound to
doubly bound motors or crosslinkers at steady state. First, the experimental
singly bound motor number per concentration of motors in solution is related to
our parameters by $n_s = K_a^i \epsilon L$. 
_In vitro_ experiments are
capable of determining the number of doubly bound crosslinks (CITE) that we
would calculate to be $n_d = \epsilon^2 K_a^i K_e^i q$, where

$$q = \int_{L}\exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(\ell_f(s)-\ell_0\right)^2\right] ds.$$

The ratio between singly and doubly bound motors then be expressed as $\rho = \frac{n_d}{n_s} = \frac{\epsilon K_e^i q}{L}$ (See "S$\rightarrow$D probability integral" for further details). Our effective association constant is then

$$K_e^i = \frac{\rho V_b L}{\epsilon q}.$$

For parallel or antiparallel configurations, $q$ is directly proportional to the length of the overlap between filaments. 
If $L \gg \sqrt{2/\beta \kappa}$ and $\ell_o$ is set to $-D$, the gaussian integral $q \approx L\sqrt{2\pi/\beta \kappa}e^{-\beta \kappa r_{\perp}^2}$ where $r_{\perp}$ is the center to center separation between filaments. 
In this simplified scenario, $K_e^i$ no longer depends on the the length of the filaments used.

$$K_e^i = \frac{\rho V_b}{\epsilon} \sqrt{\frac{\beta \kappa}{2\pi}}e^{\beta \kappa r_{\perp}^2}$$

Note: This derivation only holds true for proteins with identical binding ends. (i.e Kinesin-5, PRC1)


## The (dimensionless) S$\rightarrow$D probability integral

$$\int_{s_A}^{s_B}\exp\left[-\frac{(1-\lambda)\kappa\beta}{2}\left(\ell_f(s)-\ell_0\right)^2\right] \frac{1}{D} ds$$

The bound head $\alpha$ of the protein is already fixed at a point somewhere in space (another MT).
We first compute the perpendicular foot from head $\alpha$ to the infinite centerline, where the distance is denoted as $\ell_{m}$, meaning the 'minimal distance' from $\alpha$ to the infinite centerline.
$d$ corresponds to the variable `distPerp` in the code.

We then put an axis towards the plus direction of the MT on the MT centerline, where $s=0$ marks the perpendicular foot corresponding to $d$.
The tubule is marked by bounds $s_A$ and $s_B$ on its infinite centerline.
$s_B-s_A=L$, the MT length.
Depending on the geometry, there are in total three cases:

- $0<s_A<s_B$
- $s_A<0<s_B$
- $s_A<s_B<0$

The integral is then non-dimensionalized with MT diameter $D$:

- $s' = s/D$
- $\ell' = \ell/D$

The reason why choosing $D$ as the length scale is because

- Different proteins with different $\ell_0$ exists in the system. Choosing $D$ gives a consistent basis for all proteins.
- Actual stretching of proteins happens on the length scale of $D$ (collision, alignment, etc.), instead of $\ell_0$.

The integral becomes

$$\int_{s_{A}'}^{s_{D}'}\exp\left[-M\left(\ell_f(s)'-\ell_0'\right)^2\right] ds'$$

A dimensionless parameter $M$ can be identified, as the binding free energy factor:

$$M=\frac{(1-\lambda)\kappa\beta D^2}{2}$$

In general, $M\in(0.1,10)$ for proteins used in this model.

This nondimensionalization avoids numerical integration of a sharp peak over length scale $L$.
Instead, the integration of length scale $D$ is soft and smooth, straightforward to be integrated by Gauss-Kronrod or other 1D quadrature schemes.

The behavior of the integrand depends on the minimal distance $d'$:

- $d' > 1+\ell_0'$: the integrand is single peaked at $s'=0$. The protein is always stretched.
- $d' < 1+\ell_0'$: the integrand is double peaked at $\ell_f'(\pm s')= \ell_0$. The protein can be compressed or stretched depending on $s$.

Due to symmetry about $s'=0$, the integral is tabulated as follows:

$$\int_0^{s_{\rm bound}'}\exp\left[-M\left(\ell_f'(s')-\ell_0'\right)^2\right] ds'$$
Where $M$ and $\ell_0'$ are fixed for each species. The integral is tabulated over $d'$ and $s_{\rm bound}'$ as a 2D regular table. The tabulation is truncated at the integrand is smaller than $eps$. This gives the bounds $s_{\rm bound}'$.
The table can be constructed as a $256\times256$ matrix. $256$ grids along $d'$ and $256$ grids along $s_{\rm bound}'$.
Lookup entries lower than the bounds generate an error, while values higher than the bounds return zero.

**Note** All non-dimensionalization is handled within the lookup table class. Values returned from lookup table are dimensional and need not be scaled again.

## Dimensionless groups for binding/unbinding process

Identifying the dimensionless groups helps understanding the effects of the numerical parameter values.
| Name | Kinesin-5 | Explanation |
| ----------------------------------------------------------------------------------------------------------------------------- | --------- | ------------------------------------------------------------------------------------------ |
| $\dfrac{3 K_a}{4\pi r_c^3}$ | ? | Dimensionless association constant in the singly bound stage |
| $2\epsilon r_c$ | ? | Total number of binding sites in capture radius |
| $\dfrac{L_c}{2r_c}$ | $[0,1]$ | Geometric capturing factor depending on the configuration of MP and MTs |
| $\left(\dfrac{3 K_a}{4\pi r_c^3} \right) \left( 2 \epsilon r_c \right)$ | ? | Ratio between on and off rate |
| $\left( \dfrac{L_c}{2 r_c} \right) \left(\dfrac{3 K_a^i }{4\pi r_c^3} \right) \left( 2r_c\epsilon\right) k_{o}^{s,i}\Delta t$ | $[0,1]$ | Scale of binding probability per timestep. Should be $<1$, otherwise $\Delta t$ too large. |
| $D\epsilon$ | ? | Dimensionless binding site density |
| $M=\frac{(1-\lambda)\kappa\beta D^2}{2}$ | ? | Scale of doubly bound unbinding rate. |

## Model time step limitations

This model assumes that it is unlikely for any two binding/unbinding processes to take place within a single time step \Delta t.

# The struct `ProteinData`

## Definition

This struct describes one protein, for all different types.
This is the data structure used in actual dynamic simulations.
Each protein has one unique non-negative integer global ID `gid`, a `ProteinProperty` struct, a `ProteinBindStatus` struct, a `forceBind[2][3]`, and a `torqueBind[2][3]`.

## Member functions

The member functions (besides simple `set` and `get` functions) are grouped into two categories:

- `updateXXX()` functions update the internal member variables. No return values (`void updateXXX()`).
- `calcXXX() const` functions calculate something with member variables and some given non-member information. Internel variables are kept constant, and the calculated values are returned.

## Initialization

If the initial protein configuration is not provided as an input file `ProteinInitial.txt`, initial configurations will be generated according to the configuration file `ProteinConfig.yaml`.
Each protein will be sequentially numbered and generated on mpi rank 0 after the maximum `gid` of sylinders.
The position of each protein, `pos[3]`, is either generated from uniform random numbers inside the initBox, or set according to tubule configuration if `fixedEnd0: true` is set in `ProteinConfig.yaml`.

If the protein initial configuration file `ProteinInitial.txt` is provided, either from a history file of the previous simulation, or hand-written, the `gid` will be read from the configuration file.
_No check_ on this `gid` is performed.
It is the user's job to provide meaningful data.

## Protein position in simulation

There are two scenarios where the position of a protein is read or set: (1) inside FDPS routines and (2) outside FDPS routines.

Inside FDPS routines, the `ProteinData` struct appears in two cases: as a `FDPS::FullParticle` or a `FDPS::EssentialParticleI`.

- When `ProteinData` is used as a `FDPS::FullParticle`, `void setPos(const PS::F64vec3 &)` is called **ONLY** inside the `adjustPositionIntoRootDomain()` function, to set the position according to the periodic boundary condition along a certain direction of a box.
- When `ProteinData` is used as a `FDPS::EssentialParticleI`, `void setPos(const PS::F64vec3 &)` is **NEVER** called. This conclusion also applies to the case when `ProteinData` is used as `EPI` in `MixedPairInteration`. In this case, `PS::F64vec3 getPos() const` is called, and must return a meaningful position in the (possibly periodic) simulation box `RootDomain`, otherwise `FDPS` may work abnormally and is hard to detect at runtime.

Therefore, when `void setPos(const PS::F64vec3 &newPos)` is called, the jump (`newPos-pos`) must represent periodic jumps across periodic images of the simulation box.
The `double centerBind[2][3]`, if not `NAND`, should also be moved as the same jump `newPos-pos`.
The implicit assumption is that there is always a MT across the periodic image jump there, so the `double centerBind[2][3]` remains always valid after the jump.
This is necessary to make FDPS work.

Outside FDPS routines, the location is accessed by overloaded functions `const double * getPosPtr() const` and `void setPos(const pos[3])`.
These functions get and set the `pos` of the protein, with the same behavior as the FDPS interface versions.

## Protein cycle within one timestep

At each time step, the following events happen **in the given order**.

| Protein                    | MT                   | Comment                                                              |
| -------------------------- | -------------------- | -------------------------------------------------------------------- |
|                            | `prepareStep()`      | partition, exchange, update map, etc                                 |
| apply PBC                  |                      | `FDPS::adjustPositionIntoRootDomain()`                               |
| mpi exchange               |                      | `FDPS::exchangeParticle()` accoridng to `domainInfo` of MT partition |
| run KMC step               |                      | `calcBindInteraction()`, which calls `MixedPairInteraction`.         |
| `updateGeometryWithBind()` |                      |                                                                      |
| walk or diffuse            |                      | if walk, check `walkOff`                                             |
| `updateForceTorqueBind()`  |                      |                                                                      |
| save protein info to file  |                      | save to XML vtk file                                                 |
| apply force & torque to MT | `setForceNonBrown()` | `ZDD` and `CommMPI` classes are involved                             |
|                            | `runStep()`          | MT history file is saved inside `runStep()` before MT movement       |

**Note**

- When `walkOff == true`, protein unbinds when it walks past the ends of bound MT. If two ends simultaneously walk past the minus ends, one end is chosen (with equal probability for both ends) to walk off and the other end is clamped to the end and is still bound. This protein then continues the normal KMC cycle in the next timestep

<!-- ### Pre-KMC
The property `fixedEnd0` is handled here. If `true`, the `changeBind[0]=false` is set to bypass the KMC step for this end.
In this case, the protein switches between singly and doubly bound states only.
In the `CalcProteinBind()` step, the information of the fixed MT is updated because MTs are moving.
### KMC (near-neighbor)
There are two cases: `changeBind[2]={false, true}` or `changeBind[2]={true,true}`. In the first case, `KMC_1_u` or `KMC_2_1` is executed depending on the binding status. In the second case, `KMC_0_1`, `KMC_1_u`, or `KMC_2_1` is executed depending on the binding status.
### Post-KMC
After this check, the protein motion is moved according to given diffusivity and walking velocity.
Finally, the protein force is calculated and added to the tubules by using the `ZDD` and `CommMPI` classes. The protein binding stress is also calculated here. -->

# The struct `ProteinConfig`

This is a class reading in a `yaml` file to specify protein properties and initial configuration.
An example config `yaml` file is included, named `ProteinConfig_example.yaml`.

The file specified all protein types as a `yaml` list:

```yaml
proteins:
  - ....# specify protein type 0
  - ....# specify protein type 1
  - ....# specify protein type 2
  - ....
```

Currently there is no limit of how many types of proteins are specified.

Each type of protein includes a property section:

```yaml
tag: 0 # Type 0, some randomly set meaningful parameters
#properties:
walkOff: true
bindAntiParallel: false
fixedEnd0: false
freeLength: 0.08 # um
kappa: 1.0 # pN/um
fstall: 1.0 # pN
lambda: 0.5 # dimensionless
diffUnbound: 5.0 # um^2/s
vmax: [1.0, 1.0] # um/s, positive towards plus end
vmaxAP: [1.0, 1.0] # um/s, positive towards plus end
# KMC parameters
useBindVol: true
eps: 1625 # 13 protofilaments at 8nm per block
Ka: [53.05, 10.0] # (uM)^{-1}
ko_s: [0.003916, 0.001958] # 1/s
Ke: [53.05, 10.0] # (uM)^{-1}, same as Ka
ko_d: [0.003916, 0.001958] # 1/s
```

and an initial configuration section:

```yaml
freeNumber: 0
fixedLocationPerMT: [-2, 0, 1] # given in [-1,1]. otherwise random
```

The order of properties and configurations specified in this file does not matter.
Most properties are self-explainatory.
A single parameter refers to the property per protein, and an array of two parameters refers to the properties of each head (end).
These are read and stored in the `ProteinProperty` struct and every protein has an independent copy of all the properties.

The special property `tag` is an arbitrary integer given by the user to mark the types. The tags do not have to be sequentially ordered.
For example, the tags can be specified like this:

```yaml
proteins:
  - tag: 5 # for Kinesin-5
    .....  # other Kinesin-5 properties
  - tag: 1 # for Kinesin-1
    .....  # other Kinesin-1 properties
```

The initial configuration section specifies specifies the free and bound proteins in the beginning.
`freeNumber` denotes the number of proteins that are uniformly and randomly distributed in the simulation box, regardless of whether periodic boundary conditions are set.
`fixedLocationPerMT` is an array and denotes the number and locations of proteins singly bound to each MT in the beginning.
An entry in this array specifies the location.
A number $\in[-1,1]$ specifies the location, from minus end (-1) to positive end (1).
A number out of that range (e.g. 2, -2, 3, etc) specifies that a random location is chosen.
For example:

- `fixedLocationPerMT: []` means no initially bound protein.
- `fixedLocationPerMT: [-3,-2,0,1]` means 4 proteins (of this type) are bound to each MT in the beginning.
  One at the center (0), one at the positive end (1), two are randomly chosen (-3,-2).

# Protein IO

## Data output

Output frequency is controlled by

```yaml
dt: 0.001 # s
timeTotal: 0.1 # s
timeSnap: 0.001 # s
```

in `RunConfig.yaml` file.
This example will generate 1 (`timeSnap` divide by `dt`) snapshots per timestep, in total 100 (`timeTotal` divide by `timeSnap`) snapshots.

Each snapshot includes one ASCII file (`.dat`) and a group of binary (base64) XML VTK file (`.pvtp` and `.vtp`).

### ProteinAscii_NUMBER.dat

This is the ASCII file to record protein geomery and bound MT ID only.
`NUMBER` is sequentially ordered from 0 to mark the sequence of files.
All proteins from all mpi ranks are combined into a single file.
Each protein takes a line, and is exported as

```c
        if (bind.idBind[0] != ID_UB && bind.idBind[1] != ID_UB) {
            // protein has finite length
            // this should NOT out put nan
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.posEndBind[0][0], bind.posEndBind[0][1],
                    bind.posEndBind[0][2], //
                    bind.posEndBind[1][0], bind.posEndBind[1][1],
                    bind.posEndBind[1][2], //
                    bind.idBind[0], bind.idBind[1]);
        } else {
            // protein has zero length
            fprintf(fptr, "P %d %d %.6g %.6g %.6g %.6g %.6g %.6g %d %d\n", //
                    gid, property.tag,                                     //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.pos[0], bind.pos[1], bind.pos[2],                 //
                    bind.idBind[0], bind.idBind[1]);
        }
```

**NOTE** The order of protein `gid` appears in this dat file is **NOT** specified and is **NOT** fixed.

### Protein_NUMBER.pvtp and Protein_rX_NUMBER.vtp

Each mpi rank writes its own partition of protein data. `rX` marks the mpi rank. The data can be read and visualized by `Paraview` using the included `TubuleVis.pvsm` file for `Paraview` higher than version 5.6.

The explanation of binary (base64) XML VTK file cna be found on VTK official website.

## Data input

The file `ProteinAscii_NUMBER.dat` can be copied to `ProteinInitial.dat` to be read by the main executable program to set the initial configuration and binding status of proteins.
Binding information will be reconstructed with proper `TubuleInitial.dat`, `RunConfig.yaml`, and `ProteinConfig.yaml`.

**Note**

- The protein tags in `ProteinConfig.yaml` must include all tags appearing in `ProteinInitial.dat`. Properties of the same tag can be different from a previous run.
- The `TubuleInitial.dat` must match the `ProteinInitial.dat` if any protein is in S or D bound stage.
- The simbox geometry in `RunConfig.yaml` must match if any protein is in S or D bound stage.
