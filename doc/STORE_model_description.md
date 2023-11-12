# Introduction
This document describes the model. We are modeling  the interaction between *Prochlorococcus*, an autotrophic phytoplankton, and co-occurring marine heterotrophic bacterias.

The model represents the interaction between the autotroph and heterotroph through the lens of nutrien
# Model variables
All variables are in $\mu M$. 

*i* is in *Prochlorococcus* and Heterotroph.

* $B_i$ - biomass
* $N_i$ - N store
* $C_i$ - C store
* $DON$ - dissolved organic N
* $RDON$ - recalcitrant dissolved organic N
* $DIN$ - dissolved inorganic N
* $DOC$ - dissolved organic C
* $RDOC$ - recalcitrant dissolved organic C
* $DIC$ - dissolved inorganic C
* $ROS$ - reactive oxygen species

# Model equations
## Death and losses
Loss processes (death and leakage) is modeled linearly as a precentage of the biomass/stores.
We use exponential decay – in ISMEJ we show that other formulations are better for co-cultures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.

$deathB_i = M_i * B_i$

$deathN_i = M_i * N_i$

$deathC_i = M_i * C_i$

$B_i = B_i - deathB_i$

$C_i = N_i - deathN_i$

$N_i = N_i - deathC_i$


## Gross Uptake
### Monod
We are using monod limits to model the uptake dynamics based on resource availability

$limit_{ij} = R_j / (R_j + Kn_{ij})$

Where $R_j$ is one of the resources ($DOC$, $DON$, $DIC$, $DIN$),
and $Kn_{ij}$ is the $Kn$ affinity for of organism *i* for resource *j*.


### Droop

Uptake is regulated by not letting the stores grow too large. Using Droop-like equations.


$QC_i = (C_i + B_i * C2N_i) / (N_i + B_i)$

$regC_i = 1 - (QC_i / QCmax_i)$

$regN_i = 1 - (QCmin_i / QC_i)$

Where $C2N_i$ is a parameter describing the C/N ratio in the biomass of organism *i*.
$QCmax_i$ and $QCmin_i$ are the minimum and maximum C/N ratios for organism *i*.
Both $regC_i$ and $regN_i$ are clipped to be between 0 and 1.
$QC$ is the C/N ratio of organism *i*.

**Question: is DIC uptake by photosynthesis regulated? **


### ROS
ROS penalty is imposed in the ROS model. 

$ROSpenalty_i = e^{- \omega_{i}*ROS}$


### Gross uptake 
$grossUptake_{ij} = Vmax_{ij} * B_{i} * limit_{ij}  * reg_{ij} * ROSpenalty_{i}$

$uptakeN_i = grossUptake_{i,DIN} +  grossUptake_{i,DON}$

$uptakeC_i = grossUptake_{i,DIC} +  grossUptake_{i,DOC}$


## Net uptake 

$biosynthesisN_i = min(N_i + uptakeN_i, (C_i + uptakeC_i) / C2N_i) * Kmtb_i$


Respiration – growth associated bp/bh and maintenance associated r0p/r0h.
Modeled as: b * growth + r0 * biomass

$respirationC_i = (b_i * biosynthesisN_i + B_i * r0_i) * C2N_i$


If C store is not big enough, break some of the biomass into the stores. Make sure $C_i$ is not negative

$biomassbreakdownC_i = max(0, respirationC_i +  biosynthesisN_i * C2N_i - C_i - uptakeC_i)$


The store change is modeled as uptake minus biosynthesis and respiration.

$netDeltaN_i = uptakeN_i + biomassbreakdownC_i / C2N_i - biosynthesisN_i $

$netDeltaC_i = uptakeC_i + biomassbreakdownC_i - biosynthesisN_i * C2N_i  - respirationC_i$


## Overflow

In overflow mode, the uptake that is over the C2N ratio is released as overflow.
In photosynthetic organisms this includes photosyntheate that cannnot be utilized due to N limitation.
In the heterotrophs, under C limitation, carbon is extracted from organic compouds and the N containing waste (e.g. NH4) is released.
Make the store maintain the C:N ratio and exude the rest.

$storekeepN = min(netDeltaN, netDeltaC_i / C2N_i) $

$overflowN_i = netDeltaN_i - storekeepN_i$

$overflowC_i = netDeltaC_i - storekeepN_i * C2N_i$

## ROS
ROS is modeled in the ROS model. ROS is produced by both organisms and is toxic, limiting growth.
ROS is an unstable compounds and decays over time.

$ROSdecay = ROS * ROSdecayRate$

$netROS = ROS - ROSdecay$

ROS production depends on the biomass

$ROSrelease_i = E_{i,ROS} * B_i$

ROS breakdown is a common good that may be part of the positive interaction between *Prochlorococcus* and the heterotrophs.
*Prochlorococcus* cannot break down ROS, and may benefit from the breakdown performed by the heterotroph.

$ROSbreakdown_H = Vmax_{H,ROS} * B_H * netROS / (netROS + Kn_{H,ROS})$

$dROS/dt = \sum{ROSrelease_i} - ROSdecay - ROSbreakdown_H$

## DON breakdown due to exoenzymes. 
In the exoenzyme model, DON is degraded to DIN by exoenzymes released by the bacteria. This DIN is then available to both bacteria for uptake.
We did not model spontaneous breakdown of DON, only exoenzyme mediated. Also there is no cost for the bacteria to release these enzymes.

$DON2DIN_i = DON2DINrate_i * B_i * DON$


## Air water exchange
$DICAirWaterExchange   = - (DIC - csat) / airWaterExchangeConstant$

Where $csat$ is the saturated DIC concentration. and $airWaterExchangeConstant$ is a constant computed as:

$airWaterExchangeConstant = (h / (Kg * B * 0.01))$

Where $h$ is the height of the media in meters,  $Kg$ is the exchange rate in m sec-1 
and $B$ is the Revelle buffer factor. 

## C differential equations

$totaldeathC_i = deathB_i * C2N_i + deathC_i$

$dC_i/dt = netDeltaC_i - overflowC_i - deathC_i$

$dDIC/dt = DICAirWaterExchange + \sum{(respirationC_i - grossUptake_{i,DIC})}$

$dDOC/dt = \sum{(totaldeathC_i * \gamma_{i} + overflowC_i - grossUptake_{i,DOC})}$

$dRDOC/dt = \sum{totaldeathC_i * (1 - \gamma_{i})}$

Where $\gamma_{i}$ is the portion of the losses that is recalcitrant. 


## N differential equations

$totaldeathN_i = deathB_i + deathN_i$

$dB_i/dt = biosynthesisN_i - biomassbreakdownC_i / C2N_i - deathB_i  $

$dN_i/dt = netDeltaN_i - overflowN_i - deathN_i$

$dDIN/dt = \sum{(overflowN_i + DON2DIN_i - grossUptake_{i,DIN})}$

$dDON/dt = \sum{(totaldeathN_i * \gamma_{i} - grossUptake_{i,DON} - DON2DIN_i)}$

$dRDON/dt = \sum{totaldeathN_i * (1 - \gamma_{i})}$


Assuming that recalcitrant DON is released only during mortality/leakiness.
Assuming RDON/RDOC is recalcitrant to both organisms.






