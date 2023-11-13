# Introduction
This document describes the model. We are modeling  the interaction between *Prochlorococcus*, an autotrophic phytoplankton, and co-occurring marine heterotrophic bacterias.

The model represents the interaction between the autotroph and heterotroph through the lens of nutrients exchanges and detoxification.

![image](https://github.com/wosnat/recycle_model/assets/22752755/cf814c98-3190-48fc-801d-5aba284048e7)

# Model variables
All variables are in $\mu M$. 

In the equation below, *i* is the subscript for population, and can be either *Prochlorococcus* or Heterotroph.


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

# Main model overview
In our model, C and N are taken up from dissolved organic and/or inorganic forms into C and N stores ($C_i$ and $N_i$ respectively), which are then combined through biosynthesis into functional biomass  ($B_i$, figure 1). Functional biomass is either degraded back into stores to support respiration or lost through mortality (Eq 1)

The functional biomass has a fixed C/N ratio ($C2N_i$). $B_i$, the biomass variable is maintained in $\mu M  N$, representing also the corresponding C biomass of $C2N_i * B_i$.

The stores are utilized for biosynthesis, respiration (C store only) and are degraded into dessolved organic matter upon death (Eq2, Eq3). When the stores are imbalanced, some of them may be exuded as overflow (see below).

Eq1:    $dB_i/dt = biosynthesisN_i - biomassbreakdownC_i / C2N_i - deathB_i  $

Eq2:    $dN_i/dt = uptakeN_i + biomassbreakdownC_i / C2N_i - biosynthesisN_i - overflowN_i - deathN_i$

Eq3:    $dC_i/dt = uptakeC_i + biomassbreakdownC_i - biosynthesisN_i * C2N_i  - respirationC_i - overflowC_i - deathC_i$

<br />

The synthesis of functional biomass ($\mu M  N/sec$)) is determined by the limiting store (C or N in $\mu M$) times a biosynthesis rate ($1/sec$) (Eq4). The C and N stores available for biosynthesis include storage and uptake.

Eq4:   $biosynthesisN_i = min(N_i + uptakeN_i, (C_i + uptakeC_i) / C2N_i) * Kmtb_i$


Respiration is comprised of growth associated respiration controlled by $b_i$ and maintenance associated respiration controlled by $r0$ (Eq5).

Eq5:   $respirationC_i = (b_i * biosynthesisN_i + B_i * r0_i) * C2N_i$


The breakdown of biomass does not need occur if the cell is in balanced growth (i.e. can be equal to zero). However, biomass may be broken down to support respiration. Respiration is assumed to occur from the C storage (e.g. carbohydrates or lipids), but if there is not enough C storage biomass will be degraded instead (Eq7). Note that equation 6 is calculated after the calculation of biosynthesis but before the calculation of the C store, and hence the needs to incorporate also the C used for biosynthesis (Eq6,7).

Eq6: $requiredCStore_i = C_i + uptakeC_i - biosynthesisN_i * C2N_i $

Eq7: $biomassbreakdownC_i = max(0, respirationC_i - requiredCStore_i)$

# Model variants

There is a baseline model and additional 4 model variants encoding different interaction mechanisms between the 2 populations (Exoenzymes, Excess Overflow, Mixotrophy, ROS detoxyfication). 

### Baseline model
In the baseline model, the heterotroph consumes organic N and C as well is inorganic N. While the autotroph consumes inorganic C (photosynthesis) and inorganic N. The only processes in the base model are nutrients uptake, biosynthesis, respiration and death.

### Excess Overflow

In the overflow model, the uptake that is over the C2N ratio is released as overflow.
In photosynthetic organisms this includes photosyntheate that cannnot be utilized due to N limitation.
In the heterotrophs, under C limitation, carbon is extracted from organic compouds and the N containing waste (e.g. NH4) is released.
The model makes the store maintain the C:N ratio and exude the rest.

### Mixotrophy
In the mixotrophy model, *Prochlorococcus* can feed on organic compounds as well as inorganic compounds (a mixotroph). This is implemented by a non-zero $Vmax_{i,DOC}$, $Vmax_{i,DON}$ for *Prochlorococcus*. 

### Exoenzymes
In the exoenzyme model, DON is degraded to DIN by exoenzymes released by the bacteria. This DIN is then available to both bacteria for uptake.
We did not model spontaneous breakdown of DON, only exoenzyme mediated. Also there is no cost for the bacteria to release these enzymes.

### ROS detoxyfication
ROS (reactive oxygen species) is a natural byproduct of core metabolic processes as well as photosynthesis.
It is toxic to living organism. *Prochlorococcus* lost the cabability of to degrade ROS, and relies on the heterotrophs to degrade it.
In the ROS model, ROS is produced by both organisms and is toxic, limiting growth. The heterotroph degrades it.
ROS is an unstable compounds and also decays over time.

# Carbon Resources

Carbon resources are split into organic (DOC), inorganic (DIC) and recalcitrant, non-labile organic (RDOC). The carbon system is not closed, as new DIC is absorbed from the air via air-water exchange of $CO_2$.

Eq8: $totaldeathC_i = deathB_i * C2N_i + deathC_i$

Eq9: $dDIC/dt = DICAirWaterExchange + \sum{(respirationC_i - grossUptake_{i,DIC})}$

Eq10: $dDOC/dt = \sum{(totaldeathC_i * \gamma_{i} + overflowC_i - grossUptake_{i,DOC})}$

Eq11: $dRDOC/dt = \sum{totaldeathC_i * (1 - \gamma_{i})}$

Where $\gamma_{i}$ is the portion of the losses (death and leakiness) that is recalcitrant. 

# Nitrogen Resources

The system is closed to nitrogen. The nitrogen budget consists of the initial nitrogen in organic (DON), inroganic (DIN) and  recalcitrant, non-labile organic (RDON), as well as the nitrogen in the  cells biomass and nitrogen store. Thus the system as a whole is nitrogen limited.

Eq12: $totaldeathN_i = deathB_i + deathN_i$

Eq13: $dDIN/dt = \sum{(overflowN_i + DON2DIN_i - grossUptake_{i,DIN})}$

Eq14: $dDON/dt = \sum{(totaldeathN_i * \gamma_{i} - grossUptake_{i,DON} - DON2DIN_i)}$

Eq15: $dRDON/dt = \sum{totaldeathN_i * (1 - \gamma_{i})}$


Assuming that recalcitrant DON is released only during mortality/leakiness.
Assuming RDON/RDOC is recalcitrant to both organisms.



# Additional Model equations
## Death and losses
Loss processes (death and leakage) is modeled linearly as a precentage of the biomass/stores.
We use exponential decay – in ISMEJ we show that other formulations are better for co-cultures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.

Eq16: $deathB_i = M_i * B_i$

Eq17: $deathN_i = M_i * N_i$

Eq18: $deathC_i = M_i * C_i$

## Uptake

Eq19: $grossUptake_{ij} = Vmax_{ij} * B_{i} * limit_{ij}  * reg_{ij} * ROSpenalty_{i}$

Eq20: $uptakeN_i = grossUptake_{i,DIN} +  grossUptake_{i,DON}$

Eq21: $uptakeC_i = grossUptake_{i,DIC} +  grossUptake_{i,DOC}$



We are using monod limits to model the uptake dynamics based on resource availability

Eq22: $limit_{ij} = R_j / (R_j + Kn_{ij})$

Where $R_j$ is one of the resources ($DOC$, $DON$, $DIC$, $DIN$),
and $Kn_{ij}$ is the $Kn$ affinity for of organism *i* for resource *j*.


Uptake is regulated by not letting the stores grow too large. Using Droop-like equations.

Eq23: $QC_i = (C_i + B_i * C2N_i) / (N_i + B_i)$

Eq24: $regC_i = 1 - (QC_i / QCmax_i)$

Eq25: $regN_i = 1 - (QCmin_i / QC_i)$

Where $C2N_i$ is a parameter describing the C/N ratio in the biomass of organism *i*.
$QCmax_i$ and $QCmin_i$ are the minimum and maximum C/N ratios for organism *i*.
Both $regC_i$ and $regN_i$ are clipped to be between 0 and 1.
$QC$ is the C/N ratio of organism *i*.

**Question: is DIC uptake by photosynthesis regulated? **

ROS penalty is imposed in the ROS model. 

Eq26: $ROSpenalty_i = e^{- \omega_{i}*ROS}$

## ROS
ROS is modeled in the ROS model. ROS is produced by both organisms and is toxic, limiting growth.
ROS is an unstable compounds and decays over time.

Eq27: $ROSdecay = ROS * ROSdecayRate$

Eq28: $netROS = ROS - ROSdecay$

ROS production depends on the biomass

Eq29: $ROSrelease_i = E_{i,ROS} * B_i$

ROS breakdown is a common good that may be part of the positive interaction between *Prochlorococcus* and the heterotrophs.
*Prochlorococcus* cannot break down ROS, and may benefit from the breakdown performed by the heterotroph.

Eq30: $ROSbreakdown_H = Vmax_{H,ROS} * B_H * netROS / (netROS + Kn_{H,ROS})$

Eq31: $dROS/dt = \sum{ROSrelease_i} - ROSdecay - ROSbreakdown_H$

## Overflow

In overflow mode, the uptake that is over the C2N ratio is released as overflow.
In photosynthetic organisms this includes photosyntheate that cannnot be utilized due to N limitation.
In the heterotrophs, under C limitation, carbon is extracted from organic compouds and the N containing waste (e.g. NH4) is released.
Make the store maintain the C:N ratio and exude the rest.

Eq32: $netDeltaN_i = uptakeN_i + biomassbreakdownC_i / C2N_i - biosynthesisN_i $

Eq33: $netDeltaC_i = uptakeC_i + biomassbreakdownC_i - biosynthesisN_i * C2N_i  - respirationC_i$

Eq34: $storekeepN = min(netDeltaN, netDeltaC_i / C2N_i) $

Eq35: $overflowN_i = netDeltaN_i - storekeepN_i$

Eq36: $overflowC_i = netDeltaC_i - storekeepN_i * C2N_i$

## DON breakdown due to exoenzymes. 
In the exoenzyme model, DON is degraded to DIN by exoenzymes released by the bacteria. This DIN is then available to both bacteria for uptake.
We did not model spontaneous breakdown of DON, only exoenzyme mediated. Also there is no cost for the bacteria to release these enzymes.

Eq37: $DON2DIN_i = DON2DINrate_i * B_i * DON$


## Air water exchange
Eq38: $DICAirWaterExchange   = - (DIC - csat) / airWaterExchangeConstant$

Where $csat$ is the saturated DIC concentration. and $airWaterExchangeConstant$ is a constant computed as:

Eq39: $airWaterExchangeConstant = (h / (Kg * B * 0.01))$

Where $h$ is the height of the media in meters,  $Kg$ is the exchange rate in m sec-1 
and $B$ is the Revelle buffer factor. 


