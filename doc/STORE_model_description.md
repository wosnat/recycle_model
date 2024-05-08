# Introduction
This document describes the model. We are modeling the interaction between *Prochlorococcus*, an autotrophic phytoplankton, and co-occurring marine heterotrophic bacteria.

The model represents the interaction between the autotroph and heterotroph through the lens of nutrients exchanges and detoxification.


# Model variables

|Variable | Unit         | Description                          |
|---------|--------------|--------------------------------------|
| $B_i$   | $\mu M N$    | Functional biomass in N              |
| $N_i$   | $\mu M N$    |  N store                             |
| $C_i$   | $\mu M C$    |  C store                             |
| $DON$   | $\mu M N$    |  Dissolved Organic N                 |
| $RDON$  | $\mu M N$    |  Recalcitrant Dissolved Organic N    |
| $DIN$   | $\mu M N$    |  Dissolved Inorganic N               |
| $DOC$   | $\mu M C$    |  Dissolved Organic C                 |
| $RDOC$  | $\mu M C$    |  Recalcitrant Dissolved Organic C    |
| $DIC$   | $\mu M C$    |  Dissolved Inorganic C               |
| $ROS$   | $\mu M ROS$  |  Reactive Oxygen Species             |


$i$ is the subscript for population and can be either *Prochlorococcus* or Heterotroph.

# Main model overview
![image](https://github.com/wosnat/recycle_model/assets/22752755/93bd90ac-fb83-4eb8-be1f-3a8d5a945e48)

In our model, C and N are taken up from dissolved organic and/or inorganic forms into C and N stores ($C_i$ and $N_i$ respectively), which are then combined through biosynthesis into functional biomass ($B_i$, figure 1). Functional biomass is either degraded back into stores to support respiration or lost through mortality (Eq 1)

The functional biomass has a fixed C/N ratio ($C2N_i$). $B_i$, the biomass variable is maintained in $\mu M  N$, representing also the corresponding C biomass of $C2N_i  \cdot  B_i$.

The stores are utilized for biosynthesis, respiration (C store only) and are degraded into dissolved organic matter upon death (Eq2, Eq3). When the stores are imbalanced, some of them may be exuded as overflow (see below).

EQ1:    $dB_i/dt = biosynthesisN_i - biomassbreakdownC_i / C2N_i - deathB_i  $

EQ2:    $dN_i/dt = uptakeN_i + biomassbreakdownC_i / C2N_i - biosynthesisN_i - overflowN_i - deathN_i$

EQ3:    $dC_i/dt = uptakeC_i + biomassbreakdownC_i - biosynthesisN_i  \cdot  C2N_i  - respirationC_i - overflowC_i - deathC_i$

<br />

The synthesis of functional biomass ($\mu M N sec-1$)) is determined by the limiting store (C or N in $\mu M$) times a biosynthesis rate ($sec-1$) (EQ4). The C and N stores available for biosynthesis include storage and uptake.

EQ4:   $biosynthesisN_i = \min(N_i + uptakeN_i, \frac{C_i + uptakeC_i}{C2N_i})  \cdot  Kmtb_i$


Respiration is comprised of growth associated respiration controlled by $b_i$ and maintenance associated respiration controlled by $r0$ (Eq5).

EQ5:   $respirationC_i = (b_i  \cdot  biosynthesisN_i + B_i  \cdot  r0_i)  \cdot  C2N_i$


The breakdown of biomass does not need occur if the cell is in balanced growth (i.e. can be equal to zero). However, biomass may be broken down to support respiration. Respiration is assumed to occur from the C storage (e.g. carbohydrates or lipids), but if there is not enough C storage biomass will be degraded instead (EQ7). Note that equation 6 is calculated after the calculation of biosynthesis but before the calculation of the C store, and hence the needs to incorporate also the C used for biosynthesis (EQ6,7).

EQ6: $requiredCStore_i = C_i + uptakeC_i - biosynthesisN_i  \cdot  C2N_i $

EQ7: $biomassbreakdownC_i = \max(0, respirationC_i - requiredCStore_i)$

# Model variants

There is a baseline model and additional 4 model variants encoding different interaction mechanisms between the 2 populations (Exoenzymes, Excess Overflow, Mixotrophy, ROS detoxification). 

### Baseline model
In the baseline model, the heterotroph consumes organic N and C as well is inorganic N. While the autotroph consumes inorganic C (photosynthesis) and inorganic N. The only processes in the base model are nutrients uptake, biosynthesis, respiration, and death. In addition, there is a very slow abiotic degradation of DON to DIN. All other models are a superset of the baseline model, adding additional interaction-related processes.

### Excess Overflow

In the overflow model, the uptake that is over the C2N ratio is released as overflow.
In photosynthetic organisms this includes photosynthate that cannot be utilized due to N limitation.
In the heterotrophs, under C limitation, carbon is extracted from organic compounds and the N containing waste (e.g. NH4) is released.
The model makes the store maintain the C:N ratio and exude the rest.

### Mixotrophy
In the mixotrophy model, *Prochlorococcus* can feed on organic compounds as well as inorganic compounds (a mixotroph). This is implemented by a non-zero $Vmax_{i,DOC}$, $Vmax_{i,DON}$ for *Prochlorococcus*. 

### Exoenzymes
In the exoenzyme model, DON is degraded to DIN by exoenzymes released by the heterotrophic bacteria. This DIN is then available to both bacteria for uptake. There is no cost for the bacteria to release these enzymes.

### ROS detoxyfication
ROS (reactive oxygen species) is a natural byproduct of core metabolic processes as well as photosynthesis.
It is toxic to living organism. *Prochlorococcus* lost the capability of to degrade ROS and relies on the heterotrophs to degrade it.
In the ROS model, ROS is produced by both organisms and is toxic, increasing the loss rate. Both organisms can also degrade it. 
ROS is an unstable compound and also decays over time.

# Carbon Resources

Carbon resources are split into organic (DOC), inorganic (DIC) and recalcitrant, non-labile organic (RDOC). The carbon system is not closed, as new DIC is absorbed and lost from the air via air-water exchange of $CO_2$.

Eq8: $totaldeathC_i = deathB_i  \cdot  C2N_i + deathC_i$

Eq9: $dDIC/dt = DICAirWaterExchange + \sum{(respirationC_i - grossUptake_{i,DIC})}$

Eq10: $dDOC/dt = \sum{(totaldeathC_i  \cdot  \gamma_{i} + overflowC_i - grossUptake_{i,DOC})}$

Eq11: $dRDOC/dt = \sum{totaldeathC_i  \cdot  (1 - \gamma_{i})}$

Where $\gamma_{i}$ is the portion of the losses (death and leakiness) that is labile (avaliable for consumption). 

# Nitrogen Resources

The system is closed to nitrogen. The nitrogen budget consists of the initial nitrogen in organic (DON), inorganic (DIN) and recalcitrant, non-labile organic (RDON), as well as the nitrogen in the  cells biomass and nitrogen store. Thus, the system as a whole is nitrogen limited.
DON is constantly degraded abiotically into DIN (in all models). In the EXOENZYME model it is also degraded by exoenzymes released be the heterotroph (see DON breakdown due to exoenzymes)

Eq12: $totaldeathN_i = deathB_i + deathN_i$

Eq13: $dDIN/dt = \sum{(overflowN_i + DON2DIN_i - grossUptake_{i,DIN})} + globalDON2DIN$

Eq14: $dDON/dt = \sum{(totaldeathN_i  \cdot  \gamma_{i} - grossUptake_{i,DON} - DON2DINexo_i)} - globalDON2DIN $
Eq35: $DON2DINexo_i = KprodEXO_i  \cdot  B_i  \cdot  DON$


Eq15: $dRDON/dt = \sum{totaldeathN_i  \cdot  (1 - \gamma_{i})}$


Assuming that recalcitrant DON is released only during mortality/leakiness.
Assuming RDON/RDOC is recalcitrant to both organisms.



# Additional Model equations
## Uptake

Eq20: $uptakeN_i = grossUptake_{i,DIN} +  grossUptake_{i,DON}$

Eq21: $uptakeC_i = grossUptake_{i,DIC} +  grossUptake_{i,DOC}$

We are using monod limits to model the uptake dynamics based on resource availability.
Where $R_j$ is one of the resources ($DOC$, $DON$, $DIC$, $DIN$),
and $Kn_{ij}$ is the $Kn$ affinity for of organism $i$ for resource $j$ (EQ19).

Eq19: $grossUptake_{ij} = Vmax_{ij}  \cdot  B_{i}  \cdot  \frac{R_j}{R_j + Kn_{ij}}   \cdot  reg_{ij}$

Uptake is regulated by not letting the stores grow too large. Using Droop-like equations.

Eq22: $Q_{i,C} = \frac{C_i + B_i  \cdot  C2N_i}{N_i + B_i}$

Eq23: $Q_{i,N} =  \frac{N_i + B_i}{C_i + B_i  \cdot  C2N_i}$

Eq24: $reg_{ij} = \frac{Qmax_{ij} - Q_{ij}}{Qmax_{ij} - Qmin_{ij}}$


Where $C2N_i$ is a parameter describing the C/N ratio in the biomass of organism $i$.
$Qmax_i$ and $Qmin_i$ are the minimum and maximum C/N ratios for organism $i$.
$reg_{ij}$ is clipped to be between 0 and 1.
$Q_{i,C}$ is the C/N ratio of organism $i$, $Q_{i,N}$ is the N/C ratio of organism $i$.

**Question: is DIC uptake by photosynthesis regulated? **

## Death and losses
Loss processes (death and leakage) is modeled linearly as a precentage of the biomass/stores.
We use exponential decay – in ISMEJ we show that other formulations are better for co-cultures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.

In the ROS model, ROS is a toxin, killing some of the bacteria, and ROS influence is modeled as additional mortality. The additional mortality is capped by $ROSmaxD$, the saturated maximum toxicity of ROS. (see ROS section for additional details).


Eq12: $additionalLossRate_i = \min(ROS * \omega_i, ROSmaxD)$

Eq12: $lossRate_i = M_i + additionalLossRate_i$ 

Eq16: $deathB_i = lossRate_i   \cdot  B_i$

Eq17: $deathN_i = lossRate_i   \cdot  N_i$

Eq18: $deathC_i = lossRate_i  \cdot  C_i$


## ROS
ROS is modeled in the ROS model. ROS is produced by both organisms and is toxic, limiting growth. Both organisms can also degrade ROS. 
ROS is an unstable compounds and decays over time.

Eq29: $ROSdecay = ROSdecayRate  \cdot  ROS $

ROS production depends on the biomass

Eq29: $ROSrelease_i = KprodROS_{i}  \cdot  B_i$


ROS breakdown is a common good that may be part of the positive interaction between *Prochlorococcus* and the heterotrophs.
*Prochlorococcus* has very limited ability to break down ROS, and may benefit from the breakdown performed by the heterotroph.

Eq30: $ROSbreakdown_i = KlossROS_{i}  \cdot  B_{i}  \cdot  ROS$

Eq31: $dROS/dt = \sum{(ROSrelease_i - ROSbreakdown_i)} - ROSdecay $

ROS penalty is imposed in the ROS model. It is capped by a max toxicity rate.  
Eq26: ROSDeathRate_i = \max(\omega_{i} \cdot ROS, maxROSDeathRate_i)

## Overflow

In overflow mode, the uptake that is over the C2N ratio is released as overflow.
In photosynthetic organisms this includes photosyntheate that cannnot be utilized due to N limitation.
In the heterotrophs, under C limitation, carbon is extracted from organic compouds and the N containing waste (e.g. NH4) is released.
Make the store maintain the C:N ratio and exude the rest.

Eq32: $delta_{i,N} = N_i - C_i  / C2N_i $

Eq33: $delta_{i,C} = C_i - N_i  \cdot  C2N_i  $

Eq24: $Oreg_{ij} = [\frac{Q_{ij}}{Qmax_{ij}}]^4$

Eq34: $overflow_{ij} = max(0, delta_{ij})  \cdot  Oreg_{ij}  \cdot  Koverflow_{ij}$

Where $Oreg_{ij}$ is the regulation of the overflow release (increases as the nutrient imbalance increases), see Uptake section for description of $Q_{ij}$ and $Qmax_{ij}$ and $Koverflow_{i}$ is the overflow rate (sec-1).


## DON breakdown due to exoenzymes
In the exoenzyme model, DON is degraded to DIN by exoenzymes released by the bacteria. This DIN is then available to both bacteria for uptake.
We did not model spontaneous breakdown of DON, only exoenzyme mediated. Also there is no cost for the bacteria to release these enzymes.

Eq35: $DON2DINexo_i = KprodEXO_i  \cdot  B_i  \cdot  DON$

Eq36: $globalDON2DIN = KdecayDON  \cdot  DON$

## Air water exchange
Eq37: $DICAirWaterExchange   = - \frac{DIC - csat}{airWaterExchangeConstant}$

Where $csat$ is the saturated DIC concentration and $airWaterExchangeConstant$ is a constant computed as:

Eq38: $airWaterExchangeConstant = \frac{h}{Kg  \cdot  B  \cdot  0.01}$

Where $h$ is the height of the media in meters,  $Kg$ is the exchange rate in m sec-1 
and $B$ is the Revelle buffer factor. 


