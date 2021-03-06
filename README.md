A Joint Bayesian Space-Time Model to Integrate Spatially Misaligned Air Pollution Data in R-INLA
------------------------------------------------------------------------------------------------

This repository contains the code to reproduce the results presented in
the paper Forlani C., Bhatt S., Cameletti M., Krainski E., Blangiardo
M., *A Joint Bayesian Space-Time Model to Integrate Spatially Misaligned
Air Pollution Data in R-INLA*, Environmetrics.
2020;e2644.<DOI:https://doi.org/10.1002/env.2644>

<!-- The ArXiv preprint can be found at http://arxiv.org/abs/2006.08988. -->
We are trying to improve the predictions of NO2 concentration
integrating data from two deterministic models called AQUM and PCM. The
observations from the monitoring stations represent our response
variable.

We want to model these data jointly in order to take into account the
uncertainty associated with the deterministic model outputs.

### Data

The data workspace can be requested directly to the main author.

The workspace contains the following objects:

-   `aqum`: daily output of the AQUM model for 495 locations and 1823
    days
-   `pcm`: annual output of the PCM model for 44117 locations and 5
    years
-   `final_dataset`: daily observations of NO2 concentration for 126
    monitors and 1826 days with associated site type
-   `boundary`: boundary of the study area
-   `london.shape`: boundary of Greater London
-   `england.shape`: boundary of England
-   `shape`: combines boundary of study area and geographical boundary
    of England
-   `coordinates.aqum`: coordinates of the 495 AQUM locations
-   `coordinates.pcm`: coordinates of the 44117 PCM locations
-   `monitors`: coordinates of the 126 monitoring stations
-   `monitors_val`: list of the 6 sets of monitors used for validation
-   `pred.grid`: locations where we want to make predictions and
    corresponding site type
    <!-- * `pred.grid.aqum`: values of AQUM interpolated at the prediction grid nodes, for each day (not used for the joint model predictions) -->
    <!-- * `pred.grid.pcm`: values of PCM interpolated at the prediction grid nodes, for each year (not used for the joint model predictions) -->
-   `roads_major`: shapefile of motorways over the study area

Note that all coordinates are in UTM projection so that distances are
computed in meters. Results do not change when rescaling UTM coordinates
to Kilometers.

Please see the paper for data description.

### How to run the code

To run the joint model presented in section 3.2 and produce the plots of
section 4.2 you just need to run the script `joint_model.R` which
contains instructions for data praparation and sources all the other
necessary scripts. The script should be ideally called from an array job
passing DATA\_ID as a variable taking values 1 to 6, in order to run the
model for the 6 validation sets in parallel.

The script `separate_models.R` contains the code to run the models for
AQUM and PCM described in section 3.1.

The script `competitor_models.R` contains the code to run the models
used for comparison described in section 3.3. This script should be
ideally called from an array job passing DATA\_ID as a variable taking
values 1 to 6, in order to run the model for the 6 validation sets in
parallel, and MODEL\_ID taking values 1 to 7 for the 7 competitor
models.

<!-- These can also be run manually setting `data_id` (1 to 6) and `formula_id` (1 to 9) once the data are loaded. -->
<!-- ### Notes on the notation -->
<!-- Please see the paper for details on the model specification. -->
<!-- The following are differences in the notation from the paper to the code: -->
<!-- * $y_1$ is referred to as `pcm` -->
<!-- * $y_2$ is referred to as `aqum` -->
<!-- * $y_3$ is referred to as `y` -->
<!-- * $z_1$ is referred to as `z1`, the copy for y as `z13`, the copy for AQUM as `z12` -->
<!-- * $z_2$ is referred to as `z2`, and the copy for y as `z23` -->
<!-- * $z_3$ is referred to as `z3` -->
<!-- * $\lambda_{1,2}$ is automatically called `Beta for z12` and its prior is called `lambda12`  -->
<!-- * $\lambda_{1,3}$ is automatically called `Beta for z13` and its prior is called `lambda13` -->
<!-- * $\lambda_{2,3}$ is automatically called `Beta for z23` and its prior is called `lambda23` -->
<!-- * The intercepts are `alpha1`, `alpha2` and `alpha3` -->
<!-- * $\beta_{k_s}$ are `betaURB` and `betaRKS` respectively -->
<!---

# Model specification


### Data 

* PCM covariate ($Y_{1}$): annual ($t_1=5$) NO2 concentration output from deterministic model at $s_1=44117$ grid nodes
* AQUM covariate ($Y_{2}$): daily ($t_2=1826$) NO2 concentration output from deterministic model at $s_2=495$ grid nodes
* Response variable ($Y_{3}$): daily NO2 concentration at $s_3=126$ monitoring stations for $t=1,\dots,t_2$, $t_2=1826$ consecutive days (time series are not all complete)
* Site type covariate: categorical variable associated with the monitoring stations ($k=0$: rural (reference), $k=1$: urban and $k=2$: road-kerb side)

### Assumptions
* AQUM provides spatial and temporal information 
* PCM provides spatial-only information because temporal resolution is very low
* Underlying spatial field ($z_1(s)$) is the same for PCM, AQUM and monitor observations (they all measure NO$_2$ concentration)
* Underlying temporal field ($z_2(t)$) is the same for AQUM and monitor observations (they both measure NO$_2$ concentration)

### Level 1:

Let $y_i(s,t)$ denote the PCM ($i=1$) and AQUM ($i=2$) data and the observed NO$_2$ concentration ($i=3$) at the generic time point $t$ and site $s$, on the logarithmic scale. These are assumed to be normally distributed, with mean $\eta_i(s,t)$ and measurement error variance $\sigma^2_{\epsilon_i}$: 

$y_1(s,t) \sim N(\eta_1(s), \sigma^{2}_{\epsilon_1}) \qquad \text{(PCM)}$

$y_2(s,t) \sim N(\eta_2(t), \sigma^{2}_{\epsilon_2}) \qquad \text{(AQUM)}$

$y_3(s,t) \sim N(\eta_3(s,t), \sigma^{2}_{\epsilon_3}) \qquad \text{(Ground observations)}$

The three corresponding linear predictors are:

$\eta_1(s) = \alpha_1 + z_1(s) \qquad \text{(PCM)}$

$\eta_2(t) = \alpha_2 + z_2(t) \qquad \text{(AQUM)}$

$\eta_3(s,t) = \alpha_1 + \alpha_2 + \alpha_3  + \beta_{k_s} + \lambda_{1,3} z_1(s) + \lambda_{2,3} z_2(t) + z_3(t,k{_s})$

where 

* $\alpha_i$ are the intercepts
* $\lambda_{i,j}$ are optional scaling parameters
* $\beta_{k_s}$ is the fixed effects for the site type as categorical variable ($k=0$: rural (reference), $k=1$: urban and $k=2$: road-kerb side)
* $z_1$ and $z_2$ are the shared random effects
* $z_3(t,k{_s})$ is an interaction term which allows for a different residual temporal trend for each site type

### Level 2

* $\boldsymbol{z}_1 \sim MVN(\textbf{0}, \sigma^2_{z_1}\boldsymbol{\Sigma})$ is the common spatial latent field, with $\boldsymbol{\Sigma}$ being the correlation matrix defined by the Mat\'ern stationary and isotropic covariance function. Note that $\boldsymbol{z}_1$ is then calibrated on the monitor observations through $\lambda_{1,3}$ and $\lambda_{1,3}$. 
* $z_2(t)$ is the $t$-th element of the temporal latent field $\boldsymbol{z_2}$, and is modelled as a random walk: $z_2(t) \sim N(z_2(t-1), \sigma^2_{z_2})$. %This is defined on AQUM and copied in $\boldsymbol{\eta_3}$ with a calibration coefficient.
Similarly to $\boldsymbol{z}_1$, $\boldsymbol{z}_2$ is calibrated on the monitor observations through $\lambda_{2,3}$. 
* $z_3(t,k{_s})$ is the residual temporal trend assumed to be different for each site type (rural, urban, road-kerb side), and modelled as first order autoregressive $z_3(t,k{_s}) \sim N(\rho z_3(t-1,k{_s}), \sigma^2_{z_3})$. In other words, we assume conditionally independent replications of the same latent field for each site type, with shared hyperparameters.



### Level 3

* $log(1/\sigma^2_{z_1}) \sim logGamma(1, 5e-05)$
 
* $log(1/\sigma^2_{z_3}) \sim logGamma(1, 5e-05)$
 
* $P(\sigma_{z_2}>SD(AQUM))=0.01$

* $\rho \sim N(0.3,0.5)$ 
 
* $P(r<r_0)=0.95$, where $r_0 = 1/5$ of the domain size
 
* $P(\sigma_{z_2}>\sigma_0)=0.5$, where $\sigma_0 = 100$.

* $log(1/\sigma^2_{\epsilon_i}) \sim logGamma(1, 5e-05)$, $i=1,2,3$
 
* $\alpha_i, \beta_{k{_s}} \sim N(0,1000)$
 
* $\lambda_{1,3}, \lambda_{2,3} \sim N(1, 0.01)$


-->
