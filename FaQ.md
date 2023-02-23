![logo def-github-04](https://user-images.githubusercontent.com/30837036/91882995-13c90200-ec84-11ea-9643-0191dfbca995.jpg)
#  :computer: Frequently Asked Questions 

**As fingerprinting model developers and model users many times ourselves, we often see many questions that surface repeatedly. This FAQ section attempts to gather some of those and provide some answers for researchers and other users!**

-------------

<details open>

<summary><strong>Table of Contents (FaQs)</strong></summary>

<!-- toc -->

- [Is the tracer selection step needed in fingerprinting studies?](#Is-the-tracer-selection-step-needed-in-fingerprinting-studies?)
- [How many methods for tracer selection are out there? ](#How-many-methods-for-tracer-selection-are-out-there?)
- [What are the advantages of the new CI, CR and CTS?](#What-are-the-advantages-of-the-new-CI,-CR-and-CTS?)
- [Can the CI, CR and CTS be implemented before using Bayesian models?](#Can-the-CI,-CR-and-CTS-be-implemented-before-using-Bayesian-models?)
- [Are tracer selection methods only needed in Frequentist models?](#Are-tracer-selection-methods-only-needed-in-Frequentist-models?)
- [References](#References)

<!-- tocstop -->

</details>


### Is the tracer selection step needed in fingerprinting studies?
As each dataset is unique and has tracers of different natures, such as elemental compositions, radionuclides, stable isotopes, magnetic susceptibility… a unilateral response is not possible. However, an apparent agreement exists for real studies regarding the need to select an optimal set of tracers for every mixture. Thus, based on most of the previously published research and our experience, we would suggest implementing it.

### How many methods for tracer selection are out there? 
The sediment fingerprinting technique dates back to the 70s. Thus, several tracer selection methods have been proposed and implemented since then. Included in the package, you can find the three-step tracer selection method proposed by Collins and Walling (2002) together with the new CI, CR and CTS methods recently created by Lizaga et al. (2020) and Latorre et al. (2021). Furthermore, additional techniques, such as biplots suggested by different authors (e.g., Pulley et al., 2015), can be extracted from the correlationPlot() displaying similar information to the mixing polygon (Phillips and Gregg 2003), together with the sources and tracers correlation.

### What are the advantages of the new CI, CR and CTS?
Overall, the CI, CR and CTS methods do not produce a fixed tracer selection but supply the user with valuable information to support a justified selection or to understand the datasets. These three methods provide different information:
- The CI is a non-parametric test that analyses tracer conservativeness using the mixture's and source's data to create an index highlighting how conservative a tracer is. Combined with this method are the ternary plots representing all possible predictions that each individual tracer introduces into the models.
- The CR creates a ranking based on the consensus between the different tracers in your dataset. Consensus-reaching processes proceed in a convergent multistage way, where experts present their opinions and discuss and negotiate to bring positions closer by modifying their initial opinions (Pérez et al., 2018). In this case, these experts are the tracers in the dataset.
- The CTS, similar to the DFA, identifies the most discriminant tracers, but it also analyses their mathematical properties to avoid inconsistency in overdetermined systems (e.g., 3 sources and 3 or more tracers). One handy feature of this method is that it displays all possible combinations of tracers and orders them in terms of their discriminant capacity. Besides, by using the CTS method, you can check out the presence of multiple solutions in your dataset and use expert knowledge to understand their origin and select between them if the case.

### Are tracer selection methods only needed in Frequentist models?
As proven in several studies, although different models are not equally affected, including non-conservative or erroneous traces affects them all (Cooper & Krueger, 2017).

### Can the CI, CR and CTS be implemented before using Bayesian models?
Of course, these methods extract essential information from your dataset and are all implemented before modelling. As seen in Lizaga et al. (2020) & Latorre et al. (2021), implementing these methods improves the model performance in both Frequentist and Bayesian models.

### References
- Collins, A.L.,Walling, D.E., 2002. Selecting fingerprint properties for discriminating potential suspended sediment sources in river basins. J. Hydrol. 261, 218–244. https://doi.org/10.1016/S0022-1694(02)00011-2.
- Cooper, R.J., Krueger, T., Hiscock, K.M., Rawlins, B.G., 2014. Sensitivity of fluvial sediment source apportionment to mixing model assumptions: A Bayesian model comparison. Water Resour. Res. 50, 9031–9047. https://doi.org/10.1002/2014WR016194
- Latorre, B., Lizaga, I., Gaspar, L., Navas, A., 2021. A novel method for analysing consistency and unravelling multiple solutions in sediment fingerprinting. Sci. Total Environ. 789, 147804. https://doi.org/10.1016/j.scitotenv.2021.147804
- Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2022. Combined use of geochemistry and compound-specific stable isotopes for sediment fingerprinting and tracing. Sci. Total Environ. 154834. https://doi.org/10.1016/j.scitotenv.2022.154834 
- Pérez, I.J., Cabrerizo, F.J., Alonso, S., Dong, Y.C., Chiclana, F., Herrera-Viedma, E., 2018. On dynamic consensus processes in group decision making problems. Inf. Sci. 459, 20–35. https://doi.org/10.1016/j.ins.2018.05.017.
- Phillips, D.L., Gregg, J.W., 2003. Source partitioning using stable isotopes: coping with too many sources. Oecologia 136, 261–269. https://doi.org/10.1007/s00442-003-1218-3.
- Pulley, S., Rowntree, K., Foster, I., 2015. Conservatism of mineral magnetic signatures in farm dam sediments in the South African Karoo: the potential effects of particle size and post-depositional diagenesis. J Soils Sediments 15, 2387–2397. https://doi.org/10.1007/s11368-015-1265-5
--------------------------------------------------------------------------------