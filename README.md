# Identification of Bouc–Wen type models using the transitional Markov chain Monte Carlo method

Source code associated to the paper:

Gilberto A. Ortiz, Diego A. Alvarez, Daniel Bedoya-Ruíz (2015). *Identification of Bouc–Wen type models using the Transitional Markov Chain Monte Carlo method*. Computers & Structures, Volume 146, Pages 252-269. https://doi.org/10.1016/j.compstruc.2014.10.012.

> **Abstract:** Bayesian model updating techniques are becoming the standard tool for the identification of nonlinear dynamical systems, because unlike other identification schemes which compute maximum likelihood values of parameters, Bayesian techniques provide probabilistic information of the estimates, which can be useful at the moment of making decisions with respect to the selection of parameters and/or mathematical models that simulate the nonlinear behavior experienced by the system. The aim of this paper is to provide an overview of the application of the Transitional Markov Chain Monte Carlo (TMCMC) method to the identification of the parameters of Bouc–Wen type models. The TMCMC method is a Bayesian model updating technique which not only finds the most plausible model parameters but also estimates the probability distribution of those parameters given the data measured at the laboratory. The TMCMC method identifies the structural system and allows the observation of multi-modality of the Bouc–Wen–Baber–Noori (BWBN) model of hysteresis. The performance of the algorithm is assessed using simulated and real data.

```
@article{Ortiz2015,
    title = {Identification of Bouc–Wen type models using the Transitional 
    Markov Chain Monte Carlo method},
    author = {Gilberto A. Ortiz and Diego A. Alvarez and Daniel Bedoya-Ruíz},    
    journal = {Computers \& Structures},
    volume = {146},
    pages = {252--269},
    year = {2015},
    doi = {https://doi.org/10.1016/j.compstruc.2014.10.012},
    keywords = {Bouc–Wen–Baber–Noori model, Hysteresis, Bayesian analysis, 
                System identification, Transitional Markov Chain Monte Carlo 
                (TMCMC)},
}
```

Run the files:
* `BW_exp_data.m`       Fitting a BWBW to experimental data
* `example_bouc_wen.m`  Fitting a BWBW to simulated data

