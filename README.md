Here I have uploaded the codes of the Compartmental Models and Meta Analysis.
 The codes for the Compartmental Models defines using the Stan Modelling language.The model includes six compartments: susceptible(s);infected with sensitive strain (I_s),infected with resistant strain (I_R),recoverd (R),and two compartments for deaths(D! and D2),The model also includes several parameters ,including transmission rates,recovery rates,death rates.
 The code compiles the Stan Model and fits it to the data using (MCMC) sampling.
 The Code extracts the posterior distributions for several parameters ,including the transmission rates,recovery rates and death rates.
 It also calculates the recovery times for the sensitive and resistant strains.
For the Meta Analysis part;we are aiming to calculate the effect of hand hygiene compliance on the reduction of nosocomial/healthcare associated infection.We are calculating the pooled RR(Relevant Risk) and then the PAF(POpulation Attributable Fraction).
