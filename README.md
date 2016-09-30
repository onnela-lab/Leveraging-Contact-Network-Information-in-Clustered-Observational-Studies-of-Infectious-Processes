# Leveraging-Contact-Network-Information-in-Clustered-Observational-Studies-of-Infectious-Processes

This README describes how to replicate the analysis in our manuscript.

**Part A: Simulation of Infectious Process on Contact Networks.**

   1. Use R to run `simulation.R` For values `1-2000` of `SEED`.  (The code is arranged this way for easy paralellization.)
        Note that this code calls `data_code.py`, which generates the synthetic contact networks and performs the epidemic process.       
   2. Use R to run `analysis.R` for  values `1-2000` of `SEED`.  This code is also designed for paralellization.  This creates all relevant simulation output.
        Note that this will require downloading the `OBSgeeDR` package for observational augmented GEE analysis.        
   3. Use R to run "global analysis.R" for final formatting.

**Part B: Application to Microfinance Data.**

   1. The microfinance data used in this paper is publically available, referenced here:  http://economics.mit.edu/faculty/eduflo/social
   2. use Python to run `data generation.py` to compile the microfinance data into a single dataset.   
   3. Use R to run `Data Completion.R` to concatenate the data with the defined exposures.   
   4. Use R to run `Empirical Data Analysis.R` to create the tables in the paper and supplement.

For any questions, please contact the corresponding author (Patrick Staples) at *patrickstaples* at *fas* dot *harvard* dot *edu*.
