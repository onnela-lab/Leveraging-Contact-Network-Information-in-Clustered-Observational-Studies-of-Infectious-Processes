# Leveraging-Contact-Network-Information-in-Clustered-Observational-Studies-of-Infectious-Processes

This README describes how to replicate the analysis in our paper.

Part A: Simulation of Infectious Process on Contact Networks.
   `A1.)` Use R to run "Simulation.R" For values `1-2000` of `SEED`.  (The code is arranged this way for easy paralellization.)
        Note that this code calls `data_code.py`, which generates the synthetic contact networks and performs the epidemic process.
        
   `A2.)` Use R to run "Analysis.R" for  values `1-2000` of `SEED`.  This code is also designed for paralellization.  This creates all relevant simulation output.
        Note that this will require downloading the "OBSgeeDR" package for observational augmented GEE analysis.
        
   `A3.)` Use R to run "global analysis.R" for final formatting.

Part B: Application to Microfinance Data.

   `B1.)` The microfinance data used in this paper is publically available, referenced here:  http://economics.mit.edu/faculty/eduflo/social
   
   `B2.)` use Python to run `data generation.py` to compile the microfinance data into a single dataset.
   
   `B3.)` Use R to run `Data Completion.R` to concatenate the data with the defined exposures.
   
   `B4.)` Use R to run `Empirical Data Analysis.R` to create the tables in the paper and supplement.

For any questions, please contact the corresponding author (Patrick Staples) at `patrickstaples@fas.harvard.edu`.
