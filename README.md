# DITRAS

DITRAS (DIary-based TRAjectory Simulator) is a framework to simulate the spatio-temporal patterns of human mobility in a realistic way. DITRAS operates in two steps: the generation of a mobility diary and the translation of the mobility diary into a mobility trajectory. The mobility diary is constructed by a Markov model which captures the tendency of individuals to follow or break their routine. The mobility trajectory is produced by a model based on the concept of preferential exploration and preferential return. DITRAS can reproduce the statistical properties of real human mobility trajectories in an accurate way.

How to use DITRAS
-----------------
To run the code of DITRAS you must put on the terminal the following command format:   
`python DITRAS.py n_agents time_slots spatial_tessellation od_matrix diary_generator_1hour.pkl`, where:
- **n_agents** is an integer indicating the number of agents to simulate
- **time_slots** is the length of the time period, in hours
- **spatial_tessellation** is the of a file where the spatial tessellation is store. The format of the spatial_tessellation file must be the following: `latitude,longitude,relevance`. Latitude and longitude are float numbers indicating the geographic position of a location and relevance is an integer indicating the importance of the location (e.g., the number of residents in the locations or the number of calls made in that location during a given period). This an example of the file:   
`42.546,2.435,125`   
`36.435,12.324,20`  
`43.435,2.121,120`  
`...,...,...`  
- **od_matrix** is the name of a file where the od_matrix based on the spatial_tessellation is stored. If the file specified does not exist, the file will be automatically computed and stored in a file with the name specified in input.
- **diary_generator_1hour.pkl** the diary generator we provide stored as a pickle file.




