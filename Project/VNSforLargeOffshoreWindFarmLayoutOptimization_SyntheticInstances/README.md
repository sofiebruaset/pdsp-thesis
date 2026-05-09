# Synthetic set of instances

## Layouts
10 layout instances have been generated and used in the tests of the paper "Variable Neighborhood Search for Large Offshore Wind Farm Layout Optimization" by Davide Cazzaro and David Pisinger. We make available this set of data with the purpose to foster research on realistic layouts, which are seldomly used in literature, because the real data are often highly confidential. We call an instance the data that describe a wind farm layout, in which we want to place a certain number of turbines to maximize its power production and minimize the wake effects. We call these instances "synthetic" because they try to be realistic cases. If a technique can solve these instances close to optimality, the same technique will be ready to optimize real wind farms.

The following table reports the details of each instance. The number of available positions is the number of positions where turbines can be placed in the layout. Zones report how many sub areas are in the layout (each zone has a number of turbines to place in it). Number of Fixed turbines reports how many fixed turbines are in the layout. These turbines are not subject to the optimization (they cannot be moved), but their wake effects influence the optimization. In the runs we did in our paper, the fixed turbines always are 10MW turbines, for which we provide the data below. Finally, we report the number of turbines of 15MW that we place in each zone respectively, in the case with a high power density of 8MW/km2. We used this high density in the tests reported in the paper, because this makes the optimization more challenging when considering the minimum distance.

| Synthetic instance | Avaialable Positions | Zones | Number of Fixed turbines (10MW) | Number of Turbines (15MW) per zone (density 8MW/km2) |
|--------------------|----------------------|-------|---------------------------------|------------------------------------------------------|
| Instance A         | 3238                 | 2     | 42                              | 26, 8                                                |
| Instance B         | 6989                 | 1     | 15                              | 99                                                   |
| Instance C         | 7098                 | 2     | 8                               | 60, 30                                               |
| Instance D         | 10443                | 1     | 45                              | 170                                                  |
| Instance E         | 11518                | 3     | 40                              | 7, 94, 36                                            |
| Instance F         | 11548                | 2     | 12                              | 132, 26                                              |
| Instance G         | 14637                | 1     | 35                              | 140                                                  |
| Instance H         | 19498                | 2     | 40                              | 158, 30                                              |
| Instance I         | 20247                | 1     | 36                              | 313                                                  |
| Instance J         | 21709                | 3     | 75                              | 136, 74, 25                                          |

In the ```site/``` folder there is the list of all 10 instances, one per sub-folder. Each subfolder contains 5 files.

The most important file is ```availablePositions.txt``` which contains the coordinates of the positions in which the turbines can be placed.
The format of each row of the file is:
```x_coordinate y_coordinate z_coordinate foundation_cost zone_number```
Values are separated by whitespace. The x and y coordinates are the (arbitrary) latitude and longitude of the position of turbine (since the area covered by a wind farm is relatively small, we can assume these are coordinates in a flat 2D surface) in meters. The distance between two positions is just the euclidean distance between the two points. z_coordinate reports, in meters, the water depth of the position. The fourth column reports the estimated foundation cost, of building a turbine in that position, in Euros.

In each instance folder there is also the file ```fixed_wf.txt```, which contains the locations of the fixed turbines. The format is the same as for the available positions file: ```x_coordinate y_coordinate z_coordinate foundation_cost zone_number```, but in this case the ```z_coordinate``` and ```foundation_cost``` parameters have dummy values and must not be used (since these turbines are part of an existing wind farm, their foundation cost is not relevant because we do not optimize these positions).

The third file is ```geometry.txt``` file. This file is not needed to optimize the layout. It reports the coordinates of the boundaries and obstacles in the layout (that we use for plotting) and the list connection edges.

The fourth file is ```info.json``` file. This file reports additional informations about the generated layout, but it is not needed for the optimization. It reports the number of points (locations) for turbines in the layout, an internal seed used to generate the layout itself, and the list of zones. For each zone, it reports the area of the zone, the number of obstacles inside this area (called holes), the number of available positions in this area, and for each power density of 4, 6, and 8 MW/km2 the corresponding number of turbines to place if we use 10MW or 15MW turbines. In our optimization runs we use the value for 15MW turbines with the highest density, as we reported in the table above.

Finally, each instance folder contains ```plot.png``` file, which is the plot of the layout. Here we can see the geometry of the zones in magenta (borders and obstacles), a color map with the depth of each available position (in most plot it seems a continuous surface because of the density of the points, but it is a discrete set). In orange we report the positions of the fixed turbines (existing turbines) in the area.

## Wind Turbine Generators

The data of the characteristic curves of 2 WTGs are available, which are publicly available from NREL, in the ```wtg/``` folder.
The performance characteristics of wind turbines are given by the power coefficient and by the thrust coefficient. We report in addition the hub height and the rotor diameter of each turbine below.


Data format:
```wind_speed    power_coefficient    thrust_coefficient```
Where ```wind_speed``` is the wind velocity in m/s, ```power``` is the power produced by the turbine at this wind speed in MW, and ```thrust_coefficient``` is the thrust coefficient of the turbine at this wind speed.

### 10 MW turbine
File: ```NREL-10-179.txt```
Rotor diameter: 179 m.
Hub height: 119 m.

https://github.com/IEAWindTask37/IEA-10.0-198-RWT
@techreport{RWT,
Author = {Pietro Bortolotti and Helena Canet Tarres and Katherine Dykes and Karl Merz and Latha Sethuraman and David Verelst and Frederik Zahle},
Howpublished = {NREL/TP-73492},
Institution = {International Energy Agency},
Title = {IEA Wind Task 37 on Systems Engineering in Wind Energy -- WP2.1 Reference Wind Turbines},
Url ={https://www.nrel.gov/docs/fy19osti/73492.pdf},
Year = {2019}}

### 15 MW turbine:
File: ```NREL-15-240.txt```
Rotor diameter: 240 m.
Hub height: 150 m.

https://github.com/IEAWindTask37/IEA-15-240-RWT
@techreport{IEA15MW_ORWT,
author = {Evan Gaertner and Jennifer Rinker and Latha Sethuraman and Frederik Zahle and Benjamin Anderson and Garrett Barter and Nikhar Abbas and Fanzhong Meng and Pietro Bortolotti and Witold Skrzypinski and George Scott and Roland Feil and Henrik Bredmose and Katherine Dykes and Matt Sheilds and Christopher Allen and Anthony Viselli},
Howpublished = {NREL/TP-75698},
institution = {International Energy Agency},
title = {{Definition of the IEA 15 MW Offshore Reference Wind Turbine}},
URL = {https://www.nrel.gov/docs/fy20osti/75698.pdf}
}


## Wind scenarios
Wind data are taken from RVO website (Netherlands ministry) as of July 2020: https://offshorewind.rvo.nl/
In the folder ```wind/```, 5 files provide wind scenarios, for the sites Hollandse Kust (Noord) ```RVO_HKN.txt```, Hollandse Kust (West) ```RVO_HKW.txt```, Hollandse Kust (Zuid) ```RVO_HKZ.txt```, IJmuiden Ver ```RVO_IJV.txt```, and Ten Noorden van de Waddeneilanden ```RVO_TNW.txt```.
For HKW, IJV, and TNW the wind data have been adapted from Metocean report, because wind report was not available yet.
Note that in the benchmark reported in the paper we used the wind scenario of Ten Noorden van de Waddeneilanden ```RVO_TNW.txt```.

Data format for each row:
```wind_direction wind_speed probability```
The first column contains the direction from which the wind is blowing, in degrees. The data we report are quite sparse, with a 30 degrees of resolution. We interpolate these data in our optimization, and recommend to reach a much higher resolution to have a proper modelling of the wakes. The second column contains the wind speed in m/s of the wind. The third columns reports the probability of seeing this wind direction with this wind speed over the year. The probabilities sum to 1.
