# Efficient Local Map Search Algorithms for the Placement of Flying Relays

These codes generate simulation results depicted in the paper titled


J. Chen and D. Gesbert, "Efficient Local Map Search Algorithms for the
Placement of Flying Relays", IEEE Trans.  Wireless Commun., 2019.

## How we modeled the city:
We modeled the city as a set of buildings, where each building is modeled as a set of cubes. You can imagine that a cylindrical building can be geometrically constructed by a large number of long cubes. Each cube can be described by 4 line segments and a height. 

Whether a UAV-user link is blocked is determined by whether the line segment joining the UAV and the user penetrate one of these cubes. 

We simulate a city by randomly dropping thousands of cubes according to some distributions that were manually designed to mimic a Manhattan-like city. We also drop tens of thousands users on the streets, where the streets are predefined strip-like locations and no cubes are dropped at the streets. The dataset is not uploaded here due to its huge volume. 

## Equivalent dataset
We upload an equivalent dataset topologyK2.mat to describe the equivalent user-city-UAV topology. The idea is the following. 

Despite there are tens of thousands of buildings (or, in our terminology, cubes), for each user, there are only a few of them being essential. It can be proven that, a UAV-user link is blocked, if and only if, the link is blocked by those "essential cubes". 

The variable "Topology" in topologyK2.mat contains 10,000 cells, where each cell represents the equivalent city topology given a specific user position. For example,

>> Topology{253}

ans = 

  struct with fields:

         Blds: {22×2 cell}
      cAngles: [1×314 double]
        PosUE: [395.1281 524.9686]
     BldLines: {22×1 cell}
     BldTypes: [22×1 double]
    BldHeight: [22×1 double]

The above shows that for the 253th user, there are 22 essential cubes. Each cube can be described by a rectangle stored in a cell of "BldLines" and a height stored in an element in "BldHeight". The user locates at "PosUE". 

- Author: Junting CHEN
- Email: juntingc@cuhk.edu.cn
- Date: Jan 2, 2020
