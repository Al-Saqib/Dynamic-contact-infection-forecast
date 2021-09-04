# Predicting hospital aquired infections using dynamic networks of contact

## Paper overview

**Citation**

Myall, A. C., Peach, R. L., Wan, Y., Mookerjee, S., Jauneikaite, E., Bolt, F., Price, J. R., Davies, F., Wiesse, A. Y., Holmes, A., & Barahona, M. (2021). Characterising contact in disease outbreaks via a network model of spatial-temporal proximity. *MedRxiv*, 2021.04.07.21254497. https://doi.org/10.1101/2021.04.07.21254497.

This article introduced, tested, and validated an infection forecasting tool for primary use in hospitals. The tool incorporates dynamically changing variables to make predictions over each patient for the risk of acquiring a disease. Specifically, our framework accounts for contact patterns, a significant driver of many infectious diseases, to capture disease aquistion risk. In the article, we deployed the forecasting tool onto hospital-onset COVID-19 infections, demonstrating both its high efficacy and generalisability.

<p align="center">
  <img src="images/fig_method_overview.png" width="400">
</p>

<sub>*Forecasting framework overview. Patient pathways are extracted from electronic health records which specify the locations each patient has visited over the duration of their hospital stay. These pathways are then overlaid with COVID-19 testing results, capturing the space-time positions of patients that tested positive for COVID-19. For predictions we define a time window and for each patient extract clinical variables (fixed) and hospital contextual variables (dynamic), and their position in a contact network (dynamic). This process of variable extraction is then repeated for multiple time windows, with total information used for model training and predictions.*</sub>


## Repo overview

This repo provides an implementable example of the model proposed in *Myall et al. 2021.*  `R/` is the folder for scripts that contain R functions. All functions are documented with [roxygen2](https://roxygen2.r-lib.org/) syntax.

### Example data

There are 2 datasets required for StEP. Firstly, is the total mobility patterns over which spatial-temporal contact is measured. Secondly are the trajectories of infected individuals under investigation.

Total mobility patterns capture how diseases will spread since individuals and their person-person are a primary vector in transmission. Routes of disease spread are often dominated by a set of most probable trajectories (*Brockmann and Helbing 2013*). However, as we previously showed, heterogeneous structure existed amongst movements patterns (*Myall et al. 2020*). Probable trajectories then must be derived based on mobility which captures general patterns of movement. Specifically, we capture these most probable trajectories, `D`, which contains the shortest paths between location nodes in terms of effective distance (Brockmann and Helbing 2013). Here we can call the function `eff_dist()`, which computes matrix `D` when provided a data frame containing the weighted directed edges of total mobility patterns.

```R
D = eff_dist(read_csv("data/background_movement.csv")) 
```

The dataset `data/background_movement.csv` is a simple example with 12 locations (wards) which the subsequent `trajectories` are recorded over:

<p align="center">
  <img src="images/background_movement_example.png" width="400">
</p>


<sub>*Network visualisation of example background mobility data*</sub>

The function `example_trajectories()` generates a long dataframe with 12 pre-set trajectories. Each trajectory is a series of locations (rows in the dataframe) with both a location (wards in our paper and example, but are network nodes) and a time component. For later computation, we also split the `trajectories` dataframe into a list of dataframes `traj.l` based on the `trajectories$patient.ID` column.


```R
trajectories = example_trajectories()

traj.l <- split(trajectories , f = trajectories$patient.ID)
```

### Computing proximity between trajectories

The function `getSpatialTempProx()` computes total spatial-temporal proximities using our kernel function and returns a weighted undirected dataframe of `edges`. We measure spatial-temporal proximity between ward-time locations for every pair of individuals using our spatial-temporal distance kernel found in *Myall et al. 2021.*. The parameter beta represents a propagation speed across background movement and is specific to the pathogen under study.


```R
edges = getSpatialTempProx(traj.l,   # list of trajectories
                 D,                  # efffective distance matrix 
                 beta = 0.6)         # paramter for speed of propergation
```
For more individuals or larger pathways the computation will take increasingly longer:

```
[1] "Computing proximities"
  |=============================================================================================           |  89%
```

### Graph construction with Continuous *k*-nearest neighbors (Cknn)

Incorporating graph structure between data points can aid classification through the emphasis of strong relationships. Hence, to reveal stronger contacts capturing transmission, we remove weak connections by sparsifying the corresponding fully connected graph made from the `edges` in `getSpatialTempProx()`. 

Although several alternative graph construction methods exist, we focus on Continuous *k*-nearest neighbors (Cknn) an extension to *k*-nearest neighbors (Berry and Timothy 2016). Implemented in the function `cknneighbors_graph`, it takes provide the fully connected `edges`, the paramter `k` for number of nearest neighbors to include, and an optinal paramter `lambda` (set by defaut to `lambda=1` if not provided), and returns a sparsified edges `edges_cknn` according to the Cknn algorithm.

```R
edges_cknn = cknneighbors_graph(k=3,             # parameter for k-nearest neighbors
                                #lambda = 1,     # data point density
                                edges = edges)   # fully connected edges
```


### Visualising final contact network

```R
# Preprocess data
netDat = preproNet(trajectories,edges_cknn)

# Visualise network
visNetwork(nodes = netDat$nodes,edges = netDat$edges)
```
<p align="center">
  <img src="images/contact_network_vis.png" width="600">
</p>



## References

<div id="refs" class="references">

<div id="ref-brockmann_2013">

Brockmann, Dirk, and Helbing, Dirk. 2013. “The hidden geometry of complex, network-driven contagion phenomena.” science 342.6164 (2013): 1337-1342. https://science.sciencemag.org/content/342/6164/1337

</div>

<div id="ref-myall_2021">

Myall, Ashleigh and Peach, Robert and Wan, Yu and Mookerjee,Siddharth and Jauneikaite, Elita and Bolt, Frankie and Price, James and Davies, Frances and Wiesse, Andrea and Holmes, Alison and Barahona, Mauricio. 2021. "Characterising contact in disease outbreaks via a network model of spatial-temporal proximity." medRxiv preprint. https://doi.org/10.1101/2021.04.07.21254497

</div>


<div id="ref-myall_2020">

Myall, Ashleigh and Peach, Robert and Weiße, Andrea and Davies, Frances and Mookerjee, Siddharth and Holmes, Alison and Barahona, Mauricio. 2020. "Network memory in the movement of hospital patients carrying antimicrobial-resistant bacteria." arXiv preprint. https://arxiv.org/abs/2009.14480v2

</div>


<div id="ref-myall_2020">

Berry, Tyrus, and Timothy Sauer. 2016. "Consistent manifold representation for topological data analysis." Foundations of Data Science. http://aimsciences.org//article/id/2556e6c9-b4b9-455a-9d9e-886ef0cd166f

</div>
