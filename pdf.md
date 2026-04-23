# Design-Based Inference for Spatial Experiments under Unknown Interference

Ye Wang, Cyrus Samii, Haoge Chang, and P. M. Aronow∗. August 6, 2024 

∗Wang is Assistant Professor, Department of Political Science, University of North Carolina, Chapel Hill, USA (Email: yewang@unc.edu). Samii (contact author) is Associate Professor, Department of Politics, New York University, New York, USA (Email: cds2083@nyu.edu). Chang is Assistant Professor, Department of Economics, Columbia University, New York, USA (Email: hc3615@columbia.edu). Aronow is Associate Professor, Departments of Statistics and Data Science, Political Science, Biostatistics and Economics, Yale University, New Haven, USA (Email: p.aronow@yale.edu). For their comments and suggestions, we thank Kirill Borusyak, Stephen Cole, Alexander Demin, Naoki Egami, Jiawei Fu, Michael Hudgens, Peter Hull, Molly Roberts, Fredrik S¨avje, Davide Viviano, and seminar participants at Harris School at University of Chicago, New York University Abu Dhabi, Princeton University, Rochester University, Texas A&M., UNC, and UCSD 

# Design-Based Inference for Spatial Experiments under Unknown Interference

# Abstract

We consider design-based causal inference for spatial experiments in which treatments may have effects that bleed out and feed back in complex ways. Such spatial spillover effects violate the standard “no interference” assumption for standard causal inference methods. The complexity of spatial spillover effects also raises the risk of misspecification and bias in model-based analyses. We offer an approach for robust inference in such settings without having to specify a parametric outcome model. We define a spatial “average marginalized effect” (AME) that characterizes how, in expectation, units of observation that are a specified distance from an intervention location are affected by treatment at that location, averaging over effects emanating from other intervention nodes. We show that randomization is sufficient for non-parametric identification of the AME even if the nature of interference is unknown. Under mild restrictions on the extent of interference, we establish asymptotic distributions of estimators and provide methods for both sample-theoretic and randomization-based inference. We show conditions under which the AME recovers a structural effect. We illustrate our approach with a simulation study. Then we re-analyze a randomized field experiment and a quasi-experiment on forest conservation, showing how our approach offers robust inference on policy-relevant spillover effects. 

Keywords: causal inference, design-based inference, experiments, interference. 

# 1 Introduction

Consider a spatial experiment where an intervention is randomly assigned to specific locations in a geographic space. Then, we observe how outcomes are distributed over this geography. Figure 1 illustrates the generic structure of such experiments. The left panel presents a hypothetical point-intervention experiment. An experimental design could treat half of such points to receive the intervention (gray shaded points), with the rest of the points remaining in a control condition without intervention (unshaded points). The shading in the background raster indicates outcome values. The right panel illustrates a similar experiment, but with interventions assigned to polygons instead of points. 

A recent example of a spatial experiment is from Jayachandran et al. (2017), who study a forest conservation intervention by comparing forest cover outcomes in villages that hosted the intervention to those that did not. (In this case, the villages are the polygons.) For this application, and others like it, a major concern is the possibility of detrimental spillover effects: while the intervention may reduce deforestation within a village in which it is applied, it may simply push the deforestation into other nearby areas. The dashed lines in Figure 1 display possible zones into which effects may bleed out in our hypothetical examples. The extent and manner in which such spillover effects occur could be shaped by how the intervention is distributed over the entire space. The spillover effects might be smaller when there is a high saturation of treatment rather than a low saturation. The situation at hand is one in which the potential outcomes at a point in space depend not only on that point’s treatment status but rather on how treatments are distributed elsewhere—what is referred to as “interference” in the causal inference literature (Cox, 1958, p. 19). 

As Reich et al. (2021) discuss, in the current literature on spatial causal inference, established ways to approach spillover effects expand the set of potential outcomes to account for different types of direct and indirect exposure to a treatment. As an alternative to paramet-


Point intervention


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/bb7a4a9776a03bfaa4c60afdd7ba5a65537a746ec1ce6507044b03c5d1f1b25f.jpg)



Polygon intervention


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/9334a06743144ce90cd541004e6b86f5905c546581ceaf139c38681f01774ec6.jpg)



Figure 1: Illustrations of hypothetical spatial experiments in which interventions are applied to points (left) or polygons (right). The background raster captures the geographic outcome data. Effects may bleed out in space, as illustrated by the concentric dashed lines.


ric model-based approaches, these potential outcomes-based approaches are design-based, meaning that they aim for robustness by leveraging what is known about the experimental design while remaining agnostic about features of the outcome distributions.1 One such approach is the “exposure mapping” method of Aronow and Samii (2017), in which one defines discrete potential outcomes for points in the geography depending on whether the points are within a treated polygon and whether the points are within some set distance from any other treated points. In circumstances where the geography is partitioned into distinct areas within which spillover effects are contained, an alternative and somewhat more agnostic approach focuses on marginalized “direct” and “indirect” effects as defined by Hudgens and Halloran (2008) or Zigler and Papadogeorgou (2018). 

The challenge in applied settings is that the appropriate exposure mapping may be unknown and spillovers may not be neatly contained within discrete areas. Spillover effects may bleed out in more complex ways and interact with each other. These challenges apply to the interference problem in non-spatial settings as well, e.g., in settings of network in-

terference, leading to recent theoretical work to consider what kinds of inference might be achievable under “unknown interference” (S¨avje et al., 2021; Li and Wager, 2022; Hu et al., 2022). We bring these ideas to the spatial context in order to approach problems such as the concern about detrimental spillovers in the forest conservation application. We show that when interference is present, rather than estimating the ATE, contrasts between treated and control areas capture a more nuanced, but nonetheless, policy-relevant quantity that we call the “average marginalized effect” (AME). We can define AMEs at different distance intervals from where the intervention is applied, thereby capturing how effects transmit in space in potentially non-monotonic ways. The AME is similar to the average direct effect of Halloran and Struchiner (1995) (see also Hudgens and Halloran (2008), VanderWeele and Tchetgen (2011), S¨avje et al. (2021), Li and Wager (2022), and Hu et al. (2022)), except that it is indexed by distance for application to the spatial case. The AME measures how, on average, outcomes within the specified distance interval from an intervention node are affected by activating a treatment at that node, taking into account ambient effects emanating from treatments at other intervention nodes. Taking into account the ambient effects also allows us to capture the consequences of more or less uniform, and more or less dense, distributions of treatments. There is a direct mapping from the AME to effects that are assumed by parametric models of spatial effects: when effects emanating from different intervention nodes are additive, the AME recovers the average of these additive effects. If spatial effects are not simply additive but exhibit complex interactions or feedback, the AME still yields an interpretable and policy-relevant quantity. This interpretation is robust and does not depend on, e.g., a parametric spatial lag, autoregressive, or weighting structure (Golgher and Voss, 2016). 

We also develop inferential methods to account for the dependencies that interference creates. We work with a Horvitz-Thompson estimator and a Hajek estimator for the AME. We show that these estimators are consistent and asymptotically normal under weak re-

strictions on the degree of interdependence induced by interference. We further prove that the commonly-used spatial heteroscedasticity and autocorrelation consistent (HAC) variance estimator of Conley (1999) provides conservative estimates for the true variance of these estimators under conditions that are often satisfied in practice. 

Our analysis is related to a few streams of current methodological research on spatial causal inference and interference. First, our approach draws most directly on recent designbased analyses of causal effects under interference that consider estimands marginalized over the randomization distribution, as in Hudgens and Halloran (2008), S¨avje et al. (2021), Papadogeorgou et al. (2020), Li and Wager (2022), and Hu et al. (2022). As in these approaches, our inference does not require that we specify the functional form of all interference-induced potential outcomes precisely (i.e., the “exposure mapping” of Aronow and Samii (2017)). It thus skirts the issue of non-overlap caused by a potentially high-dimensional set of potential outcomes (Leung, 2022). Second, our analysis is related to recent work on non-parametric estimation of spatial effects, including work on “bipartite causal inference” by Zigler and Papadogeorgou (2018) and on cluster-randomized designs by Leung (2022). These works focus on cases where points of intervention are far enough apart to yield disjoint clusters that interfere with each other minimally. Such designs are appealing, but they are not always feasible. Similar to $m$ -dependence for a time series, we assume hard limits on the extent of interference, but we do not assume that the set of units can be partitioned into a set of disjoint clusters. Third, our inferential results rely on the contributions of Ogburn et al. (2020). We also draw connections to inferential results in the spatial regression literature (Arbia, 2006; Jenish, 2016; Kelejian and Piras, 2017). We justify regression estimators on spatial data from the design-based perspective and provide causal interpretations for the coefficients. We clarify the connection to Conley (1999)’s spatial HAC variance estimator. 

We begin by developing the formal inferential setting and main theoretical results, using a toy example to illustrate concepts. We then consider extensions and refinements. We 

provide simulation evidence of the performance of our proposed estimators and then turn to the forest conservation application that we used above to motivate the analysis. 

# 2 Setting

Suppose a set of intervention nodes $\cal { S } = \{ 1 , . . . , N \}$ . Each node $i \in S$ can be either a point or a collection of points (e.g., a polygon) that resides in a two-dimensional set $\mathcal { X }$ indexed by $\boldsymbol { x } = ( x _ { 1 } , x _ { 2 } )$ (e.g., latitude and longitude). An experimental design assigns a binary treatment $Z _ { i } ~ \in ~ \{ 0 , 1 \}$ to each intervention node. The ordered vector of experimental assignment variables is $\mathbf { Z } \equiv ( Z _ { 1 } , . . . , Z _ { N } )$ , and the ex post realized assignment from the experiment is given by ${ \bf z } \equiv ( z _ { 1 } , . . . , z _ { \mathrm { N } } ) \in \{ 0 , 1 \} ^ { N }$ . The experimental design fixes the set of possible assignment vectors as well as a probability distribution over that set. Our analysis considers the case of Bernoulli randomization for each $Z _ { i }$ , i.e., (possibly weighted) coin flips to determine treatment status at each node.2 

Potential outcomes at any point $x \in \mathcal { X }$ are defined for each value of $\mathbf { z }$ , $( Y _ { x } ( \mathbf { z } ) ) _ { \mathbf { z } \in \{ 0 , 1 \} ^ { N } }$ . Given a realized treatment assignment $\mathbf { z }$ , we observe the corresponding potential outcome at $x$ : 

$$
Y _ {x} = \sum_ {\mathbf {z} \in \{0, 1 \} ^ {N}} Y _ {x} (\mathbf {z}) I (\mathbf {Z} = \mathbf {z}), \tag {1}
$$

where $I ( \cdot )$ is the indicator function and $I ( \mathbf { Z } = \mathbf { z } )$ takes on the value 1 if the random assignments $\mathbf { Z }$ realize to be $\mathbf { z }$ , and 0 otherwise. Data for points in $\mathcal { X }$ may come in various formats including raster data or data on a discrete set of points in $\mathcal { X }$ . Let $\mathbf { Y } ( \mathbf { z } ) = ( Y _ { x } ( \mathbf { z } ) ) _ { x \in \mathcal { X } }$ denote the full set of potential outcomes when $\mathbf { Z } = \mathbf { z }$ and $\mathbf { Y } = ( Y _ { x } ) _ { x \in \mathcal { X } }$ denote the full set of realized outcomes. 

We map the full set of potential outcomes $\mathbf { Y } ( \mathbf { z } )$ for all points in the outcome space $\mathcal { X }$ 

back to the intervention nodes by defining the “circle average” function. For $i \in S$ , we define 

$$
\mu_ {i} (\mathbf {Y} (\mathbf {z}); \Omega_ {d}) = \frac {\int_ {x : d _ {i} (x) \in \Omega_ {d}} Y _ {x} (\mathbf {z}) \mathrm {d} \zeta}{\int_ {x : d _ {i} (x) \in \Omega_ {d}} \mathrm {d} \zeta}. \tag {2}
$$

In the expression above, $d _ { i } ( x )$ measures the distance between point $x$ and intervention node $i$ . When $i$ is a point located at $x ( i )$ , $d _ { i } ( x ) = \gamma ( x ( i ) , x )$ , where $\gamma ( \cdot , \cdot )$ is a well-defined metric (e.g., Euclidean, geodesic, or a least-cost distance) that satisfies triangular inequality. If $i$ is a collection of points, then $\begin{array} { r } { d _ { i } ( x ) = \operatorname* { m i n } _ { x ^ { \prime } \in i } | | x ^ { \prime } - x | | } \end{array}$ , the minimal distance between $x$ and points belonging to $i$ . $\Omega _ { d }$ is a set of distance values and $\zeta$ is a suitable measure on $\mathcal { X }$ . Therefore, $\mu _ { i } ( \mathbf { Y } ( \mathbf { z } ) ; \Omega _ { d } )$ is the average outcome across points whose distance to $i$ falls in $\Omega _ { d }$ . 

If the points in $\mathcal { X }$ are dense and spaced evenly, then $\Omega _ { d }$ could be a singleton, $\Omega _ { d } = \{ d \}$ . The circle average amounts to taking the average across points along the edge of a circle of radius $d$ around $i$ . If the points are spaced such that there are few or no points precisely at the edge of the circle, $\Omega _ { d }$ could be a donut where $\Omega _ { d } = \left\lfloor d - \kappa , d \right\rfloor$ (with $\kappa$ a user-chosen constant dictating the donut’s thickness), or a disk where $\Omega _ { d } = [ 0 , d ]$ . By considering a collection of sets across different d values, $\{ \Omega _ { d } \} _ { d \in \mathcal { D } }$ , we will be able to examine how the circle average’s value varies over the geography. When it does not cause confusion, we write $\mu _ { i } ( \mathbf { Y } ( \mathbf { z } ) ; \Omega _ { d } )$ simply as $\mu _ { i } ( \mathbf { Y } ( \mathbf { z } ) ; d )$ . Similarly, the realized circle average for intervention node $i$ at $d$ is 

$$
\mu_ {i} (\mathbf {Y}; d) = \sum_ {\mathbf {z} \in \{0, 1 \} ^ {N}} \mu_ {i} (\mathbf {Y} (\mathbf {z}); d) I (\mathbf {Z} = \mathbf {z}). \tag {3}
$$

This representation allows us to see how an experiment is a process of sampling potential circle averages for intervention nodes, and therefore allows us to apply sample theoretic results in our analysis below. 

The left plot in Figure 2 illustrates a toy example for a point intervention ( $N = 4$ ) and raster outcome data. The plot shows a “null raster” for which none of the intervention nodes has been assigned to treatment and so ${ \bf z } = ( 0 , 0 , 0 , 0 )$ . The outcomes are $Y _ { x } ( 0 , 0 , 0 , 0 )$ for all 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/d391900899abe2e64430d791b36eeb38315c09ad9f332a629ec8a0c59db4ef12.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/991bcae174557c02605e2028721fbe2d582b33edbe3767f02e65baa084bbb8b0.jpg)



Figure 2: Left: Illustration of a “null raster,” with $N = 4$ intervention nodes (points), none of which are assigned to treatment. Raster cells are colored according to outcome levels. White circles around the nodes are where circle averages are computed. Lighter colors represent larger outcome values. Right: a possible effect function that is non-monotonic in distance.


$x$ in the space. As we can see, outcomes are defined for any point $x$ in the space, although outcomes are constant within raster cells. This is a feature of the raster data. Other types of data may exhibit finer levels of granularity — e.g., data produced from kriging interpolation that varies smoothly in space. We take these outcome data, and any coarsening or smoothing operations that they incorporate, as fixed. For our design-based inference, the only source of stochastic variation is from $\mathbf { Z }$ . For data that are smoothed using kriging, one could tune smoothing parameters on auxiliary data so that they are fixed with respect to $\mathbf { Z }$ . White circles around the intervention nodes demonstrate one possible way to construct the circle average. We use the Euclidean distance and take averages across all the raster cells passed by the edge of the circles. Note that we do not prohibit circles around different nodes to intersect with each other. 

Spatial effects can exhibit considerable complexity. For the sake of illustration, our toy example supposes that treatments transmit effects that are non-monotonic in distance and that effects from different intervention nodes accumulate in an additive manner. The right plot in Figure 2 illustrates such an effect function. Then, the net result would depend on how treatments are distributed over the intervention points. Figure 3 illustrates how 

outcomes would be affected over different allocations of the treatment given that effects take the form as in Figure 2. In the analysis below, we do not assume that effects are additive or homogeneous in form—this is done here merely to provide a simple illustration. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/1e60b42f908308d685c3e6dfd6a35d3ff1c1fbb8f5f6f67e25f3a8aa2daf7664.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ee55041382e0f521030f5ab5f56a10f85450988f59e844f851483bf4269e16fe.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/23e043a0b3f619dcb0a2821b34bbe76185e20d8d12324de8a95258a069909504.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/5d223a4d69a21c3dc9dd6ea9cb1b20d229725ed20475ebc3d0af68124c3f5acc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/3e3e577adfaacea2efdfe677b77fcf8f7e95921d98bfaa43759c2946cd05ed86.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/1387b6fe6c0303053304105e2ad3e27f8d1fbdeef6d364fd78faf2c1ee164fb5.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/a9e9dd03fbb40ec41a3ac98e247a1e0375fadfc4c54962309c1da81cd948c848.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/86e55d0240ffe3b8f3d3afa6d976ae3c2d9f1dba85cf953328c48de8b2cafc99.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/775e9d905332d62e3ef641a1202b43501100c843b48c2cf8c6184c6b51eda71c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/75ca1230102931ec6c6a08c51e4b958f0a9f35e00d17108bb31db6757169399b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b3365648f3721c8b017c6fb18bdaaf34703d14182ec4fae7327a1f854d39172e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/0e9e2a365ef2c155232cd2b531e4785682744b5f1d86b4a7ae532017e59fc22d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/474916b5d68c1ccac92e7db2c95e2cf3e5f5863da61499d95285cb304e8dd77a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b2e237f6e6e6195430a989a94cb1ff9d9123c8d41843ddbeb9ce5f1b6f41db42.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/a2bee42ee9c8197f04d030f0c73e20973b708d0d58f988154fa4dfc21304ecb4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/a7d1f87150b5c320c8c016037b89ea465c9fd568fabc802ffdad75883ee9afe8.jpg)



Figure 3: Illustration of how outcomes are affected given different treatment allocations given the effect function in Figure 2. Treated intervention points are white, while nontreated intervention points are black.


# 3 Defining a marginal spatial effect

As the potential outcome notation for $Y _ { x } ( \mathbf { z } )$ indicates, the outcome at any point may depend on the full vector of realized treatment assignments $\mathbf { z }$ . Similarly, the potential outcome notation for the circle average, $\mu _ { i } ( \mathbf { Y } ( \mathbf { z } ) ; d )$ , shows that the realized circle average for node $i$ may depend on treatment assignments for nodes other than $i$ . As such, the circle averages are potentially subject to causal “interference.” 

We now define a spatial effect that we call the “average marginalized effect” (AME). The AME is a marginal effect that accounts for interference. It can be defined for any distance value $d$ and enables researchers to examine how the impacts of the intervention nodes on the outcome of their neighboring points vary in space. The usual definition of a unit-level treatment effect takes the difference between a unit’s potential outcomes under one treatment condition versus under another treatment condition. A unit-level marginal effect is different because it takes the difference between the average of a unit’s potential outcomes over a set of potential outcomes versus the average over another set. We apply this idea to the spatial setting. In doing so, we consider effects that may bleed out in ways that are not necessarily contained within pre-defined strata, as in Hudgens and Halloran (2008), or summarized by a simple statistic, as in Aronow and Samii (2017). 

To define the spatial AME, let us first rewrite the potential outcome at point $x$ as $Y _ { x } ( z _ { i } , { \bf z } _ { - i } )$ , where $\mathbf { z } _ { - i }$ is a vector equaling $\mathbf { z }$ except that the value for intervention node $i$ is omitted. This allows us to pay special attention to how variation in treatment at node $i$ relates to potential outcomes at point $x$ , given the variation in treatment values in $\mathbf { z } _ { - i }$ . We can marginalize over variation in $\mathbf { z } _ { - i }$ to define an “individualistic” average of potential 

outcomes for point $x$ , holding the treatment at intervention node $i$ to treatment value $z$ : 

$$
Y _ {i x} (z _ {i}; \eta) = \operatorname {E} _ {\mathbf {Z} _ {- i}} [ Y _ {x} (z _ {i}, \mathbf {Z} _ {- i}) ] = \sum_ {\mathbf {z} _ {- i} \in \{0, 1 \} ^ {N - 1}} Y _ {x} (z _ {i}, \mathbf {z} _ {- i}) \Pr (\mathbf {Z} _ {- i} = \mathbf {z} _ {- i}; \eta), \tag {4}
$$

where $\eta$ is an experimental design parameter that is an index for the distribution of $\mathbf { Z }$ (that is, the probability of treatment assignments). This is the individualistic marginal potential outcome at point $x$ given that node $i$ is assigned to treatment condition $z$ , marginalizing over possible assignments to other nodes. We can use Figure 3 to illustrate. To construct $Y _ { 1 x } ( 0 ; \eta )$ , one would take a weighted average of the potential outcomes at point $x$ under assignments labeled in the figure as Z1, Z3, Z4, Z5, Z9, Z10, Z11, and Z15, where the weights would be proportional to the probability of each assignment. 

We can define a similar marginal quantity at the level of the circle averages: 

$$
\mu_ {i} (z _ {i}; d, \eta) = \operatorname {E} _ {\mathbf {Z} _ {- i}} [ \mu_ {i} (\mathbf {Y} (z _ {i}, \mathbf {Z} _ {- i}); d) ] = \sum_ {\mathbf {z} _ {- i} \in \{0, 1 \} ^ {N - 1}} \mu_ {i} (\mathbf {Y} (z _ {i}, \mathbf {z} _ {- i}); d) \Pr (\mathbf {Z} _ {- i} = \mathbf {z} _ {- i}; \eta), \quad (5)
$$

where we use $\mathbf { Y } ( z _ { i } , \mathbf { z } _ { - i } )$ to denote the vector of potential outcomes over points in $\mathcal { X }$ that obtain under treatment assignment $\left( z _ { i } , \mathbf { z } _ { - i } \right)$ . This is the potential circle average at distance $d$ around node $i$ , given that $i$ is assigned to treatment condition $z _ { i }$ , marginalizing over possible assignments to other nodes. 

We can now define an individual marginalized effect at point $x$ of intervening on node $i$ , allowing other nodes to vary as they otherwise would under $\eta$ : 

$$
\tau_ {i x} (\eta) = Y _ {i x} (1; \eta) - Y _ {i x} (0; \eta). \tag {6}
$$

This defines the response at point $x$ of switching node $i$ from no treatment to active treatment, averaging over possible treatment assignments to nodes other than $i$ . At the level 

of circle averages, we can define 

$$
\tau_ {i} (d; \eta) = \mu_ {i} (1; d, \eta) - \mu_ {i} (0; d, \eta), \tag {7}
$$

which is the average of individual responses for points along the circle at distance $d$ around node $i$ . Using Figure 3 to illustrate, one would construct $\tau _ { 1 } ( d ; \eta )$ by working with the $d$ - radius circle averages around intervention node 1, taking the difference between the mean of the circle averages under assignments Z2, Z6, Z7, Z8, Z12, Z13, Z14, and Z16 minus the mean of circle averages under assignments Z1, Z3, Z4, Z5, Z9, Z10, Z11, and Z15. 

Finally, define the average marginalized effect (AME) for distance $d$ by taking the mean over the intervention nodes: 

$$
\operatorname {A M E} (d; \eta) = \frac {1}{N} \sum_ {i = 1} ^ {N} \tau_ {i} (d; \eta) \tag {8}
$$

The interpretation of the AME for distance $d$ is the average effect of switching a node $i \in S$ to treatment on points at distance $d$ from that node, marginalized over possible realizations of treatment statuses in other intervention nodes. The distribution of these possible realizations of treatment statuses depends on the experimental design. When $d = 0$ , the AME captures the direct effect generated by the treatment at the location of intervention, in a way similar to the “expected average treatment effect” in S¨avje et al. (2021). For $d > 0$ , the AME resembles the “average indirect causal effect” in Hu et al. (2022) but is defined for specific distance values. In practice, researchers may select a series of distance values, $\{ d _ { l } \} _ { l = 1 } ^ { L }$ , based on the resolution of $\mathcal { X }$ and the potential magnitude of spillover effects. The resulting collection of AMEs demonstrates how effects vary over the distance from an intervention node. 

Before ending this section, we note that our analysis focuses on experimental designs with Bernoulli assignment, in which case the possible assignments consists of the $2 ^ { N }$ possible 

vectors that could be obtained from $N$ (possibly differentially weighted) coin flips. This allows for a clean definition of causal effects (S¨avje et al., 2021). This is because Bernoulli assignment ensures that $( 1 , { \bf z } _ { - i } )$ and $( 0 , { \bf z } _ { - i } )$ each has positive probability of occurring. In this case, the marginal quantities $Y _ { i x } ( 1 ; \eta )$ and $Y _ { i x } ( 0 ; \eta )$ are defined by marginalizing over the same sets of $\mathbf { z } _ { - i }$ values, and the individualistic response has a clear ceteris paribus interpretation. Things are different under completely randomized assignment, where a fixed number $N _ { 1 }$ of nodes are assigned to treatment. Then, for $Y _ { i x } ( 1 ; \eta )$ , one marginalizes over assignments with $N _ { 1 } - 1$ units assigned to treatment, while for $Y _ { i x } ( 0 ; \eta )$ , one marginalizes over assignments with $N _ { 1 }$ units assigned to treatments. As $N$ grows, differences between AMEs in Bernoulli and complete random assignment typically become negligible when interference is local, as shown in S¨avje et al. (2021). 

# 4 Inferential assumptions

In this section, we lay out assumptions on the experimental design and potential outcomes, including restrictions on the extent of interference for the inferential results in Section 5. Our asymptotic analysis considers a sequence of sets indexed by the sample size $N$ . The set of intervention nodes is denoted as $S _ { N }$ and the set of outcome points $\mathcal { X }$ . Note that assumptions below are assumed to hold uniformly for all large sample sizes. We begin with the following assumptions: 

C 1. (Bernoulli design) $\left( Z _ { 1 } , . . . , Z _ { N } \right)$ is a vector of independent Bernoulli(p) draws. 

C 2. (Bounded potential outcomes) $| Y _ { x } ( \mathbf { z } ) | < b$ for some finite constant b and all $x \in \mathcal { X }$ and $\mathbf { z } \in \{ 0 , 1 \} ^ { N }$ . 

Assumption C1 defines the experimental design. As discussed above, condition C1 ensures that individualistic responses are ceteris paribus for variation in treatment assignment 

at a given node. We work with the assumption that the assignment probability, $p$ , is constant over intervention nodes and discuss the extension to cases where assignment probabilities vary in Section 6.4 below. Assumption C2 is a common regularity condition on the potential outcomes. It ensures the boundedness of higher-order moments for the distribution of functions of the potential outcomes. 

Our next assumption follows S¨avje et al. (2021) by using a dependency graph to characterize interference-induced dependencies among the circle averages defined at a specific distance value $d$ . Let $I _ { i j } ( d )$ be an indicator for whether assignment at intervention node $j$ interferes with the $d$ -radius circle average at node $i$ : 

$$
I _ {i j} (d) = \left\{ \begin{array}{l l} 1 & \text {i f} \mu_ {i} (\mathbf {Y} (\mathbf {z}); d) \neq \mu_ {i} (\mathbf {Y} (\mathbf {z} ^ {\prime}); d) \text {f o r s o m e} \mathbf {z}, \mathbf {z} ^ {\prime} \in \{0, 1 \} ^ {N} \text {s u c h t h a t} \mathbf {z} _ {- j} = \mathbf {z} _ {- j} ^ {\prime} \\ 1 & \text {i f} i = j, \\ 0 & \text {o t h e r w i s e .} \end{array} \right. \tag {9}
$$

Then, let $s _ { i j } ( d )$ be an indicator for whether $d$ -radius circle averages at $i$ and $j$ are subject to interference from treatment at some common intervention node $\ell$ (which could be $i$ , $j$ , or some other third intervention node): 

$$
s _ {i j} (d) = \left\{ \begin{array}{l l} 1 & \text {i f} I _ {i \ell} (d) I _ {j \ell} (d) = 1 \text {f o r s o m e} \ell \in \mathcal {S} _ {N}, \\ 0 & \text {o t h e r w i s e .} \end{array} \right. \tag {10}
$$

If $s _ { i j } ( d ) = 1$ , then circle averages at $i$ and $j$ will vary together whenever there is variation in treatment values at the relevant ℓs, meaning non-independence over possible values of $\mathbf { Z }$ . 

Using this dependency graph, our third assumption is a restriction on the extent of interference dependencies for circle averages at $d$ . Let us denote the distance between two intervention nodes $i$ and $j$ as $d _ { i j }$ .3 Then we have: 

C 3. (Local interference.) Let $h : [ 0 , \infty ) \to [ 0 , \infty )$ be a function independent of sample sizes. For each d and all large sample sizes $N$ , and all pairs of intervention nodes $i$ and $j$ in $S _ { N }$ , there exists a constant $h ( d )$ such that if $d _ { i j } > h ( d )$ , then $s _ { i j } ( d ) = 0$ . 

Assumption C3 means that there are hard limits to the spatial extent of the interference: nodes that are beyond some distance from each other have no interference-induced dependencies. Assumption C3 is an assumption on the possible extent of spillovers. For intervention node $j$ to satisfy C3 with respect to $i$ , it would require that the outcomes $Y _ { x } ( \mathbf { z } )$ used to construct the circle average of $i$ at distance $d$ is unaffected by not only $j$ ’s treatment value but also the treatment values at any intervention nodes that affect the circle average of $j$ at distance $d$ . For example, consider a simple case where the outcomes used to construct the circle average of $j$ at distance $d$ depend on $Z _ { j }$ and another assignment variable $Z _ { k }$ . If $s _ { i j } ( d ) = 0$ , this implies that the outcomes used to construct the circle average of $i$ at distance $d$ cannot depend on either $Z _ { j }$ or $Z _ { k }$ . In other words, $i$ and $j$ cannot share any sources of variation in their circle averages at distance $d$ . This assumption hence implies that for $i , j$ with $d _ { i j } > h ( d )$ , $\mu _ { i } ( Y ( z , { \bf Z } _ { - i } ) ; d )$ and $\mu _ { j } ( Y ( z ^ { \prime } , \mathbf { Z } _ { - j } ) ; d )$ are independent with each other. 

To give a more interpretable sufficient condition, suppose that the intervention at every node has no influence on the outcome points more than $d$ away from it. Then, C3 would be satisfied if the distance between the boundary of the radius- $d$ circle around node $i$ and the radius-d circle around node $j$ is larger than $2 d$ . In this case, $h ( d )$ in C3 can be set to $2 d + 2 d$ . When nodes $i$ and $j$ are $2 d + 2 d$ apart, they neither interfere with each other, nor do they have any common neighbors that influence both nodes. To further illustrate the idea, a simple example is given in Figure 4. In this example with three intervention nodes, intervention node 2 is on the radius- $d$ circle centering node 3. Because nodes 1 and 3 are sufficiently far apart, neither node 3 and node 2 affect the d-circle average of node 1 and 

$$
\min  _ {x \in i, x ^ {\prime} \in j} \gamma (x, x ^ {\prime}).
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/c759d4eda4c5de78d2c56e1568acbd00496a1c48c629044a4b935f5a018249be.jpg)



Figure 4: Illustration of condition C3. White circles around the nodes are circle averages defined at the distance value $d$ . Black circles with a radius of $d + d$ depict the maximal range that interference can happen for the $d$ -circle averages. As the distance between node 1 and node 3 is larger than $2 \bar { d } + 2 d$ , the circle averages of the two nodes do not depend on each other, even when a third node lies between them.


vice versa. 

A final assumption defines an increasing-domain asymptotic growth process in which the number of independent pairs of intervention nodes increases. Define $S _ { N } ( i ; d ) \equiv \{ j \in { \mathcal { S } } _ { N } :$ $d _ { i j } \leq d \}$ , the set of intervention nodes whose distance to node $i$ is less than $d$ . We have the following asymptotic restriction on the spacing of the intervention nodes: 

C 4. (Intervention node spacing) Let $b : [ 0 , \infty ) \to [ 0 , \infty )$ be a function independent of sample sizes. For each d and all large sample sizes $N$ , the sequence of intervention nodes satisfies $\begin{array} { r } { \operatorname* { s u p } _ { i \in { \cal S } _ { N } } | { \cal S } _ { N } ( i ; d ) | \le b ( d ) . } \end{array}$ 

C4 ensures that as the size of intervention node grows, the number of intervention nodes that reside within a given distance of a node is bounded. It is satisfied when the intervention nodes are deliberately chosen such that they are adequately spaced out geographically. For 

example, for point-intervention experiments, we can choose the set of nodes $S _ { N }$ from a meshgrid where the distance between any two points $i$ and $j$ on the grid is bounded from below, i.e. $d _ { i j } \geq d _ { 0 }$ . In practice, researchers can first divide the space into disjoint areas and select one intervention node from each area to make C4 plausible.4 For polygon-intervention experiments, one may require that the size of each polygon is larger than a threshold value (thus ensuring adequate spacing between non-adjacent polygons).5 

We note that C3 and C4 should not be confused with each other. C3 is an assumption on the extent of spillover effects. C4 is an assumption on the spacings of the intervention nodes which prevents nodes from concentrating densely in a particular region. To see how the two assumptions fit together, we complete our specification of the extent of interference. Let us define the neighborhood $B ( i ; d )$ that includes all the nodes whose circle averages at $d$ may interfere with that of node $i$ : 

$$
\mathcal {B} _ {N} (i; d) = \{j \in \mathcal {S} _ {N}: d _ {i j} \leq h (d) \}. \tag {11}
$$

From C3, we know that $s _ { i j } ( d ) = 0$ for $j \notin B _ { N } ( i ; d )$ . From C4, we know that $| B _ { N } ( i ; d ) | \le$ $b ( h ( d ) )$ . Define $c _ { i } ( d ) = | B _ { N } ( i ; d ) |$ and $c _ { N } ( d ) = \operatorname* { m a x } _ { i \in { \cal S } _ { N } } c _ { i } ( d )$ . Conditions C3 and C4 imply the following condition: 

C 4a. (Limited local interference with respect to $d$ ) For each $d$ , $c _ { N } ( d ) = { \cal O } ( 1 )$ . 

We note that this asymptotic implication is only pointwise in $d$ : for each $d$ , we have $c _ { N } ( d ) =$ $O ( 1 )$ but it is not true that $\begin{array} { r } { \operatorname* { s u p } _ { d \in \mathbb { R } } c _ { N } ( d ) = O ( 1 ) } \end{array}$ . In practice, researchers should focus on 

distance values within a moderate range, such that the area covered by any circle is not excessive relative to the whole geography. In our analyses below, we suppress the subscripts for asymptotic sequences unless they are needed to add clarity. 

# 5 Estimation and inference

In this section, we study the Horvitz-Thompson (HT) and Hajek estimators for the AME and establish their consistency and asymptotic normality. We recommend the use of the Hajek estimator, as it is usually more efficient in practice compared with the HT estimator.6 We express the Hajek estimator as a regression estimator and propose a variance estimator based on the spatial heteroscedasticity and autocorrelation consistent (spatial HAC) estimator of Conley (1999). Hence, our proposed Hajek estimator for the AME is equivalent to a regression of the circle average on a constant term and the treatment indicator. The inference is carried out using a spatial HAC standard error estimator and normal approximation. All proofs are contained in the appendix. 

Consider the following Horvitz-Thompson (HT) estimator: 

$$
\widehat {\tau} _ {\mathrm {H T}} (d) = \frac {1}{N p} \sum_ {i = 1} ^ {N} Z _ {i} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) \mu_ {i} (\mathbf {Y}; d). \tag {12}
$$

The terms on the right-hand side consist of design parameters $N$ and $p$ , assignment indicators $\{ Z _ { i } \} _ { i = 1 } ^ { N }$ , and observed circle average $\{ \mu _ { i } ( \mathbf { Y } ; d ) \} _ { i = 1 } ^ { N }$ as defined in (2). Hence the quantity is computable from observed data alone. 

Our first two results show that ${ \widehat { \tau } } _ { \mathrm { H T } } ( d )$ is unbiased for the AME at distance $d$ under C1, and is consistent and asymptotically normal under C1-C4. 

Proposition 1 (Unbiasedness). Under C1, 

$$
\operatorname {E} _ {\mathbf {Z}} \left[ \widehat {\tau} _ {\mathrm {H T}} (d) \right] = \operatorname {A M E} (d; \eta), \tag {13}
$$

where the expectation is taken over the random assignment variables. 

Let $N ( 0 , 1 )$ denote the standard Gaussian distribution with mean 0 and variance 1. 

Proposition 2 (Asymptotic Distribution for the Horvitz-Thompson estimator). Under C1- C4 and if $N \times \mathrm { V a r } \left( \widehat { \tau } _ { \mathrm { H T } } ( d ) \right)$ is uniformly bounded below for all large $N$ , then, as $N  \infty$ ， 

$$
\frac {\widehat {\tau} _ {\mathrm {H T}} (d) - \operatorname {A M E} (d ; \eta)}{\sqrt {\operatorname {V a r} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right)}} \xrightarrow {d} N (0, 1), ^ {7} \tag {14}
$$

The Hajek estimator is a refinement to the Horvitz-Thompson estimator and, in this setting, is equivalent to a difference-in-means estimator: instead of using $p N$ and $( 1 - p ) N$ as the denominator, the Hajek estimator replaces them with $\begin{array} { r } { N _ { 1 } = \sum _ { i = 1 } ^ { n } Z _ { i } } \end{array}$ and $N _ { 0 } = N - N _ { 1 }$ , respectively: 

$$
\widehat {\tau} _ {\mathrm {H A}} (d) = \frac {1}{N _ {1}} \sum_ {i = 1} ^ {N} Z _ {i} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N _ {0}} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) \mu_ {i} (\mathbf {Y}; d) \tag {15}
$$

The Hajek estimator is usually more efficient in pracitce, in terms of variances , to a Horvitz-Thompson estimator. 8 

Proposition 3 (Asymptotic Distribution for the Hajek estimator). Let AVar ( $\widehat { \tau } _ { \mathrm { H A } } ( d ) )$ denote the asymptotic variance of the Hajek estimator.9 Under C1-C4 and if $N \times \mathrm { A V a r } ( \widehat { \tau } _ { \mathrm { H A } } ( d ) )$ is 

uniformly bounded below for all large $N$ , then, as $N \to \infty$ 

$$
\frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \operatorname {A M E} (d ; \eta)}{\sqrt {\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)}} \xrightarrow {d} N (0, 1). \tag {16}
$$

As a simple difference in means, the Hajek estimator is algebraically equivalent to a least square regression of the circle averages on an intercept and treatment indicators of the intervention nodes: 

$$
\left( \begin{array}{l} \widehat {\mu} _ {0} (d) \\ \widehat {\tau} _ {\mathrm {H A}} (d) \end{array} \right) = \arg \min  _ {(\mu_ {0}, \tau)} \sum_ {i = 1} ^ {N} \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {0} - \tau Z _ {i}\right) ^ {2}. \tag {17}
$$

Our approach to variance estimation borrows from the spatial econometrics literature and works with the spatial heteroskedasticity and autocorrelation consistent (spatial HAC) variance estimator of Conley (1999). This estimator takes the form, 

$$
\widehat {\Sigma} _ {\mathrm {H A C}} (d) = \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1} \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \mathbf {X} _ {i} ^ {\prime} \mathbf {X} _ {j} \hat {e} _ {i} \hat {e} _ {j} K \left(\frac {d _ {i j}}{\tilde {d}}\right)\right) \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1}, \tag {18}
$$

where $\mathbf { X } = \left( \begin{array} { l l l } { 1 } & { 1 } & { \ldots , 1 } \\ { } & { } & { } \\ { Z _ { 1 } } & { Z _ { 2 } } & { \ldots Z _ { N } } \end{array} \right) ^ { \prime } \in \mathbb { R } ^ { N \times 2 }$ and $\mathbf { X } _ { i }$ denotes the $i$ th row of the matrix. The ${ \hat { e } } _ { i }$ ’s are the residuals from the regression, where $\hat { e } _ { i } = \mu _ { i } ( \mathbf { Y } ; d ) - \widehat { \mu } _ { 0 } ( d ) - \widehat { \tau } _ { \mathrm { H A } } ( d ) Z _ { i }$ . $K ( \cdot )$ is a kernel function. $\hat { d }$ is a cutoff value and our setup suggests setting it at $\tilde { d } = h ( d )$ . The (2,2)-entry of the estimator $\widehat { \Sigma } _ { \mathrm { H A C } } ( d )$ is our estimator for the variance of $\widehat { \tau } _ { \mathrm { H A } } ( d )$ and we denote it as $\widehat { \mathrm { V } } _ { \mathrm { H A C } } ( d )$ . In practice, $h ( d )$ is unknown and we can examine the robustness of the results by varying $\hat { d }$ . In the appendix, we show that the regression estimator combined with the spatial HAC variance estimator with the uniform kernel provides asymptotically valid inference for the AME under an extra assumption: 

C 5. (Homophily in treatment effects) $\begin{array} { r } { \frac { 1 } { N } \sum _ { i = 1 } ^ { N } ( \tau _ { i } ( d ; \eta ) - \mathrm { A M E } ( d ; \eta ) ) \sum _ { j \in \{ i \} \cup \mathcal { B } ( i ; d ) } ( \tau _ { j } ( d ; \eta ) - } \end{array}$ AME $( d ; \eta ) ) \ge 0$ for each value $d \geq 0$ .10 

The assumption is that the expected treatment effect generated by node $i$ at distance $d$ is positively correlated with that generated by its neighbors in $B ( i ; d )$ and itself. In other words, there is homophily in treatment effects in space: nodes that generate larger-thanaverage effects reside close to each other. We consider that a positive spatial correlation assumption reasonable in many applied settings so we recommend the use of HAC variance estimator in practice. 

When researchers have concern over C5, one can use an alternative variance-bound estimator proposed in S¨avje et al. (2021): 

$$
\widehat {\mathrm {V}} _ {\mathrm {S A H}} (d) = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {Z _ {i} c _ {i} (d) \widehat {e} _ {i} ^ {2}}{p ^ {2}} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {(1 - Z _ {i}) c _ {i} (d) \widehat {e} _ {i} ^ {2}}{(1 - p) ^ {2}}. \tag {19}
$$

The validity of this variance estimator does not depend on C5, although the estimates tend to be overly conservative in realistic dataset. We provide a comparison using the simulation exercises below. 

We summarize all the inferential results in the following proposition. For some $\alpha < 1$ , let $z _ { \frac { \alpha } { 2 } }$ and $z _ { 1 - \frac { \alpha } { 2 } }$ be the $\scriptstyle { \frac { \alpha } { 2 } } \mathrm { t h }$ and $\scriptstyle ( 1 - { \frac { \alpha } { 2 } } ) { \mathrm { t h } }$ quantiles of the standard normal distribution, respectively. 

Proposition 4. Let $K ( \cdot )$ in (18) be the uniform kernel.11 Under C1-C4 and if $N \times$ AVar $( \widehat { \tau } _ { \mathrm { H A } } ( d ) )$ is uniformly bounded below for all large $N$ , we have, for each $\alpha < 1$ , 

$$
(i) \lim _ {N \to \infty} \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\widehat {\mathrm {V}} _ {\mathrm {S A H}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) \geq 1 - \alpha ;
$$

$$
\mathrm {(i i)} \mathrm {a d d i t i o n a l l y u n d e r C 5 ,} \lim _ {N \to \infty} \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) \geq 1 - \alpha .
$$

One possible concern for the HAC variance estimator is the possibility of negative estimates. This occurs in our simulation when estimating AMEs with large d values. For regular grids like $\mathbb { Z } ^ { 2 }$ and Euclidean distance metric, one can design positive semidefinite HAC variance estimators (Newey and West, 1987; Conley, 1999). However, for irregular grids and arbitrary distance metric, we are not aware of a general method for creating an exact positive semidefinite HAC variance estimator. Nonetheless, it is possible to design a biased-upward positive definite estimator.12 The procedure is outlined below. 

Note that the HAC variance estimator can be expressed as 

$$
\widehat {\Sigma} _ {\mathrm {H A C}} (d) = \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1} \left(\mathbf {X} ^ {\prime} \left(\left(\widehat {e} \widehat {e} ^ {\prime}\right) \circ \mathcal {K}\right) \mathbf {X}\right) \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1}, \tag {20}
$$

where ${ \widehat { e } } = ( { \widehat { e } } _ { 1 } , . . . . , { \widehat { e } } _ { n } )$ , $\mathcal { K }$ is a $n$ -by- $n$ symmetric matrix with $\begin{array} { r } { \mathcal { K } _ { i j } = K ( \frac { d _ { i j } } { \tilde { d } } ) } \end{array}$ , and $\cup$ denotes the pointwise (Hadamard) matrix product. Denote the eigenvalue decomposition of the matrix $\mathcal { K }$ as $\begin{array} { r } { K = \sum _ { i = 1 } ^ { n } \lambda _ { i } v _ { i } v _ { i } ^ { \prime } } \end{array}$ , where $\{ \lambda _ { i } \} _ { i = 1 } ^ { n }$ are the eigenvalues and $\{ v _ { i } \} _ { i = 1 } ^ { n }$ are the eigenvectors. We define $\begin{array} { r } { \mathcal { K } ^ { P D } = \sum _ { i = 1 } ^ { n } \operatorname* { m a x } \{ \lambda _ { i } , 0 \} v _ { i } v _ { i } ^ { \prime } } \end{array}$ , and the corresponding positive semidefinite variance estimator 

$$
\widehat {\Sigma} _ {\mathrm {H A C}} ^ {\mathrm {P D}} (d) = \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1} \left(\mathbf {X} ^ {\prime} \left(\left(\widehat {e} \widehat {e} ^ {\prime}\right) \circ \mathcal {K} ^ {\mathrm {P D}}\right) \mathbf {X}\right) \left(\mathbf {X} ^ {\prime} \mathbf {X}\right) ^ {- 1}, \tag {21}
$$

We note that $( \widetilde { e e } ) \circ \mathcal { K } ^ { \mathrm { P D } }$ is a positive semidefinite matrix, and hence so is the estimator $\widehat { \Sigma } _ { \mathrm { H A C } } ^ { \mathrm { P D } } ( d )$ .13 In addition, $\widehat { \Sigma } _ { \mathrm { H A C } } ^ { \mathrm { P D } } ( d ) \geq \widehat { \Sigma } _ { \mathrm { H A C } } ( d )$ .14 We shall refer to this positive semidefinite 

variance estimator as the HAC-PD estimator. We investigate its performance in the simulation section. In general, the HAC-PD estimator returns nonnegative variance estimates with little loss in efficiency. 

In addition to various variance estimators, we include a discussion of empirical degree of freedom (edof) adjustment in Appendix A.7. We find in our simulations that the edof adjustment is important for improving the finite sample performance of the confidence intervals. Since the derivation is fairly standard (Imbens and Kolesar, 2012; Bell and McCaffrey, 2002; Young, 2015), we leave the derivation in the appendix. 

# 6 Extensions

# 6.1 Structural Interpretation of the AME

Recall that the AME can be interpreted as the average effect of switching an intervention node from control to treatment, given ambient interference emanating from other intervention nodes. The degree of such ambient interference is dictated by the experimental design and in particular the level of treatment saturation ( $p$ ). Generally speaking, the AME is not invariant with respect to the experimental design. Here we show that the AME can have a structural interpretation (i.e., invariant over designs) if spatial effects are additive. This particular case aligns with standard model-based spatial analyses (Darmofal, 2015). 

Suppose that for each outcome node $x$ , its potential outcome value is generated additively: 

$$
Y _ {x} (\mathbf {Z}) = \sum_ {i = 1} ^ {N} Z _ {i} g _ {i} (x) + f (x), \tag {22}
$$

where $f ( x )$ captures spatial trends in the absence of any intervention, and then $g _ { i } ( x )$ captures effects that emanate, perhaps idiosyncratically, from each of the intervention nodes. This model covers a wide variety of more restrictive models of homogeneous spatial effects. 

Under this restriction on the potential outcomes, we have that the effect of assigning treatment to an intervention node $i$ shifts outcomes at point $x$ by $g _ { i } ( x )$ : 

$$
\begin{array}{l} \tau_ {i x} (\eta) = \mathrm {E} _ {\mathbf {Z} _ {- \mathrm {i}}} \left[ Y _ {x} (1, \mathbf {Z} _ {- i}) \right] - \mathrm {E} _ {\mathbf {Z} _ {- \mathrm {i}}} \left[ Y _ {x} (0, \mathbf {Z} _ {- i}) \right] \\ = \operatorname {E} _ {\mathbf {Z} _ {- \mathrm {i}}} \left[ g _ {i} (x) + \sum_ {j \neq i} ^ {N} Z _ {j} g _ {j} (x) + f (x) \right] - \operatorname {E} _ {\mathbf {Z} _ {- \mathrm {i}}} \left[ \sum_ {j \neq i} ^ {N} Z _ {j} g _ {j} (x) + f (x) \right] \tag {23} \\ = g _ {i} (x). \\ \end{array}
$$

The circle-average at distance $d$ from intervention node $i$ would be equal to the added effect that emanates from node $i$ : 

$$
\tau_ {i} (d; \eta) = \frac {\int_ {x : d _ {i} (x) \in \Omega_ {d}} g _ {i} (x) \mathrm {d} \zeta}{\int_ {x : d _ {i} (x) \in \Omega_ {d}} \mathrm {d} \zeta}. \tag {24}
$$

Unlike the general case, this does not depend on the distribution of treatments over intervention nodes other than $i$ . The AME for distance $d$ is then the average of the ways that each intervention point individually affects outcomes at distance $d$ , regardless of the treatment assignment. Thus, we can interpret the AME as a structural quantity with such an additive potential outcome model. 

# 6.2 Smoothing

The AME at a particular distance $d$ with $\Omega _ { d } = \{ d \}$ , as defined in (2), may be a noisy quantity to estimate well in practice. When the AME curve as a function of $d$ is considered to be smooth, it may be a good idea to estimate an alternative quantity, the smoothed AME at d, defined as: 

$$
\operatorname {s A M E} (d; \eta) = \int_ {\mathbb {R}} \operatorname {A M E} (t; \eta) K _ {h} \left(\frac {d - t}{h}\right) d t, \tag {25}
$$

where $K : \mathbb { R }  \mathbb { R } ^ { + }$ is a nonnegative kernel and $h$ is a user-specified bandwidth.15 The integral can be similarly defined for measures with discrete supports. This quantity is a smoothed version of the AME function, defined with respect to a chosen kernel function and the bandwidth. The donut and disk AMEs defined in Section 3 are special cases with (properly-normalized) uniform kernels. 

We can define a similar smoothing operation on observed circle means: 

$$
\mu_ {i} ^ {\mathrm {s m}} (\mathbf {Y}; d) = \int_ {\mathbb {R}} \mu_ {i} (\mathbf {Y}; t) K _ {h} \left(\frac {d - t}{h}\right) d t. \tag {26}
$$

With this definition, we can estimate the smoothed-AME using the same regression approach studied in Section 5. It should be clear that the statistical results developed in Section 5 remain valid provided that the new quantities satisfies the assumptions with proper parameters (e.g., with a larger interference neighborhood).16 

# 6.3 Randomization Tests

In our analysis above, we discussed how to construct pointwise confidence intervals for the AME values at different distance values. An alternative for testing is the randomization test, albeit under a stronger sharp null hypothesis. In addition to pointwise tests, the 

and/or consider a more refined estimation strategy: 

$$
(\widehat {\mu}, \widehat {\tau}, \widehat {\beta}, \widehat {\delta}) = \arg \min _ {(\mu , \tau , \beta , \delta)} \sum_ {i = 1} ^ {N} \sum_ {d ^ {\prime} \in \mathcal {D}} (\mu_ {i} (\mathbf {Y} (\mathbf {z}); d ^ {\prime}) - \mu - \tau Z _ {i} - \beta (d ^ {\prime} - d) - \delta Z _ {i} (d ^ {\prime} - d)) ^ {2} K _ {h} \left(\frac {d - d ^ {\prime}}{h}\right). \quad (2 8)
$$

The estimation theory for (27), based on our setup, is similar to that of (26) and we omit here. We leave the development of the estimation theory of (28) for a future study. 

randomization test can be flexibly adapted to test other type of hypothese, for example, researchers may be interested in whether effects are statistically significant on a particular interval rather than at some point. For this purpose, one can use a randomization test with test statistics $\mathrm { m a x } _ { d \in [ d _ { 1 } , d _ { 2 } ] } \widehat { \tau } _ { \mathrm { H A } } ( d )$ . 

Under the sharp null hypothesis that $Y _ { x } ( \mathbf { z } ) = Y _ { x } ( \mathbf { 0 } )$ for any $\mathbf { z }$ , we know the full distribution of potential outcomes. Denote the statistic of interest as $T ( \mathbf { Y } , \mathbf { Z } )$ . Examples include estimates of the AMEs at each distance value or the average/maximum of such estimates on an interval $[ d _ { 1 } , d _ { 2 } ]$ . As all the potential outcomes are known under the sharp null, we can resample the assignment $\mathbf { z }$ for $P$ times and calculate the corresponding $T ( \mathbf { Y } , \mathbf { Z } _ { p } )$ , $p = 1 , . . . P$ . The resampling distribution of $\{ T ( \mathbf { Y } , \mathbf { Z } _ { p } ) \} _ { p = 1 , \ldots P }$ will approximate the distribution of $T ( \mathbf { Y } , \mathbf { Z } )$ under the sharp null. As a result, rejecting the null if $\begin{array} { r } { \frac { 1 } { P } \sum _ { p } \mathbf { 1 } \{ | T ( \mathbf { Y } , \mathbf { Z } _ { p } ) | \geq T ( \mathbf { Y } , \mathbf { Z } ) \} \leq \alpha } \end{array}$ for some fixed large enough $P$ gives an $\alpha$ -level test of the sharp null hypothesis. We include simulation results for a class of randomization tests in Section 6.3. 

# 6.4 Observational studies

Our framework can be generalized to observational studies. We first comment on this extension. In many observational studies, potential outcomes, treatment status, and confounders are assumed to be drawn from a superpopulation and uncertainties are assessed with respect to all three components. In our generalization, we treat the potential outcome and confounders as fixed and only study the uncertainty arising from the treatment assignments. In this way, we focus on estimating a sample AME instead of a population AME. Estimating a population AME in our setting will require more assumptions on the data generating process of the potential outcomes. 

Recall our notations for the evaluation set $\mathcal { X }$ and potential outcomes $Y _ { x } ( \mathbf { z } )$ in Section 2. Additionally, we denote the collection of confounders for intervention node $i$ as $O _ { i }$ , and make the following independent assignment assumption. 

C 6. (Probablistic Assignment) For all sample sizes $N$ , 

(i) The random assignment variables $\{ Z _ { i } \} _ { i = 1 } ^ { N }$ are jointly independent. 

(ii) There exists a treatment probability model $p ( \cdot ) : \mathcal { O }  [ 0 , 1 ]$ such that for each $i \in S _ { N }$ 

$$
\operatorname {P r o b} \left(Z _ {i} = 1\right) = p \left(O _ {i}\right) \tag {29}
$$

(iii) For an $\epsilon \in ( 0 , \frac { 1 } { 2 } )$ and for all $i \in S _ { N }$ , Prob $\left( Z _ { i } = 1 \right) \in \left( \epsilon , 1 - \epsilon \right)$ . 

C6-(i) assumes that the treatments are assigned independently. C6-(ii) states that the $i$ th node’s treatment probability can be fully described by the confounders $O _ { i }$ , independent of the potential outcomes. The implication that the assignments are independent of the potential outcomes conditioning on the confounders is the same as what the unconfoundness assumption establishes for observational studies under a super-population assumption. C6- (ii) implies that we can model the treatment probability solely as a function of $O _ { i }$ and that the probability of being treated at node $i$ is not predicted by the potential outcomes at any evaluation point $x$ . Under C6, one can show that the Horvitz-Thompson estimator with known propensity scores, defined as17 

$$
\widehat {\tau} _ {\mathrm {H T}} ^ {\mathrm {o b s}} (d) = \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i}}{p (O _ {i})} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {1 - Z _ {i}}{p (O _ {i})} \mu_ {i} (\mathbf {Y}; d).
$$

is unbiased for AME $( d ; \eta )$ . 18 

In practice, $p ( O _ { i } )$ is unknown to researchers. In such cases, parametric methods, such as logistic regression, and nonparametric methods, such as the sieve estimator in Hirano 

et al. (2003), can be used to estimate the propensity score. In Appendix A.6, we develop a complete inferential theory for the inverse probability weighted (IPW) estimator where propensity scores are modeled using a logistic model. Results there include: i) additional assumptions; 2) asymptotic linear expansion and asymptotic distribution characterization; 3) variance estimation and inference. 

# 6.5 Weaker assumptions on the extent of interference

The limited interference assumption C3 can be relaxed to accommodate the cases where the potential outcome of one intervention node is affected by all intervention nodes but the effect decreases as the distance between the nodes increases. In Section C.1, we extend inferential results on the Hajek estimator by relaxing the limited interference assumption C3. We follow the literature on spatial near-epoch dependence (Jenish and Prucha, 2012) and provide results on root-N consistency, asymptotic normality and HAC variance estimation.19 

# 7 Simulation

In this section, we use simulated datasets to illustrate propositions introduced in the previous sections and examine the performance of inferential methods based on our analytical results. In the main text, we present simulation results of a point intervention with two different effect functions.20 For each simulation design, we run simulations with three sample sizes 64, 100 and 144. 

For the first simulation scenario, the effect function is non-monotonic and additive. Let $Y _ { x } ( 0 )$ be the control outcome at an outcome point $x$ , the outcomes are generated as: 

$$
Y _ {x} (\mathbf {Z}) = Y _ {x} (0) + \sum_ {i = 1} ^ {n} f _ {x} \left(d _ {i x}\right) Z _ {i} \tag {30}
$$

where $n$ denotes the number of intervention nodes, $Z _ { i }$ is the treatment status of $i$ th intervention node, and $d _ { i x }$ is the distance from the outcome point $x$ to the intervention node $i$ . $f _ { x } ( \cdot )$ is an effect function and is constructed by mixing the density of two gamma-distributions. 

For the second simulation scenario, the effect is interactive. The outcome at the outcome point $x$ is generated by 

$$
Y _ {x} (\mathbf {Z}) = Y _ {x} (0) + \sum_ {i = 1} ^ {n} f _ {x} \left(d _ {i x}\right) Z _ {i} + \sum_ {i = 1} ^ {n} g _ {x} \left(d _ {i x}\right) Z _ {i} Z _ {\mathcal {N} (i)}, \tag {31}
$$

where $\mathcal { N } ( i )$ denotes the intervention node that is closest to the intervention node $i$ and $g _ { x } ( \cdot )$ is an additional effect function. This design reflects the story that the treatment effect may be stronger when two nearby nodes are treated. 

The AME curves for both cases are shown in Figure 5. The left figure displays the additive-effect case and the right figure the interactive-effect case. When effects are additive, the effect curve (i.e., $f _ { x } ( \cdot )$ in (30)) and the AME curve are the same as expected. This follows from our analysis of the structural interpretation of the AME above. Th interactive effect function emanating from a treated intervention node has the same shape as the additive one only when its nearest neighbor is not treated. Otherwise, it is monotonic. Therefore, the AME curve looks like the average of two effect functions. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/de1e7f1a720fd26f1aa9a92e3346f868a0817a6fd9363fa0f1bb6bf70641dcc2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/004fbfe579a3addb5b8a3c509fc6a65e8e818960d79f952a19057d5737b2b93c.jpg)



Figure 5: This figure displays the AME curves in solid lines for both the additive case (30) and the interactive case (31). For the interactive case, the dashed lines are the marginalized effect curves when the nearest neighbor is treated and when it is not.


Figure 6 shows the Mean Squared Errors of the Hajek estimator for both additive-effect (30) and interactive-effect cases (31), with intervention-node sample sizes of 64, 100, and 144. In both cases the MSEs decrease as sample sizes increase, as predicted by our theory. 

Figure 7 and Figure 8 report coverage rates and median half-lengthes for the Hajek estimator with different confidence interval constructions in the additive-effect case and the interactive-effect case, respectively. For brevity, we only display the case with a sample size of 144. The results are illustrated for AMEs at different distance values. 

We highlight two observations from the these figures. Firstly, for small distance values, all confidence intervals have proper coverage rates. However, for large distance values, the empirical degree of correction is important to improve finite-sample performance. This happens because for large distance values, the effective sample size becomes small and some finite sample adjustment is necessary to better reflect the randomness with a small sample. Secondly, the confidence interval procedure based on the SAH variance estimator tends to be overly conservative. The confidence interval generated by the positive semidefinite HAC variance estimator is not significantly longer than the one with the HAC variance estimator. Based on these observations, we recommend researchers to use the positive semidefinite HAC variance estimator with the empirical degree of adjustment when applying our method. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b610b30e84164c4af43d13bcd7c31b946e6e8aafbff869cbc40c2fa8ddf86716.jpg)



(a) MSE for the additive case (30)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ed61264c9bbdb7049266bddffdd5f53c66cc0662d065d1112fe0fab8cd886890.jpg)



(b) MSE for the interactive case (31)



Figure 6: The left and right figures report the Mean Squared Errors of the Hajek estimator in the additive-effect case and the interactive-effect case, respectively.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/4468d4650525ba9451990140fc1530a28458bca86e3437978e866679f3170610.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/f3cd5426b9740d6b972de8bf12ce19de69337fd625041173ca87d18e88fa8cf5.jpg)



Figure 7: Point-intervention simulation results on the coverage rates and half lengthes of twosided 95% confidence intervals with the Hajek estimator and different variance estimators in the additive effect case (30). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/f4546ba9c04241da605c8e892d61ce00b3e1955ceef6d2dae8a2ae4d2e33807d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/5d6eb7a6ce4a59536c72a0b433f3cb42cb127107ee13b8cb3c09ec9a0a48e444.jpg)



Figure 8: Point-intervention simulation results on the coverage rates and half lengthes of twosided 95% confidence intervals with the Hajek estimator and different variance estimators in the interactive effect case (31). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


We also evaluate the performance of the randomization tests for testing the sharp null hypothesis, as discussed in Section 6.3. We use the AME estimator at each distance value $d$ as the test statistics and the size of our test is 5%. Figure 9 reports the rejection probability of the tests (pointwise-in-d) in the null-effect scenario,21, additive-effect scenario (30), and interactive-effect scenario (31). It can be seen that under the null-effect case, the rejection probability is around 5%. In the additive-effect and interactive-effects, the rejection probabilities are high at some locations (i.e., locations with anon-null effect) and approaches 1 as the sample size gets larger for these locations. We also evaluate the performance of the randomization tests based on the statistics $\mathrm { s u p } _ { d \in { \mathcal { D } } } \left| \widehat { \tau } _ { \mathrm { H A } } ( d ) \right|$ .22 With a sample size of 64, the rejection probability is 0.056 in the null-effect case, and 1 in the additive-effect and interactive-effect case. The rejection probabilities remain similar for larger sample sizes 100 

and 144. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/6036021037852b7a5e71c5e4c67fc31ceb687f040d2c515b83db3c6f05655293.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b35107eae33222a4710e1825abfbbffbb4df113e477e203b3136ae0987be1227.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/679c3ed28aa124163d33a16f588a1cc5e863ad17e236988b58c746911a549363.jpg)



Figure 9: The figures display rejection probabilities of the pointwise randomization tests of the sharp null hypothesis for three simulation scenarios: null effect, additive effect (30), and interactive effect (31).


# 8 Application

# 8.1 A forest conservation experiment (Jayachandran et al., 2017)

We now return to the forest conservation experiment from Jayachandran et al. (2017). The authors evaluate the effects of a “payments for ecosystems services” (PES) program based on 121 villages in Hoima and northern Kibaale districts of Uganda. 60 villages were randomly assigned to the treatment group. Private forest owners in these villages were paid to reduce deforestation on their own land over a course of two years, from 2011 to 2013. Figure 10 shows the location of each village in the experiment and its treatment status. 

A primary concern in both academic and policy discussions about PES programs is what conservation scientists refer to as “leakage,” which in forest conservation contexts refers to 

negative spillover effects such that interventions reduce deforestation in targeted locations only to lead to its increase in others (Wunder, 2008; Alix-Garcia et al., 2012; Samii et al., 2014). In the area of Uganda in which Jayachandran et al. (2017) were researching, private forest owners cleared forest for either agricultural land or for timber sales into local markets. As such, the concern would be that forest conservation in targeted areas would cause the private forest owners to shift to clearing in other nearby forests. 

Jayachandran et al. (2017) originally estimated effects assuming no interference between villages, and only measure the outcome variable (forest cover) within the sampled village boundaries. They assessed the potential for leakage by studying whether the beneficial effects were larger in areas near forest reserves, with the idea being that these would be areas in which farmers would more easily shift forest clearing from their farmland to forest reserve land. They did not find such a pattern. Jayachandran et al. also examined whether deforestation was higher in control villages that were near treated villages, and found no such pattern. This is essentially an “exposure mapping” approach, and its validity depends on proper specification of indirect exposure. 

We use the methods above to conduct another analysis that also accounts for possible leakage into areas outside the sampled village boundaries. To do this, we construct a deforestation outcome variable using forest cover data from Hansen et al. (2013) for years 2012 and 2013. We code a pixel as deforested if a pixel goes from forest coverage rate greater than 25% in 2012 to one that is below 25% in 2013. In Figure 10, the dark spots indicate where deforestation happened. To construct the circle averages, we generate buffers around each of the village polygons. The distance range for estimating the AME is set to run from 0 km to 15 km. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/c8a876c03586ee46877d9aaaed8d4a34aaa7dc5fba2cc3ba2847053cadf7f439.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/adc9e5a8fb7ed71ce74e56c9d2156f8249fd127b98ccdab771d37a6c2fd15443.jpg)



Figure 10: Both plots show the boundary of the 121 villages in Jayachandran et al. (2017). Villages with a golden boundary are treated and those with a turquoise boundary are under control. Dark spots on the map represent deforestation during the experiment. The left plot shows buffers around each treated village and the right one shows buffers around each untreated village.


Figure 11 displays results from the Hajek estimator and the smoothed Hajek estimator.23 The point estimates show decreased deforestation within the treated village boundaries (Distance = 0km), similar to the authors’ original analysis. Then, the estimates for distance values greater than 0 km capture the spatial spillovers. The issue that we seek to address is whether there is any indication of leakage—i.e., increases in deforestation within the vicinity of treated villages. Given that we do not have a strong substantive basis to select a cutoff value for the spatial HAC variance estimator, we assess the robustness of inferences by considering a range of values (2km, 5km and 10km).24 We display two-sided 90% confidence intervals. The results do not offer any indication of substantial leakage and the point estimates actually suggest some beneficial spillovers. 

Alternatively, we use a randomization test to evaluate the cumulative effects between 0 km and 5 km. The test statistic is the sum of the point Hajek estimates between 0 and 5 km. The estimated value is -0.034 (a 3.4 percentage point decline in deforested area). For the randomization test, we use 10000 random draws and the p-value is around 0.065. Thus we would reject with 90% confidence of the sharp null hypothesis. These results indicate that the net gains from the intervention were indeed beneficial, with no indication of substantial leakage. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/34ce595672397daea72dee07a6d477885c3a686ef3a7e8da2c2ab9e25c82dc6b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/f0776985fe56515a7647abd7c75460e1351591e79df2fe01fa807d9b1441e47a.jpg)



Figure 11: Two plots present results based on (Jayachandran et al., 2017). The top plot presents results for the Hajek estimator. The bottom plot presents results for the smoothed Hajek estimator with a triangular kernel and a bandwidth of 5km. We plot pointwise twosided 90% intervals based on the positive-semidefinite HAC variance estimators and empirical degree of adjustments with three different cutoffs (2km, 5km, and 10km). (U) indicates the upper end of the interval and (L) indicates the lower end of the interval.


# 8.2 A forest conservation observational study (Ferraro et al., 2011)

In this section, we detail results based on Ferraro et al. (2011), an observational study that scrutinized the efficacy of Costa Rica’s protected areas in forest conservation. The outcome is measured at the level of parcels with a size of 3 hectares. Each parcel has a value of 1 if deforestation occurred inside it before 1980. The parcels are incorporated in a raster object for analysis. Originally, the study compared deforestation in parcels within protected areas against those external, after matching them on observable characteristics. 

This approach neglects potential spillover effects on outer parcels proximate to protected areas. Our method helps to alleviate this concern. 

We construct the intervention nodes by converting the map of Costa Rica into a separate raster object. Tiles in this raster represent potential intervention nodes and are larger in size compared to parcels in which we measure the outcome. We first find all the tiles within the protected areas and those that are no more than 5 km away from the boundary of these areas. The former set of nodes are defined as treated and the latter as untreated. Then, we remove all the treated nodes that are not adjacent to untreated ones to ensure that the two sets are more comparable. We end up with 233 treated nodes and 522 untreated ones. The geographic distribution of outcomes and intervention node placement is visualized in the left panel of Figure 12. We estimate the probability for each node to be treated by running a logistic regression model of the treatment indicator on three covariates aggregated to the level of intervention nodes: soil quality, distance to the nearest road, and distance to the nearest large city. The distribution of the the propensity score estimates is presented in the right panel of Figure 12. 

In the analysis, we generate buffers around each intervention node and construct circle averages accordingly. The distance range is set to run from 0 km and 20 km. The cutoff value is set at 2km, 5km or 10km. We present results from the Hajek estimator and the smoothed Hajek estimator in Figure 13. Similar to what we observe from the replication of Jayachandran et al. (2017), deforestation activities diminished within the protected areas and nearby areas. Notably, spillover effects wane with increasing distance from the boundaries of the protected areas, ceasing to be significant past 5 km. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b2c3f415a393377d88ac16db6765df005048e951b7650f194b71b5ea3dfff510.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b5dc688643c58c6f79975a313ae5a6369371f3ea1bd22bdc5ba2abdbb69fe547.jpg)



Figure 12: The left plot shows the boundary of all the intervention nodes in our replication of Ferraro et al. (2011). Tiles with a golden boundary are treated and those with a turquoise boundary are under control. Dark spots on the map represent the occurrence of deforestation in the outcome parcels. The right plot shows the distribution of the propensity score estimates based on logistic regression.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/4ce60cf1aa312b0755e4ae52cf6fc98116b38d35cae764b27a5796403f8f1d8f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/d398497a0ca3ec93d2a2a7b5945568941858c2525f1b14db8e88695e0aca1d16.jpg)



Figure 13: Two plots present results based on (Ferraro et al., 2011). The top plot presents results from the Hajek estimator. The bottom plot presents results from the smoothed Hajek estimator with a triangular kernel and a bandwidth of 5km. We plot pointwise two-sided 90% intervals based on the positive semidefinite HAC variance estimators and empirical degree of adjustments with three different cutoffs (2km, 5km, and 10km). (U) indicates the upper end of the interval and (L) indicates the lower end of the interval.


# 9 Conclusion

When treatments are applied at locations in space, the effects may bleed out and feed back in complex ways. As a result, outcomes at any point can depend on the distribution of treatments over the space, rather than on the treatment status of, e.g., the nearest intervention site. Such effects have important implications for policy. For example, in approaches to forest conservation such “payments for ecosystems services” (PES) interventions, a primary concern is that gains in targeted areas are negated by “leakage” of negative spillover effects into non-targeted areas. To capture such effects, one needs to account for spatial interference. Standard approaches, which ignore such interference, yield conclusions about average policy impacts that may be unwarranted. 

This paper explains how one can account for such interference in a randomized spatial experiment in which available information or knowledge is limited, and so we cannot confidently specify a parametric outcome model or non-parametric “exposure mapping,” nor can we be confident that interference is neatly contained within discrete geographical regions. We show that even in this situation of “unknown interference,” we can still estimate a meaningful spatial effect—what we call the “average marginalized effect” (AME). The AME tells us what would happen, on average, if we switch an intervention node at a given distance into treatment, averaging over ambient effects emanating from other intervention nodes. We can construct AME estimates for different distances, yielding a spatial effect curve. The AME is identified under random assignment as a simple contrast. Under restrictions on the spatial extent of interference, we can estimate the AME consistently and perform accurate inference using simple difference-in-means estimators and readily-available spatial standard error estimators. 

We also develop extensions. This includes specifying conditions under which the AME can be interpreted as a structural quantity that does not depend on the experimental design. 

We offer an approach for smoothing over distance, and explain how to test hypotheses on joint effects using Fisher-style randomization test under the sharp null. 

We illustrate our approach using simulation and applications to two real-world studies on forest conservation. The examples show the soundness of our proposed methods but also point to areas for further research. These include introducing methods to increase precision, for example, through covariate adjustment, variance estimation under less restrictive conditions, and inference for joint hypotheses. 

Disclosure: The authors report there are no competing interests to declare. 

# References



Alix-Garcia, J. M., E. N. Shapiro, and K. R. Sims (2012). Forest conservation and slippage: Evidence from mexico’s national payments for ecosystem services program. Land Economics 88 (4), 613–638. 





Andrews, D. W. (1984). Non-strong mixing autoregressive processes. Journal of Applied Probability 21 (4), 930–934. 





Arbia, G. (2006). Spatial econometrics: statistical foundations and applications to regional convergence. New York: Springer. 





Aronow, P. M. and C. Samii (2017). Estimating average causal effects under general interference, with application to a social network experiment. Annals of Applied Statistics 11 (4), 1912–1947. 





Bell, R. M. and D. F. McCaffrey (2002). Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology ${ \it 2 8 ( 2 ) }$ , 169–181. 





Chang, H. (2023). Design-based estimation theory for complex experiments. arXiv preprint arXiv:2311.06891 . 





Chen, L. H., L. Goldstein, and Q.-M. Shao (2010). Normal approximation by Stein’s method. Springer Science & Business Media. 





Conley, T. G. (1999). GMM estimation with cross sectional dependence. Journal of Econometrics 92, 1–45. 





Cox, D. R. (1958). Planning of Experiments. Wiley. 





Darmofal, D. (2015). Spatial Analysis for the Social Sciences. Cambridge: Cambridge University Press. 





Davidson, J. (1994). Stochastic limit theory: An introduction for econometricians. OUP Oxford. 





Davidson, J. (2020). A new consistency proof for hac variance estimators. Economics Letters 186, 108811. 





Doukhan, P. and G. Lang (2002). Rates in the empirical central limit theorem for stationary weakly dependent random fields. Statistical inference for stochastic processes 5, 199–228. 





Ferraro, P. J., M. M. Hanauer, and K. R. Sims (2011). Conditions associated with protected area success in conservation and poverty reduction. Proceedings of the National Academy of Sciences 108 (34), 13913–13918. 





Gao, M. and P. Ding (2023). Causal inference in network experiments: regression-based analysis and design-based properties. arXiv preprint arXiv:2309.07476 . 





Gibbons, S. and H. G. Overman (2012). Mostly pointless spatial econometrics? Journal of regional Science 52 (2), 172–191. 





Golgher, A. B. and P. R. Voss (2016). How to interpret the coefficients of spatial models: Spillovers, direct and indirect effects. Spatial Demography 4, 175–205. 





Halloran, M. E. and C. J. Struchiner (1995). Causal inference in infectious diseases. Epidemiology, 142–151. 





Hansen, M. C., P. V. Potapov, R. Moore, M. Hancher, S. A. Turubanova, A. Tyukavina, D. Thau, S. V. Stehman, S. J. Goetz, T. R. Loveland, et al. (2013). High-resolution global maps of 21st-century forest cover change. science 342 (6160), 850–853. 





Hirano, K., G. W. Imbens, and G. Ridder (2003). Efficient estimation of average treatment effects using the estimated propensity score. Econometrica 71 (4), 1161–1189. 





Horn, R. A. and C. R. Johnson (2012). Matrix analysis. Cambridge university press. 





Hu, Y., S. Li, and S. Wager (2022). Average direct and indirect causal effects under interference. Biometrika. 





Hudgens, M. G. and M. E. Halloran (2008). Toward causal inference with interference. Journal of the American Statistical Association 103 (482), 832–842. 





Imbens, G. W. and M. Kolesar (2012). Robust standard errors in small samples: Some practical advice. NBER Working Paper Series 14726. 





Imbens, G. W. and D. B. Rubin (2015). Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction. Cambridge: Cambridge University Press Press. 





Jayachandran, S., J. de Laat, E. F. Lambin, C. Y. Stanton, R. Audy, and N. E. Thomas (2017). Cash for carbon: A randomized trial of payments for ecosystem services to reduce deforestation. Science 357 (6348), 267–273. 





Jenish, N. (2016). Spatial semiparametric model with endogenous regressors. Econometric Theory 32 (3), 714–739. 





Jenish, N. and I. R. Prucha (2009). Central limit theorems and uniform laws of large numbers for arrays of random fields. Journal of Econometrics 150 (1), 86–98. 





Jenish, N. and I. R. Prucha (2012). On spatial processes and asymptotic inference under near-epoch dependence. Journal of econometrics 170 (1), 178–190. 





Kelejian, H. and G. Piras (2017). Spatial Econometrics. New York: Elsevier. 





Kelejian, H. H. and I. R. Prucha (2007). Hac estimation in a spatial framework. Journal of Econometrics 140 (1), 131–154. 





Leung, M. P. (2022). Rate-optimal cluster-randomized designs for spatial interference. The Annals of Statistics 50 (5), 3064–3087. 





Leung, M. P. (2023). Design of cluster-randomized trials with cross-cluster interference. arXiv preprint arXiv:2310.18836 . 





Li, S. and S. Wager (2022). Random graph asymptotics for treatment effect estimation under network interference. The Annals of Statistics 50 (4), 2334–2358. 





Newey, W. K. and D. McFadden (1994). Large sample estimation and hypothesis testing. Handbook of econometrics 4, 2111–2245. 





Newey, W. K. and K. D. West (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. Econometrica 55 (3), 703–708. 





Ogburn, E. L., O. Sofrygin, I. Diaz, and M. J. van der Laan (2020). Causal inference for social network data. arXiv preprint arXiv:1705.08527 . 





Papadogeorgou, G., K. Imai, J. Lyall, and F. Li (2020). Causal inference with spatiotemporal data: estimating the effects of airstrikes on insurgent violence in iraq. arXiv preprint arXiv:2003.13555 . 





Reich, B. J., S. Yang, Y. Guan, A. B. Giffin, M. J. Miller, and A. Rappold (2021). A review of spatial causal inference methods for environmental and epidemiological applications. International Statistical Review 89 (3), 605–634. 





Ross, N. et al. (2011). Fundamentals of Stein’s method. Probability Surveys 8, 210–293. 





Samii, C., M. Lisiecki, P. Kulkarni, L. Paler, L. Chavis, B. Snilstveit, M. Vojtkova, and E. Gallagher (2014). Effects of payment for environmental services (pes) on deforestation and poverty in low and middle income countries: a systematic review. Campbell Systematic Reviews 10 (1), 1–95. 





S¨avje, F., P. M. Aronow, and M. G. Hudgens (2018). Average treatment effects in the presence of unknown interference. arXiv:1711.06399 [math.ST]. 





S¨avje, F., P. M. Aronow, and M. G. Hudgens (2021). Average treatment effects in the presence of unknown interference. The Annals of Statistics 49 (2), 673–701. 





VanderWeele, T. J. and E. J. T. Tchetgen (2011). Effect partitioning under interference in two-stage randomized vaccine trials. Statistics & probability letters 81 (7), 861–869. 





Wunder, S. (2008). How do we deal with leakage. Moving ahead with REDD: issues, options and implications 1, 65–75. 





Young, A. (2015). Improved, nearly exact, statistical inference with robust and clustered covariance matrices using effective degrees of freedom corrections. Unpublished Manuscript, London School of Economics. 



Zigler, C. M. and G. Papadogeorgou (2018). Bipartite causal inference with interference. arXiv:1807.08660 [stat.ME]. 

# Appendix to Design-Based Inference for Spatial Experiments under Unknown Interference

# Contents

A.1 Results on Horvitz-Thompson Estimator 47 

A.1.1 Proof for Proposition 1 48 

A.1.2 Variance Characterization 48 

A.1.3 Unbiasedness of the HT estimator in observaitional studies 51 

A.2 Results on Hajek Estimator 51 

A.2.1 Linearization and Asymptotic Variance Characterization 51 

A.3 Results on Asymptotic Distribution 53 

A.3.1 Proofs for Propositions 2 and 3 54 

A.4 Variance Estimation 55 

A.4.1 HAC Variance Estimator 55 

A.4.2 SAH Variance Estimator 59 

A.4.3 Proof of Proposition 4 61 

A.5 Efficiency Comparison between the Hajek and HT estimators . 62 

A.6 Theoretical Results on Observational Studies 66 

A.7 Effective Degree of Freedom Adjustment 71 

B.1 Simulation Designs 73 

B.2 More Point-Intervention Simulation Results in Section 7 75 

B.3 Results for a Polygon-Intervention Simulation 75 

B.4 Additional Results for Empirical Applications 80 

C.1 Weaker Assumptions on the Extent of Interference . 82 

D.1 Notation 105 

# A Analytical Results

This section contains technical results discussed in the paper. It includes the proof of Propositions 1, 2, 3, and 4. The proof of Proposition 1 is in Section A.1.1, proofs of 2 and 3 are in Section A.3.1, and results on HAC variance estimator are included in Section A.4. 

We first prove the unbiasedness of the Horvitz-Thompson estimator for the AME. We then characterize the variance of the Horvitz-Thompson estimator and the asymptotic variance of the Hajek estimator. Asymptotic normality follows from Lemma 1 and Lemma 2 in Ogburn et al. (2020). Finally, we show that the spatial HAC standard errors estimator estimates a quantity that is probably larger than the true asymptomatic variance under C5, enabling conservative Wald-inference. 

To reduce the complexity of notations, we make the following simplifications that will be used throughout the appendices: 

• Without noted otherwise, the expectation is always taken over random assignments $\mathbf { Z }$ 

• We write $\begin{array} { r } { \mathrm { E } \left[ \mu _ { i } \left( \mathbf { Y } \left( \mathbf { Z } \right) ; d \right) \big | Z _ { i } \ = \ 1 \right] \ = \ \mathrm { E } \left[ \mu _ { i } ( \mathbf { 1 } ; d ) \right] \ \mathrm { a n d } \ \mathrm { E } \left[ \mu _ { i } \left( \mathbf { Y } \left( \mathbf { Z } \right) ; d \right) \big | Z _ { i } \ = \ 0 \right] \ = \ 0 . } \end{array}$ $\operatorname { E } \left[ \mu _ { i } \left( \mathbf { Y } \left( \mathbf { Z } \right) ; d \right) \big | Z _ { i } \ = \ 1 \right] \ = \ \operatorname { E } \left[ \mu _ { i } ( \mathbf { 1 } ; d ) \right]$ $\operatorname { E } \left[ \mu _ { i } ( \mathbf { 0 } ; d ) \right]$ . 

• We write $\operatorname { E } [ \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) \big | Z _ { i } = a , Z _ { j } = b \big ] = \operatorname { E } [ \mu _ { i } ( \mathbf { a _ { i } } , \mathbf { b _ { j } } ; d ) ]$ for $a , b \in \{ 0 , 1 \}$ . 

• We write PNi=1 Pj∈B(i;d) as Pi;j∈B(i;d). $\textstyle \sum _ { i = 1 } ^ { N } \sum _ { j \in B ( i ; d ) }$ $\textstyle \sum _ { i ; j \in B ( i ; d ) }$ 

# A.1 Results on Horvitz-Thompson Estimator

The following lemma is useful. 

Lemma A.1. For any function $f : \{ 0 , 1 \} ^ { N } \to { \mathbb { R } }$ and assuming C1, we have: 

C1.1 E $\left[ Z _ { i } ^ { k } f ( \mathbf { Z } ) \right] = p \mathrm { E } \left[ f ( 1 , \mathbf { Z } _ { - i } ) \right]$ for any positive integer $k$ . 

C1.2 E $\left[ Z _ { i } ^ { k } Z _ { j } ^ { l } f ( { \bf Z } ) \right] = p ^ { 2 } \mathrm { E } \left[ f ( 1 , 1 , { \bf Z } _ { - ( i , j ) } ) \right]$ for any postive integers k and l. 

Proof. By Law of Iterated Expectations. 

# A.1.1 Proof for Proposition 1

Proof. The Horvitz-Thompson estimator is defined in (12) and the AME is defined in (8). We have: 

$$
\begin{array}{l} \mathrm {A M E} (d; \eta) = \frac {1}{N} \sum_ {i = 1} ^ {N} \mu_ {i} (1; d, \eta) - \frac {1}{N} \sum_ {i = 1} ^ {N} \mu_ {i} (0; d, \eta) \\ = \frac {1}{N p} \sum_ {i = 1} ^ {N} p \mathrm {E} [ \mu_ {i} (\mathbf {Y}; d) | Z _ {i} = 1 ] - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {N} (1 - p) \mathrm {E} [ \mu_ {i} (\mathbf {Y}; d) | Z _ {i} = 0 ] \\ = \frac {1}{N p} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ Z _ {i} \mu_ {i} (\mathbf {Y}; d) \right] - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ (1 - Z _ {i}) \mu_ {i} (\mathbf {Y}; d) \right] \\ = \operatorname {E} \left[ \frac {1}{N p} \sum_ {i = 1} ^ {N} Z _ {i} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) \mu_ {i} (\mathbf {Y}; d) \right] = \operatorname {E} \left[ \widehat {\tau} _ {\mathrm {H T}} (d) \right], \\ \end{array}
$$

where the second equality uses the definition, and the third equality follows from Lemma A.1. 

# A.1.2 Variance Characterization

We now characterize the variance of the Horvitz-Thompson estimator. 

Lemma A.2. Under conditions C1-C3, the variance of estimator ${ \widehat { \tau } } _ { \mathrm { H T } } ( d )$ is bounded as follows: 

$$
\begin{array}{l} \mathrm {V a r} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right) \leq \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \mu_ {i} (\mathbf {1}; d) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \mu_ {i} (\mathbf {0}; d) ^ {2} \right] \\ + \frac {1}{N ^ {2}} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \Bigl \{\mathrm {E} [ \mu_ {i} (\mathbf {a _ {i}}, \mathbf {b _ {j}}; d) \mu_ {j} (\mathbf {a _ {i}}, \mathbf {b _ {j}}; d) ] - \mathrm {E} [ \mu_ {i} (\mathbf {a}; d) ] \mathrm {E} [ \mu_ {j} (\mathbf {b}; d) ] \Bigr \}, \\ \end{array}
$$

and, in addition under C4, we have that 

$$
\operatorname {V a r} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right) = O \left(\frac {1}{N}\right).
$$

Proof. Using the expression of the Horvitz-Thompson estimator, we have: 

$$
\begin{array}{l} \mathrm {V a r} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right) = \frac {1}{N ^ {2}} \mathrm {V a r} \left[ \sum_ {i = 1} ^ {N} \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d) \right] \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \mathrm {V a r} \left[ \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d) \right] \\ + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \mathrm {C o v} \left[ \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d), \left(\frac {Z _ {j}}{p} - \frac {1 - Z _ {j}}{1 - p}\right) \mu_ {j} (\mathbf {Y}; d) \right] \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \left(\left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d)\right) ^ {2} \right] - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\operatorname {E} \left[ \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d) \right]\right) ^ {2} \\ + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \operatorname {C o v} \left[ \frac {Z _ {i}}{p} \mu_ {i} (\mathbf {Y}; d), \frac {Z _ {j}}{p} \mu_ {j} (\mathbf {Y}; d) \right] \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \mathrm {C o v} \left[ \frac {Z _ {i}}{p} \mu_ {i} (\mathbf {Y}; d), \frac {1 - Z _ {j}}{1 - p} \mu_ {j} (\mathbf {Y}; d) \right] \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \mathrm {C o v} \left[ \frac {1 - Z _ {i}}{1 - p} \mu_ {i} (\mathbf {Y}; d), \frac {Z _ {j}}{p} \mu_ {j} (\mathbf {Y}; d) \right] \\ + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \operatorname {C o v} \left[ \frac {1 - Z _ {i}}{1 - p} \mu_ {i} (\mathbf {Y}; d), \frac {1 - Z _ {j}}{1 - p} \mu_ {j} (\mathbf {Y}; d) \right]. \\ \end{array}
$$

We further expand the first two terms in the above expression: 

$$
\begin{array}{l} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \left(\left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d)\right) ^ {2} \right] - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d) \right] ^ {2} \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \frac {Z _ {i} ^ {2}}{p ^ {2}} \mu_ {i} ^ {2} (\mathbf {Y}; d) \right] + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \frac {(1 - Z _ {i}) ^ {2}}{(1 - p) ^ {2}} \mu_ {i} ^ {2} (\mathbf {Y}); d) \right] \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \mathrm {E} ^ {2} \left[ \left(\frac {Z _ {i}}{p} - \frac {1 - Z _ {i}}{1 - p}\right) \mu_ {i} (\mathbf {Y}; d) \right] \\ = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) ^ {2} \right] - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) \right]\right) ^ {2} \\ \leq \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) ^ {2} \right]. \\ \end{array}
$$

We also have: 

$$
\frac {1}{N p} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) ^ {2} \right] - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) \right]\right) ^ {2} = O (1),
$$

since C2 implies that all the moments are bounded. Next, we examine the first covariance term, which equals 

$$
\begin{array}{l} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \mathrm {C o v} \left[ \frac {Z _ {i}}{p} \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d), \frac {Z _ {j}}{p} \mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) \right] = \frac {1}{N ^ {2}} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \mathrm {C o v} \left[ \frac {Z _ {i}}{p} \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d), \frac {Z _ {j}}{p} \mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) \right] \\ = \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \mathrm {E} [ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) \mu_ {j} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) ] - \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \mathrm {E} [ \mu_ {i} (\mathbf {1}; d) ] \mathrm {E} [ \mu_ {j} (\mathbf {1}; d) ]. \\ \end{array}
$$

The first equality holds because of assumption C3 on local interference. Moreover, 

$$
\frac {1}{N ^ {2}} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \mathrm {E} [ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) \mu_ {j} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) ] - \frac {1}{N ^ {2}} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \mathrm {E} [ \mu_ {i} (\mathbf {1}; d) ] \mathrm {E} [ \mu_ {j} (\mathbf {1}; d) ] = O (1),
$$

since $| B ( i ; d ) |$ is bounded by (C4a). Other covariance terms have similar forms. We obtain 

the bound of the variance and its convergence rate by combining these terms together. 

# A.1.3 Unbiasedness of the HT estimator in observaitional studies

Proof. Note that in our setup the only randomness comes from the random assignment. Hence functions of the covariates are considered fixed. In the scenario of observational studies and known propensity scores, we have 

$$
\begin{array}{l} \mathrm {A M E} (d; \eta) = \frac {1}{N} \sum_ {i = 1} ^ {N} \mu_ {i} (1; d, \eta) - \frac {1}{N} \sum_ {i = 1} ^ {N} \mu_ {i} (0; d, \eta) \\ = \frac {1}{N} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \frac {p (O _ {i})}{p (O _ {i})} \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) | Z _ {i} = 1 \right] - \frac {1}{N} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \frac {1 - p (O _ {i})}{1 - p (O _ {i})} \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) | Z _ {i} = 0 \right] \\ = \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {p (O _ {i})}{p (O _ {i})} \mathrm {E} [ Z _ {i} \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) | Z _ {i} = 1 ] - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {1 - p (O _ {i})}{1 - p (O _ {i})} \mathrm {E} [ (1 - Z _ {i}) \mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) | Z _ {i} = 0 ] \\ = \operatorname {E} \left[ \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i}}{p \left(O _ {i}\right)} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {1 - Z _ {i}}{p \left(O _ {i}\right)} \mu_ {i} (\mathbf {Y}; d) \right] \\ \end{array}
$$

The fourth equality uses C6-(ii) and the law of total expectations. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/c5b57ec4cf6d73d4c8e8e2bccf54e92e69efffc127fa29c146dbdafbd74e394d.jpg)


# A.2 Results on Hajek Estimator

# A.2.1 Linearization and Asymptotic Variance Characterization

We derive the limiting variance of the Hajek estimator using linearization technique. 

Lemma A.3. Consider the estimator $\widehat { \tau } _ { \mathrm { H A } } ( d )$ defined in (15). It has the following asymptotic linear expansion: 

$$
\widehat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d) = \mathrm {A M E} (d; \eta) + \frac {1}{N p} \sum_ {i = 1} ^ {N} Z _ {i} (\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {1} (d)) - \frac {1}{N p} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) (\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {0} (d)).
$$

Such an expansion satisfies $\sqrt { N } ( \widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) - \widehat { \tau } _ { \mathrm { H A } } ( d ) ) = o _ { p } ( 1 ) .$ 

Proof. Denote $\begin{array} { r } { \widehat { \mu } ^ { 1 } ( d ) \ = \ \frac { 1 } { N p } \sum _ { i = 1 } ^ { N } Z _ { i } \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) } \end{array}$ $\begin{array} { r } { \widehat { \mu } ^ { 1 } ( d ) \ = \ \frac { 1 } { N p } \sum _ { i = 1 } ^ { N } Z _ { i } \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) , \ \widehat { \mu } ^ { 0 } ( d ) \ = \ \frac { 1 } { N ( 1 - p ) } \sum _ { i = 1 } ^ { N } ( 1 - Z _ { i } ) \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) . } \end{array}$ , $\begin{array} { r } { \widehat { N } _ { 1 } = \frac { \sum _ { i = 1 } ^ { N } Z _ { i } } { N p } } \end{array}$ , $\begin{array} { r } { \widehat { N } _ { 0 } = \frac { \sum _ { i = 1 } ^ { N } \left( 1 - Z _ { i } \right) } { N \left( 1 - p \right) } } \end{array}$ , and ${ \bf W } = ( \widehat { \mu } ^ { 1 } ( d ) , \widehat { \mu } ^ { 0 } ( d ) , \widehat { N } _ { 1 } , \widehat { N } _ { 0 } )$ . Further define $\bar { \mu } ^ { 1 } ( d ) =$ $\begin{array} { r } { \frac { 1 } { N } \sum _ { i = 1 } ^ { N } \mathrm { ~ E ~ } [ \mu _ { i } ( \mathbf { 1 } ; d ) ] } \end{array}$ and $\begin{array} { r } { \bar { \mu } ^ { 0 } ( d ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } \mathrm { ~ E ~ } [ \mu _ { i } ( \mathbf { 0 } ; d ) ] } \end{array}$ . 

We know that E $[ \widehat { \mu } ^ { 1 } ( d ) ] = \bar { \mu } ^ { 1 } ( d )$ , E $[ \widehat { \mu } ^ { 0 } ( d ) ] = \bar { \mu } ^ { 0 } ( d )$ , E [N1] = E [N0] = 1. Thus, $\operatorname { E } \left[ \mathbf { W } \right] =$ $( \bar { \mu } ^ { 1 } ( d ) , \bar { \mu } ^ { 0 } ( d ) ) , 1 , 1 )$ . Define $\begin{array} { r } { f ( w ) = f ( a , b , c , d ) = \frac { a } { c } - \frac { b } { d } } \end{array}$ . Then the Hajek estimator can be written as $f ( { \bf W } ) = f ( \widehat { \mu } ^ { 1 } ( d ) , \widehat { \mu } ^ { 0 } ( d ) , \widehat { N } _ { 1 } , \widehat { N } _ { 0 } )$ . 

With probability approaching 1, we have the following Taylor expansion of the Hajek estimator:25 

$$
\widehat {\tau} _ {\mathrm {H A}} (d) = f (\mathbf {W}) = f (\operatorname {E} [ \mathbf {W} ]) + \left(\nabla f (\operatorname {E} [ \mathbf {W} ])\right) ^ {T} (\mathbf {W} - \operatorname {E} [ \mathbf {W} ]) + O _ {P} (| | \mathbf {W} - \operatorname {E} [ \mathbf {W} ] | | _ {2} ^ {2}).
$$

Following the same argument as in Lemma A.2, we have that $N | | \mathbf { W } - \mathrm { E } \left[ \mathbf { W } \right] | | _ { 2 } ^ { 2 } = O _ { p } ( 1 )$ and $O _ { P } ( \sqrt { N } | | \mathbf { W } - \mathrm { E } \left[ \mathbf { W } \right] | | _ { 2 } ^ { 2 } ) = o _ { P } ( 1 )$ . It is easy to see that $\nabla f ( \mathrm { E } \left[ { \mathbf W } \right] ) = ( 1 , - 1 , - \bar { \mu } ^ { 1 } ( d ) , \bar { \mu } ^ { 0 } ( d ) ) ^ { \prime }$ . Some algebraic manipulations prove that the first two terms simplify to the expressions in $\widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d )$ . □ 

Lemma A.4. The variance of the linearized Hajek estimator can be expressed as 

$$
\begin{array}{l} \operatorname {V a r} (d) \left(\widehat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d)\right) (32) \\ = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \left(\mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right] (33) \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \mathrm {E} ^ {2} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) - \left(\bar {\mu} ^ {1} (d) - \bar {\mu} ^ {0} (d)\right) \right] (34) \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \left(\mu_ {i} \left(\mathbf {a _ {i}}, \mathbf {b _ {j}}; d\right) - \bar {\mu} ^ {a} (d)\right) \left(\mu_ {j} \left(\mathbf {a _ {i}}, \mathbf {b _ {j}}; d\right) - \bar {\mu} ^ {b} (d)\right) \right] (35) \\ - \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} [ \mu_ {i} (\mathbf {a}; d) - \bar {\mu} ^ {a} (d) ] E [ \mu_ {j} (\mathbf {b}; d) - \bar {\mu} ^ {b} (d) ], (36) \\ \end{array}
$$

Under conditions C1-C5, we have the following variance bound for Var $\left( \widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) \right)$ 

$$
\begin{array}{l} \tilde {\mathrm {V}} _ {\mathrm {H A}} = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {0}; d); d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right] (37) \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a, b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \right.\left(\mu_ {i} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \bar {\mu} ^ {a} (d)\right)\left(\mu_ {j} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right); d\right) - \bar {\mu} ^ {b} (d)\left. \right)\left. \right]. (38) \\ \end{array}
$$

Proof. The characterization of the asymptotic variance is similar to that Lemma A.2. We omit the details. Note that for the terms in lines (34) and (36), we have, up to a minus sign, that 

$$
\begin{array}{l} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \mathrm {E} ^ {2} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) - (\bar {\mu} ^ {1} (d) - \bar {\mu} ^ {0} (d)) \right] \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a, b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \mu_ {i} (\mathbf {a}; d) - \bar {\mu} ^ {a} (d) \right] \mathrm {E} \left[ \mu_ {j} (\mathbf {b}; d) - \bar {\mu} ^ {b} (d) \right] \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\tau_ {i} (d; \eta) - \operatorname {A M E} (d; \eta)\right) ^ {2} \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \left\{\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d) \right] \times \operatorname {E} \left[ \mu_ {j} (\mathbf {1}; d) - \mu_ {j} (\mathbf {0}; d) - \operatorname {A M E} (d; \eta)) \right] \right\} \\ - \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \left\{\operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d) \right] \times \operatorname {E} \left[ \mu_ {j} (\mathbf {1}; d) - \mu_ {j} (\mathbf {0}; d) - \operatorname {A M E} (d; \eta) \right] \right\} \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\tau_ {i} (d; \eta) - \operatorname {A M E} (d; \eta)\right) \sum_ {j \in \mathcal {B} (i; d)} \left(\tau_ {j} (d; \eta) - \operatorname {A M E} (d; \eta)\right) \geq 0. \\ \end{array}
$$

by C5. Hence the term combining expressions in (33) and (36) is non-positive in the limit under C5, which proves the lemma. □ 

# A.3 Results on Asymptotic Distribution

The consistency of the proposed estimator follows from the fact that both $\mathrm { V a r } \left( \widehat { \tau } _ { \mathrm { H T } } ( d ) \right)$ and $\mathrm { V a r } \left( \widehat { \tau } _ { \mathrm { H A } } ( d ) \right)$ converge to zero as $N \to \infty$ under conditions C1-C4. The asymptotic normality 

of the Horvitz-Thompson estimator can be derived using classic central limit theorems for finitely dependent random variables based on the Stein’s method (Chen et al., 2010; Ross et al., 2011; Ogburn et al., 2020). The Hajek estimator’s asymptotic distribution can be then obtained via the linear expansion in Proposition A.4 and a similar application of the CLT result. For this purpose, we adpat the results in Ogburn et al. (2020) using the terms defined in our paper. 

Lemma A.5. (Ogburn et al. (2020), Lemma 1 and 2) Consider a set of $N$ units. Let $U _ { 1 } , \dots , U _ { N }$ be bounded mean-zero random variables with finite fourth moments and dependency neighborhoods $B ( i ; d )$ . If $c _ { i } ( d ) \leq \tilde { c }$ for all $i$ and $\tilde { c } ^ { 2 } / N \to 0$ , then 

$$
\frac {\sum_ {i = 1} ^ {N} U _ {i}}{\sqrt {\operatorname {V a r} \left(\sum_ {i = 1} ^ {N} U _ {i}\right)}} \to N (0, 1).
$$

# A.3.1 Proofs for Propositions 2 and 3

Proof. We first prove the case for the HT estimator. Define $U _ { i }$ as 

$$
U _ {i} = \frac {1}{\sqrt {\mathrm {V a r} (\widehat {\tau} _ {\mathrm {H T}} (d)}} \left(\frac {Z _ {i} \mu_ {i} ({\bf 1} ; d)}{N p} - \frac {(1 - Z _ {i}) \mu_ {i} (\mu_ {i} ({\bf 0} ; d)}{N (1 - p)} - \frac {\mathrm {E} [ \mu_ {i} ({\bf 1} ; d) ] - \mathrm {E} [ \mu_ {i} ({\bf 0} ; d) ]}{N}\right).
$$

We have $\begin{array} { r l } { \quad } & { { } \frac { \widehat \tau _ { \mathrm { H T } } ( d ) - \mathrm { A M E } ( d ; \eta ) } { \sqrt { \mathrm { V a r } \left( \widehat \tau _ { \mathrm { H T } } ( d ) \right) } } = \sum _ { i = 1 } ^ { N } U _ { i } } \end{array}$ . Obviously $\operatorname { E } \left[ U _ { i } \right] = 0$ and $\begin{array} { r } { \mathrm { V a r } \left( \sum _ { i = 1 } ^ { N } U _ { i } \right) = 1 } \end{array}$ . Under the premises of Proposition 2 and by C1 and C2, we have that the fourth moment of $U _ { i }$ is bounded for all $i$ . By condition C 4a, $c _ { i } ( d ) ~ \leq ~ c _ { N } ( d )$ in our case and $c _ { N } ^ { 2 } ( d ) / N \to 0$ . From Lemma A.5, we know that $\begin{array} { r } { \frac { \widehat \tau _ { \mathrm { H T } } ( d ) - \mathrm { A M E } ( d ; \eta ) } { \sqrt { \mathrm { V a r } ( \widehat \tau _ { \mathrm { H T } } ( d ) ) } } \to N ( 0 , 1 ) } \end{array}$ . Similarly by Lemma A.3, we have $\sqrt { N } ( \widehat { \tau } _ { \mathrm { H A } } ( d ) - \mathrm { A M E } ( d ; \eta ) ) = \sqrt { N } ( \widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) - \mathrm { A M E } ( d ; \eta ) ) + o _ { p } ( 1 )$ . Under the premises in Proposition 3, a similar argument proves the normality for the Hajek estimator. □ 

# A.4 Variance Estimation

# A.4.1 HAC Variance Estimator

In Section 5, we showed that the Hajek estimator can be interpreted as an OLS estimator that regresses circle averages on a constant term and the treatment indicator. We further suggested the use of a spatial HAC estimator for quantifying uncertanties. This section studies the behavior of the spatial HAC estimator. Importantly, we show that the HAC variance estimator is consistent for the variance bound defined in Lemma A.3 and lines (37)-(38). This result suggests that the standard Wald type inference is valid for AMEs. 

With a uniform kernel and $\tilde { d } = h ( d )$ , expanding the expression in (18) we have: 

$$
\begin{array}{l} \widehat {\Sigma} _ {\mathrm {H A C}} \left( \begin{array}{c} \hat {\mu} _ {0} (d) \\ \widehat {\tau} _ {\mathrm {H A}} (d) \end{array} \right) = \left( \begin{array}{c} N, N _ {1} \\ N _ {1}, N _ {1} \end{array} \right) ^ {- 1} \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \mathbf {X} _ {i} \mathbf {X} _ {j} ^ {\prime} \hat {e} _ {i} (d) \hat {e} _ {j} (d) \mathbf {1} \{j \in \mathcal {B} (i; d) \}\right) \left( \begin{array}{c} N, N _ {1} \\ N _ {1}, N _ {1} \end{array} \right) ^ {- 1} \\ = \frac {1}{N _ {1} ^ {2} N _ {0} ^ {2}} \left( \begin{array}{c} N _ {1}, - N _ {1} \\ - N _ {1}, N \end{array} \right) \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \left( \begin{array}{c} 1, Z _ {j} \\ Z _ {i}, Z _ {i} Z _ {j} \end{array} \right) \hat {e} _ {i} (d) \hat {e} _ {j} (d) {\bf 1} \{j \in {\cal B} (i; d) \}\right) \left( \begin{array}{c} N _ {1}, - N _ {1} \\ - N _ {1}, N \end{array} \right), \\ \end{array}
$$

$( 2 , 2 )$ entry of $\begin{array} { r } { \left( N _ { 1 } , - N _ { 1 } \right) \left( \begin{array} { c } { 1 , Z _ { j } } \\ { Z _ { i } , Z _ { i } Z _ { j } } \end{array} \right) \left( \begin{array} { c } { N _ { 1 } , - N _ { 1 } } \\ { - N _ { 1 } , N } \end{array} \right) \mathrm { ~ e q u a l s ~ t o ~ } N _ { 1 } ^ { 2 } - N _ { 1 } N Z _ { i } - } \\ { - N _ { 1 } , N } \end{array}$ $N _ { 1 } ^ { 2 } - N _ { 1 } N Z _ { i } -$ $N _ { 1 } N Z _ { j } + N ^ { 2 } Z _ { i } Z _ { j } .$ . 

Reindex the sample such that treated observations lie before observations under control 

and plug in the expression of $\hat { e } _ { i } ( d )$ , we can see that: 

$$
\begin{array}{l} \widehat {\mathrm {V}} _ {\mathrm {H A C}} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right) (39) \\ = \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \hat {e} _ {i} ^ {2} (d) + \frac {1}{N _ {0} ^ {2}} \sum_ {i = N _ {1} + 1} ^ {N} \hat {e} _ {i} ^ {2} (d) + \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 1} \hat {e} _ {i} (d) \hat {e} _ {j} (d) (40) \\ - \frac {1}{N _ {1} N _ {0}} \sum_ {i = 1} ^ {N _ {1}} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 0} \hat {e} _ {i} (d) \hat {e} _ {j} (d) - \frac {1}{N _ {1} N _ {0}} \sum_ {i = N _ {1} + 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 1} \hat {e} _ {i} (d) \hat {e} _ {j} (d) (41) \\ + \frac {1}{N _ {0} ^ {2}} \sum_ {i = N _ {1} + 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 0} \hat {e} _ {i} (d) \hat {e} _ {j} (d) (42) \\ = \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) ^ {2} + \frac {1}{N _ {0} ^ {2}} \sum_ {i = N _ {1} + 1} ^ {N} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {0} (d)\right) ^ {2} (43) \\ + \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 1} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) \left(\mu_ {j} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) (44) \\ - \frac {1}{N _ {1} N _ {0}} \sum_ {i = 1} ^ {N _ {1}} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 0} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) \left(\mu_ {j} (d) - \widehat {\bar {\mu}} ^ {0} (d)\right) (45) \\ - \frac {1}{N _ {1} N _ {0}} \sum_ {i = N _ {1} + 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 1} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {0} (d)\right) \left(\mu_ {j} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) (46) \\ + \frac {1}{N _ {0} ^ {2}} \sum_ {i = N _ {1} + 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 0} \left(\mu_ {i} (d) - \widehat {\bar {\mu}} ^ {0} (d)\right) \left(\mu_ {j} (d) - \widehat {\bar {\mu}} ^ {0} (d)\right), (47) \\ \end{array}
$$

where $\begin{array} { r } { \widehat { \bar { \mu } } ^ { 1 } ( d ) = \frac { \sum _ { i = 1 } ^ { N } Z _ { i } \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) } { \sum _ { i = 1 } Z _ { i } } } \end{array}$ d) and µ¯0(d) = PNi=1 (1−Zi)µi(Y(Z);d) . $\begin{array} { r } { \widehat { \bar { \mu } } ^ { 0 } ( d ) = \frac { \sum _ { i = 1 } ^ { N } ( 1 - Z _ { i } ) \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) } { \sum _ { i = 1 } ( 1 - Z _ { i } ) } } \end{array}$ 

We now show that the variance estimate $\widehat { \mathrm { V } } _ { \mathrm { H A C } } ( \widehat { \tau } _ { \mathrm { H A } } ( d ) )$ is consistent for the rescaled variance bound defined in Lemma A.4, $\tilde { \mathrm { V } } _ { \mathrm { H A } } ( d )$ . Note that $\tilde { \mathrm { V _ { H A } } } ( d )$ is provably larger than the asymptotic variance of the Hajek estimator. The result below thus suggests that the normal confidence interval with the HAC variance provides conservative coverage for the Hajek estimator asymptotically. 

Proposition 5. If $N \times \tilde { \mathrm { V } } _ { \mathrm { H A } } ( d )$ is uniformly bounded below for large $N$ , then we have 

$$
\frac {\widehat {\mathrm {V}} _ {\mathrm {H A C}} (\widehat {\tau} _ {\mathrm {H A}} (d)) - \tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)} \xrightarrow {p} 0
$$

Proof. We show below that $N \times ( \hat { \mathrm { V } } _ { \mathrm { H A C } } ( \hat { \tau } _ { \mathrm { H A } } ( d ) ) - \tilde { \mathrm { V } } _ { \mathrm { H A } } ( d ) ) \stackrel { p } {  } 0$ . This, together with the premise in the lemma, leads to the claim that $\frac { \widehat { \mathrm { V } } _ { \mathrm { H A C } } ( \widehat { \tau } _ { \mathrm { H A } } ( d ) ) - \widetilde { \mathrm { V } } _ { \mathrm { H A } } } { \widetilde { \mathrm { V } } _ { \mathrm { H A } } } \ \overset { p } { \longrightarrow } \ 0$ . We first study terms in (43). For the treated group, we have 

$$
\begin{array}{l} N \left(\frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \widehat {\bar {\mu}} ^ {1} (d)\right) ^ {2} - \frac {1}{p} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} E [ (\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)) ^ {2} ]\right) \\ = N \left(\frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \mu_ {i} ^ {2} (\mathbf {Y} (\mathbf {Z}); d) - \frac {1}{N _ {1}} \left(\widehat {\bar {\mu}} ^ {1} (d)\right) ^ {2} - \frac {1}{p} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} E [ \mu_ {i} (\mathbf {1}; d) ^ {2} ] + \frac {1}{p N} \left(\bar {\mu} ^ {1} (d)\right) ^ {2}\right) \\ = \left(p \frac {N ^ {2}}{N _ {1} ^ {2}} \frac {1}{N p} \sum_ {i = 1} ^ {N _ {1}} \mu_ {i} ^ {2} (\mathbf {Y} (\mathbf {Z}); d) - \frac {1}{p} \frac {1}{N} \sum_ {i = 1} ^ {N} E [ \mu_ {i} (\mathbf {1}; d) ^ {2} ]\right) - \left(\frac {N}{N _ {1}} \left(\widehat {\bar {\mu}} ^ {1} (d)\right) ^ {2} - \frac {1}{p} \left(\bar {\mu} ^ {1} (d)\right) ^ {2}\right) \\ \stackrel {p} {\rightarrow} 0 \\ \end{array}
$$

The convergence in probability is justified by noting $\begin{array} { r } { \frac { 1 } { N p } \sum _ { i = 1 } ^ { N _ { 1 } } \mu _ { i } ^ { 2 } ( { \bf Y } ( { \bf Z } ) ; d ) } \end{array}$ $\begin{array} { r } { \frac { 1 } { N p } \sum _ { i = 1 } ^ { N _ { 1 } } \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) ; d ) } \end{array}$ and $\textstyle { \frac { N _ { 1 } } { N } }$ are all Horvitz-Thompson estimators, and under C1-C4 they converge to their mean in probability. For the control group, we can similarly show 

$$
N \left(\frac {1}{N _ {0} ^ {2}} \sum_ {i = 1} ^ {N _ {0}} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \widehat {\bar {\mu}} ^ {0} (d)\right) ^ {2} - \frac {1}{1 - p} \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} E [ (\mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d)) ^ {2} ]\right) \xrightarrow {p} 0
$$

Now we consider terms (44)-(47). All terms are similar in stuctural so for simplicity we only include the calculation for (44). We first have the algebraic identity: 

$$
\begin{array}{l} \frac {N}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N _ {1}} \sum_ {j \in \mathcal {B} (i; d), j \neq i, Z _ {j} = 1} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \widehat {\bar {\mu}} ^ {1} (d)\right) \left(\mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) - \widehat {\bar {\mu}} ^ {1} (d)\right) (48) \\ = \frac {N}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i,} Z _ {i} Z _ {j} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {1} (d)\right) \left(\mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {1} (d)\right) (49) \\ + \frac {N}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} Z _ {i} Z _ {j} \left(\bar {\mu} ^ {1} (d) - \hat {\bar {\mu}} ^ {1} (d)\right) \left(\mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) - \hat {\bar {\mu}} ^ {1} (d)\right) (50) \\ + \frac {N}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} Z _ {i} Z _ {j} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}; d) - \bar {\mu} ^ {1} (d)) \left(\bar {\mu} ^ {1} (d) - \hat {\bar {\mu}} ^ {1} (d)\right) \right. (51) \\ = \frac {N ^ {2}}{N _ {1} ^ {2}} \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} Z _ {i} Z _ {j} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {1} (d)\right) \left(\mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) - \bar {\mu} ^ {1} (d)\right) (52) \\ + \frac {N ^ {2}}{N _ {1} ^ {2}} \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} Z _ {i} Z _ {j} \left(\bar {\mu} ^ {1} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) \left(\mu_ {j} (\mathbf {Y} (\mathbf {Z}); d) - \widehat {\bar {\mu}} ^ {1} (d)\right) (53) \\ + \frac {N ^ {2}}{N _ {1} ^ {2}} \frac {1}{N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} Z _ {i} Z _ {j} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}; d) - \bar {\mu} ^ {1} (d)) \left(\bar {\mu} ^ {1} (d) - \widehat {\bar {\mu}} ^ {1} (d)\right) \right. (54) \\ \end{array}
$$

# Notice

1. Terms in (53) and (54) are of order $o _ { p } ( 1 )$ . For (53) , for example, by C1-C4, 

$$
\begin{array}{l} | \frac {N ^ {2}}{N _ {1} ^ {2}} (\bar {\mu} ^ {1} (d) - \widehat {\bar {\mu}} ^ {1} (d)) \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d)} Z _ {i} Z _ {j} (\mu_ {i} (\mathbf {Y} (\mathbf {Z}) - \widehat {\bar {\mu}} ^ {1} (d)) | \\ \leq \frac {N ^ {2}}{N _ {1} ^ {2}} | \bar {\mu} ^ {1} (d) - \widehat {\bar {\mu}} ^ {1} (d)) | \times \frac {1}{N} \sum_ {i = 1} ^ {N} c _ {i} (d) | \mu_ {i} (\mathbf {Y} (\mathbf {Z}) - \widehat {\bar {\mu}} ^ {1} (d) | \\ = o _ {p} (1) \times O _ {p} (1) = o _ {p} (1) \\ \end{array}
$$

where $c _ { i } ( d )$ is defined in Section 4.26 The argument is the same for the third term. 

2. For first term, we have, under C1-C4, 

$$
(5 2) - \frac {1}{N} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} E \left[ \left(\mu_ {i} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) - \bar {\mu} ^ {1} (d))\right) \left(\mu_ {j} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) - \bar {\mu} ^ {1} (d))\right) \right] = o _ {p} (1).
$$

The proof is similar as in Proposition 6.2 in Aronow and Samii (2017). 

Collecting all terms, we have proved that $\begin{array} { r } { N \times \left( \hat { \mathrm { V } } _ { \mathrm { H A C } } ( \hat { \tau } _ { \mathrm { H A } } ( d ) ) - \tilde { \mathrm { V } } _ { \mathrm { H A } } ( d ) \right) \xrightarrow { p } 0 . } \end{array}$ . 

# A.4.2 SAH Variance Estimator

In the previous section, we discussed the inference procedure for the Hajek estimator under C5. We now provide an alternative approach, based on a proposal in S¨avje et al. (2021). We have the following lemma: 

Proposition 6. Under C1-C4, we have, $\cdot$ 

$$
\begin{array}{l} \mathrm {V a r} \left(\hat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d)\right) \leq \bar {\mathrm {V}} _ {\mathrm {H A}} (d) \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} c _ {i} (d) \frac {\mathrm {E} [ (\mu_ {i} ({\bf 1} ; d) - \bar {\mu} ^ {1} (d)) ^ {2} ]}{p} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} c _ {i} (d) \frac {\mathrm {E} [ (\mu_ {i} ({\bf 1} ; d) - \bar {\mu} ^ {0} (d)) ^ {2} ] ^ {2}}{1 - p}. \\ \end{array}
$$

Define an estimator $\widehat { \mathrm { V } } _ { \mathrm { S A H } } ( d )$ 

$$
\widehat {\mathrm {V}} _ {\mathrm {S A H}} (d) = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} c _ {i} (d) \frac {(\mu_ {i} (\mathbf {Z} , d) - \hat {\mu} ^ {1} (d)) ^ {2}}{p ^ {2}} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) c _ {i} (d) \frac {(\mu_ {i} (\mathbf {Z} , d) - \hat {\mu} ^ {0} (d)) ^ {2}}{(1 - p) ^ {2}}. (5 5)
$$

Provided that $N \times \mathrm { V _ { H A } } ( d )$ is uniformly bounded below for large $N$ , $\widehat { \mathrm { V } } _ { \mathrm { S A H } } ( d )$ is consistent for $\mathrm { V _ { H A } } ( d )$ : 

$$
\frac {\widehat {\mathrm {V}} _ {\mathrm {S A H}} (d) - \bar {\mathrm {V}} _ {\mathrm {H A}} (d)}{\bar {\mathrm {V}} _ {\mathrm {H A}} (d)} \xrightarrow {p} 0.
$$

Proof. We first prove the upper bound. We have the following identity: 

$$
\begin{array}{l} \text {V a r} \left(\frac {Z _ {i} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {1} (d)\right)}{p} - \frac {\left(1 - Z _ {i}\right) \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {0} (d)\right)}{1 - p}\right) \\ = \operatorname {V a r} \left(\frac {Z _ {i} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {1} (d)\right)}{p}\right) + \operatorname {V a r} \left(\frac {\left(1 - Z _ {i}\right) \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {0} (d)\right)}{1 - p}\right) \\ - 2 \operatorname {C o v} \left(\frac {Z _ {i} \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {1} (d)\right)}{p}, \frac {\left(1 - Z _ {i}\right) \left(\mu_ {i} (\mathbf {Y} (\mathbf {Z}) , d) - \bar {\mu} ^ {0} (d)\right)}{1 - p}\right) \\ = \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {1} ; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right]}{p} - \left(\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right]\right) ^ {2} \\ + \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {0} ; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right]}{1 - p} - \left(\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right]\right) ^ {2} \\ + 2 \mathrm {E} [ \mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d) ] \times \mathrm {E} [ \mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d) ] \\ = \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {1} ; d) - \mu^ {1} (d)\right) ^ {2} \right]}{p} + \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {0} ; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right]}{1 - p} \\ \left. - \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d) \right] - \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d) \right]\right) ^ {2} \right. \\ \leq \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {1} ; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right]}{p} + \frac {\operatorname {E} \left[ \left(\mu_ {i} (\mathbf {0} ; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right]}{1 - p} \\ \end{array}
$$

Let’s define $\begin{array} { r } { A _ { i } = \frac { Z _ { i } \left( \mu _ { i } \left( \mathbf { Y } \left( \mathbf { Z } \right) , d \right) - \bar { \mu } ^ { 1 } \left( d \right) \right) } { p } } \end{array}$ and $\begin{array} { r } { B _ { i } = \frac { ( 1 - Z _ { i } ) ( \mu _ { i } ( \mathbf { Y } ( \mathbf { Z } ) , d ) - \bar { \mu } ^ { \mathrm { { U } } } ( d ) ) } { 1 - p } } \end{array}$ , then 

$$
\begin{array}{l} \operatorname {V a r} \left[ \hat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d) \right] = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {V a r} \left[ A _ {i} - B _ {i} \right] + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \operatorname {C o v} \left[ A _ {i} - B _ {i}, A _ {j} - B _ {j} \right] \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {V a r} [ A _ {i} - B _ {i} ] + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} \operatorname {C o v} [ A _ {i} - B _ {i}, A _ {j} - B _ {j} ] \\ \leq \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {V a r} \left[ A _ {i} - B _ {i} \right] + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; d), j \neq i} \frac {\operatorname {V a r} \left[ A _ {i} - B _ {i} \right] + \operatorname {V a r} \left[ A _ {j} - B _ {j} \right]}{2} \\ = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \operatorname {V a r} [ A _ {i} - B _ {i} ] + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} (c _ {i} (d) - 1) \operatorname {V a r} [ A _ {i} - B _ {i} ] \\ \leq \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} c _ {i} (d) \left(\frac {\operatorname {E} \left[ (\mu_ {i} (\mathbf {1} ; d) - \bar {\mu} ^ {1} (d)) ^ {2} \right]}{p} + \frac {\operatorname {E} \left[ (\mu_ {i} (\mathbf {0} ; d) - \bar {\mu} ^ {0} (d)) ^ {2} \right]}{1 - p}\right). \\ \end{array}
$$

Note that we used the fact that, by definition, $j \in B ( i ; d )$ if and only if $i \in B ( j ; d )$ . This 

proves the upper bound for the variance. The consistency of the variance estimator is analogous to that in Proposition 5 and we omit for brevity. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/607864af6d1824a7fc3a0cc72ac440146055a87ad03a41fa97bb9dc7dec4ec58.jpg)


# A.4.3 Proof of Proposition 4

We use AVar ( $\widehat { \tau } _ { \mathrm { H A } } ( d ) )$ to denote $\mathrm { V a r } \left( \tau _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) \right)$ . By Proposition 3 and Proposition 5, we have: 

$$
\frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} = \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\mathrm {A V a r} (\widehat {\tau} _ {\mathrm {H A}} (d))}} \times \sqrt {\frac {\mathrm {A V a r} (\widehat {\tau} _ {\mathrm {H A}} (d))}{\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}} \times \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}}.
$$

Thus we have, for each $\alpha < \textstyle { \frac { 1 } { 2 } }$ 

$$
\begin{array}{l} \textbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) \\ = \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\operatorname {A V a r} (\widehat {\tau} _ {H A} (d))}} \times \sqrt {\frac {\operatorname {A V a r} (\widehat {\tau} _ {\mathrm {H A}} (d))}{\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}} \times \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) \\ = \mathbf {P r o b} \left(\sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)}} z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \operatorname {A M E} (d ; \eta)}{\sqrt {\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)}} \times \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}} \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)}}\right) \\ \geq \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\mathrm {A V a r} (\widehat {\tau} _ {\mathrm {H A}} (d))}} \times \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right), \\ \end{array}
$$

where the last line follows because $\begin{array} { r } { \sqrt { \frac { \tilde { \mathrm { V } } _ { \mathrm { H A } } ( d ) } { \mathrm { A V a r } ( \widehat { \tau } _ { H A } ( d ) ) } } \geq 1 } \end{array}$ . Hence we have: 

$$
\begin{array}{l} \lim _ {N \to \infty} \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d)}{\sqrt {\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) \\ \geq \lim _ {N \to \infty} \mathbf {P r o b} \left(z _ {\frac {\alpha}{2}} \leq \frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \mathrm {A M E} (d)}{\sqrt {\mathrm {A V a r} (\widehat {\tau} _ {H A} (d))}} \times \sqrt {\frac {\tilde {\mathrm {V}} _ {\mathrm {H A}} (d)}{\widehat {\mathrm {V}} _ {\mathrm {H A C}} (d)}} \leq z _ {1 - \frac {\alpha}{2}}\right) = 1 - \alpha , \\ \end{array}
$$

because $\begin{array} { r } { \frac { \widehat \tau _ { \mathrm { H A } } ( d ) - \mathrm { A M E } ( d ; \eta ) } { \sqrt { \mathrm { A V a r } \left( \widehat \tau _ { \mathrm { H A } } ( d ) \right) } } \times \sqrt { \frac { \widetilde \mathrm { V } _ { \mathrm { H A } } ( d ) } { \widehat \mathrm { V } _ { \mathrm { H A C } } ( d ) } } \ \overset { d } { \to } \ N ( 0 , 1 ) } \end{array}$ V˜ HA(d)V (d) d→ N (0, 1). A similar calculation follows for the case using the variance estimator $\widehat { \mathrm { V } } _ { \mathrm { S A H } } ( d )$ . 

# A.5 Efficiency Comparison between the Hajek and HT estimators

This section compares the estimation efficiency, in terms of asymptotic variances, between the Hajek and HT estimators. Although it is not true that the Hajek estimator is not always more efficient than the HT estimator in our setting, the Hajek estimator has some efficiency properties that make it appealing in practice: 

1. When there is no interactive effect, $^ { 2 8 }$ the Hajek estimator is optimal among the estimators that use treatment-arm-specific intercepts. 

2. When the interactive effect size is small, the Hajek estimator can be expected to be more efficient than the HT estimator. 

3. When an interactive effect quantity ( $\Delta ( 1 )$ below), average treated and control outcomes have the same sign, the Hajek estimator is more efficient than the HT estimator. 

For other cases, there are potential outcomes that make Hajek estimator more efficient than the HT estimator and vice versa. 

Consider the following class of estimators, which adjust the treated and control outcomes 

Note, by C1, this also implies that 

If we assume that the potential outcome has the form $\begin{array} { r } { \mu _ { i } ( \mathbf { Z } ) = \beta _ { i } \mathbf { Z } _ { i } + \sum _ { j \in \mathcal { B } ( i : d ) } \beta _ { i , j } \mathbf { Z } _ { j } + \sum _ { k , l \in \mathcal { B } ( i : d ) } \beta _ { i , k l } \mathbf { Z } _ { k } \mathbf { Z } _ { l } } \end{array}$ , setting all coefficients $\beta _ { i , k l } = 0$ will rule out the interactive effect. 

with a treatment-arm-specific intercepts. 

$$
\widehat {\tau} \left(\mu_ {1}, \mu_ {0}; d\right) = \mu_ {1} - \mu_ {0} + \left(\frac {1}{N p} \sum_ {i = 1} ^ {n} \mathbf {Z} _ {i} \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {1}\right) - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {n} \left(1 - \mathbf {Z} _ {i}\right) \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {0}\right)\right) \tag {58}
$$

Both the HT estimator and the linearized version of the Hajek estimator have this form: $\widehat { \mu } ( d ) = \widehat { \tau } ( 0 , 0 )$ and $\widehat { \mu } _ { \mathrm { H A } } ( d ) = \widehat { \tau } ( \bar { \mu } ^ { 1 } ( d ) , \bar { \mu } ^ { 0 } ( d ) )$ . It is in this class of the estimators we discuss the efficiency property of the Hajek estimator. We define the interactive effect quantity mentioned above 

$$
\Delta (1) = \frac {1}{N} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) \right] - \operatorname {E} \left[ \mu_ {i} (\mathbf {0 _ {i}}, \mathbf {1 _ {j}}; d) \right] - \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) \right] - \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) \right]\right)\right). ^ {2 9}
$$

$^ { 2 9 }$ One can also define a similar interactive effect quantity 

$$
\Delta (0) = \frac {1}{N} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {0 _ {i}}, \mathbf {0 _ {j}}; d) \right] - \operatorname {E} \left[ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {0 _ {j}}; d) \right] - \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) \right] - \operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) \right]\right)\right).
$$

These two quantities are related by the identity $p \Delta ( 1 ) = ( 1 - p ) \Delta ( 0 )$ . This follows by a calculation 

$$
\begin{array}{l} N \times (p \Delta (1) - (1 - p) \Delta (0)) = \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} (p \mathrm {E} [ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {1 _ {j}}; d) ] + (1 - p) \mathrm {E} [ \mu_ {i} (\mathbf {1 _ {i}}, \mathbf {0 _ {j}}; d) ]) \\ \left. - \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \left(p \operatorname {E} \left[ \mu_ {i} \left(\mathbf {0} _ {\mathbf {i}}, \mathbf {1} _ {\mathbf {j}}; d\right) \right] + (1 - p) \operatorname {E} \left[ \mu_ {i} \left(\mathbf {0} _ {\mathbf {i}}, \mathbf {0} _ {\mathbf {j}}; d\right) \right]\right) - \left(\sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) \right] - \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \operatorname {E} \left[ \mu_ {i} (\mathbf {0}; d) \right]\right)\right) \\ = 0. \\ \end{array}
$$

Proposition 7. Under C1-C4, we have, for any scalars $\mu _ { 1 }$ and $\mu _ { 0 }$ 

$$
\begin{array}{l} \mathrm {V} \left(\widehat {\tau} \left(\mu_ {1}, \mu_ {0}; d\right)\right) = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \mu_ {1}\right) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {0}; d) - \mu_ {0}\right) ^ {2} \right] (59) \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) - \left(\mu_ {1} - \mu_ {0}\right) \right]\right) ^ {2} (60) \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \left(\mu_ {i} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \mu_ {a}\right) \left(\mu_ {j} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \mu_ {b}\right) \right] (61) \\ - \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} [ \mu_ {i} (\mathbf {a}; d) - \mu_ {a} ] E [ \mu_ {j} (\mathbf {b}; d) - \mu_ {b} ]. (62) \\ \end{array}
$$

In particular, we have the following algebraic identity: 

$$
\mathrm {V} (\widehat {\tau} \left(\mu_ {1}, \mu_ {0}; d)\right) = \frac {(1 - p) p}{N} \left(\frac {\mu_ {1}}{p} + \frac {\mu_ {0}}{1 - p} - \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {p \Delta (1)}{p (1 - p)}\right)\right) ^ {2} + C, \tag {63}
$$

where $C$ is a constant independent of $\mu _ { 1 }$ and $\mu _ { 0 }$ 

Proof. The characterization of the variance is similar to that in Proposition A.4. Expanding 

the variance expression we have the following identity: 

$$
\begin{array}{l} \mathrm {V} \left(\widehat {\tau} \left(\mu_ {1}, \mu_ {0}; d\right)\right) = \left(\frac {1}{N p} - \frac {1}{N}\right) \mu_ {1} ^ {2} + \left(\frac {1}{N (1 - p)} - \frac {1}{N}\right) \mu_ {0} ^ {2} + \frac {2}{N} \mu_ {1} \mu_ {0} \\ + \left(\left(- \frac {2}{N p} + \frac {2}{N}\right) \bar {\mu} _ {1} (d) - \frac {2}{N} \bar {\mu} _ {0} (d) - \frac {2}{N} \Delta (1)\right) \mu_ {1} \\ + \left(\left(- \frac {2}{N (1 - p)} + \frac {2}{N}\right) \bar {\mu} _ {0} (d) - \frac {2}{N} \bar {\mu} _ {1} (d) - \frac {2}{N} \Delta (0)\right) \mu_ {0} \\ + \mathrm {V} (\widehat {\tau} (0, 0; d)) \\ = \frac {1}{N} \left(\frac {1 - p}{p} \mu_ {1} ^ {2} + \frac {p}{1 - p} \mu_ {0} ^ {2} + 2 \mu_ {1} \mu_ {0}\right) - \frac {2}{N} \left(\frac {1 - p}{p} \bar {\mu} _ {1} (d) + \bar {\mu} _ {0} (d) + \Delta (1)\right) \mu_ {1} \\ - \frac {2}{N} \left(\frac {p}{1 - p} \bar {\mu} _ {0} (d) + \bar {\mu} _ {1} (d) + \Delta (0)\right) \mu_ {0} + \mathrm {V} (\widehat {\tau} (0, 0; d)) \\ = \frac {(1 - p) p}{N} \left(\frac {\mu_ {1} ^ {2}}{p ^ {2}} + \frac {\mu_ {0} ^ {2}}{(1 - p) ^ {2}} + 2 \frac {\mu_ {1} \mu_ {0}}{p (1 - p)}\right) - \frac {2 (1 - p) p}{N} \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {p \Delta (1)}{p (1 - p)}\right) \frac {\mu_ {1}}{p} \\ - \frac {2 (1 - p) p}{N} \left(\frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {\bar {\mu} _ {1} (d)}{p} + \frac {(1 - p) \Delta (0)}{p (1 - p)}\right) \frac {\mu_ {0}}{1 - p} + \mathrm {V} (\widehat {\tau} (0, 0)) \\ = \frac {(1 - p) p}{N} \left(\left(\frac {\mu_ {1}}{p} + \frac {\mu_ {0}}{1 - p}\right) ^ {2} - 2 \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {p \Delta (1)}{p (1 - p)}\right) \left(\frac {\mu_ {1}}{p} + \frac {\mu_ {0}}{1 - p}\right)\right) + \mathrm {V} (\widehat {\tau} (0, 0; d)) \\ = \frac {(1 - p) p}{N} \left(\frac {\mu_ {1}}{p} + \frac {\mu_ {0}}{1 - p} - \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {p \Delta (1)}{p (1 - p)}\right)\right) ^ {2} + C \\ \end{array}
$$

where $C$ is independent of $\mu _ { 1 }$ and $\mu _ { 0 }$ . In the fourth equality, we used the fact that $p \Delta ( 1 ) =$ $( 1 - p ) \Delta ( 0 )$ . □ 

Expression (63) immediately leads to the conclusions at the beginning of this section. When the interactive effect is small and hence $\Delta ( 1 )$ is very close to zero, we have 

$$
\mathrm {V} \left(\widehat {\tau} _ {H A} ^ {\mathrm {L}} (d)\right) = \mathrm {V} \left(\widehat {\tau} \left(\bar {\mu} _ {1} (d), \bar {\mu} _ {1} (d)\right)\right) \approx 0 + C, \tag {64}
$$

which is close to minimizing the variance within the class of estimators. We can further 

compare the variance of the HT and HA estimator: 

$$
\begin{array}{l} \mathrm {V} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right) - \mathrm {V} \left(\widehat {\tau} _ {H A} ^ {\mathrm {L}} (d)\right) \propto \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {p \Delta (1)}{p (1 - p)}\right) ^ {2} - \left(\frac {p \Delta (1)}{p (1 - p)}\right) ^ {2} (65) \\ = \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p}\right) \left(\frac {\bar {\mu} _ {1} (d)}{p} + \frac {\bar {\mu} _ {0} (d)}{1 - p} + \frac {2 p \Delta (1)}{p (1 - p)}\right). (66) \\ \end{array}
$$

Thus when the average treated and control outcomes, and the interactive effect is of the same sign, the quantity is nonnegative and the Hajek estimator is weakly more efficient than the HT estimator. 

# A.6 Theoretical Results on Observational Studies

This section contains results for inference on the AME in observational setting. In particular, we consider the case where assignment probabilities are modeled by a logistic model and derive asymptotic linearization and variance estimator. 

Given the setup in Section 6.4. We make the further parametric assumptions on the assignment model. 

C 7. The following properties hold for all sample size $N$ and for all $i \in S _ { N }$ : 

(i) (Fixed Dimensionality)The confounder $O _ { i } = ( o _ { 1 i } , . . . , o _ { k i } ) ^ { \prime }$ is of dimension $k$ . 

(ii) (Bounded Confounders) Confouders are uniformly bounded: there exists a constant $U$ such that $\mathrm { m a x } _ { i \in { \mathcal { S } } _ { N } } | | O _ { i } | | _ { \infty } \leq U$ . 30 

(iii) (Logistic Model) The assignment probabilites follow a logistic model: for all $i \in S _ { N }$ there exists a $\theta _ { 0 } \in \mathbb { R } ^ { k }$ such that 

$$
P (Z _ {i} = 1) = p (O _ {i} | \theta) = \frac {\exp (O _ {i} ^ {\prime} \theta_ {0})}{1 + \exp (O _ {i} ^ {\prime} \theta_ {0})}
$$

(iv) Compact Parameter Space: $\theta _ { 0 } \in \Theta _ { 0 }$ where $\Theta _ { 0 } \subset \mathbb { R } ^ { k }$ is a compact set. 

(v) (Nonsingularity) The smallest eigenvalue of $\begin{array} { r } { \frac { 1 } { N } \sum _ { i = 1 } ^ { N } O _ { i } O _ { i } ^ { \prime } } \end{array}$ is uniformly bounded below. 

(vi) The MLE estimator is used to estimate the coefficient vector $\theta$ : 

$$
\widehat {\theta} _ {\mathrm {M L E}} = \arg \max _ {\theta \in \Theta_ {0}} \sum_ {i = 1} ^ {N} \big (Z _ {i} O _ {i} ^ {\prime} \theta - \log (1 + \exp (O _ {i} ^ {\prime} \theta)) \big).
$$

Lemma A.6. Under C6 and C7, $\widehat { \theta } _ { M L E } - \theta _ { 0 } = o _ { p } ( 1 )$ and 

$$
\widehat {\theta} _ {\mathrm {M L E}} - \theta_ {0} = \left(\frac {1}{N} \sum_ {i = 1} ^ {N} p (O _ {i} | \theta_ {0}) (1 - p (O _ {i} | \theta_ {0})) O _ {i} O _ {i} ^ {\prime}\right) ^ {- 1} \frac {1}{N} \sum_ {i = 1} ^ {N} ((Z _ {i} - p (O _ {i} | \theta_ {0})) O _ {i}) + o _ {p} (N ^ {- \frac {1}{2}}). ^ {3 1}
$$

In particular, $\hat { \theta _ { \mathrm { M L E } } } - \theta _ { 0 } = O _ { p } ( N ^ { - \frac { 1 } { 2 } } )$ . 

Proof. Identification is established by the standard KL-divergence argument and C7-(v). The rest of the proof is standard by Taylor expansions. See for example Newey and McFadden (1994) and Chang (2023). □ 

We defined the following IPW estimator: 

$$
\widehat {\tau} _ {\mathrm {I P W}} (d) = \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i}}{p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} \mu_ {i} (\mathbf {Y}; d) - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {1 - Z _ {i}}{1 - p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} \mu_ {i} (\mathbf {Y}; d).
$$

For brevity, we denote $p \left( O _ { i } | \theta _ { 0 } \right) = p _ { 0 } ( O _ { i } )$ . We define the following coefficients: 

$$
\beta_ {1, N} = \left(\frac {1}{N} \sum_ {i = 1} ^ {N} p _ {0} (O _ {i}) \left(1 - p _ {0} (O _ {i})\right) O _ {i} O _ {i} ^ {\prime}\right) ^ {- 1} \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (1; d) p _ {0} (O _ {i}) \left(1 - p _ {0} (O _ {i})\right) O _ {i} \right]\right),
$$

and, 

$$
\beta_ {0, N} = \left(\frac {1}{N} \sum_ {i = 1} ^ {N} p _ {0} (O _ {i}) \left(1 - p _ {0} (O _ {i})\right) O _ {i} O _ {i} ^ {\prime}\right) ^ {- 1} \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \operatorname {E} \left[ \mu_ {i} (0; d) p _ {0} (O _ {i}) \left(1 - p _ {0} (O _ {i})\right) O _ {i} \right]\right). ^ {3 2}
$$

A similar homophily condition is needed for the HAC variance estimation. 

C 8. For all sample size $N$ , 

$$
\frac {1}{N} \sum_ {i = 1} ^ {N} (\tau_ {i} (d; \eta) - (O _ {i} ^ {\prime} \beta_ {1, N} - O _ {i} ^ {\prime} \beta_ {0, N})) \sum_ {j \in \mathcal {B} (i; d)} (\tau_ {j} (d; \eta) - (O _ {i} ^ {\prime} \beta_ {1, N} - O _ {i} ^ {\prime} \beta_ {0, N})) \geq 0.
$$

We can similarly define a HAC type estimator as in (6). Let ${ \bf X } = \left( \begin{array} { l l l } { { 1 } } & { { 1 } } & { { \ldots , 1 } } \\ { { } } & { { } } & { { } } \\ { { Z _ { 1 } } } & { { Z _ { 2 } } } & { { \ldots Z _ { N } } } \end{array} \right) ^ { \prime } \in \nonumber$ RN×2. $\mathbb { R } ^ { N \times 2 }$ 

$$
\widehat {\mathrm {V}} _ {\mathrm {H A C}} ^ {\mathrm {o b s}} (d) = (\mathbf {X} ^ {\prime} \mathbf {X}) ^ {- 1} \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \mathbf {X} _ {i} \mathbf {X} _ {j} ^ {\prime} \hat {\epsilon} _ {i} \hat {\epsilon} _ {j} \mathbf {1} \{j \in \mathcal {B} (i; d) \}\right) (\mathbf {X} ^ {\prime} \mathbf {X}) ^ {- 1},
$$

where $\begin{array} { r } { \widehat { \epsilon } _ { i } = \frac { N _ { 1 } } { N } \frac { Z _ { i } } { p ( O _ { j } | \widehat { \theta } _ { \mathrm { M L E } } ) } \left( \mu _ { i } ( \mathbf { Y } ; d ) - O _ { i } ^ { \prime } \widehat { \beta } _ { 1 , N } \right) + \frac { N _ { 0 } } { N } \frac { 1 - Z _ { i } } { 1 - p ( O _ { j } | \widehat { \theta } _ { \mathrm { M L E } } ) } \left( \mu _ { i } ( \mathbf { Y } ; d ) - O _ { i } ^ { \prime } \widehat { \beta } _ { 0 , N } \right) } \end{array}$ .33 The formula is similar to (6), except now that the residuals are weighted with the inverse propensity scores. Note that the covariates $O _ { i }$ do not appear in the matrix $\mathbf { X }$ . Similarly, we can similarly define a S¨avje et al. (2018) type variance estimator as in (19): 

$$
\widehat {V} _ {\mathrm {S A H}} ^ {\mathrm {o b s}} (d) = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {Z _ {i} c _ {i} (d) \widehat {\epsilon} _ {i} ^ {2}}{p ^ {2} (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {(1 - Z _ {i}) c _ {i} (d) \widehat {\epsilon} _ {i} ^ {2}}{(1 - p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})) ^ {2}},
$$

We have the following linearization results. 

Proposition 8. Define: 

$$
\begin{array}{l} \widehat {\tau} _ {\mathrm {I P W}} ^ {\mathrm {L}} (d) = \frac {1}{N} \sum_ {i = 1} ^ {N} O _ {i} ^ {\prime} \beta_ {1, N} - \frac {1}{N} \sum_ {i = 1} ^ {N} O _ {i} ^ {\prime} \beta_ {0, N} \\ + \frac {1}{N} \sum_ {i = 1} ^ {n} \frac {Z _ {i} (\mu_ {i} (\mathbf {Y} ; d) - O _ {i} ^ {\prime} \beta_ {1 , N})}{p (O _ {i} | \theta_ {0})} - \frac {1}{N} \sum_ {i = 1} ^ {n} \frac {(1 - Z _ {i}) (\mu_ {i} (\mathbf {Y} ; d) - O _ {i} ^ {\prime} \beta_ {0 , N})}{1 - p (O _ {i} | \theta_ {0})}. \\ \end{array}
$$

Under C2, C3, C4, C6, and C7, 

$$
\widehat {\tau} _ {\mathrm {I P W}} ^ {\mathrm {L}} (d) - \widehat {\tau} _ {\mathrm {I P W}} (d) = o _ {p} (N ^ {- \frac {1}{2}}). \tag {67}
$$

Define 

$$
\begin{array}{l} \mathrm {V} \left(\widehat {\tau} _ {\mathrm {I P W}} ^ {\mathrm {L}} (d)\right) = \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {\mathrm {E} \left[ \left(\mu_ {i} (\mathbf {1} ; d) - O _ {i} ^ {\prime} \beta_ {1 , N}\right) ^ {2} \right]}{p _ {0} \left(O _ {i}\right)} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {\mathrm {E} \left[ \left(\mu_ {i} (\mathbf {0} ; d) - O _ {i} ^ {\prime} \beta_ {0 , N}\right) ^ {2} \right]}{1 - p _ {0} \left(O _ {i}\right)} (68) \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \left(\operatorname {E} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) - \left(O _ {i} ^ {\prime} \beta_ {1, N} - O _ {i} ^ {\prime} \beta_ {0, N}\right) \right]\right) ^ {2} (69) \\ + \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} \mathrm {E} \left[ \left(\mu_ {i} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - O _ {i} ^ {\prime} \beta_ {a, N}\right) \left(\mu_ {j} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - O _ {i} ^ {\prime} \beta_ {b, N}\right) \right] (70) \\ - \frac {1}{N ^ {2}} \sum_ {i; j \in \mathcal {B} (i; d), j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \mu_ {i} (\mathbf {a}; d) - O _ {i} ^ {\prime} \beta_ {a, N} \right] E \left[ \mu_ {j} (\mathbf {b}; d) - O _ {i} ^ {\prime} \beta_ {b, N} \right]. (71) \\ \end{array}
$$

In addition under C8, 

$$
\begin{array}{l} \mathrm {V} \left(\widehat {\tau} _ {\mathrm {I P W}} ^ {\mathrm {L}} (d)\right) \leq \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {\mathrm {E} \left[ (\mu_ {i} (\mathbf {1} ; d) - O _ {i} ^ {\prime} \beta_ {1 , N}) ^ {2} \right]}{p _ {0} (O _ {i})} + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \frac {\mathrm {E} \left[ (\mu_ {i} (\mathbf {0} ; d) - O _ {i} ^ {\prime} \beta_ {0 , N}) ^ {2} \right]}{1 - p _ {0} (O _ {i})} \\ + \frac {1}{N ^ {2}} \sum_ {{i; j \in \mathcal {B} (i; d), j \neq i}} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} \mathrm {E} \left[ (\mu_ {i} (\mathbf {a _ {i}}, \mathbf {b _ {j}}; d) - O _ {i} ^ {\prime} \beta_ {a, N}) (\mu_ {j} (\mathbf {a _ {i}}, \mathbf {b _ {j}}; d) - O _ {i} ^ {\prime} \beta_ {b, N}) \right] \\ \end{array}
$$

Provided that $N \times \mathrm { V } \left( \widehat { \tau } _ { \mathrm { I P W } } ^ { \mathrm { L } } ( d ) \right)$ is uniformly bounded below for all large N , ${ \widehat { \tau } } _ { \mathrm { I P W } } ( d )$ follows the asymptotic distribution: 

$$
\frac {\widehat {\tau} _ {\mathrm {I P W}} (d) - \mathrm {A M E} (d ; \eta)}{\sqrt {\mathrm {V} (\widehat {\tau} _ {\mathrm {I P W}} ^ {\mathrm {L}} (d))}} \xrightarrow {d} N (0, 1).
$$

Further more, for each $\alpha < \textstyle { \frac { 1 } { 2 } }$ , 

$\begin{array} { r } { ( i ) \ \operatorname* { l i m } _ { N  \infty } \mathbf { P r o b } ( z _ { 2 } \leq | \frac { \widehat \tau _ { \mathrm { H A } } ( d ) - \mathrm { A M E } ( d ; \eta ) } { \sqrt { \widehat \mathrm { V } _ { \mathrm { S A H } } ^ { \mathrm { o b s } } ( d ) } } | \leq z _ { 1 - \frac { \alpha } { 2 } } ) \geq 1 - \alpha ; } \end{array}$ 

(ii) additionally under C8, limN→∞ Prob(z α2 ≤ τbHAq(d)−AME(d;η)VobsHAC(d) ≤ z 1 − α2 ) ≥ 1 − α . 

Proof. We only show the calculation for the treated group. Calculation for the control group is similar. First note we have the following identitity: for any $\theta _ { 1 }$ and $\theta _ { 2 }$ 

$$
\left| p \left(O _ {i} \mid \theta_ {1}\right) - p \left(O _ {i} \mid \theta_ {2}\right) \right| \leq \left\| O _ {i} \right\| _ {2} \times \left\| \theta_ {1} - \theta_ {2} \right\| _ {2}, ^ {3 4} \tag {72}
$$

by C7-(iii) and a Taylor expansion argument. Note this implies $p ( O _ { i } | \theta _ { 2 } ) \in [ p ( O _ { i } | \theta _ { 1 } ) - | | O _ { i } | | _ { 2 } \times$ $| | \theta _ { 1 } - \theta _ { 2 } | | _ { 2 } , p ( O _ { i } | \theta _ { 1 } ) + | | O _ { i } | | _ { 2 } \times | | \theta _ { 1 } - \theta _ { 2 } | | _ { 2 } ] .$ 

$\begin{array} { r } { \frac { 1 } { N } \sum _ { i = 1 } ^ { N } \frac { Z _ { i } } { p ( O _ { i } | \widehat { \theta } _ { \mathrm { M L E } } ) } \mu _ { i } ( \mathbf { Y } ; d ) } \end{array}$ $\theta _ { 0 }$ 

$$
\begin{array}{l} \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} = \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \theta_ {0})} - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \theta_ {0})} (1 - p (O _ {i} | \theta_ {0})) O _ {i} ^ {\prime} (\widehat {\theta} _ {\mathrm {M L E}} - \theta_ {0}) \\ + \left(\widehat {\theta} _ {\mathrm {M L E}} - \theta_ {0}\right) ^ {\prime} \left(\frac {1}{2 N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p \left(O _ {i} \mid \tilde {\theta}\right)} \left(1 - p \left(O _ {i} \mid \tilde {\theta}\right)\right) O _ {i} O _ {i} ^ {\prime}\right) \left(\widehat {\theta} _ {\mathrm {M L E}} - \theta_ {0}\right), \tag {73} \\ \end{array}
$$

where $\tilde { \theta }$ is a point between $\theta _ { 0 }$ and $\widehat { \theta } _ { \mathrm { M L E } }$ . By (72) and $| | \widehat { \theta } _ { \mathrm { M L E } } - \theta _ { 0 } | | _ { 2 } = O _ { p } ( N ^ { - \frac { 1 } { 2 } } )$ , we have $\vert \vert \tilde { \theta } - \theta _ { 0 } \vert \vert ~ = ~ O _ { p } ( N ^ { - \frac { 1 } { 2 } } )$ , and $\operatorname { i n f } _ { i \in S _ { N } } p ( O _ { i } | { \tilde { \theta } } )$ is bounded away from 0 with probability one. 

Together with C2, C7-(ii) and Lemma A.6, we have that the second order term is of order $O _ { p } ( N ^ { - 1 } )$ . 

Let $\begin{array} { r } { H = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } p ( O _ { i } | \theta _ { 0 } ) ( 1 - p ( O _ { i } | \theta _ { 0 } ) O _ { i } O _ { i } ^ { \prime } } \end{array}$ . It’s clear that by C2, C6-(iii), C7-(ii) and (v), and an argument similar to A.2 : 

$$
H ^ {- 1} \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \theta_ {0})} p \left(O _ {i} \mid \theta_ {0}\right) \left(1 - p \left(O _ {i} \mid \theta_ {0}\right)\right) O _ {i}\right) - \beta_ {1, N} = o _ {p} (1). \tag {74}
$$

Hence (73) becomes 

$$
\frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} = \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \theta_ {0})} - \frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} - p (O _ {i} | \theta_ {0})}{p (O _ {i} | \theta_ {0})} O _ {i} ^ {\prime} \beta_ {1, N} + o _ {p} (N ^ {- \frac {1}{2}}).
$$

In particular, 

$$
\frac {1}{N} \sum_ {i = 1} ^ {N} \frac {Z _ {i} \mu_ {i} (\mathbf {Y} ; d)}{p (O _ {i} | \widehat {\theta} _ {\mathrm {M L E}})} - \bar {\mu} _ {1} (d) = \frac {1}{N} \sum_ {i = 1} ^ {N} O _ {i} ^ {\prime} \beta_ {1, N} - \bar {\mu} _ {1} (d) + \frac {1}{N} \sum_ {i = 1} ^ {n} \frac {Z _ {i} (\mu_ {i} (\mathbf {Y} ; d) - O _ {i} ^ {\prime} \beta_ {1 , N})}{p (O _ {i} | \theta_ {0})} + o _ {p} (N ^ {- \frac {1}{2}})
$$

This proves (67). The calculation of the asymptotic variance, variance upper bound is similar to that in Lemma A.4. The asymptotic distribution result follows similarly from Lemma A.5. 

For the variance estimation results, the proof is similar to that of Proposition 5. The only difference is to establish the equivalence between the variance estimator with estimated coefficient $\widehat { \theta } _ { \mathrm { M L E } }$ and the variance estimator with true coefficient $\theta _ { 0 }$ . This again can be shown to hold using a Taylor expansion argument together with C2, C3, C4 C6, C7-(ii) and (v). We omit the details here for brevity. □ 

# A.7 Effective Degree of Freedom Adjustment

When the distance of the AME is significant relative to the dataset’s total spatial coverage, the reliability of confidence intervals based on a normal approximation may diminish due to 

small effective sample sizes. This situation is similar to the cluster-robust inference settings with a small number of clusters. To improve the finite-sample coverage of the confidence intervals, we derive the effective degree of freedom adjustment as in (Imbens and Kolesar, 2012; Bell and McCaffrey, 2002; Young, 2015). 

Recall the regression interpretation of the Hajek estimator discussed in (5) and the HAC variance estimator $\widehat { \Sigma } _ { \mathrm { H A C } } ( d )$ in (18). Define $\mathbf { w } = ( 0 , 1 ) ^ { \prime }$ . Suppose we want to test the null hypothesis AME $( d ) = \tau _ { 0 }$ , the t-statistic under the null can be written as 

$$
\frac {\widehat {\tau} (d) - \tau_ {0}}{\sqrt {\mathbf {w} ^ {\prime} \widehat {\Sigma} _ {\mathrm {H A C}} \mathbf {w}}} = \frac {\frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \tau_ {0}}{\sqrt {\mathbf {w} ^ {\prime} (\mathbf {X} ^ {\prime} \mathbf {X}) ^ {- 1} \mathbf {w}}}}{\sqrt {\frac {\mathbf {w} ^ {\prime} \widehat {\Sigma} _ {\mathrm {H A C}} \mathbf {w}}{\mathbf {w} ^ {\prime} (\mathbf {X} ^ {\prime} \mathbf {X}) ^ {- 1} \mathbf {w}}}}.
$$

Define $\lambda = \mathbf { w } ^ { \prime } ( \mathbf { X } ^ { \prime } \mathbf { X } ) ^ { - 1 } \mathbf { X } \in \mathbb { R } ^ { 1 \times N }$ , $\mathbf { M } = \mathbf { I } - \mathbf { X } ^ { \prime } ( \mathbf { X } ^ { \prime } \mathbf { X } ) ^ { - 1 } \mathbf { X } \ \in \ \mathbb { R } ^ { N \times N }$ , and the weight matrix $\Omega \in \mathbb { R } ^ { N \times N }$ where $\Omega _ { i j } = \mathbf { 1 } \{ j \in B ( i ; d ) \}$ . Young (2015)’s effective degree of freedom adjustment in our setting is calculated as: 

$$
\mu = \mathbf {T r a c e} \left(\frac {N _ {1} N _ {0}}{N} \mathbf {M} (\Omega \circ (\lambda^ {\prime} \lambda)) \mathbf {M}\right), \tag {75}
$$

where $\cup$ denotes the pointwise (Hadamard) matrix product. The estimated variance is inflated by $\textstyle { \frac { 1 } { \mu } }$ , becoming w′ΣbHACw . We further define the quantity: $\frac { \mathbf { w } ^ { \prime } \overset { \sim } { \Sigma } _ { \mathbf { H A C W } } } { \mu }$ 

$$
\nu = 2 * \mathbf {T r a c e} \left(\left(\frac {N _ {1} N _ {0}}{N} \mathbf {M} (\Omega \circ (\lambda^ {\prime} \lambda)) \mathbf {M}\right) \left(\frac {N _ {1} N _ {0}}{N} \mathbf {M} (\Omega \circ (\lambda^ {\prime} \lambda)) \mathbf {M}\right)\right).
$$

We use the $\frac { \alpha } { 2 }$ -quantile and $\left( 1 - { \frac { \alpha } { 2 } } \right)$ -quantile of the t-distribution with $\frac { 2 \mu ^ { 2 } } { v }$ degree of freedom to construct the confidence interval. 

# B Simulation Results

# B.1 Simulation Designs

# Outcome Points and Intervention Nodes

Let $S = \{ 8 0 , 1 0 0 , 1 2 0 \}$ . We first generate a raster with $S \times S$ tiles, each of which is an outcome point. The side length of each tile is 1 generic unit. The untreated potential outcome for the outcome point $x$ , $Y _ { x } ( 0 )$ , is randomly drawn from the standard normal distribution. 

For point interventions, we coarsen the raster into $\frac { S ^ { 2 } } { 1 0 0 }$ tiles and random sample half of the tiles. For each sampled tile, we add perturbations to the center of the tile twice to create a pair of intervention points. 

For polygon interventions, we subsample $\frac { S ^ { 2 } } { 1 0 }$ tiles, use the sampled tiles and the Voronoi tessellation to generate polygons. To construct the set of intervention nodes, we randomly sample $\frac { S ^ { 2 } } { 2 0 0 }$ polygons and, for each sampled polygons, we randomly sample an adjacent polygon. This gives a total of $\frac { S ^ { 2 } } { 1 0 0 }$ polygons as intervention nodes. 

# Data Generating Process

Let $\Gamma ( d , a , b )$ denote the density of a gamma distribution at value d with shape parameter a and scale parameter b. For an outcome point $x$ , we define the effect function: 

$$
f _ {x} (d) = 3 \alpha_ {x} \times (\Gamma (d; 1, 1) - \Gamma (d; 5, 0. 5)) \times \max  \left(\left(1 - \frac {d ^ {2}}{3 6}\right), 0\right), \tag {77}
$$

where $\alpha _ { x }$ captures the heterogeneity in treatment effects among outcome points. $\alpha _ { x }$ are generated from kriging interpolation of some randomly generated values on a coarsened raster. The term $\textstyle \operatorname* { m a x } \left( \left( 1 - { \frac { d ^ { 2 } } { 3 6 } } \right) , 0 \right)$ is used to guarantee that there is no treatment effect beyond 6 units of distance. 

For the additive-effect case, the outcome at point $x$ , $Y _ { x } ( \mathbf { Z } )$ , is generated by the formula 

$$
Y _ {x} (\mathbf {Z}) = Y _ {x} (0) + \sum_ {i = 1} ^ {n} f _ {x} \left(d _ {i x}\right) Z _ {i}, \tag {78}
$$

where $n$ denotes the number of intervention nodes, and $d _ { i x }$ is the distance from the outcome point $x$ to the intervention node $i$ . 

For the interactive-effect case, we define another effect function: 

$$
g _ {x} (d) = 3 \alpha_ {x} \times (\Gamma (d; 5, 0. 5)) \times \max  \left(\left(1 - \frac {d ^ {2}}{3 6}\right), 0\right). \tag {79}
$$

The outcome at point $x$ , $Y _ { x } ( \mathbf { Z } )$ , is generated by the formula 

$$
Y _ {x} (\mathbf {Z}) = Y _ {x} (0) + \sum_ {i = 1} ^ {n} f _ {x} \left(d _ {i x}\right) Z _ {i} + \sum_ {i = 1} ^ {n} g _ {x} \left(d _ {i x}\right) Z _ {i} Z _ {\mathcal {N} (i)}, \tag {80}
$$

where $\mathcal { N } ( i )$ denotes the intervention node that is closest to the intervention node $i$ . When there are multiple polygons that are adjacent to the $i$ th polygon, we pick the polygon with the smallest assigned index to be $\mathcal { N } ( i )$ . 

# Experimental Designs

For our simulation, we use a Bernoulli design where $Z _ { i } ~ = ~ 1$ with probability 0.5 for all intervention nodes. We run the simulation for 2000 times. Denote the $p$ th assignment of $\mathbf { Z }$ as $\mathbf { Z } _ { p }$ and the value of $Z _ { i }$ under the assignment as $Z _ { p i }$ . To obtain AME, we first calculate 

$$
\tau_ {i x} (\eta) \approx \frac {\sum_ {p = 1} ^ {2 0 0 0} Z _ {p i} Y _ {x} \left(\mathbf {Z} _ {p}\right)}{\sum_ {p = 1} ^ {2 0 0 0} Z _ {p i}} - \frac {\sum_ {p = 1} ^ {2 0 0 0} \left(1 - Z _ {p i}\right) Y _ {x} \left(\mathbf {Z} _ {p}\right)}{\sum_ {p = 1} ^ {2 0 0 0} \left(1 - Z _ {p i}\right)}. \tag {81}
$$

Then $\tau _ { i } ( d ; \eta )$ and the AME can be constructed following their definitions. 

# B.2 More Point-Intervention Simulation Results in Section 7

We report simulation results for the Smoothed Hajek estimator with a triangular kernel and a bandwidth 1. Figure 14 reports results on MSE, Figure 15 on coverage and half length for the additive case, and Figure 16 on coverage and half length for the interactive case 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ffc3e44c7c9cd5addb882726c8308dd4857f30a5ab71dbc8b6c5b175c72ba1e7.jpg)



(a) MSE for the additive case (30)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/93bbe9f934505700a6572f393191220081de5a3fe5455a9bdd43618cd8889acf.jpg)



(b) MSE for the interactive case (31)



Figure 14: The left and right figures report the Mean Squared Errors of the Smoothed Hajek estimator in the additive-effect case and the interactive-effect case, respectively.


# B.3 Results for a Polygon-Intervention Simulation

We report simulation results for the Hajek estimator with a polygon intervention simulation. Figure 17 reports results on MSE, Figure 18 on coverage and half length for the additive case, and Figure 19 on coverage and half length for the interactive case. 

We report simulation results for the Smoothed Hajek estimator with a polygon intervention simulation. We use a triangular kernel and bandwidth 1. Figure 20 reports results on MSE, Figure 21 reports results on coverage and half length for the additive case, and Figure 22 on coverage and half length for the interactive case. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/0e7c8c47f62fd48f74dc2ee4a782d40ec0f10e1701aa77197e0aaf9ff6e05c90.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/aa387362f86d98bf3f249e83563f2a85cab9b5a20b980b12390adf3c1b2dc978.jpg)



Figure 15: Point-intervention simulation results on the coverage rates and half lengthes of two-sided 95% confidence intervals with the Smoothed Hajek estimator and different variance estimators in the additive effect case (30). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/70d3998d04574db9046f54383d09899a76f6187226831341177736c7a3ceac83.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/591422434fb5bacd17dd02dc6a04d95ab356ccb02da0e87c6d1856c027a1bab0.jpg)



Figure 16: Point-intervention simulation results on the coverage rates and half lengthes of two-sided 95% confidence intervals with the Smoothed Hajek estimator and different variance estimators in the interactive effect case (31). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ea6507b38a5b0bc806d1c13a13810e000677a4a1459b2fa9081bad796018dfda.jpg)



(a) MSE for the additive case (30)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/51bcba26f9b20004d9820bad3fb78efac14c85e818f447059873d1a72915fb2c.jpg)



(b) MSE for the interactive case (31)



Figure 17: The left and right figures report the Mean Squared Errors of the Hajek estimator in the additive-effect case and the interactive-effect case, respectively.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/603b4d360218dc95e10007eda7e2ca616ce3413d7eab4410cdf75cfd5be4c870.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ba10287bf511839b256f2041fff5324be9bb8a2c7a3c2540536eb5c6af99a4c2.jpg)



Figure 18: Polygon-intervention simulation results on the coverage rates and half lengthes of two-sided 95% confidence intervals with the Hajek estimator and different variance estimators in the additive effect case (30). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/f825c8629b8d8e17b8f194acca0d5c0b315857af3d954a3337db20bd3088af18.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/3a4feb0bd476ff7c39a2fd627dd81cf261aae85b4cea55bcfb9e8d1cfbdbefd6.jpg)



Figure 19: Polygon-intervention simulation results on the coverage rates and half lengths of two-sided 95% confidence intervals with the Hajek estimator and different variance estimators in the interactive effect case (31). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/9d1dda55b4a9a276a2c1cf86b6427561571b515efe01b66a09db500049b92637.jpg)



(a) MSE for the additive case (30)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/4b2cab9e5de992937d0cdd4acb4d6df212b52fbe5301d809fa6e7a4133cf326a.jpg)



(b) MSE for the interactive case (31)



Figure 20: The left and right figures report the Mean Squared Errors of the Smoothed Hajek estimator in the additive-effect case and the interactive-effect case, respectively.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/9e33ebde8a3b48d57cbbcb670e08e2c95ae2ecfe0998ac894ce5c46fcbf66fd1.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/3ca9a3263e35f5e72f3bda23943bd554434a2c8c2be23ac54dc2c02306d620b8.jpg)



Figure 21: Polygon-intervention simulation results on the coverage rates and half lengthes of two-sided 95% confidence intervals with the Smoothed Hajek estimator and different variance estimators in the additive effect case (30). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/b88f61f43b4fd60262a00a8e0e835f585ec6e6dc75894305f02c41e027b491cf.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/0941d610f98a1fac596991bb1c64bb7bf5d8b382df2462d0bbd05f077e62cfc3.jpg)



Figure 22: Polygon-intervention simulation results on the coverage rates and half lengthes of two-sided 95% confidence intervals with the Smoothed Hajek estimator and different variance estimators in the interactive effect case (31). The sample size is 144. HAC refers to the CI with the HAC variance estimator in (18) and a normal critical value. The length and coverage of the HAC CI is assesed with respect to the cases where HAC estimator returns a nonnegative value. HAC PD refers to the CI with positive-semidefinite HAC variance estimator in (21) and a normal critical value. HAC (edof) refers to the CI with HAC variance estimator in (18) and empirical degree of freedom adjustment. HAC PD (edof) refers to the CI with HAC variance estimator in (21) and empirical degree of freedom adjustment. SAH refers to the CI with SAH variance estimator (19).


# B.4 Additional Results for Empirical Applications

Figure 23 provides more details on the construction of the new outcome variable for use in our re-analysis of the Jayachandran et al. (2017) experiment. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/7d2b34d1675ae41c36c8fee692e2b81ed57ca5be66c4709daec94c0983acc34f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/139456953ef662ecb156fccc8a1564a5b24767c2b4e2cfb7123235e8dfe2feed.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/3e644ecad5bf8db3b937c0f9f6be766a3ff50f050cc7b5d1e2b26d8b7d4287e2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/f184ed5d0b0de923691ac6da539ca3fdb9f5b90c0fd1d286da23048fa0eac987.jpg)



- - 	--



- 	-
 



Figure 23: Top plot: Study area of randomized control trial for a PES program in Hoima and Kibaale district in Uganda, from Jayachandran et al. (2017). Boundaries of treatment (60) and control (61) villages were digitized using publicly available data and published maps. Bottom plots: The Global Forest Cover (GFC) dataset over a subset of the study area showing forest cover for 2012, 2013, and forest loss in 2012 (Hansen et al., 2013).


# C Weaker assumptions on the extent of interference

# C.1 Weaker Assumptions on the Extent of Interference

This section extends our inferential results on the Hajek estimator by relaxing the local interference assumption C3 such that spatial dependence does not have to be contained within a strict distance cutoff. We allow for more general spatial “near-epoch dependence” (Jenish and Prucha, 2012) and provide results on root-N consistency, asymptotic normality, and HAC variance estimation.35 We do not pursue the most general results (e.g. the most relaxed conditions on data moments and correlation structures). Rather, we impose assumptions that are commensurable with the assumptions in the main text. Proofs of root-N consistency and asymptotic normality are standard and follow as special cases of the results in Jenish and Prucha (2012). The proof of the consistency of the HAC variance estimator is relatively new, as far as we know.36 

We introduce the near-epoch dependence concept. Our version is a simplified version from Jenish and Prucha (2012). Denote the sample size by $N$ . Let $S _ { N }$ be a set of intervention nodes equipped with a distance metric $\gamma$ that satisfies positivity and the triangle inequality. Let $H _ { N } = \{ H _ { i , N } , i \in S _ { N } \}$ and $V _ { N } = \{ V _ { i , N } , i \in S _ { N } \}$ be two arbitrary sets of random variables. Define $\mathcal { F } _ { i } ( s ) \ : = \ : \sigma \left( V _ { j , N } , j \in S _ { n } , \gamma ( i , j ) \le s \right)$ , the $\sigma$ -field generated by the random variables $\{ V _ { j , N } \}$ located in the s-neighborhood of the intervention node $i$ . Denote the $L _ { p }$ -norm of a random variable as $\vert \vert X \vert \vert _ { p } = ( E [ \vert X \vert ^ { p } ] ) ^ { \frac { 1 } { p } }$ . 

Definition 1 ( $L _ { 2 }$ -NED). Let $H _ { N } = \{ H _ { i , N } , i \in S _ { N } , N \geq 1 \}$ be a random field with $| | H _ { i , N } | | _ { 2 } <$ $\infty$ . Let $V _ { N } = \{ V _ { i , N } , i \in \mathcal { S } _ { N } , N \ge 1 \}$ be a random field, where $| S _ { N } | = N$ , and let $g _ { N } =$ $\{ g _ { i , N } , i \in S _ { N } , N \ge 1 \}$ be a set of positive constants. Then the random field $H _ { N }$ is said to be $L _ { 2 } ( g _ { N } )$ -near-epoch dependent on the random field $V _ { N }$ if, for all $i$ , 

$$
\left| \left| H _ {i, N} - \operatorname {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} \leq g _ {i, N} \psi (s), \tag {82}
$$

for a non-increasing sequence $\psi ( s ) \geq 0$ with $\begin{array} { r } { \operatorname* { l i m } _ { s \to \infty } \psi ( s ) = 0 } \end{array}$ . $H _ { N }$ is said to be $L _ { \mathrm { 2 } } - N E D$ on $V _ { N }$ of size −λ if $\psi ( s ) = O ( s ^ { - \mu } )$ for some $\mu > \lambda > 0$ . Further more, if $\mathrm { s u p } _ { N } \mathrm { s u p } _ { i \in { \mathcal { S } } _ { N } } g _ { i , N } < \infty$ , then $H _ { N }$ is said to be uniformly $L _ { 2 }$ -NED on $V _ { N }$ . 

We establish several implications of the $L _ { 2 }$ -NED property. The following lemma is inspired by Theorem 17.5 in Davidson (1994). 

Let $\mathcal { F } _ { i } ^ { c } ( s ) = \sigma \left( V _ { j , N } , j \in S _ { n } , \gamma ( i , j ) > s \right)$ , the $\sigma$ -field generated by the random variables $\{ V _ { j , N } \}$ that are more than $s$ distance away from the intervention node $i$ . The first lemma below states that $L _ { 2 }$ -NED random variables are nearly constant when conditioned on random variables $\{ V _ { j , N } \}$ that are far away. This is an adaption of the mixingale property from the time series literature to our spatial process setting.37 

Lemma C.1. Let $V _ { N }$ be a set of independent random variables, and $H _ { N }$ be a set of zero-mean $L _ { 2 } ( g _ { N } )$ -NED random variables on $V _ { N }$ , then $| | \mathrm { E } \left[ H _ { i , N } | \mathcal { F } _ { i } ^ { c } ( s ) \right] | | _ { 2 } \leq g _ { i , N } \psi ( s )$ . 

Proof. We have the following chain of inequalities: 

$$
\left. \left| \left| \mathrm {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} ^ {c} (s) \right] \right| \right| _ {2} \right. \tag {83}
$$

$$
\leq \left\| \operatorname {E} \left[ \left(H _ {i, N} - \operatorname {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} (s) \right]\right) \mid \mathcal {F} _ {i} ^ {c} (s) \right] \right\| _ {2} + \left\| \operatorname {E} \left[ \left(\operatorname {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} (s) \right]\right) \mid \mathcal {F} _ {i} ^ {c} (s) \right] \right\| _ {2} \tag {84}
$$

$$
\leq \left| \left| H _ {i, N} - \mathrm {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} + \left| \left| \mathrm {E} \left[ \left(\mathrm {E} \left[ H _ {i, N} \mid \mathcal {F} _ {i} (s) \right]\right) \mid \mathcal {F} _ {i} ^ {c} (s) \right] \right| \right| _ {2} \tag {85}
$$

$$
\leq g _ {i, N} \psi (s) + 0. \tag {86}
$$

(84) is by the triangle inequality and (85) is by Jensen’s inequality. (86) follows by definition and that $\mathrm { E } \left[ \left( \mathrm { E } \left[ H _ { i , N } \vert \mathcal { F } _ { i } ( s ) \right] \right) \vert \mathcal { F } _ { i } ^ { c } ( s ) \right] = 0$ . Note that $\mathrm { E } \left[ \left( \mathrm { E } \left[ H _ { i , N } \vert \mathcal { F } _ { i } ( s ) \right] \right) \vert \mathcal { F } _ { i } ^ { c } ( s ) \right] = 0$ by Theorem 10.22 in Davidson (1994) following the facts that $H _ { i , N }$ is mean zero and that $\mathcal { F } _ { i } ( s )$ and $\mathcal { F } _ { i } ^ { c } ( s )$ are independent $\sigma$ -fields. □ 

We shall inherit C1 and C2. We relax C3 with the following weak-dependence assumption. 

C 9. For any d-circle average outcomes and for all sample sizes, $\{ \mu _ { i } ( \mathbf { Y } , d ) , i \ \in \ S _ { N } \}$ is uniform $L _ { 2 }$ -NED of size -4 on the random field $\{ Z _ { i } , i \in S _ { N } \}$ . 38 

Denote the corresponding NED-coefficient function as $\psi _ { d } ( \ v r ) : \mathbb { R } _ { + } \to \mathbb { R } _ { + }$ .39 Now recall the class of estimators defined in Section A.5: 

$$
\widehat {\tau} \left(\mu_ {1}, \mu_ {0}; d\right) = \mu_ {1} - \mu_ {0} + \left(\frac {1}{N p} \sum_ {i = 1} ^ {n} \mathbf {Z} _ {i} \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {1}\right) - \frac {1}{N (1 - p)} \sum_ {i = 1} ^ {n} \left(1 - \mathbf {Z} _ {i}\right) \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {0}\right)\right). \tag {87}
$$

This class contains the HT estimator ( $\mu _ { 1 } = \mu _ { 0 } = 0$ ) and the linearized Hajek estimator ( $\mu _ { 1 } = \bar { \mu } _ { 1 } ( d )$ and $\mu _ { 0 } = \bar { \mu } _ { 0 } ( d )$ ) as special cases. The following lemma states that provided 

$$
\begin{array}{r l} & {\{\mu_ {i} (\mathbf {Y}; d), i \in \mathcal {S} _ {N} \} \mathrm {i s u n i f o r m L 2 - N E D w i t h s i z e - 4 , s o i s t h e t r a n s f o r m a t i o n} \widehat {\tau} _ {i, N} (\mu_ {1}, \mu_ {0}; d) =} \\ & {\left(\frac {\mathbf {Z} _ {i}}{N p} - \frac {(1 - \mathbf {Z} _ {i})}{N (1 - p)}\right) \mu_ {i} (\mathbf {Y}; d) + \frac {1}{N} (\mu_ {1} - \mu_ {0}) - \left(\frac {\mathbf {Z} _ {i} \mu_ {1}}{N p} - \frac {(1 - \mathbf {Z} _ {i}) \mu_ {0}}{N (1 - p)}\right) - \tau_ {i} (d; \eta). \mathrm {N o t e E} [ \widehat {\tau} _ {i, N} (\mu_ {1}, \mu_ {0}; d) ] = 0.} \end{array}
$$

Lemma C.2. Under C1 and C9, for each pair $\mu _ { 1 } , \mu _ { 0 } \in \mathbb { R }$ , $\left\{ \widehat { \tau } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d ) , i \in \mathcal { S } _ { N } \right\}$ is uniform $L _ { 2 } ( g _ { N } ^ { \tau } )$ -NED with size −4 on the random field $\{ Z _ { i } , i \in S _ { N } \}$ . We have $\begin{array} { r } { g _ { i , N } ^ { \tau } \le \frac { C } { N p ( 1 - p ) } } \end{array}$ , where $C$ is independent of $N$ , $\mu _ { 1 }$ and $\mu _ { 0 }$ . 

Proof. Note $\mathrm { E } \left[ \widehat { \tau } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d ) | \mathcal { F } _ { i } ( s ) \right]$ is equal to: 

$$
\left(\frac {\mathbf {Z} _ {i}}{N p} - \frac {(1 - \mathbf {Z} _ {i})}{N (1 - p)}\right) \operatorname {E} \left[ \mu_ {i} (\mathbf {Y}; d) \mid \mathcal {F} _ {i} (s) \right] + \frac {1}{N} \left(\mu_ {1} - \mu_ {0}\right) - \left(\frac {\mathbf {Z} _ {i} \mu_ {1}}{N p} - \frac {(1 - \mathbf {Z} _ {i}) \mu_ {0}}{N (1 - p)}\right) - \tau_ {i} (d), \tag {88}
$$

because the σ field generated by Zi is a subset of Fi(s). Noting that | ZiN p − (1−Zi)N(1−p)| $\sigma$ $Z _ { i }$ $\mathcal { F } _ { i } ( s )$ $\begin{array} { r } { | \frac { { \bf { Z } } _ { i } } { N p } - \frac { ( 1 - { \bf { Z } } _ { i } ) } { N ( 1 - p ) } | \leq \frac { 1 } { N p ( 1 - p ) } } \end{array}$ (1-Zi) ≤ 1N p(1−p) , we have 

$$
| | \widehat {\tau} _ {i, N} (d) - \operatorname {E} \left[ \widehat {\tau} _ {i, N} (d) \mid \mathcal {F} _ {i} (s) \right] | | _ {2} \leq \frac {1}{N p (1 - p)} | | \mu_ {i} (\mathbf {Y}, d) - \operatorname {E} \left[ \mu_ {i} (\mathbf {Y}, d) \mid \mathcal {F} _ {i} (s) \right] | | _ {2} \tag {89}
$$

$$
\lesssim \frac {1}{N p (1 - p)} \psi_ {d} (s), \tag {90}
$$

where by using the symbol $\lesssim$ we suppress a constant independent of $N$ . Importantly, note that the constant is independent of $\mu _ { 1 }$ and $\mu _ { 0 }$ . □ 

The next lemma is a covariance inequality which, together with a later assumption on the density of intervention nodes, implies a $O _ { p } \big ( \frac { 1 } { \sqrt { N } } \big )$ convergence rate for the HT and Hajek estimators.40 

Lemma C.3. Let $d _ { i j } = \gamma ( i , j )$ denote the distance between intervention node i and intervention node j. Under C1, C2 and $\textit { C }$ , for any $\mu _ { 1 }$ and $\mu _ { 0 }$ such that $| \mu _ { 1 } | < B$ and $| \mu _ { 0 } | < B$ 

for some positive constant $B$ , 

$$
\left| \operatorname {C o v} \left(\widehat {\tau} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right), \widehat {\tau} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right)\right) \right| \leq \frac {C (B , p)}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \tag {91}
$$

where $C ( B , p )$ denotes a generic constant dependent on $B$ and $p$ , and independent of $N$ . 

Proof. We write $\widehat { \tau } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ as $\widehat { \tau } _ { i } ( d )$ , and $\widehat { \tau } _ { j , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ as $\widehat { \tau } _ { j } ( d )$ for brevity. Note $E [ \widehat { \tau } _ { i } ( d ) ] = 0$ and $E [ \widehat { \tau } _ { j } ( d ) ] = 0$ . we have: 

$$
\left| \operatorname {C o v} \left(\widehat {\tau} _ {i} (d), \widehat {\tau} _ {j} (d)\right) \right| = \left| \mathrm {E} \left[ \widehat {\tau} _ {i} (d) \widehat {\tau} _ {j} (d) \right] \right| \tag {92}
$$

$$
\leq \left| \operatorname {E} \left[ \left(\widehat {\tau} _ {i} (d) - \operatorname {E} \left[ \widehat {\tau} _ {i} (d) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right]\right) \widehat {\tau} _ {j} (d) \right] \right| + \left| \operatorname {E} \left[ \operatorname {E} \left[ \widehat {\tau} _ {i} (d) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right] \widehat {\tau} _ {j} (d) \right] \right| \tag {93}
$$

$$
\leq \left| \left| \widehat {\tau} _ {i} (d) - \mathrm {E} \left[ \widehat {\tau} _ {i} (d) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right] \right| \right| _ {2} \times \left| \left| \widehat {\tau} _ {j} (d) \right| \right| _ {2} \tag {94}
$$

$$
+ \left| \mathrm {E} \left[ \mathrm {E} \left[ \widehat {\tau} _ {i} (d) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right] \mathrm {E} \left[ \widehat {\tau} _ {j} (d) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right] \right] \right| \tag {95}
$$

$$
\lesssim \frac {1}{N ^ {2}} C _ {1} (B, p) \psi_ {d} \left(\frac {d _ {i j}}{3}\right) + | | \mathrm {E} [ \widehat {\tau} _ {i} (d) | \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) ] | | _ {2} \times | | \mathrm {E} [ \widehat {\tau} _ {j} (d) | \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) ] | | _ {2} \tag {96}
$$

$$
\lesssim \frac {C _ {1} (B , p)}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) + \frac {1}{N} C _ {2} (B, p) \times | | \mathrm {E} [ \widehat {\tau} _ {j} (d) | \mathcal {F} _ {j} ^ {c} \left(\frac {d _ {i j}}{3}\right) ] | | _ {2} \tag {97}
$$

$$
\lesssim \frac {C _ {1} (B , p)}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) + \frac {C _ {2} (B , p)}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \leq \frac {C (B , p)}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \tag {98}
$$

where $C _ { 1 } ( B , p )$ and $C _ { 2 } ( B , p )$ denote generic constants depending on $B$ and $p$ , and by using the symbol $\lesssim$ we suppress constants that are independent of $N$ . (93) follows by the triangle inequality. (94) follows by the Cauchy-Schwarz inequality. (96) follows by bounding the size of $\widehat { \tau } _ { j } ( d )$ using C2 and by using (90). (97) follows by bounding the size of $\widehat { \tau } _ { i } ( d )$ using C2 and by the Jensen’s inequality and the fact that $\begin{array} { r } { \mathcal { F } _ { i } \big ( \frac { d _ { i j } } { 3 } \big ) \subset \mathcal { F } _ { j } ^ { c } \big ( \frac { d _ { i j } } { 3 } \big ) } \end{array}$ . (98) follows by using Lemma C.1 and reading off the size of $g _ { i , N }$ from (90). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ccf737df83ce22269152ec488f861c15792b3d63de3cc8437d02a1e9bc2bbd51.jpg)


We need a stronger assumption to control the concentration of the intervention nodes. 

The following assumption prevents nodes from concentrating around a few locations. 

C 10. For each radius $s \geq 0$ and all large sample sizes $N$ , there exists a $C _ { b } > 0$ such that $\begin{array} { r } { \operatorname* { s u p } _ { i \in { \cal S } _ { N } } | B _ { N } ( i ; s ) | \le b ( s ) } \end{array}$ , with $b ( s ) = C _ { b } s ^ { 2 }$ . 

The following lemma establishes the root-N convergence rate of the HT estimator and the asymptotic equivalence of the Hajek estimator. We note that all lemmas and propositions below are stated for each distance value $d$ . 

Lemma C.4. Under C1, C2, C9, and C10, $\mathrm { E } \left[ \hat { \tau } _ { \mathrm { H T } } ( d ) \right] = \mathrm { A M E } ( d ; \eta )$ , and Var $( \hat { \tau } _ { \mathrm { H T } } ( d ) ) =$ $O ( \textstyle { \frac { 1 } { N } } )$ . Consider the estimator $\widehat { \tau } _ { \mathrm { H A } } ( d )$ defined in (15) and its asymptotic linear expansion $\widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d )$ defined in Lemma A.3, we have $\sqrt { N } ( \widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) - \widehat { \tau } _ { \mathrm { H A } } ( d ) ) = o _ { p } ( 1 )$ . 

Proof. Proof of unbiasedness of the HT estimator is by Proposition 1. To bound the variance, note 

$$
\mathrm {V a r} \left(\widehat {\tau} _ {\mathrm {H T}} (d)\right) = \mathrm {V a r} \left(\sum_ {i = 1} ^ {N} \widehat {\tau} _ {i, N} (0, 0; d)\right) = \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \mathrm {C o v} \left(\widehat {\tau} _ {i, N} (0, 0; d), \widehat {\tau} _ {j, N} (0, 0; d)\right)
$$

Write $\widehat { \tau } _ { i , N } ( 0 , 0 ; d )$ as $\widehat { \tau } _ { i } ( d )$ and $\widehat { \tau } _ { j , N } ( 0 , 0 ; d )$ as $\widehat { \tau } _ { j } ( d )$ . We have, for each $i$ , 

$$
\left| \sum_ {j = 1} ^ {N} \operatorname {C o v} \left(\widehat {\tau} _ {i} (d), \widehat {\tau} _ {j} (d)\right) \right| \leq \sum_ {j = 1} ^ {N} | \operatorname {C o v} \left(\widehat {\tau} _ {i} (d), \widehat {\tau} _ {j} (d)\right) | \lesssim \frac {1}{N ^ {2}} \sum_ {j = 1} ^ {N} \psi_ {d} \left(\frac {d _ {i j}}{3}\right), \tag {99}
$$

where by using $\lesssim$ we suppress a constant independent of $N$ . To prove the lemma, we only need to show the quantity $\textstyle \sum _ { j = 1 } ^ { N } \psi _ { d } ( { \frac { d _ { i j } } { 3 } } )$ is uniformly bounded for all $i$ and all sample sizes $N$ . 

By C9, there exists a $C _ { 2 }$ such that $\begin{array} { r } { \psi _ { d } \left( \frac { k } { 3 } \right) \le C _ { 2 } \operatorname* { m i n } \{ 1 , \frac { 1 } { k ^ { 2 + \epsilon } } \} } \end{array}$ .41 We have the following 

inequality: 

$$
\begin{array}{l} \sum_ {j = 1} ^ {N} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \leq \lim  _ {K \rightarrow \infty} \left(\sum_ {k = 0} ^ {K} (b (k + 1) - b (k)) \psi_ {d} (\frac {k}{3})\right) (100) \\ \leq C _ {2} \lim  _ {K \rightarrow \infty} \left(\sum_ {k = 0} ^ {K} (b (k + 1) - b (k)) \min  \{1, \frac {1}{k ^ {2 + \epsilon}} \}\right) (101) \\ = C _ {2} \lim  _ {K \rightarrow \infty} \left(- b (0) 1 + b (1) (1 - 1) + \sum_ {k = 2} ^ {K} b (k) \left(\frac {1}{(k - 1) ^ {2 + \epsilon}} - \frac {1}{k ^ {2 + \epsilon}}\right) + b (K + 1) \times \frac {1}{K ^ {2 + \epsilon}}\right) (102) \\ \lesssim - b (0) + \lim  _ {K \rightarrow \infty} \left(\sum_ {k = 2} ^ {K} b (k) \left(\frac {1}{(k - 1) ^ {2 + \epsilon}} - \frac {1}{k ^ {2 + \epsilon}}\right)\right) + \lim  _ {K \rightarrow \infty} \left(b (K + 1) \times \frac {1}{K ^ {2 + \epsilon}}\right) (103) \\ \lesssim - b (0) + C _ {b} \lim  _ {K \rightarrow \infty} \left(\sum_ {k = 2} ^ {K} k ^ {2} \left(\frac {1}{(k - 1) ^ {2 + \epsilon}} - \frac {1}{k ^ {2 + \epsilon}}\right)\right) + \lim  _ {K \rightarrow \infty} \left(\frac {C _ {b} (K + 1) ^ {2}}{K ^ {2 + \epsilon}}\right) (104) \\ <   \infty , (105) \\ \end{array}
$$

where $b ( \cdot ) : \mathbb { R } _ { + } \to \mathbb { R } _ { + }$ is the neighborhood size function defined in C10. (100) is by an inclusion criteria and the fact that $\psi$ is a non-increasing function. (101) follows because $b$ by definition is a non-decreasing function. (102) and (103) are purely algebraic. (104) follows by C10 and the fact that $\begin{array} { r } { \frac { 1 } { ( k - 1 ) ^ { 2 + \epsilon } } - \frac { 1 } { k ^ { 2 + \epsilon } } \geq 0 } \end{array}$ for all $k \geq 2$ . (105) follows because $\begin{array} { r } { \sum _ { k = 2 } ^ { \infty } k ^ { 2 } \left( \frac { 1 } { ( k - 1 ) ^ { 2 + \epsilon } } - \frac { 1 } { k ^ { 2 + \epsilon } } \right) < \infty } \end{array}$ by calculus.42 Note the bound is uniform in $i$ and sample sizes $N$ . 

Linearization proof for the Hajek estimator is identical to the one in Lemma A.3. □ 

The following lemma is needed to bound the bias of the HAC variance estimator later. 

Lemma C.5. Under C9 and C10, for any $\{ b _ { N } \} _ { N = 1 } ^ { \infty }$ such that $b _ { N } \to \infty$ , we have 

$$
\lim  _ {N \rightarrow \infty} \sup  _ {i \in \mathcal {S} _ {N}} \sum_ {\{j: d _ {i j} > b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) = 0. \tag {106}
$$

Proof. For each $i$ , similar to the proof in Lemma C.4 we have 

$$
\begin{array}{l} 0 \leq \sum_ {\{j: d _ {i j} > b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \leq \left(\sum_ {k = \lfloor b _ {N} \rfloor} ^ {\infty} (b (k + 1) - b (k)) \psi_ {d} (\frac {k}{3})\right) (107) \\ \leq b \left(\left\lfloor b _ {N} \right\rfloor\right) \psi_ {d} \left(\frac {\left\lfloor b _ {N} \right\rfloor}{3}\right) + \sum_ {k = \left\lfloor b _ {N} \right\rfloor + 1} ^ {\infty} b (k) \left(\frac {1}{(k - 1) ^ {2 + \epsilon}} - \frac {1}{k ^ {2 + \epsilon}}\right). (108) \\ \end{array}
$$

Note the right-hand side is independent of $i$ . Thus we have: 

$$
\begin{array}{l} \lim  _ {N \rightarrow \infty} \sup  _ {i \in \mathcal {S} _ {N}} \sum_ {\{j: d _ {i j} > b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) (109) \\ \leq \lim  _ {N \rightarrow \infty} \left(b (\lfloor b _ {N} \rfloor) \psi_ {d} (\frac {\lfloor b _ {N} \rfloor}{3}) + \sum_ {k = \lfloor b _ {N} \rfloor + 1} ^ {\infty} b (k) \left(\frac {1}{(k - 1) ^ {2 + \epsilon}} - \frac {1}{k ^ {2 + \epsilon}}\right)\right) = 0, (110) \\ \end{array}
$$

since under C9 and C10 we have $\begin{array} { r } { \sum _ { k = 2 } ^ { \infty } b ( k ) \left( \frac { 1 } { ( k - 1 ) ^ { 2 + \epsilon } } - \frac { 1 } { k ^ { 2 + \epsilon } } \right) < \infty } \end{array}$ as shown in Lemma C.4 and $b ( \lfloor b _ { N } \rfloor ) \psi _ { d } ( \frac { \lfloor b _ { N } \rfloor } { 3 } ) \to 0$ . 

Proposition 9. Under C1, C2, C9, and C10 and if $N \times \mathrm { A V a r } \left( \widehat { \tau } _ { \mathrm { H A } } ( d ) \right)$ is uniformly bounded below for all large $N$ , we have, 

$$
\frac {\widehat {\tau} _ {\mathrm {H A}} (d) - \operatorname {A M E} (d ; \eta)}{\sqrt {\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)}} \xrightarrow {d} N (0, 1). \tag {111}
$$

Proof. The proof of $\sqrt { N } \left( \widehat { \tau } _ { \mathrm { H A } } ( d ) - \mathrm { A M E } ( d ; \eta ) \right) = \sqrt { N } \left( \widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d ) - \mathrm { A M E } ( d ; \eta ) \right) + o _ { p } ( 1 )$ is similar as in Section A.3.1. We only need to establish the asymptotic distribution for $\widehat { \tau } _ { \mathrm { H A } } ^ { \mathrm { L } } ( d )$ . 

The proof follows from Theorem 2 in Jenish and Prucha (2012) (JP hereafter) with a slight modification of the proof. Assumption 1 in Jenish and Prucha (2012) is replaced by C10. Assumption 1 in JP is used to prove the covariance inequality in Lemma A.3 in JP, which can be replaced by the inequality in the proof of Lemma C.4. Assumption 1 is also used in established Theorem A.1 in JP for the $\alpha$ -mixing CLT. Because our assignment variables $\{ Z _ { i , N } , i \in S _ { N } \}$ are independent by C1 and by C10, we can replace the $\alpha -$ mixing CLT with a CLT with bounded-degree dependency graph. Assumption 2(a) and Assumption 4(a) in JP are trivially satisfied by C2 and choosing $\begin{array} { r } { c _ { i , n } = \frac { c } { N } } \end{array}$ for all $i$ by some constant $c > 0$ independent of $N$ . Assumption 2(b) and Assumption 3 in JP are trivially satisfied by C1. Assumption 4(b) is satisfied by the premise of this lemma. Assumption 4(c) is satisfied by let $\begin{array} { r } { d _ { i n } = \frac { d } { N } } \end{array}$ for all $i$ by some constant $d > 0$ independent of $N$ . Then the proposition follows. □ 

We study the asymptotic variance of the Hajek estimator under the current set of assumptions: 

Lemma C.6. Under C1, the variance of the linearized Hajek estimator can be expressed as 

$$
\begin{array}{l} \operatorname {V a r} \left(\widehat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d)\right) (112) \\ = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right] (113) \\ + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \left(\mu_ {i} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \bar {\mu} ^ {a} (d)\right) \left(\mu_ {j} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \bar {\mu} ^ {b} (d)\right) \right] (114) \\ \end{array}
$$

Proof. We have the expansion, similar to the derivation in Lemma A.4, 

$$
\begin{array}{l} \operatorname {V a r} \left(\widehat {\tau} _ {\mathrm {H A}} ^ {\mathrm {L}} (d)\right) (115) \\ = \frac {1}{N ^ {2} p} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {1}; d) - \bar {\mu} ^ {1} (d)\right) ^ {2} \right] + \frac {1}{N ^ {2} (1 - p)} \sum_ {i = 1} ^ {N} \mathrm {E} \left[ \left(\mu_ {i} (\mathbf {0}; d) - \bar {\mu} ^ {0} (d)\right) ^ {2} \right] (116) \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \mathrm {E} ^ {2} \left[ \mu_ {i} (\mathbf {1}; d) - \mu_ {i} (\mathbf {0}; d) - \left(\bar {\mu} ^ {1} (d) - \bar {\mu} ^ {0} (d)\right) \right] (117) \\ + \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} \left[ \left(\mu_ {i} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \bar {\mu} ^ {a} (d)\right) \left(\mu_ {j} \left(\mathbf {a} _ {\mathbf {i}}, \mathbf {b} _ {\mathbf {j}}; d\right) - \bar {\mu} ^ {b} (d)\right) \right] (118) \\ - \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \sum_ {a = 0} ^ {1} \sum_ {b = 0} ^ {1} (- 1) ^ {a + b} \mathrm {E} [ \mu_ {i} (\mathbf {a}; d) - \bar {\mu} ^ {a} (d) ] E [ \mu_ {j} (\mathbf {b}; d) - \bar {\mu} ^ {b} (d) ], (119) \\ = (1 1 6) + (1 1 8) + 0. (120) \\ \end{array}
$$

where we use the fact that $\begin{array} { r } { \sum _ { i = 1 } ^ { N } \operatorname { E } \left[ \mu _ { i } ( \mathbf { a } ; d ) - \bar { \mu } ^ { a } ( d ) \right] = 0 } \end{array}$ by definition for $a \in \{ 0 , 1 \}$ . □ 

Now we provide results on HAC variance estimation. We first establish several lemmas. 

Lemma C.7. Under C1, C2, C9, for any pair of i and $j$ and for any $\mu _ { 1 }$ and $\mu _ { 0 }$ such that $| \mu _ { 1 } | < B$ and $| \mu _ { 0 } | < B$ for some positive constant $B$ , we have 

$$
\left| \left| \widehat {\tau} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) - \mathrm {E} \left[ \widehat {\tau} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} \lesssim \left\{ \begin{array}{l l} \frac {1}{N} & \text {i f} s \leq d _ {i j} \\ \frac {1}{N} \psi_ {d} \left(s - d _ {i j}\right) & \text {i f} s > d _ {i j} \end{array} , \right. \tag {121}
$$

where by using the symbol ≲ we suppress a constant dependent of B but independent of $N$ 

Further, 

$$
\left. \left| \left| \widehat {\tau} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right) \widehat {\tau} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) - \operatorname {E} \left[ \widehat {\tau} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right) \widehat {\tau} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} \right.
$$

$$
\lesssim \left\{ \begin{array}{l l} \frac {1}{N ^ {2}} & \text {i f} s \leq d _ {i j} \\ \frac {1}{N ^ {2}} \psi_ {d} (s - d _ {i j}) & \text {i f} s > d _ {i j}, \end{array} \right. \tag {123}
$$

where by using the symbol $\lesssim$ we suppress a constant dependent of $B$ but independent of $N$ 

Proof. For brevity we write $\widehat { \tau } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ as $\widehat { \tau } _ { i , N } ( d )$ and $\widehat { \tau } _ { j } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ as $\widehat { \tau } _ { j } ( d )$ . 

For the case where $s \leq d _ { i j }$ the statements are proved by the triangle inequality and C2. For the case where $s > d _ { i j }$ , we have: 

$$
\left\| \widehat {\tau} _ {j, N} (d) - \operatorname {E} \left[ \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right\| _ {2} \tag {124}
$$

$$
\leq | | \widehat {\tau} _ {j, N} (d) - \operatorname {E} \left[ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {j} (s - d _ {i j}) \right] | | _ {2} + | | \operatorname {E} \left[ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {j} (s - d _ {i j}) \right] - \operatorname {E} \left[ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {i} (s) \right] | | _ {2} \quad (1 2 5)
$$

$$
\lesssim \frac {1}{N} \psi_ {d} (s - d _ {i j}) + | | \mathrm {E} [ \widehat {\tau} _ {j, N} (d) - \mathrm {E} [ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {j} (s - d _ {i j}) ] | \mathcal {F} _ {i} (s) ] | | _ {2} \tag {126}
$$

$$
\lesssim \frac {1}{N} \psi_ {d} (s - d _ {i j}) + | | \widehat {\tau} _ {j, N} (d) - \mathrm {E} [ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {j} (s - d _ {i j}) ] | | _ {2} \tag {127}
$$

$$
\lesssim \frac {1}{N} \psi_ {d} \left(s - d _ {i j}\right). \tag {128}
$$

(125) is by triangle inequality. (126) follows the estimates in (90) and the fact that ${ \mathcal { F } } _ { j } ( s -$ $d _ { i j } ) \subset F _ { i } ( s )$ . (127) follows by the Jensen’s inequality, and (128) follows by (90). 

$$
\begin{array}{l} \left. \left| \left| \widehat {\tau} _ {i, N} (d) \widehat {\tau} _ {j, N} (d) - \mathrm {E} \left[ \widehat {\tau} _ {i, N} (d) \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right] \right| \right| _ {2} (129) \\ \leq \left\| \widehat {\tau} _ {i, N} (d) \widehat {\tau} _ {j, N} (d) - \widehat {\tau} _ {i, N} (d) \mathrm {E} \left[ \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right\| _ {2} (130) \\ + \left\| \widehat {\tau} _ {i, N} (d) \mathrm {E} \left[ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {i} (s) \right] - \mathrm {E} \left[ \widehat {\tau} _ {i, N} (d) | \mathcal {F} _ {i} (s) \right] \mathrm {E} \left[ \widehat {\tau} _ {j, N} (d) | \mathcal {F} _ {i} (s) \right] \right\| _ {2} \\ + \left\| \operatorname {E} \left[ \widehat {\tau} _ {i, N} (d) \mid \mathcal {F} _ {i} (s) \right] \operatorname {E} \left[ \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] - \operatorname {E} \left[ \widehat {\tau} _ {i, N} (d) \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right\| _ {2} \\ \lesssim \frac {1}{N} \left| \left| \widehat {\tau} _ {j, N} (d) - \mathrm {E} \left[ \widehat {\tau} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} + \frac {1}{N} \left| \left| \widehat {\tau} _ {i, N} (d) - \mathrm {E} \left[ \widehat {\tau} _ {i, N} (d) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} (131) \\ + \frac {1}{N} | | \widehat {\tau} _ {i, N} (d) - \operatorname {E} \left[ \widehat {\tau} _ {i, N} (d) \mid \mathcal {F} _ {i} (s) \right] | | _ {2}, \\ \end{array}
$$

where by using the symbol $\lesssim$ we suppress a constant dependent of B but independent of $N$ (131) is bounded by, up to a constant independent of $N$ , $\begin{array} { r } { \frac { 1 } { N ^ { 2 } } \psi _ { d } ( s - d _ { i j } ) + \frac { 1 } { N ^ { 2 } } \psi _ { d } ( s ) } \end{array}$ if $s > d _ { i j }$ and $\textstyle { \frac { 1 } { N ^ { 2 } } }$ if $s \leq d _ { i j }$ . (130) is by triangle inequality. (131) is by C2. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/a1a8cef81b07f5a90468a1a432cd48f450e7bdb7c75843a9cb683b6b21c179b0.jpg)


Let $K : \mathbb { R } _ { + } \to \mathbb { R }$ be the kernel function used for the HAC variance estimator. We make the following assumption on $K$ . 

C 11. The kernel function $K : \mathbb { R } _ { + } \to \mathbb { R }$ has the following properties: 

1. $K ( 0 ) = 1$ 

2. It is supported on [0, 1]. 43 

3. It is uniformly bounded: there exists a $K _ { m a x } < \infty$ such that $\begin{array} { r } { \operatorname* { s u p } _ { x \in \mathbb { R } _ { + } } | K ( x ) | < K _ { m a x } } \end{array}$ 

4. It is locally Lipschitz at 0: there exists a $\epsilon \in ( 0 , 1 )$ and a positive constant $C > 0$ such that 

$$
| K (x) - 1 | \leq C | x |, \tag {132}
$$

for any $x \in [ 0 , \epsilon )$ . 

Lemmas below provide preliminary estimates for some random quantities. We define the following quantity: 

$$
\widehat {\epsilon} _ {i, N} (\mu_ {1}, \mu_ {0}; d) = \frac {\mathbf {Z} _ {i}}{N p} \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {1}\right) - \frac {1 - \mathbf {Z} _ {i}}{N (1 - p)} \left(\mu_ {i} (\mathbf {Y}; d) - \mu_ {0}\right). \tag {133}
$$

In terms of dependency, $\widehat { \epsilon } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ behaves similarly as $\widehat { \tau } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ . We summarize them in the lemma below: 

Lemma C.8. Under C1, $\it C 2$ , C9 and C10. For any pair i and $j$ and for any $\mu _ { 1 }$ and $\mu _ { 0 }$ such that $| \mu _ { 1 } | < B$ and $| \mu _ { 0 } | < B$ for some positive constant $B$ , 

$$
\left| \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right), \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right)\right) \right| \lesssim \frac {1}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) \tag {134}
$$

where by using the symbol ≲ we suppress a constant dependent of B but independent of $N$ . Moreover, 

$$
| | \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) - \mathrm {E} \left[ \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) \mid \mathcal {F} _ {i} (s) \right] | | _ {2} \lesssim \left\{ \begin{array}{l l} \frac {1}{N} & \text {i f} s \leq d _ {i j} \\ \frac {1}{N} \psi_ {d} \left(s - d _ {i j}\right) & \text {i f} s > d _ {i j} \end{array} , \right. \tag {135}
$$

where by using the symbol ≲ we suppress a constant dependent of B but independent of $N$ . Moreover, 

$$
\left. \left| \left| \widehat {\epsilon} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right) \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) - \operatorname {E} \left[ \widehat {\epsilon} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right) \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} \right.
$$

$$
\lesssim \left\{ \begin{array}{l l} \frac {1}{N ^ {2}} & \text {i f} s \leq d _ {i j} \\ \frac {1}{N ^ {2}} \psi_ {d} (s - d _ {i j}) & \text {i f} s > d _ {i j}, \end{array} \right. \tag {137}
$$

where by using the symbol ≲ we suppress a constant dependent of B but independent of $N$ . 

Proof. The first statement follows from Lemma C.3 upon noting that $\widehat { \epsilon } _ { i , N } ( \mu _ { 1 } . \mu _ { 0 } ; d )$ differs from $\widehat { \tau } _ { i , N } ( \mu _ { 1 } . \mu _ { 0 } ; d )$ by a constant. The rest statements follow the same calculation as in Lemma C.7. □ 

Define the quantity 

$$
\widehat {v} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d, b _ {N}\right) = \sum_ {j \in \mathcal {B} (i; b _ {N})} \widehat {\epsilon} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d\right) \widehat {\epsilon} _ {j, N} \left(\mu_ {1}, \mu_ {0}; d\right) K \left(\frac {d _ {i j}}{b _ {N}}\right), \tag {138}
$$

and $v _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d , b _ { N } ) = \operatorname { E } \left[ \widehat { v } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d , b _ { N } ) \right]$ 

Lemma C.9. Under C1, C2, C9, C10, and C11, for all $i \in S _ { N }$ and for any $\mu _ { 1 }$ and $\mu _ { 0 }$ such that $| \mu _ { 1 } | < B$ and $| \mu _ { 0 } | < B$ for some positive constant $B$ , we have, uniformly for all i and all sample sizes $N$ , 

(i) $\begin{array} { r } { | \widehat { v } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d , b _ { N } ) | \lesssim \frac { b _ { N } ^ { 2 } } { N ^ { 2 } } } \end{array}$ and $\begin{array} { r } { | v _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d , b _ { N } ) | \lesssim \frac { b _ { N } ^ { 2 } } { N ^ { 2 } } } \end{array}$ , 

(ii) If $s > b _ { N }$ 

$$
\left. \right.\left|\left| \widehat {v} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d, b _ {N}\right) - \mathrm {E} \left[ \widehat {v} _ {i, N} \left(\mu_ {1}, \mu_ {0}; d, b _ {N}\right) \mid \mathcal {F} _ {i} (s) \right]\right|\right| _ {2} \lesssim \frac {1}{N ^ {2}} \sum_ {j \in \mathcal {B} (i; b _ {N})} \psi_ {d} (s - d _ {i j}), \tag {139}
$$

(iii) If $s > 3 b _ { N }$ 

$$
\left| \operatorname {C o v} \left(\widehat {v} _ {i, N} (d, b _ {N}), \widehat {v} _ {j, N} (d, b _ {N})\right) \right| \lesssim \frac {b _ {N} ^ {2}}{N ^ {4}} \times \max  \left\{\sum_ {\{k: d _ {i k} \leq b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3} - d _ {i k}\right), \sum_ {\{l: d _ {j l} \leq b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3} - d _ {j l}\right) \right\}, \tag {140}
$$

where by using the symbol ≲ we suppress a constant dependent of $B$ but independent of $N$ . 

Proof. (i) follows by C1, C2, C10 and C11 and an application of the triangle inequality. We write $\widehat { v } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d , b _ { N } )$ as ${ \widehat { v _ { i , N } } } ( d , b _ { N } )$ and $\widehat { \epsilon } _ { i , N } ( \mu _ { 1 } , \mu _ { 0 } ; d )$ as $\widehat { \epsilon } _ { i , N } ( d )$ for brevity. All the 

constants omitted below depend on $B$ 

$$
\begin{array}{l} \left| \left| \widehat {v} _ {i, N} (d, b _ {N}) - \mathrm {E} \left[ \widehat {v} _ {i, N} (d, b _ {N}) \mid \mathcal {F} _ {i} (s) \right] \right| \right| _ {2} (141) \\ \leq K _ {\max } \times \sum_ {j \in \mathcal {B} (i; b _ {N})} | | \widehat {\epsilon} _ {i, N} (d) \widehat {\epsilon} _ {j, N} (d) - E [ \widehat {\epsilon} _ {i, N} (d) \widehat {\epsilon} _ {j, N} (d) \mid \mathcal {F} _ {i} (s) ] | | _ {2} (142) \\ \lesssim \frac {K _ {\max }}{N ^ {2}} \left(| j: d _ {i j} \geq s, d _ {i j} \leq b _ {N} | + \sum_ {\{j: d _ {i j} \leq s \}} \psi_ {d} (s - d _ {i j})\right) (143) \\ \end{array}
$$

(142) is by the triangle inequality and C11. (143) is by the estimates in Lemma C.8 and C2. If $s > b _ { N }$ , the term $| j : d _ { i j } \geq s , d _ { i j } \leq b _ { N } |$ vanishes, proving the result. 

As in Lemma C.1, we have, if $s > b _ { N }$ 

$$
\left| \left| \operatorname {E} \left[ \widehat {v} _ {i, N} (d, b _ {N}) - v _ {i, N} (d, b _ {N}) \mid \mathcal {F} _ {i} ^ {c} (s) \right] \right| \right| _ {2} \lesssim \frac {1}{N ^ {2}} \sum_ {j \in \mathcal {B} (i; b _ {N})} \psi_ {d} (s - d _ {i j}). \tag {144}
$$

Similarly as in Lemma C.3, if $d _ { i j } > 3 b _ { N }$ 

$$
\begin{array}{l} \left| \operatorname {C o v} \left(\widehat {v} _ {i, N} (d, b _ {N}), \widehat {v} _ {j, N} (d, b _ {N})\right) \right| \\ \leq | \operatorname {E} \left[ \left(\widehat {v} _ {i, N} (d, b _ {N}) - v _ {i, N} (d, b _ {N})\right) \left(\widehat {v} _ {j, N} (d, b _ {N}) - v _ {j, N} (d, b _ {N})\right) \right] | \\ \leq \left| \operatorname {E} \left[ \left(\widehat {v} _ {i, N} (d, b _ {N}) - v _ {i, N} (d, b _ {N}) - \operatorname {E} \left[ \left(\widehat {v} _ {i, N} (d, b _ {N}) - v _ {i, N} (d, b _ {N})\right) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right) \right]\right) \left(\widehat {v} _ {j, N} (d, b _ {N}) - v _ {j, N} (d, b _ {N})\right) \right] \right| \\ + \left| \right. \operatorname {E} \left[ \operatorname {E} \left[\left(\widehat {v} _ {i, N} (d, b _ {N}) - v _ {i, N} (d, b _ {N})\right) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right)\right] E \left[\left(\widehat {v} _ {j, N} (d, b _ {N}) - v _ {j, N} (d, b _ {N})\right) \mid \mathcal {F} _ {i} \left(\frac {d _ {i j}}{3}\right)\right]\right] \\ \end{array}
$$

$$
\lesssim \frac {b _ {N} ^ {2}}{N ^ {2}} \times \frac {1}{N ^ {2}} \max  \{\sum_ {\{k: d _ {i k} \leq b _ {N} \}} \psi_ {d} (\frac {d _ {i j}}{3} - d _ {i k}), \sum_ {\{l: d _ {j l} \leq b _ {N} \}} \psi_ {d} (\frac {d _ {i j}}{3} - d _ {j l}) \}.
$$

Recall the definition of the HAC variance estimator from (18). An calculation similar to that at the beginning of Section A.4.1 shows that the HAC variance estimator $\widehat { \mathrm { V } } _ { \mathrm { H A C } } ( d )$ with 

a kernel function $K$ can be expressed as: 

$$
\begin{array}{l} \widehat {\mathrm {V}} _ {\text {H A C}} (d; b _ {N}) (145) \\ = \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d) + \frac {1}{N _ {0} ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d) + \frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} Z _ {i} Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) (146) \\ - \frac {2}{N _ {1} N _ {0}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} (1 - Z _ {i}) Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) + \frac {1}{N _ {0} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} (1 - Z _ {i}) (1 - Z _ {j}) \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) (147) \\ \end{array}
$$

Define the following quantities 

$$
\begin{array}{l} \widehat {\mathrm {V}} _ {1} (d; b _ {N}) = \frac {1}{N ^ {2} p ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d) + \frac {1}{N ^ {2} (1 - p) ^ {2}} \sum_ {i = 1} ^ {N} (1 - Z _ {i}) \hat {e} _ {i} ^ {2} (d) (148) \\ + \frac {1}{N ^ {2} p ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} Z _ {i} Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) (149) \\ - \frac {2}{N ^ {2} p (1 - p)} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} \left(1 - Z _ {i}\right) Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) (150) \\ + \frac {1}{N ^ {2} (1 - p) ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} (1 - Z _ {i}) (1 - Z _ {j}) \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) (151) \\ \end{array}
$$

Note $\begin{array} { r } { \widehat { \mathrm { V } } _ { 1 } ( d ; b _ { N } ) = \sum _ { i = 1 } ^ { N } v _ { i } ( \widehat { \bar { \mu } } _ { 1 } ( d ) , \widehat { \bar { \mu } } _ { 0 } ( d ) ; d , b _ { N } ) } \end{array}$ . Also define 

$$
\widehat {\mathrm {V}} _ {2} (d; b _ {N}) = \sum_ {i = 1} ^ {N} v _ {i} \left(\bar {\mu} _ {1} (d), \bar {\mu} _ {0} (d); d, b _ {N}\right) \tag {152}
$$

The following assumption is similar as C5. 

C 12. Given a kernel function K and a sequence of bandwidth $\{ b _ { N } \} _ { N = 1 } ^ { \infty }$ 

$$
\lim  _ {N} \inf  _ {N} \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \left(\tau_ {i} (d) - \operatorname {A M E} (d; \eta)\right) \left(\tau_ {j} (d) - \operatorname {A M E} (d; \eta)\right)\right) \geq 0. ^ {4 4} \tag {153}
$$

Remark C.1. The following lemma shows that some limitation on the correlation of effect heterogeneity when nodes are far apart implies C12. The lemma is given with a high-level condition. It can be checked when one imposes specific low-level conditions on the majorizing function (i.e., $\phi _ { d }$ function below.) and on the spacing of intervention nodes. 

Lemma C.10. Under C2, C10 and C11, and suppose for all large N and any pair i and $j$ we have 

$$
\left| \tau_ {i} (d) - \operatorname {A M E} (d; \eta) \right| \times \left| \tau_ {j} (d) - \operatorname {A M E} (d; \eta) \right| \leq \phi_ {d} \left(d _ {i j}\right), \tag {154}
$$

and $\begin{array} { r } { \operatorname* { l i m } _ { N } \operatorname* { s u p } _ { i \in \mathcal { S } _ { N } } \sum _ { j \in \mathcal { B } ( i , b _ { N } ) } \phi _ { d } ( d _ { i j } ) \to 0 } \end{array}$ for any $b _ { N } \to \infty$ 

We have 

$$
\lim _ {N} \inf _ {N} \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K (\frac {d _ {i j}}{b _ {N}}) (\tau_ {i} (d) - \mathrm {A M E} (d; \eta)) (\tau_ {j} (d) - \mathrm {A M E} (d; \eta))\right) = 0.
$$

Proof. Note we have $\begin{array} { r } { \sum _ { i = 1 } ^ { N } \sum _ { j = 1 } ^ { N } \left( \tau _ { i } ( d ) - \mathrm { A M E } ( d ; \eta ) \right) \left( \tau _ { j } ( d ) - \mathrm { A M E } ( d ; \eta ) \right) = 0 } \end{array}$ by definition. 

Let $\tilde { b } _ { N } = o ( b _ { N } ^ { \frac 1 3 } )$ . Then, for large $N$ , 

$$
\begin{array}{l} \lim  _ {N} \left| \left(\frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \left(\tau_ {i} (d) - \operatorname {A M E} (d; \eta)\right) \left(\tau_ {j} (d) - \operatorname {A M E} (d; \eta)\right)\right) \right| \\ = \lim  _ {N} \left| \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \left(K \left(\frac {d _ {i j}}{b _ {N}}\right) - 1\right) \left(\tau_ {i} (d) - \operatorname {A M E} (d; \eta)\right) \left(\tau_ {j} (d) - \operatorname {A M E} (d; \eta)\right) \right| \\ \lesssim \lim _ {N} \left| \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i, \tilde {b} _ {N})} \frac {d _ {i j}}{b _ {N}} \right| + \left| \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \notin \mathcal {B} (i, \tilde {b} _ {N})} | (\tau_ {i} (d) - \mathrm {A M E} (d; \eta)) (\tau_ {j} (d) - \mathrm {A M E} (d; \eta)) | \right| \\ \lesssim \lim _ {N} \frac {\tilde {b} _ {N} ^ {3}}{b _ {N}} + \lim _ {N} \sup _ {i \in \mathcal {S} _ {N}} \sum_ {j \notin \mathcal {B} (i, \tilde {b} _ {N})} | (\tau_ {i} (d) - \mathrm {A M E} (d; \eta)) (\tau_ {j} (d) - \mathrm {A M E} (d; \eta)) | = o (1). \\ \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/ee403ebffc2fb45055d922fa6f64c1f9772ee9d0985f4b32146ff2413b3b6547.jpg)


We conclude the section with a proof of the consistency of the HAC variance estimator and establish the asymptotic validity of the Wald confidence interval. 

Proposition 10. Let $b _ { N } = o ( N ^ { \frac { 1 } { 6 } } )$ and $b _ { N } \to \infty$ . For $\alpha < 1$ and under C1, C2, C9, C10, C11 and C12, 

(i) $\begin{array} { r } { N \times \left( \widehat { \mathrm { V } } _ { \mathrm { H A C } } ( d ; b _ { N } ) - \widehat { \mathrm { V } } _ { 1 } ( d , b _ { N } ) \right) = o _ { p } ( 1 ) . } \end{array}$ 

(ii) $\begin{array} { r } { N \times \left( \widehat { \mathrm { V } } _ { 1 } ( d ; b _ { N } ) - \widehat { \mathrm { V } } _ { 2 } ( d , b _ { N } ) \right) = o _ { p } ( 1 ) . } \end{array}$ 

(iii) $\begin{array} { r } { N \times \left( \widehat { \mathrm { V } } _ { 2 } ( d ; b _ { N } ) - \mathrm { E } \left[ \widehat { \mathrm { V } } _ { 2 } ( d ; b _ { N } ) \right] \right) = o _ { p } ( 1 ) . } \end{array}$ 

(iv) $\begin{array} { r } { \operatorname* { l i m } \operatorname* { i n f } _ { N } \Big ( N \times \Big ( \mathrm { E } \left[ \widehat { \mathrm { V } } _ { 2 } ( d ; b _ { N } ) \right] - \mathrm { A V a r } \left( \widehat { \tau } _ { \mathrm { H A } } ( d ) \right) \Big ) \Big ) \geq 0 . } \end{array}$ 

$\begin{array} { r } { ( v ) \ \operatorname* { l i m } _ { N \to \infty } \mathbf { P r o b } \left( z _ { 2 } ^ { } \leq \frac { \widehat \gamma _ { \mathrm { H A } } ( d ) - \mathrm { A M E } ( d ; \eta ) } { \sqrt { \widehat \mathrm { V } _ { \mathrm { H A C } } ( d , b _ { N } ) } } \leq z _ { 1 - \frac { \alpha } { 2 } } \right) \geq 1 - \alpha } \end{array}$ ≤ z1− 

Proof. We first prove (i). Notice by C1, $\begin{array} { r } { \frac { N _ { 1 } } { N } - p = O _ { p } ( N ^ { - \frac { 1 } { 2 } } ) } \end{array}$ . We have 

$$
\left| \frac {N ^ {2}}{N _ {1} ^ {2}} - \frac {1}{p ^ {2}} \right| = \left| \frac {p ^ {2} - \left(\frac {N _ {1}}{N}\right) ^ {2}}{p ^ {2} \left(\frac {N _ {1}}{N}\right) ^ {2}} \right| = \left| \frac {\left(\frac {N _ {1}}{N} - p\right) \left(\frac {N _ {1}}{N} + p\right)}{p ^ {2} \left(\frac {N _ {1}}{N}\right) ^ {2}} \right| = O _ {p} \left(N ^ {- \frac {1}{2}}\right) \tag {155}
$$

We have 

$$
\begin{array}{l} N \times \left(\frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d) - \frac {1}{N ^ {2} p ^ {2}} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d)\right) (156) \\ = \left(\frac {N ^ {2}}{N _ {1} ^ {2}} - \frac {1}{p ^ {2}}\right) \frac {1}{N} \sum_ {i = 1} ^ {N} Z _ {i} \hat {e} _ {i} ^ {2} (d) = O _ {p} \left(N ^ {- \frac {1}{2}}\right), (157) \\ \end{array}
$$

where we use the fact that $\begin{array} { r } { \frac { 1 } { N } \sum _ { i = 1 } ^ { N } Z _ { i } \hat { e } _ { i } ^ { 2 } ( d ) = O ( 1 ) } \end{array}$ by C2.45 

Also, 

$$
\begin{array}{l} \left| N \times \left(\frac {1}{N _ {1} ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} Z _ {i} Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) - \frac {1}{N ^ {2} p ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} Z _ {i} Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right)\right) \right| (158) \\ \leq \left(\frac {N ^ {2}}{N _ {1} ^ {2}} - \frac {1}{p ^ {2}}\right) \left| \frac {1}{N} \sum_ {i = 1} ^ {N} \sum_ {j \neq i} Z _ {i} Z _ {j} \hat {e} _ {i} (d) \hat {e} _ {j} (d) K \left(\frac {d _ {i j}}{b _ {N}}\right) \right| (159) \\ \lesssim \left(\frac {N ^ {2}}{N _ {1} ^ {2}} - \frac {1}{p ^ {2}}\right) \left(b _ {N} ^ {2} \times K _ {\max }\right) = O _ {p} \left(N ^ {- \frac {1}{2}} b _ {N} ^ {2}\right) = o _ {p} (1) (160) \\ \end{array}
$$

where in (160) we use C2, C10 and C11. Other terms can be processed similarly. This proves (i). 

Now we prove (ii). We have: 

$$
\begin{array}{l} \widehat {\mathrm {V}} _ {1} (d; b _ {N}) - \widehat {\mathrm {V}} _ {2} (d; b _ {N}) (161) \\ = \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; b _ {N})} K \left(\frac {d _ {i j}}{b _ {N}}\right) \times (162) \\ \underbrace {\left(\widehat {\epsilon} _ {i , N} \left(\widehat {\mu} ^ {1} (d) , \widehat {\mu} ^ {0} (d) ; d\right) \widehat {\epsilon} _ {j , N} \left(\widehat {\mu} ^ {1} (d) , \widehat {\mu} ^ {0} (d) ; d\right) - \widehat {\epsilon} _ {i , N} \left(\bar {\mu} ^ {1} (d) , \bar {\mu} ^ {0} (d) ; d\right) \widehat {\epsilon} _ {j , N} \left(\bar {\mu} ^ {1} (d) , \bar {\mu} ^ {0} (d) ; d\right)\right)} _ {(u _ {i j})} (163) \\ \end{array}
$$

For $( u _ { i j } )$ we have the expressions 

$$
\begin{array}{l} u _ {i j} = \frac {Z _ {i} Z _ {j}}{N ^ {2} p ^ {2}} \left(- \mu_ {i} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {1} (d)\right) - \mu_ {j} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {1} (d)\right) + \left(\widehat {\bar {\mu}} ^ {1} (d)\right) ^ {2} - \left(\bar {\mu} ^ {1} (d)\right) ^ {2}\right) (164) \\ + \frac {\left(1 - Z _ {i}\right) \left(1 - Z _ {j}\right)}{N ^ {2} (1 - p) ^ {2}} \left(- \mu_ {i} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {0} (d) - \bar {\mu} ^ {0} (d)\right) - \mu_ {j} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {0} (d) - \bar {\mu} ^ {0} (d)\right) + \left(\widehat {\bar {\mu}} ^ {0} (d)\right) ^ {2} - \left(\bar {\mu} ^ {0} (d)\right) ^ {2}\right) (165) \\ - \frac {\left(1 - Z _ {i}\right) Z _ {j}}{N ^ {2} (1 - p) ^ {2}} \left(- \mu_ {i} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {1} (d)\right) - \mu_ {j} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {0} (d) - \bar {\mu} ^ {0} (d)\right) + \widehat {\bar {\mu}} ^ {0} (d) \widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {0} (d) \bar {\mu} ^ {1} (d)\right) (166) \\ - \frac {\left(1 - Z _ {j}\right) Z _ {i}}{N ^ {2} (1 - p) ^ {2}} \left(- \mu_ {i} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {0} (d) - \bar {\mu} ^ {0} (d)\right) - \mu_ {j} (\mathbf {Y}; d) \left(\widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {1} (d)\right) + \widehat {\bar {\mu}} ^ {0} (d) \widehat {\bar {\mu}} ^ {1} (d) - \bar {\mu} ^ {0} (d) \bar {\mu} ^ {1} (d)\right) (167) \\ \end{array}
$$

Some calculation will show that by C1, C2, C9, C10, and C11 we have 

$$
\widehat {\mathrm {V}} _ {1} (d; b _ {N}) - \widehat {\mathrm {V}} _ {2} (d; b _ {N}) = O _ {p} \left(\frac {1}{N} b _ {N} ^ {2} N ^ {- \frac {1}{2}}\right). \tag {168}
$$

We then have $N \times ( \widehat { \mathrm { V } } _ { 1 } ( d ; b _ { N } ) - \widehat { \mathrm { V } } _ { 2 } ( d ; b _ { N } ) ) = o _ { p } ( 1 )$ because $b _ { N } = o ( N ^ { \frac { 1 } { 6 } } )$ . This proves (ii). 

Now we prove (iii). Write for brevity $\widehat { v } _ { i , N } ( \bar { \mu } ^ { 1 } ( d ) , \bar { \mu } ^ { 0 } ( d ) ; d , b _ { N } ) = \widehat { v } _ { i , N } ( d ; b _ { N } )$ . Note that 

$\bar { \mu } ^ { 1 } ( d )$ and $\bar { \mu } ^ { 0 } ( d )$ are uniformly bounded by C2 for all sample sizes $N$ . 

$$
\begin{array}{l} \operatorname {V a r} \left(\widehat {\mathrm {V}} _ {2} (d; b _ {N})\right) = \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \operatorname {C o v} \left(\widehat {v} _ {i, N} (d; b _ {N}), \widehat {v} _ {j, N} (d; b _ {N})\right) (169) \\ \leq \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; 6 b _ {N})} | \operatorname {C o v} \left(\widehat {v} _ {i, N} (d; b _ {N}), \widehat {v} _ {j, N} (d; b _ {N})\right) | + \sum_ {i = 1} ^ {N} \sum_ {j \notin \mathcal {B} (i; 6 b _ {N})} | \operatorname {C o v} \left(\widehat {v} _ {i, N} (d; b _ {N}), \widehat {v} _ {j, N} (d; b _ {N})\right) | (170) \\ \lesssim \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; 6 b _ {N})} \sqrt {\operatorname {V a r} \left(\widehat {v} _ {i , N} (d ; b _ {N})\right) \times \operatorname {V a r} \left(\widehat {v} _ {j , N} (d ; b _ {N})\right)} (171) \\ + \sum_ {i = 1} ^ {N} \sum_ {j \notin \mathcal {B} (i; 6 b _ {N})} \frac {b _ {N} ^ {2}}{N ^ {4}} \times \max  \left\{\sum_ {\{k: d _ {i k} \leq b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3} - d _ {i k}\right), \sum_ {\{l: d _ {j l} \leq b _ {N} \}} \psi_ {d} \left(\frac {d _ {i j}}{3} - d _ {j l}\right) \right\} (172) \\ \lesssim N \times 3 6 b _ {N} ^ {2} \times \frac {b _ {N} ^ {4}}{N ^ {4}} + \sum_ {i = 1} ^ {N} \sum_ {j \notin \mathcal {B} (i; 6 b _ {N})} \frac {b _ {N} ^ {2}}{N ^ {4}} \times C _ {b} b _ {N} ^ {2} \times \frac {1}{b _ {N} ^ {4 + \epsilon}} (173) \\ \lesssim \frac {b _ {N} ^ {6}}{N ^ {3}} + \frac {1}{N ^ {2}} \frac {1}{b _ {N} ^ {\epsilon}} (174) \\ \end{array}
$$

Thus we have $N ^ { 2 } \times \mathrm { V a r } \left( \widehat { \mathrm { V } } _ { 2 } ( d , b _ { N } ) \right) = o ( 1 )$ because $b _ { N } = o ( N ^ { \frac { 1 } { 6 } } )$ 

Now we prove (iv). Write for brevity $\widehat { \epsilon } _ { i , N } ( \bar { \mu } ^ { 1 } ( d ) , \bar { \mu } ^ { 0 } ( d ) ; d ) = \widehat { \epsilon } _ { i , N } ( d )$ .46 First note we have: 

$$
\operatorname {E} \left[ \widehat {\mathrm {V}} _ {2} (d; b _ {N}) \right] = \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; b _ {N})} K \left(\frac {d _ {i j}}{b _ {N}}\right) \operatorname {E} \left[ \widehat {\epsilon} _ {i, N} (d) \widehat {\epsilon} _ {j, N} (d) \right], \tag {175}
$$

and, 

$$
\operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right) = \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right). \tag {176}
$$

We have the identity: 

$$
\begin{array}{l} \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \operatorname {E} \left[ \widehat {\epsilon} _ {i, N} (d) \widehat {\epsilon} _ {j, N} (d) \right] - \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) (177) \\ = \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) E [ \widehat {\epsilon} _ {i, N} (d) ] E [ \widehat {\epsilon} _ {j, N} (d) ] (178) \\ = \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \left(\tau_ {i} (d) - \operatorname {A M E} (d; \eta)\right) \left(\tau_ {j} (d) - \operatorname {A M E} (d; \eta)\right). (179) \\ \end{array}
$$

Note by Lemma C.8, 

$$
\left| \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) \right| \lesssim \frac {1}{N ^ {2}} \psi_ {d} \left(\frac {d _ {i j}}{3}\right). \tag {180}
$$

Let $\tilde { b } _ { N } = o ( b _ { N } ^ { \frac 1 3 } )$ . For large $N$ , we then have: 

$$
\begin{array}{l} \left| \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) - \operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right) \right| (181) \\ \leq \sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} \left| K \left(\frac {d _ {i j}}{b _ {N}}\right) - 1 \right| \times | \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) | (182) \\ \leq \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; \bar {b} _ {N})} \left| K \left(\frac {d _ {i j}}{b _ {N}}\right) - 1 \right| \times | \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) | (183) \\ + \left(K _ {\max } + 1\right) \sum_ {i = 1} ^ {N} \sum_ {j \notin \mathcal {B} (i; \bar {b} _ {N})} | \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) | (184) \\ \lesssim \frac {1}{N ^ {2}} \sum_ {i = 1} ^ {N} \sum_ {j \in \mathcal {B} (i; \tilde {b} _ {N})} \frac {d _ {i j}}{b _ {N}} + \frac {1}{N} \sup  _ {i \in \mathcal {S} _ {N}} \sum_ {j \notin \mathcal {B} (i; \tilde {b} _ {N})} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) (185) \\ \lesssim \frac {1}{N} \frac {\tilde {b} _ {N} ^ {3}}{b _ {N}} + \frac {1}{N} \sup  _ {i \in \mathcal {S} _ {N}} \sum_ {j \notin \mathcal {B} (i; \tilde {b} _ {N})} \psi_ {d} \left(\frac {d _ {i j}}{3}\right) = o \left(\frac {1}{N}\right) (186) \\ \end{array}
$$

(185) is by C11. (186) is by C10 and Lemma C.5. 

Thus we have $\begin{array} { r } { N \times \left( \sum _ { i = 1 } ^ { N } \sum _ { j = 1 } ^ { N } K ( \frac { d _ { i j } } { b _ { N } } ) \mathrm { C o v } \left( \widehat { \epsilon } _ { i , N } ( d ) , \widehat { \epsilon } _ { j , N } ( d ) \right) - \mathrm { A V a r } \left( \widehat { \tau } _ { \mathrm { H A } } ( d ) \right) \right) = o ( 1 ) } \end{array}$ . We 

have: 

$$
\left. \lim  _ {N} \inf  _ {N} \left(N \times \left(\mathrm {E} \left[ \widehat {\mathrm {V}} _ {2} (d, b _ {N}) \right] - \operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)\right)\right) \right. \tag {187}
$$

$$
\geq \lim  _ {N} \inf  _ {N} N \times \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) E [ \widehat {\epsilon} _ {i, N} (d) ] E [ \widehat {\epsilon} _ {j, N} (d) ]\right) \tag {188}
$$

$$
+ \lim  _ {N} \inf  _ {N} N \times \left(\sum_ {i = 1} ^ {N} \sum_ {j = 1} ^ {N} K \left(\frac {d _ {i j}}{b _ {N}}\right) \operatorname {C o v} \left(\widehat {\epsilon} _ {i, N} (d), \widehat {\epsilon} _ {j, N} (d)\right) - \operatorname {A V a r} \left(\widehat {\tau} _ {\mathrm {H A}} (d)\right)\right) \geq 0, \tag {189}
$$

where (189) is by (179) and C12, and the calculation from (181) to (186). 

This proves (iv). (v) follows similarly as in Proposition 5. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-16/0550a198-f6d4-4d7a-885d-98f44a791d2a/6a2a096c37b76be2edf2e56876e9046914e14219f624d7137436386e610b3030.jpg)


# D Notations

# D.1 Notation

<table><tr><td>Notations</td><td>Definitions</td><td>First Appear in</td></tr><tr><td>S</td><td>The set of intervention nodes</td><td>Page 6</td></tr><tr><td>X</td><td>Two-dimensional set locations for outcomes</td><td>Page 6</td></tr><tr><td>Z</td><td>Ordered vector of experimental assignment variable</td><td>Page 6</td></tr><tr><td>z</td><td>Realized assignment</td><td>Page 6</td></tr><tr><td>Yx(z)</td><td>Potential outcome when assignment is z at location x</td><td>Page 6</td></tr><tr><td>Yx</td><td>Observed outcome at location x</td><td>Page 6</td></tr><tr><td>Y(z)</td><td>Full set of potential outcomes when assignment is z</td><td>Page 6</td></tr><tr><td>Y</td><td>Full set of observed outcomes</td><td>Page 6</td></tr><tr><td>μi(Y(z); Ωd)</td><td>circle averages with Y(z) at distance d</td><td>Page 7</td></tr><tr><td>μi(Y; Ωd)</td><td>observed circle averages with Y(z) at distance d</td><td>Page 7</td></tr><tr><td>Yx(zi; η)</td><td>Marginalized potential outcome at point x</td><td>Page 11</td></tr><tr><td>μi(zi; d, η)</td><td>Marginalized circle averages</td><td>Page 11</td></tr><tr><td>τix(η)</td><td>Individual marginalized effect</td><td>Page 11</td></tr><tr><td>τi(d; η)</td><td>average of individual marginalized effect at distance d</td><td>Page 12</td></tr></table>