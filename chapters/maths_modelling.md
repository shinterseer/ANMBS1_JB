# Mathematical modelling

The mathematical modelling of the problem will follow {cite}`A7`.

Moisture flow is the sum of liquid and vapor flow:

$$
g = g_{v} + g_{l}
$$

with

$$
g_{v} = -\delta_{p} \frac{\partial p_{v}} {\partial x}
$$

and

$$
g_{v} = K_{w} \frac{\partial P_{suc}} {\partial x}
$$

with $P_{suc}$ denoting the suction pressure
With $w$ in $kg/m^3$ to denote the total moisture content we get the governing equation of the problem:

$$
\frac {\partial w}{\partial t} = \dot w =  - \frac {\partial g_{v}}{\partial x}  - \frac {\partial g_{l}}{\partial x}
$$

## Discretisation scheme
The domain is split into control volumes (also referred to as cells in this document), 
that are assumed to have a constant moisture level. 
The discretization scheme is depicted in {numref}`Figure {number} <discretization-fig>`.


```{figure} ./discretization.jpg
---
name: discretization-fig
width: 500px
align: center
---
Discretization of the domain: Capital letters denote center points of cells (W...West, P...Point, E...East), 
small letters denote midpoints between cells. (source: {cite}`A7`)
```


## Approximation of spatial derivatives

The operator $\frac {\partial}{\partial x} $ has to be apprximated twice:
* once in $\dot w = - \frac{\partial g}{\partial x}$ with $g=g_{v}+g_{l}$
* and once in $g_{v} = - \delta_{p} \frac{\partial p}{\partial x}$ and 
$g_{l} = - K \frac{\partial P_{suc}}{\partial x}$

The idea is to do the first one as the centered derivative between the points $e$ and $w$ and the second one as one sided derivatives between the points $W$ and $P$ and $P$ and $E$ respectively.

Starting with the first one:

$$
\frac{\partial g}{\partial x} \approx \frac {g_{PE} - g_{WP}} {d_{we}}
$$

this leads to:

$$
\begin{eqnarray}
\dot w &=& -\frac{\partial g}{\partial x} = \frac{\partial g-{v}}{\partial x} - \frac{\partial g_{l}}{\partial x}\\
&\approx& -\frac{g_{v,PE} - g_{v,WP}}{d_{we}} -\frac{g_{l,PE} - g_{l,WP}}{d_{we}} \\
&=& - \frac{1}{d_{we}}(g_{v,PE} - g_{v,WP} + g_{l,PE} - g_{l,WP})
\end{eqnarray}
$$