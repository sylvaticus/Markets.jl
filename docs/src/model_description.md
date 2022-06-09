# Model description 

## Explicit regional trade version

Current status: This model returns a "Too few degrees of freedom" error
Todo: (a) make it working :-); (b) add transport costs, policies, exogenous region and exogenous independent variables in the supply/demand functions; (c) estimate parameters (stochastic gradient descent with norm-1 - à la Lasso - penalty for sparce coefficients ??)

### Sets
- | $p$: products (with alias $pin$)
- | $r$: regions (with aliases $rfrom$ and $rto$)

### Variables
- | $sc_{p,r}$: Composite supply (quantity)
- | $dc_{p,r}$: Composite demand (quantity)
- | $pcs_{p,r}$: Composite supply price
- | $pcd_{p,r}$: Composite demand price
- | $s_{p,rfrom,rto}$: Supply of product $p$ from region $rfrom$ to region $rto$
- | $d_{p,rfrom,rto}$: Demmand of product $p$ coming from region $rfrom$ in region $rto$
- | $pa_{p,rfrom,rto}$: Actual market price for product $p$ coming from region $rin$ and going to region $rto$

### Parameters
- | $ϵsI_{p,r,pin}$: elasticity of the (composite) supply of $p$ with respect to the (composite) price of input product $pin$ in region $r$
- | $ϵsO_{p,r,pout}$: elasticity of the (composite) supply of $p$ with respect to the (composite) price of output product $pout$ in region $r$
- | $ϵdI_{p,r,pin}$: elasticity of the (composite) demand of $p$ with respect to the (composite) price of input product $pin$ in region $r$
- | $ϵdO_{p,r,pout}$: elasticity of the (composite) demand of $p$ with respect to the (composite) price of output product $pout$ in region $r$
- | $consts_{p,r}$: constant coefficient of the composite supply function for product $p$ in region $r$
- | $constd_{p,r}$: constant coefficient of the composite supply function for product $p$ in region $r$
- | $a_{p,r,rto}$: share parameter of the supply of product $p$ from $r$ to $rto$
- | $b_{p,rfrom,r}$: share parameter of the demand of product p in region $r$ coming from region $rfrom$
- | $ϵts_p$: trade elasticity of the supply of product $p$ 
- | $ϵtd_p$: trade elasticity of the demand of product $p$

### Constraints (equations)
- **Composite supply [p,r]** : $sc_{p,r} = consts_{p,r} * \prod_{pI} pcd_{pI,r}^{ϵsI_{p,r,pI}} * \prod_{pO} pcs_{pO,r}^{ϵsO_{p,r,pO}}$
- **Composite demand [p,r]** : $dc_{p,r} = constd_{p,r} * \prod_{pI} pcd_{pI,r}^{ϵdI_{p,r,pI}} * \prod_{pO} pcs_{pO,r}^{ϵdO_{p,r,pO}}$
- **Composite supply aggregation [p,r]** : $sc_{p,r} =  (\sum_{rto} a_{p,r,rto}* s_{p,r,rto}^{ϵts_p})^{1/ϵts_p}$
- **Composite demand aggregation [p,r]** : $dc_{p,r} =  (\sum_{rfrom} b_{p,rfrom,r}* d_{p,rfrom,r}^{ϵtd_p})^{1/ϵtd_p}$
- **Supply by destination [p,r]** : $s_{p,r,rto} =  a_{p,r,rto} * \frac{pcs_{p,r}}{pa_{p,r,rto}}^\frac{1}{1-ϵts_p}$
- **Demand from origin [p,r]** : $d_{p,rfrom,r} =  b_{p,rfrom,r} * \frac{pcd_{p,r}}{pa_{p,rfrom,r}}^\frac{1}{1-ϵtd_p}$
- **Supply monetary budget [p,r]** : $pcs_{p,r} * sc_{p,r}  =  \sum_{rto} pa_{p,r,rto} * s_{p,r,rto}$
- **Demand monetary budget [p,r]** : $pcd_{p,r} * dc_{p,r}  =  \sum_{rfrom} pa_{p,rfrom,r} * s_{p,rfrom,r}$
- **Physical balance [p, r]** :  $\sum_{rto} s_{p,r,rto} =  \sum_{rfrom} d_{p,rfrom,r}$

### Objective
$\min \sum_{p,rfrom,rto} (d_{p,rfrom,rto} - s_{p,rfrom,rto})^2$


## Implicit regional trade version

Current status: Model correctly solves (very quickly!) but it doesn't return the information on the amount of production allocated to own region or exported not the share of demand satisfied by own production or imported.

### Sets
- | $p$: products (with alias $pin$)
- | $r$: regions (with aliases $rin$ and $rto$)

### Variables
- | $s_{p,r}$: Supply (quantity)
- | $d_{p,r}$: Demand (quantity)
- | $p_{p,r}$: Price


### Parameters
- | $ϵs_{p,r,pin,rin}$: elasticity of the (total) supply of product $p$ in region $r$ with respect to the price of product $pin$ in region $rin$
- | $ϵd_{p,r,pin,rin}$: elasticity of the (total) demand of product $p$ in region $r$ with respect to the price of product $pin$ in region $rin$
- | $consts_{p,r}$: constant coefficient of the supply function for product $p$ in region $r$
- | $constd_{p,r}$: constant coefficient of the supply function for product $p$ in region $r$


### Constraints (equations)
- **Total supply [p,r]** : $s_{p,r} = consts_{p,r} * \prod_{pin,rin} p_{pin,r}^{ϵs_{p,r,pin,rin}}$
- **Total demand [p,r]** : $d_{p,r} = constd_{p,r} * \prod_{pin,rin} p_{pin,r}^{ϵd_{p,r,pin,rin}}$

### Objective
$\min \sum_{p,r} (d_{p,r} - s_{p,r})^2$