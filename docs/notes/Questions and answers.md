# Questions and answers for Claude for the market module

## Clarifications

### What is a process/nest ?

A process is the rappresentation of a transformation. It divide the inputs in nests (groups). Within an individual nest, input products are substitutable (if more than one), across nests they are required in the same rations.
A process can make multiple products:

### What are shares ? Why we need them ?

### What about sigma in between 0 and 1 ?
In the docstring of Nest it is said to keep sigma > 1, but the in the leontief function it is set to zero. Why these 2 disjoint requirements (zero or > 1)? What for processes that are only partially substitute ?
ok, it sais these are hard for the solver to go to an equilibrium

### [solved] What is vacost ?
VACost is the overall costs of the process not inputed explicitly to the provided production inputs (electricity, labour...)

### Why the producer surplus don't enter the maximisation ?
Because it is just integral of the demand minus integral of the supply: the two pq rectangle cancel them out

### Are the frigth costs per unit of value ( if so, why?) or unit of quantity ?
per unit of quantity

## Possible additions to the model

### How to integrate in a process a nest with products that are exogenous, e.g. for which we have an exogenous price (like the "resin" example) ?

### How to model in the products chains situations with multiple transformations, e.g. roundwood -> sawnwood -> pallet

### How do I model subsitution in the final demand between modelled products and exogenous products and how do I link demand to exogenous conditions (e;g. income levels, population, nex construction houses) ?

### How do I connect primary supply to other exogenous conditions (mainly the availability of a given set of available resource for a primary product)

### [done] How do I introduce imperfect substitution between regions (eg Armington) ?
Implemented: add an `Armington(product = ..., sigma = ..., shares = ...)` entry to
`MarketData.armington`. Each destination then uses a CES composite of the varieties of
all the origins it can buy from (its own included), so regions cross-haul and every
price responds to supply and demand everywhere. `sigma = Inf` — the default for any
product not listed — is exactly the homogeneous/Samuelson case, so nothing changes
unless you ask for it. `shares` are keyed `(origin, destination)` and are the shares
at equal delivered prices: calibrate them on a base-year trade matrix, as the default
of equal shares implies a very strong taste for imports. The `prices` table now has
two columns: `price` (what local users pay for the composite) and `producer_price`
(what local producers get). See the "Imperfect substitution between origins" sections
of the "Using the module" and "Modelling choices" pages.

## How to introduce substitution with exogenous products in the demand (demand depend on the ratio of the relative prices, like fuelwood depends on price of fossil fuels) ?

## How to introduce exogenous elements in the supplies? Eg hardwood supply depends on availability of forest resources or the price of fuels ? Can it be addictive or multiplicative ?

## How to introduce limits in the upper bounds of transformations or trade (capacity) ?

## References
Need to add proper references, optimally 2 reference for each concept: the theorical one - the paper that introduced the concept - and the implementaitonal one - a (possibly open source) repository where the concept is implemented in code.
So, accross the documentation, there will be references and these will end up in a "references" page (still in the documentation)

