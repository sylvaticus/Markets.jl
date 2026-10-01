# Questions and answers for Claude for the market module

## Clarifications

### What is a process/nest ?

A process is the rappresentation of a transformation. It divide the inputs in nests (groups). Within an individual nest, input products are substitutable (if more than one), across nests they are required in the same rations.
A process can make multiple products:

### What are shares ? Why we need them ?

### What about sigma in between 0 and 1 ?
In the docstring of Nest it is said to keep sigma > 1, but the in the leontief function it is set to zero. Why these 2 disjoint requirements (zero or > 1)? What for processes that are only partially substitute ?

### [solved] What is vacost ?
VACost is the overall costs of the process not inputed explicitly to the provided production inputs (electricity, labour...)

### Why the producer surplus don't enter the maximisation ?

### Are the frigth costs per unit of value ( if so, why?) or unit of quantity ?


## Possible additions to the model

### How to integrate in a process a nest with products that are exogenous, e.g. for which we have an exogenous price (like the "resin" example) ?

### How to model in the products chains situations with multiple transformations, e.g. roundwood -> sawnwood -> pallet

### How do I model subsitution in the final demand between modelled products and exogenous products and how do I link demand to exogenous conditions (e;g. income levels, population, nex construction houses) ?

### How do I connect primary supply to other exogenous conditions (mainly the availability of a given set of available resource for a primary product)

### How do I introduce imperfect substitution between regions (eg Armington) ?

