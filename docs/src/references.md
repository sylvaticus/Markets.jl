# References

Two lists, kept apart on purpose, because they answer two different questions
about any piece of a model.

**Where does this idea come from?** is answered by the theory list, cited as
**Author (year)** — the paper that introduced the concept, or the one that
settled how it is used.

**Where can I see it in working code?** is answered by the implementation list,
cited as **[C1]**, **[C2]** and so on — a repository, preferably open source,
in which the same concept is written out and runs. Reading a CES nest in an
article and reading one in a solver's input are different kinds of
understanding, and a modelling package is in a position to offer both.

So a sentence in these pages may read: *the spatial price equilibrium of
[Samuelson (1952)](@ref samuelson1952), implemented for the forest sector in
[[C3]](@ref C3)*. The first takes you to the argument, the second to the code.

Entries are listed alphabetically within each list. Every DOI and repository
below was checked against CrossRef or the GitHub API when the page was written;
where no DOI exists — papers of the 1950s mostly predate them — the journal,
volume and pages are given instead.

## Theory

#### [Armington, P. S. (1969)](@id armington1969)

A theory of demand for products distinguished by place of production.
*Staff Papers — International Monetary Fund* 16(1), 159–178.
[doi:10.2307/3866403](https://doi.org/10.2307/3866403)

> The origin of treating goods from different places as imperfect substitutes.
> Used throughout [Imperfect substitution between origins (Armington)](@ref).

#### [Arrow, K. J., Chenery, H. B., Minhas, B. S. & Solow, R. M. (1961)](@id acms1961)

Capital-labor substitution and economic efficiency.
*The Review of Economics and Statistics* 43(3), 225–250.
[doi:10.2307/1927286](https://doi.org/10.2307/1927286)

> The constant-elasticity-of-substitution function itself, which this package
> uses three times over: for process inputs, for origins and for substitutes in
> final demand.

#### [Buongiorno, J. (2015)](@id buongiorno2015)

Income and time dependence of forest product demand elasticities and
implications for forecasting. *Silva Fennica* 49(5), 1395.
[doi:10.14214/sf.1395](https://doi.org/10.14214/sf.1395)

> Evidence on the size of forest-product demand elasticities, and so on the
> inelastic range this package deliberately supports; see
> [Demand and supply curves](@ref).

#### [Buongiorno, J., Zhu, S., Zhang, D., Turner, J. & Tomberlin, D. (2003)](@id gfpm2003)

*The Global Forest Products Model: structure, estimation, and applications.*
Academic Press.
[doi:10.1016/b978-012141362-0/50004-2](https://doi.org/10.1016/b978-012141362-0/50004-2)
(the chapter stating the formulation)

> The reference forest-sector partial-equilibrium model, and the lineage this
> package's example belongs to.

#### [Enke, S. (1951)](@id enke1951)

Equilibrium among spatially separated markets: solution by electric analogue.
*Econometrica* 19(1), 40–47.
[doi:10.2307/1907907](https://doi.org/10.2307/1907907)

> The problem posed, one year before Samuelson turned it into an optimisation —
> and solved, memorably, by building the circuit.

#### [Ferris, M. C. & Pang, J.-S. (1997)](@id ferris1997)

Engineering and economic applications of complementarity problems.
*SIAM Review* 39(4), 669–713.
[doi:10.1137/S0036144595285963](https://doi.org/10.1137/S0036144595285963)

> The survey of what complementarity formulations can express that an
> optimisation cannot, which is the case for leaving the route this package
> takes; see
> [Why processes, and not supply and demand per transformation](@ref).

#### [Gerlagh, R., Lorang, E., Pommeret, A. & Pottier, A. (2026)](@id glpp2026)

Source taxes versus end-of-chain taxes in general equilibrium. Working paper,
July 2026.

> When taxing an externality where it arises and taxing it as embodied in
> final goods are equivalent, and why the second needs prices that depend on
> each good's content. See [Where along the chain to tax](@ref).

#### [Howitt, R. E. (1995)](@id howitt1995)

Positive mathematical programming.
*American Journal of Agricultural Economics* 77(2), 329–342.
[doi:10.2307/1243543](https://doi.org/10.2307/1243543)

> Calibrating a programming model so that the observed activity is an interior
> optimum rather than a corner held by a constraint — the standard remedy for
> the constant-returns behaviour discussed in
> [Why processes, and not supply and demand per transformation](@ref).

#### [Latta, G. S., Sjølie, H. K. & Solberg, B. (2013)](@id latta2013)

A review of recent developments and applications of partial equilibrium models
of the forest sector. *Journal of Forest Economics* 19(4), 350–360.
[doi:10.1016/j.jfe.2013.06.006](https://doi.org/10.1016/j.jfe.2013.06.006)

> What the forest-sector partial-equilibrium literature contains and how its
> models are built, which is the lineage this package's example belongs to.

#### [Lorang, E. & Gerlagh, R. (2026)](@id lorang2026)

Resource extraction taxes versus waste taxes in general equilibrium. Working
paper, April 2026.

> Material balances kept jointly with market balances by letting the material
> intensity of every good follow from the allocation, which is what
> [`material_accounts`](@ref) does after the solve; and the efficiency gap
> between taxing material at extraction and at disposal. See
> [Material accounts](@ref) and [Where along the chain to tax](@ref).

#### [Lubin, M., Dowson, O., Garcia, J., Huchette, J. & Legat, B. (2023)](@id jump2023)

JuMP 1.0: recent improvements to a modeling language for mathematical
optimization. *Mathematical Programming Computation* 15(3), 581–589.
[doi:10.1007/s12532-023-00239-3](https://doi.org/10.1007/s12532-023-00239-3)

> The modelling layer this package is written in; see [[C1]](@ref C1).

#### [Mathiesen, L. (1985)](@id mathiesen1985)

Computation of economic equilibria by a sequence of linear complementarity
problems. *Mathematical Programming Studies*, 144–162.
[doi:10.1007/bfb0121030](https://doi.org/10.1007/bfb0121030)

> Equilibrium as a complementarity problem — the alternative to the
> optimisation route taken here, and the one that would make ad valorem
> tariffs exact; see
> [Why ad valorem rates on trade are special](@ref).

#### [McCarl, B. A. & Spreen, T. H. (1980)](@id mccarl1980)

Price endogenous mathematical programming as a tool for sector analysis.
*American Journal of Agricultural Economics* 62(1), 87–102.
[doi:10.2307/1239475](https://doi.org/10.2307/1239475)

> The case for representing a sector as a programming problem whose prices are
> duals, including what integrability demands of it. Directly on the choice
> made in [Why processes, and not supply and demand per transformation](@ref).

#### [McDaniel, C. A. & Balistreri, E. J. (2003)](@id mcdaniel2003)

A review of Armington trade substitution elasticities.
*Économie internationale* 94–95(2), 301–313.
[doi:10.3917/ecoi.094.0301](https://doi.org/10.3917/ecoi.094.0301)

> What the Armington elasticity has actually been estimated to be, and how
> widely the estimates vary — worth reading before choosing a `sigma`.

#### [Rutherford, T. F. (1995)](@id rutherford1995)

Extension of GAMS for complementarity problems arising in applied economic
analysis. *Journal of Economic Dynamics and Control* 19(8), 1299–1324.
[doi:10.1016/0165-1889(94)00831-2](https://doi.org/10.1016/0165-1889(94)00831-2)

> MPSGE: equilibrium models written as complementarity problems, with CES and
> Armington nests as primitives. Implemented in [[C4]](@ref C4).

#### [Samuelson, P. A. (1952)](@id samuelson1952)

Spatial price equilibrium and linear programming.
*American Economic Review* 42, 283–303. (Predates DOIs; available on JSTOR.)

> The result this package rests on: that the competitive spatial equilibrium is
> the maximum of net social surplus. See
> [The equilibrium as surplus maximisation](@ref).

#### [Sato, K. (1967)](@id sato1967)

A two-level constant-elasticity-of-substitution production function.
*The Review of Economic Studies* 34(2), 201–218.
[doi:10.2307/2296809](https://doi.org/10.2307/2296809)

> Nesting one CES inside another, which is how this package lets substitution
> differ between pairs — between origins, and between products in demand —
> without an inconsistent matrix of elasticities.

#### [Takayama, T. & Judge, G. G. (1964)](@id takayama1964)

Equilibrium among spatially separated markets: a reformulation.
*Econometrica* 32(4), 510–524.
[doi:10.2307/1910175](https://doi.org/10.2307/1910175)

> The quadratic-programming formulation of the spatial equilibrium, and the
> standard reference for solving it numerically. Extended in their 1971 book
> *Spatial and Temporal Price and Allocation Models* (North-Holland).

#### [van der Mensbrugghe, D. & Peters, J. C. (2016)](@id vdmp2016)

Volume preserving CES and CET formulations. Paper presented at the 19th
Annual Conference on Global Economic Analysis, Washington DC, June 2016.
Center for Global Trade Analysis, Purdue University.

> The additive CES and CET, which keep imperfect substitution and make the
> aggregate a plain sum, and how far the standard forms drift from one. See
> [What it would take to make the composites add up](@ref).

#### [van Tongeren, F., van Meijl, H. & Surry, Y. (2001)](@id vantongeren2001)

Global models applied to agricultural and trade policies: a review and
assessment. *Agricultural Economics* 26(2), 149–172.
[doi:10.1111/j.1574-0862.2001.tb00060.x](https://doi.org/10.1111/j.1574-0862.2001.tb00060.x)

> A classification of sector models by what they assume rather than by what
> they cover, and the source for which paradigm prevails in which community.

#### [Wächter, A. & Biegler, L. T. (2006)](@id wachter2006)

On the implementation of an interior-point filter line-search algorithm for
large-scale nonlinear programming. *Mathematical Programming* 106(1), 25–57.
[doi:10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y)

> The algorithm behind Ipopt, this package's default solver — and the reason
> quantities raised to a power below one need a floor away from zero; see
> [Numerical details](@ref). Implemented in [[C2]](@ref C2).

## Implementations

#### [[C1] JuMP.jl](@id C1)

<https://github.com/jump-dev/JuMP.jl>

> The modelling language this package is built on: every variable, constraint
> and objective in `src/model.jl` is a JuMP expression, and the duals that
> become prices are read back through it. Described in
> [Lubin et al. (2023)](@ref jump2023).

#### [[C2] Ipopt](@id C2)

<https://github.com/coin-or/Ipopt> — EPL-2.0

> The interior-point solver used by default, through
> [JuMP](@ref C1). Any solver handling nonlinear constraints and returning
> duals can replace it. Described in [Wächter & Biegler (2006)](@ref wachter2006).

#### [[C3] TiMBA](@id C3)

<https://github.com/TI-Forest-Sector-Modelling/TiMBA> — AGPL-3.0 —
[doi:10.21105/joss.08034](https://doi.org/10.21105/joss.08034)

> The closest open-source relative of this package: a spatial partial
> equilibrium of the global forest sector, in Python, built on the
> [GFPM](@ref gfpm2003) formulation and solved — as here — by maximising social
> surplus period by period. Read it to see the same economics with different
> choices made about dynamics, trade and data.

#### [[C4] MPSGE.jl](@id C4)

<https://github.com/julia-mpsge/MPSGE.jl> — MIT

> The complementarity route in Julia, implementing
> [Rutherford's (1995)](@ref rutherford1995) MPSGE. Nested CES and Armington
> structures are primitives of the language rather than constraints one writes
> out, and prices are variables rather than duals — which is what lets it treat
> an ad valorem tax exactly where this package has to iterate.

#### [[C5] PATHSolver.jl](@id C5)

<https://github.com/chkwon/PATHSolver.jl> — MIT

> The mixed-complementarity solver those formulations are solved with, wrapped
> for Julia; [Complementarity.jl](https://github.com/chkwon/Complementarity.jl)
> provides the modelling layer above it.
