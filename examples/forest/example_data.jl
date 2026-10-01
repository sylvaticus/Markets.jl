# =============================================================================
#  Example economy for the Markets.jl engine: the forest-products sector.
#
#  THIS is the file you edit to amend the model.  Add a product to `products`,
#  give it a demand or supply curve, wire it into a `Process`, list it as
#  `tradable` — the engine (the Markets package) rebuilds itself automatically.
#
#  Taxonomy follows the FAOSTAT / UNECE forest-products chain, split by
#  coniferous (softwood) vs non-coniferous (hardwood) roundwood.  All numbers
#  are ILLUSTRATIVE but chosen to be the right order of magnitude:
#    * roundwood / sawnwood / panels in million m³, pulp & paper in million t;
#    * prices in USD per m³ (wood) or per t (pulp, paper);
#    * conversion coefficients near published wood-use / recovery factors.
#  Replace them with calibrated values from your own statistics when ready.
# =============================================================================

using Markets

# ---- regions ----------------------------------------------------------------
regions = [:EU, :NA, :AS]               # Europe, North America, Asia-Pacific

# ---- products ---------------------------------------------------------------
#   primary (from the forest):  swr  softwood roundwood,  hwr  hardwood roundwood
#   residue / intermediate:      chips (sawmill residues),  pulp (wood pulp)
#   final (consumed):            sawn_sw, sawn_hw, panel, paper
products = [:swr, :hwr, :chips, :pulp, :sawn_sw, :sawn_hw, :panel, :paper]

# ---- final demand (constant elasticity, η > 1) ------------------------------
# reference quantities differ by region → different regional markets.
demand = DemandSpec[
    # softwood sawnwood — construction grade
    DemandSpec(:sawn_sw, :EU; p0 = 250, q0 =  90, elasticity = 1.3),
    DemandSpec(:sawn_sw, :NA; p0 = 250, q0 = 120, elasticity = 1.3),
    DemandSpec(:sawn_sw, :AS; p0 = 250, q0 =  80, elasticity = 1.3),
    # hardwood sawnwood — furniture / appearance grade
    DemandSpec(:sawn_hw, :EU; p0 = 300, q0 =  15, elasticity = 1.3),
    DemandSpec(:sawn_hw, :NA; p0 = 300, q0 =  20, elasticity = 1.3),
    DemandSpec(:sawn_hw, :AS; p0 = 300, q0 =  60, elasticity = 1.3),
    # wood-based panels
    DemandSpec(:panel,   :EU; p0 = 350, q0 =  55, elasticity = 1.4),
    DemandSpec(:panel,   :NA; p0 = 350, q0 =  40, elasticity = 1.4),
    DemandSpec(:panel,   :AS; p0 = 350, q0 =  90, elasticity = 1.4),
    # paper & paperboard (priced higher: it sits at the end of a long, wood- and
    # conversion-intensive chain, so its equilibrium price is naturally high)
    DemandSpec(:paper,   :EU; p0 = 1000, q0 =  90, elasticity = 1.2),
    DemandSpec(:paper,   :NA; p0 = 1000, q0 =  80, elasticity = 1.2),
    DemandSpec(:paper,   :AS; p0 = 1000, q0 = 130, elasticity = 1.2),
]

# ---- primary supply (roundwood from the forest, elasticity ε > 0) -----------
# NA is softwood-rich, AS is hardwood-rich → comparative advantage drives trade.
supply = SupplySpec[
    SupplySpec(:swr, :EU; p0 = 70, q0 = 250, elasticity = 0.6),
    SupplySpec(:swr, :NA; p0 = 65, q0 = 380, elasticity = 0.6),
    SupplySpec(:swr, :AS; p0 = 80, q0 = 150, elasticity = 0.6),
    SupplySpec(:hwr, :EU; p0 = 80, q0 =  90, elasticity = 0.5),
    SupplySpec(:hwr, :NA; p0 = 80, q0 =  70, elasticity = 0.5),
    SupplySpec(:hwr, :AS; p0 = 75, q0 = 180, elasticity = 0.5),
]

# ---- processes (multi-stage transformation with residues) -------------------
# Sawmilling is fixed-proportion joint production (sawnwood + chip residues).
# Panel & pulp mills use CES nests → smooth substitution between wood inputs:
# when chips get scarce/expensive the mills shift toward roundwood, and a
# softwood/hardwood mill substitutes between fibre types along the σ elasticity.
processes = [
    Process(:sawmill_sw,
            [leontief(:swr, 1.0)],                       # 1 m³ softwood log in …
            [:sawn_sw => 0.50, :chips => 0.35];          # … → lumber + residues
            vacost = 40),
    Process(:sawmill_hw,
            [leontief(:hwr, 1.0)],
            [:sawn_hw => 0.45, :chips => 0.35];
            vacost = 45),
    Process(:panelmill,                                  # particleboard / fibreboard
            [ces(1.30, [:chips, :swr, :hwr], [0.55, 0.30, 0.15]; sigma = 2.5)],
            [:panel => 1.0];
            vacost = 120),
    Process(:pulpmill,                                   # softwood = long fibre (preferred)
            [ces(4.0, [:swr, :hwr, :chips], [0.45, 0.30, 0.25]; sigma = 1.8)],
            [:pulp => 1.0];
            vacost = 300),
    Process(:papermill,
            [leontief(:pulp, 1.10)],                     # ~1.1 t pulp per t paper
            [:paper => 1.0];
            vacost = 200),
]

# ---- which products can be traded -------------------------------------------
tradable = [:swr, :hwr, :chips, :pulp, :sawn_sw, :sawn_hw, :panel, :paper]

# ---- transport costs  (freight rate × relative distance) --------------------
# Bulky low-value goods (roundwood, chips) cost more per unit of value to ship.
freight = Dict(:swr => 12.0, :hwr => 12.0, :chips => 14.0, :pulp => 35.0,
               :sawn_sw => 30.0, :sawn_hw => 32.0, :panel => 35.0, :paper => 55.0)
distance = Dict((:EU, :NA) => 1.0, (:EU, :AS) => 1.2, (:NA, :AS) => 1.1)
reldist(a, b) = a == b ? 0.0 : haskey(distance, (a, b)) ? distance[(a, b)] : distance[(b, a)]

transport = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
for p in tradable, a in regions, b in regions
    a == b && continue
    transport[(p, a, b)] = freight[p] * reldist(a, b)
end

# ---- assemble ---------------------------------------------------------------
example_market = MarketData(; regions, products, demand, supply,
                              processes, tradable, transport)
