using Mimi
using MimiGIVE

include("new_sector_damages.jl")
include("DamageAggregator_NewSectorDamages.jl")

function get_new_sector_model()

    # Obtain MimiGIVE model
    m = MimiGIVE.get_model()

    # Add new damage sector component
    add_comp!(m, NewSectorDamages, first = 2020, after = :energy_damages)

    # Replace Damage Aggregator component with modified one
    replace!(m, :DamageAggregator => DamageAggregator_NewSectorDamages)

    # Need to set this damage aggregator to run from 2020 to 2300, currently picks up
    # 1750 to 2300 from replace!
    Mimi.set_first_last!(m, :DamageAggregator, first=2020);

    # Connections
    connect_param!(m, :NewSectorDamages => :temperature, :temperature => :T)
    connect_param!(m, :NewSectorDamages => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :DamageAggregator => :damage_new_sector, :NewSectorDamages => :damages)

    return m
end