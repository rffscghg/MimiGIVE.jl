# constants needed for SCC calculations

const _model_years = collect(1750:2300)
const _damages_years = collect(2020:2300)
const _damages_idxs = indexin(_damages_years, _model_years)

const scc_gas_molecular_conversions = Dict(:CO2 => 12/44, # C to CO2
                                            :N2O => 28/44, # N2 to N2O,
                                            :CH4 => 1., # CH4 to CH4
                                            :HFC23 => 1., # HFC23 to HFC23
                                            :HFC32 => 1., # HFC32 to HFC32
                                            :HFC43_10 => 1., # HFC43_10 to HFC43_10
                                            :HFC125 => 1., # HFC125 to HFC125
                                            :HFC134a => 1., # HFC134a to HFC134a
                                            :HFC143a => 1., # HFC143a to HFC143a
                                            :HFC227ea => 1., # HFC227ea to HFC227ea
                                            :HFC245fa => 1.) # HFC245fa to HFC245fa

const scc_gas_pulse_size_conversions = Dict(:CO2 => 1e9, # Gt to t
                                        :N2O => 1e6, # Mt to t
                                        :CH4 => 1e6, # Mt to t
                                        :HFC23 => 1e3, # kt to t
                                        :HFC32 => 1e3, # kt to t
                                        :HFC43_10 => 1e3, # kt to t
                                        :HFC125 => 1e3, # kt to t
                                        :HFC134a => 1e3, # kt to t
                                        :HFC143a => 1e3, # kt to t
                                        :HFC227ea => 1e3, # kt to t
                                        :HFC245fa => 1e3) # kt to t