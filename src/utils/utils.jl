using Interpolations, Mimi

"""
Return the name of the module of a given component with `comp_name` in model `m`.
This is a small helper function useful for internals like the mcs.
"""
function _get_module_name(m::Model, comp_name::Symbol)
    return nameof(m.md.namespace[comp_name].comp_id.module_obj)
end

"""
Return the name of the Moore agriculture GTAP damage function specification in 
model `m`. This is a small helper function useful for internals like the mcs.
"""
function _get_mooreag_gtap(m::Model)

    # model may not have been run yet, so need to get model parameter name to look 
    # up the value
    model_param_name = Mimi.get_model_param_name(m, :Agriculture, :gtap_name)
    Agriculture_gtap = Mimi.model_param(m, model_param_name).value

    return Agriculture_gtap
end
