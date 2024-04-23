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
    return m.md.namespace[:Agriculture].namespace[:gtap_label]
end
