using Interpolations, Mimi

"""
Return the name of the module of a given component with `comp_name` in model `m`.
This is a small helper function useful for internals like the mcs.
"""
function _get_module_name(m::Model, comp_name::Symbol)
    return nameof(m.md.namespace[comp_name].comp_id.module_obj)
end