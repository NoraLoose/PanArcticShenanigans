# Workaround for NumericalEarth #676 (nested_ocean_model): the coupled model applies the atmosphere's pressure to the
# ocean through `forcing_barotropic_potential(coupled_model.ocean)`, which has no method for a NestedModel, so it falls
# back to `nothing` and the pressure is never applied to a nested ocean. The child's barotropic potential then stays
# zero, and so does #676's inverse-barometer sea level (`inverse_barometer = true` changes nothing). Forward the call
# to the child, which holds the pressure forcing. To be fixed in #676; remove this file then.

NumericalEarth.Oceans.forcing_barotropic_potential(model::NumericalEarth.NestedModels.NestedModel) =
    NumericalEarth.Oceans.forcing_barotropic_potential(model.child)
