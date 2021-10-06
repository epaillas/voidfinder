from voidfinder.box import SphericalVoidFinder

parameter_file = "find_voids_periodic_box.yaml"

SVF = SphericalVoidFinder(parameter_file=parameter_file)

SVF.full_pipeline()
