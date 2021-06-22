import voidfinder.box as vf

# first step is to convert our tracer file to unformatted F90 file.
data = np.genfromtxt(input_ascii)
vf.save_as_unformatted(data, tracers_filename)


# find prospective void centres by running a Delaunay triangulation
centres = vf.delaunay_triangulation(tracers_filename, centres_filename)


# grow spheres
voids = vf.grow_spheres(tracers_filename, centres_filename, voids_filename)

# remove overlapping voids
voids = vf.overlapping_filter(voids_filename)



# alternative is to run the full pipeline
vf.full_pipeline(input_ascii)
