import numpy
from matplotlib import pyplot
from amuse.plot import xlabel, ylabel, effective_iso_potential_plot
from amuse.units import units, constants, nbody_system, quantities
from amuse.community.hermite0.interface import Hermite
from amuse.datamodel import Particles
from amuse import plot as aplot
from amuse.lab import *

def read_files(filename):
    b = read_set_from_file(filename+"_binary.amuse", "hdf5", close_file=True)
    print b
    c = read_set_from_file(filename+"_core.amuse", "hdf5", close_file=True)
    print c
    g = read_set_from_file(filename+"_gas.amuse", "hdf5", close_file=True)
    print g
    return b, c, g

def setup_gravity_code(triple):
    converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
    gravity = Hermite(converter)
    gravity.particles.add_particles(triple)
    return gravity

def make_effective_iso_potential_plot(gravity_code, lim):
    omega = (constants.G * gravity_code.particles.total_mass()
              / (1.0|units.AU**3)).sqrt()
    center_of_mass = gravity_code.particles.center_of_mass()[:2]

    """
    pyplot.rcParams.update({'font.size': 30})
    figure = pyplot.figure(figsize = (12, 12))
    ax = pyplot.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.minorticks_on() 
    """

    lim = lim.value_in(units.au)
    #current_axes = pyplot.subplot(1, 1, 1)
    #current_axes.set_aspect("equal", adjustable = "box")
    effective_iso_potential_plot(gravity_code, 
                                 omega,
                                 xlim = [-lim, lim] | units.AU,
                                 ylim = [-lim, lim] | units.AU,
                                 center_of_rotation = center_of_mass,
                                 fraction_screen_filled=0.8,
                                 resolution = [4000, 4000],
                                 number_of_contours = 20)
    
#    xlabel('x [au]')
#    ylabel('y [au]')
    """
    pyplot.text(0.6, -0.06, "$L_1$")
    pyplot.text(1.35, -0.06, "$L_2$")
    pyplot.text(-0.99, -0.06, "$L_3$")
    pyplot.text(0.40, 0.82, "$L_4$")
    pyplot.text(0.40, -0.92, "$L_5$")
    """
   

def get_densities(snapshot, N=200, plot_range=10 | units.au):
    x = quantities.linspace(-plot_range, plot_range, N)
    y = quantities.linspace(-plot_range, plot_range, N)
    X, Y = quantities.meshgrid(x, y)

    x = X.flatten()
    y = Y.flatten()
    z = [0.] | units.au

    converter=nbody_system.nbody_to_si(snapshot.mass.sum(), 1000|units.au)
    hydro = Fi(converter)
    hydro.parameters.eps_is_h_flag = True
    hydro.parameters.timestep = 1|units.minute
    hydro.gas_particles.add_particles(snapshot)
    print "particles added."

    rho = hydro.get_hydro_state_at_point(x, y, z)[0]
    rho = rho.reshape((N, N))
    #hydro.close()
    return X, Y, rho
    
def plot_gas(snapshot, N, plot_range):
    X, Y, rho = get_densities(snapshot, N=N, plot_range=plot_range)

    print "extrema:", rho.max().in_(units.g/units.cm**3), rho.min().in_(units.g/units.cm**3)
    vmax = 0.1*rho.max().in_(units.g/units.cm**3)
    vmin = 10.*rho.min().in_(units.g/units.cm**3)

    g = aplot.imshow_color_plot(X.value_in(units.au), Y.value_in(units.au), rho.value_in(units.g/units.cm**3), cmap='inferno', vlog=True, aspect=N)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "hydro_triple_0600",
                      help="input filename [%default]")
    return result

if __name__ == "__main__":
    o, arguments  = new_option_parser().parse_args()

    b, c, g = read_files(o.filename)
    b.y -= c.y
    g.y -= c.y
    c.y -= c.y
    triple = b.copy()
    c.mass += g.mass.sum()
    triple.add_particles(c)
    
    figure = pyplot.figure(figsize=(14, 14))
    current_axes = figure.gca()
    #current_axes = pyplot.subplot(1, 1, 1)
    current_axes.set_facecolor('#101010')
    #current_axes.set_aspect("equal")
    current_axes.set_aspect("equal", adjustable = "box")

    lim = 4.0 | units.au

    gravity = setup_gravity_code(triple)
    make_effective_iso_potential_plot(gravity, lim)
    gravity.stop()
    
    plot_gas(g, N=1000, plot_range=lim)
    #pyplot.tight_layout()

    pyplot.xlabel('x [au]', fontsize=20)
    pyplot.ylabel('y [au]', fontsize=20)
    
    save_file = 'lagrange_points_and_sph.png'
    pyplot.savefig(save_file)
    print "\nOutput saved in", save_file
    pyplot.show()
