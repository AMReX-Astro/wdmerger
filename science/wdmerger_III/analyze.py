import wdmerger

linestyles = ['-', '--', ':', '-.']
markers = ['o', 's', 'D', '^']

def get_rank():
    """Return MPI rank, or 0 if not using MPI."""

    try:

        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

    except:

        rank = 0

    return rank



def is_root():
    """Return whether we are the root processor (0, by convention)."""

    rank = get_rank()

    if (rank == 0):
        return True
    else:
        return False



# Gravitational wave signal

def gravitational_wave_signal(eps_filename, results_dir, max_time = 1.e20, markevery = None,
                              scale_exp = -22, vary_lines = False):
    """Plot the gravitational radiation waveform."""

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_dir):
        return

    from matplotlib import pyplot as plt
    import numpy as np

    print "Generating plot with filename " + eps_filename

    scale = 10.0**scale_exp

    if vary_lines:
        linestyles = ['-', '--', ':', '-', '--', ':']
        markers = ['', '', '', 'o', 's', 'D']
        colors = ['b', 'g', 'r', 'c', 'm', 'k']
    else:
        linestyles = ['-', '--', '-', '--', '-', '--']
        markers = ['', '', '', '', '', '']
        colors = ['b', 'g', 'b', 'g', 'b', 'g']

    diag_file = results_dir + '/output/grid_diag.out'

    time       = np.array(wdmerger.get_column('TIME', diag_file))

    strain_p_1 = np.array(wdmerger.get_column('h_+ (axis 1)', diag_file)) / scale
    strain_p_2 = np.array(wdmerger.get_column('h_+ (axis 2)', diag_file)) / scale
    strain_p_3 = np.array(wdmerger.get_column('h_+ (axis 3)', diag_file)) / scale

    strain_x_1 = np.array(wdmerger.get_column('h_x (axis 1)', diag_file)) / scale
    strain_x_2 = np.array(wdmerger.get_column('h_x (axis 2)', diag_file)) / scale
    strain_x_3 = np.array(wdmerger.get_column('h_x (axis 3)', diag_file)) / scale

    # Generate an offset for the x- and y- strains. Then round this to the nearest
    # integer for the convenience of the plot legend.

    offset = max( int(round(2.5 * (max(np.max(strain_p_3), np.max(strain_x_3)) - min(np.min(strain_p_2), np.min(strain_x_2))))),
                  int(round(2.5 * (max(np.max(strain_p_2), np.max(strain_x_2)) - min(np.min(strain_p_1), np.min(strain_x_1))))) )
    z_offset = 0
    y_offset = offset
    x_offset = 2 * offset

    strain_p_1 += x_offset
    strain_p_2 += y_offset
    strain_p_3 += z_offset

    strain_x_1 += x_offset
    strain_x_2 += y_offset
    strain_x_3 += z_offset

    z_label_offset = ''
    y_label_offset = '+ {:d}'.format(y_offset)
    x_label_offset = '+ {:d}'.format(x_offset)

    label_p_1 = r'$h^x_+{}$'.format(x_label_offset)
    label_p_2 = r'$h^y_+{}$'.format(y_label_offset)
    label_p_3 = r'$h^z_+{}$'.format(z_label_offset)

    label_x_1 = r'$h^x_\times{}$'.format(x_label_offset)
    label_x_2 = r'$h^y_\times{}$'.format(y_label_offset)
    label_x_3 = r'$h^z_\times{}$'.format(z_label_offset)

    markersize = 12

    mask = np.where(time <= max_time)

    plot_p_1, = plt.plot(time[mask], strain_p_1[mask], lw = 4.0, color = colors[0], linestyle = linestyles[0],
                         marker = markers[0], markersize = markersize, markevery = markevery, label = label_p_1)
    plot_p_2, = plt.plot(time[mask], strain_p_2[mask], lw = 4.0, color = colors[2], linestyle = linestyles[2],
                         marker = markers[2], markersize = markersize, markevery = markevery, label = label_p_2)
    plot_p_3, = plt.plot(time[mask], strain_p_3[mask], lw = 4.0, color = colors[4], linestyle = linestyles[4],
                         marker = markers[4], markersize = markersize, markevery = markevery, label = label_p_3)

    plot_x_1, = plt.plot(time[mask], strain_x_1[mask], lw = 4.0, color = colors[1], linestyle = linestyles[1],
                         marker = markers[1], markersize = markersize, markevery = markevery, label = label_x_1)
    plot_x_2, = plt.plot(time[mask], strain_x_2[mask], lw = 4.0, color = colors[3], linestyle = linestyles[3],
                         marker = markers[3], markersize = markersize, markevery = markevery, label = label_x_2)
    plot_x_3, = plt.plot(time[mask], strain_x_3[mask], lw = 4.0, color = colors[5], linestyle = linestyles[5],
                         marker = markers[5], markersize = markersize, markevery = markevery, label = label_x_3)

    plt.tick_params(labelsize=20)

    plt.xlabel(r'$t\ \mathrm{(s)}$', fontsize = 24)
    plt.ylabel(r'$h\, /\, 10^{}$'.format('{' + str(scale_exp) + '}'), fontsize = 24)

    xlim = [time[mask][0], time[mask][-1]]
    plt.xlim(xlim)
    ylim = [min(np.min(strain_p_3), np.min(strain_x_3)) - 0.05 * offset, max(np.max(strain_p_1), np.max(strain_x_1)) + 0.5 * offset]
    plt.ylim(ylim)

    # Taking a trick from http://matplotlib.org/users/legend_guide.html#multiple-legends-on-the-same-axes
    # to split the legend into three parts.

    legend_1 = plt.legend(handles = [plot_p_1, plot_x_1], loc = [0.1, 0.85], prop = {'size':20}, ncol = 2, fancybox = True)
    legend_2 = plt.legend(handles = [plot_p_2, plot_x_2], loc = [0.1, 0.55], prop = {'size':20}, ncol = 2, fancybox = True)
    legend_3 = plt.legend(handles = [plot_p_3, plot_x_3], loc = [0.1, 0.25], prop = {'size':20}, ncol = 2, fancybox = True)

    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.gca().add_artist(legend_3)

    plt.tight_layout()

    plt.savefig(eps_filename)

    wdmerger.insert_commits_into_eps(eps_filename, diag_file, 'diag')

    plt.close()



def vol_renders():
    """Volume render the simulations."""

    import os
    import yt
    yt.enable_parallelism()

    results_base = 'results/approximate/'
    plots_base = 'plots/'

    for mass_P in ['0.90']:
        for mass_S in ['0.60', '0.90']:
            for roche in ['0.90', '1.00']:
                for rot in ['0', '1']:
                    for hybrid in ['0', '1']:
                        for ncell in ['256', '512', '1024']:

                            mass_string = "_m_P_" + mass_P + "_m_S_" + mass_S

                            results_dir = 'mass_P_' + mass_P + '/mass_S_' + mass_S + '/roche' + roche + '/rot' + rot + '/hybrid' + hybrid + '/n' + ncell + '/'

                            output_dir = results_base + results_dir + '/output/'

                            if not os.path.exists(output_dir):
                                continue

                            plot_list = wdmerger.get_plotfiles(output_dir, prefix='smallplt')

                            plot_list = [output_dir + plot for plot in plot_list]

                            plot_dir = plots_base + results_dir

                            if not os.path.exists(plot_dir) and is_root():
                                os.makedirs(plot_dir)

                            ts = yt.DatasetSeries(plot_list)

                            for ds in ts.piter():

                                pltname = ds.basename

                                outfile_name = plot_dir + 'approximate' + mass_string + '_roche_' + roche + '_rot_' + rot + '_hybrid_' + hybrid + '_n_' + ncell + '_' + pltname + '.png'

                                if os.path.isfile(outfile_name):
                                    continue

                                wdmerger.vol_render_density(outfile_name, ds)



def stellar_distance_convergence(eps_filename, results_base):
    """Plot the WD distance as a function of resolution."""

    if not is_root():
        return

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_base):
        return

    import matplotlib.pyplot as plt
    import numpy as np

    res_list = wdmerger.get_parameter_list(results_base)
    res_list = sorted([int(res[1:]) for res in res_list])

    diag_list = [results_base + '/n' + str(res) + '/output/star_diag.out' for res in res_list]

    for i, diag in enumerate(diag_list):

        time = wdmerger.get_column('TIME', diag)
        dist = wdmerger.get_column('WD DISTANCE', diag)

        mask = np.where(time <= 50.0)

        plt.plot(np.array(time)[mask], np.array(dist)[mask] / np.array(dist)[0], lw = 4.0, linestyle = linestyles[i],
                 marker = markers[i], markersize = 12.0, markevery = 2000, label = 'n = ' + str(res_list[i]))

    plt.tick_params(labelsize = 20)

    plt.xlabel('Time (s)', fontsize = 24)
    plt.ylabel('WD Distance / Initial Distance', fontsize = 24)
    plt.legend(loc = 'best', prop = {'size':20}, fancybox = True, framealpha = 0.5)
    plt.tight_layout()
    plt.savefig(eps_filename)

    plt.close()



def stellar_distance(eps_filename, results_base):
    """Plot the WD distance as a function of time."""

    if not is_root():
        return

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_base):
        return

    import matplotlib.pyplot as plt
    import numpy as np

    res_list = wdmerger.get_parameter_list(results_base)
    res_list = sorted([int(res[1:]) for res in res_list])

    diag_list = [results_base + '/n' + str(res) + '/output/star_diag.out' for res in res_list]

    for i, diag in enumerate(diag_list):

        eps_filename_curr = eps_filename.replace('.eps', '_n_' + str(res_list[i]) + '.eps')

        if os.path.isfile(eps_filename_curr):
            return

        time = np.array(wdmerger.get_column('TIME', diag))
        dist = np.array(wdmerger.get_column('WD DISTANCE', diag))

        mask = np.where(dist > 0.0)

        plt.plot(time[mask], dist[mask] / dist[0], lw = 4.0, linestyle = linestyles[0])

        plt.tick_params(labelsize = 20)

        plt.xlabel('Time (s)', fontsize = 24)
        plt.ylabel('WD Distance / Initial Distance', fontsize = 24)
        plt.tight_layout()
        plt.savefig(eps_filename_curr)

        plt.close()



def angular_momentum_comparison(eps_filename, results_base, ncell, subtract_grid_losses = False):
    """Plot the system angular momentum as a function of evolution method."""

    if not is_root():
        return

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_base):
        return

    import matplotlib.pyplot as plt
    import numpy as np

    idx = 0

    for rot in ['0', '1']:
        for hybrid in ['0', '1']:

            results_dir = results_base + '/rot' + rot + '/hybrid' + hybrid + '/n' + ncell

            diag = results_dir + '/output/grid_diag.out'

            time = np.array(wdmerger.get_column('TIME', diag))
            lerr = np.array(wdmerger.get_column('ANG. MOM. Z', diag))

            if subtract_grid_losses:
                bndy_diag = results_dir + '/output/bndy_diag.out'
                bndy = np.array(wdmerger.get_column('ANG. MOM. Z LOST', bndy_diag))
                lerr = lerr + bndy

            mask = np.where(time <= 100.0)

            title = ''

            if rot == '0':
                title += 'Inertial Frame'
            else:
                title += 'Rotating Frame'

            if hybrid == '0':
                title += ', Standard Equations'
            else:
                title += ', Hybrid Equations'

            plt.plot(time[mask], abs((lerr[mask] - lerr[0]) / lerr[0]), lw = 4.0, linestyle = linestyles[idx], marker = markers[idx], markersize = 12.0, markevery = 5000, label = title)

            idx += 1

    plt.tick_params(labelsize = 20)

    plt.xlabel('Time (s)', fontsize = 24)
    plt.ylabel(r'Relative Change in $L_z$', fontsize = 24)
    plt.legend(loc = 'best', prop = {'size':16}, fancybox = True, framealpha = 0.5)
    plt.ylim([-5e-5, 1.e-3])
    plt.tight_layout()
    plt.savefig(eps_filename)

    plt.close()



def stellar_masses(eps_filename, results_base):
    """Plot the WD masses as a function of time."""

    if not is_root():
        return

    import os

    if not os.path.exists(results_base):
        return

    import matplotlib.pyplot as plt
    import numpy as np

    res_list = wdmerger.get_parameter_list(results_base)
    directory_list = [results_base + '/' + res for res in res_list]

    for i, directory in enumerate(directory_list):

        eps_filename_curr = eps_filename.replace('.eps', '_n_' + res_list[i][1:] + '.eps')

        if os.path.isfile(eps_filename_curr):
            continue

        p_diag = directory + '/output/primary_diag.out'
        s_diag = directory + '/output/secondary_diag.out'

        p_time = wdmerger.get_column('TIME', p_diag)
        p_mass = wdmerger.get_column('PRIMARY MASS', p_diag)

        s_time = wdmerger.get_column('TIME', s_diag)
        s_mass = wdmerger.get_column('SECONDARY MASS', s_diag)

        p_mask = np.where(p_time <= 100.0)
        s_mask = np.where(s_time <= 100.0)

        plt.plot(np.array(p_time)[p_mask], np.array(p_mass)[p_mask] / np.array(p_mass)[0], lw = 4.0, linestyle = linestyles[0], markersize = 12.0, markevery = 2000, label = 'Primary')
        plt.plot(np.array(s_time)[s_mask], np.array(s_mass)[s_mask] / np.array(s_mass)[0], lw = 4.0, linestyle = linestyles[1], markersize = 12.0, markevery = 2000, label = 'Secondary')

        plt.tick_params(labelsize = 20)

        plt.xlabel('Time (s)', fontsize = 24)
        plt.ylabel('WD Mass / Initial Mass', fontsize = 24)
        plt.legend(loc = 'best', prop = {'size':20}, fancybox = True, framealpha = 0.5)
        plt.tight_layout()
        plt.savefig(eps_filename_curr)

        plt.close()



def time_series_plot(eps_filename, results_dir, times, field = 'density',
                     zlim = [1e-4, 1e8], nrows = 3, ncols = 3, smallplt = False,
                     colormap = 'hot', zoom = 1.0, xticks = None, yticks = None):
    """Make a time series plot of the density."""

    if not is_root():
        return

    import os

    if os.path.isfile(eps_filename):
        return

    if not os.path.exists(results_dir):
        return

    print "Generating plot with filename " + eps_filename

    directory = results_dir + '/output/'
    pltfiles = [directory + plt for plt in wdmerger.get_plotfiles(directory)]

    wdmerger.multipanel_slice_plot(field, eps_filename,
                                   [wdmerger.get_nearest_plotfile(pltfiles, time) for time in times],
                                   nrows = nrows, ncols = ncols,
                                   annotate_time = True,
                                   zlim = zlim,
                                   zoom = zoom,
                                   xticks = xticks,
                                   yticks = yticks,
                                   colormap = colormap
                                  )



# Main execution: do all of the analysis routines.

if __name__ == "__main__":
    """Generate the plots for the 3D mergers."""

    def results_dir():

        if results_base is not (None or ''):
            results_dir = results_base
        else:
            print "Error: results_base undefined"
            exit

        if mass_P is not (None or ''):
            results_dir += "/mass_P_" + mass_P

        if mass_S is not (None or ''):
            results_dir += "/mass_S_" + mass_S

        if roche is not (None or ''):
            results_dir += "/roche" + roche

        if rot is not (None or ''):
            results_dir += "/rot" + rot

        if hybrid is not (None or ''):
            results_dir += "/hybrid" + hybrid

        if ncell is not (None or ''):
            results_dir += "/n" + ncell

        return results_dir



    def plot_string(prefix, file_format = ".eps"):

        if plots_base is not (None or ''):
            plot_string = plots_base
        else:
            print "Error: plots_base undefined"
            exit

        if prefix is not (None or ''):
            plot_string += prefix
        else:
            print "Error: no prefix given to plot_string."
            exit

        if mass_P is not (None or ''):
            plot_string += "_m_P_" + mass_P

        if mass_S is not (None or ''):
            plot_string += "_m_S_" + mass_S

        if roche is not (None or ''):
            plot_string += "_roche_" + roche

        if rot is not (None or ''):
            plot_string += "_rot_" + rot

        if hybrid is not (None or ''):
            plot_string += "_hybrid_" + hybrid

        if ncell is not (None or ''):
            plot_string += "_n_" + ncell

        plot_string += file_format

        return plot_string



    vol_renders()

    results_base = 'results/approximate/'
    plots_base = 'plots/'

    mass_P = '0.90'
    mass_S = '0.60'

    for mass_P in ['0.90']:
        for mass_S in ['0.60', '0.90']:
            for roche in ['0.90', '1.00']:
                for rot in ['0', '1']:
                    for hybrid in ['0', '1']:

                        # Plots that combine data from multiple resolutions

                        ncell = ''

                        stellar_distance(plot_string("merger_distance"), results_dir())
                        stellar_distance_convergence(plot_string("merger_convergence"), results_dir())
                        stellar_masses(plot_string("merger_masses"), results_dir())

                        # Plots that we want at each resolution

                        for ncell in ['256', '512', '1024']:

                            max_time = 1.e20

                            if mass_P == '0.90' and mass_S == '0.60':
                                max_time = 150.0

                            gravitational_wave_signal(plot_string("merger_gw_signal"), results_dir(), max_time = max_time)



    ncell = '512'
    rot = ''
    hybrid = ''

    for roche in ['0.90', '1.00']:
        angular_momentum_comparison(plot_string("merger_angular_momentum_comparison"), results_dir(), ncell)
        angular_momentum_comparison(plot_string("merger_angular_momentum_comparison_minus_grid_losses"), results_dir(), ncell, subtract_grid_losses = True)

    # Create some time series plots. We need to construct these by hand in the sense that we want
    # control over the simulation times shown.

    mass_P = '0.90'
    mass_S = '0.60'

    roche = '0.90'
    ncell = '512'
    hybrid = '1'

    zoom = 1.0 / 0.75
    xticks = [-3.0e9, 0.0, 3.0e9]
    yticks = [-3.0e9, 0.0, 3.0e9]

    rot = '0'
    times = [0.0, 25.0, 50.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0]

    time_series_plot(plot_string("merger_time_series"), results_dir(), times, smallplt = True, zoom = zoom,
                     xticks = xticks, yticks = yticks, colormap = 'inferno')

    rot = '1'

    time_series_plot(plot_string("merger_time_series"), results_dir(), times, smallplt = True, zoom = zoom,
                     xticks = xticks, yticks = yticks, colormap = 'inferno')


    mass_S = '0.90'

    rot = '0'
    times = [0.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

    time_series_plot(plot_string("merger_time_series"), results_dir(), times, smallplt = True, zoom = zoom,
                     xticks = xticks, yticks = yticks, colormap = 'inferno')

    rot = '1'

    time_series_plot(plot_string("merger_time_series"), results_dir(), times, smallplt = True, zoom = zoom,
                     xticks = xticks, yticks = yticks, colormap = 'inferno')
