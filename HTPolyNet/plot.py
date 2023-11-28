"""

.. module:: plot
   :synopsis: provides plotting functionality
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from HTPolyNet.gromacs import *
from HTPolyNet.utils import *
from HTPolyNet.banner import banner
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import logging
import networkx as nx
from datetime import datetime
import yaml

logger = logging.getLogger(__name__)

# prevents "RuntimeError: main thread is not in main loop" tk bug
plt.switch_backend('agg')


def scatter(df, xcolumn, columns=[], outfile='plot.png', **kwargs):
    """scatter generic scatter plot generator

    :param df: dataframe containing data
    :type df: pd.DataFrame
    :param xcolumn: name of column holding x-data
    :type xcolumn: str
    :param columns: list of y-value columns to be plotted vs. x, defaults to []
    :type columns: list, optional
    :param outfile: name of output image file, defaults to 'plot.png'
    :type outfile: str, optional
    """
    logging.disable(logging.DEBUG)
    cmapname = kwargs.get('colormap', 'plasma')
    size = kwargs.get('size', (8, 6))
    yunits = kwargs.get('yunits', None)
    cmap = cm.get_cmap(cmapname)
    fig, ax = plt.subplots(1, 1, figsize=size)
    ax.set_xlabel(xcolumn)
    for n in columns:
        ax.scatter(df[xcolumn], df[n], label=n)
    plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    logging.disable(logging.NOTSET)


def trace(qty, edrs, outfile='plot.png', **kwargs):
    """trace generates a plot of the energy-like quantity named by 'qty' vs time by reading data from the list of edr files named in 'edrs'

    :param qty: name of energy-like quantity; must conform to menu generated by 'gmx energy'
    :type qty: str
    :param edrs: list of names of edr files to scan, in order
    :type edrs: list
    :param outfile: name of output image file, defaults to 'plot.png'
    :type outfile: str, optional
    :return: the list of average values
    :rtype: list
    """
    # disable debug-level logging and above since matplotlib has a lot of debug statements
    logging.disable(logging.DEBUG)
    df = pd.DataFrame()
    cmapname = kwargs.get('colormap', 'plasma')
    size = kwargs.get('size', (8, 6))
    yunits = kwargs.get('yunits', None)
    avgafter = kwargs.get('avgafter', 0)
    cmap = cm.get_cmap(cmapname)
    xshift = 0.0
    chkpt = []
    for edr in edrs:
        if not df.empty:
            xshift = df.tail(1).iloc[0, 0]
        data = gmx_energy_trace(edr, [qty], xshift=xshift)
        lastchkpt = 0
        if len(chkpt) > 0:
            lastchkpt = chkpt[-1]
        chkpt.append(data.shape[0] + lastchkpt)
        df = pd.concat((df, data), ignore_index=True)
    fig, ax = plt.subplots(1, 1, figsize=size)
    nseg = len(chkpt)
    beg = 0
    avg = []
    for c in df.columns[1:]:
        for seg in range(nseg):
            ax.plot(
                df.iloc[beg:chkpt[seg], 0],
                df[c].iloc[beg:chkpt[seg]],
                label=(c if seg == 0 else None),
                color=cmap(seg / nseg)
            )
            beg = chkpt[seg]
        if 'avgafter' in kwargs:
            if avgafter > 0:
                pass
            else:
                avgafter = df['time(ps)'].iloc[-1] / 2
            sdf = df[df['time(ps)'] > avgafter]
            avg.append(sdf[c].mean())
            ax.plot(
                df.iloc[:, 0], [avg] * df.shape[0],
                'k-',
                alpha=0.3,
                label=f'{avg:0.2f}'
            )
        else:
            avg.append(df[c].mean())
    if not yunits:
        plt.ylabel(qty)
    else:
        plt.ylabel(f'{qty} ({yunits})')
    plt.xlabel('time(ps)')
    plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)
    return avg


def multi_trace(
    dfL,
    xnames,
    ynames,
    labels=[],
    xlabel='time [ps]',
    ylabel='',
    outfile='plot.png',
    **kwargs
):
    """multi_trace generates a plot of each y vs x in df

    :param df: a pandas dataframe
    :type df: pandas.DataFrame
    :param xnames: list of x-column names
    :type xnames: list
    :param ynames: list of y-column names, parallel to xnames
    :type ynames: list
    :param outfile: name of output image file, defaults to 'plot.png'
    :type outfile: str, optional
    """
    # disable debug-level logging and above since matplotlib has a lot of debug statements
    default_units = {
        'Temperature': 'K',
        'Pressure': 'bar',
        'Density': 'kg/m^3',
        'Potential': 'kJ/mol'
    }
    units = kwargs.get('units', default_units)
    logging.disable(logging.DEBUG)
    size = kwargs.get('size', (16, 4))
    legend = kwargs.get('legend', True)
    fig, ax = plt.subplots(1, 1, figsize=size)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cmapname = kwargs.get('colormap', 'plasma')
    cmap = cm.get_cmap(cmapname)

    ndatasets = len(xnames)
    assert ndatasets == len(ynames)

    for i, (df, x, y, l) in enumerate(zip(dfL, xnames, ynames, labels)):
        ax.plot(df[x], df[y], label=l, color=cmap(i / ndatasets))

    if legend:
        plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)


def global_trace(
    df,
    names,
    outfile='plot.png',
    transition_times=[],
    markers=[],
    interval_labels=[],
    y2names=[],
    **kwargs
):
    """global_trace generates custom-formatted multiplots of energy-like quantities named in 'names' in the input dataframe df

    :param df: pandas dataframe containing all data
    :type df: pd.DataFrame
    :param names: list of quantity names (Density, Temperature, etc)
    :type names: list
    :param outfile: name of output image file, defaults to 'plot.png'
    :type outfile: str, optional
    :param transition_times: time values at which vertical lines are drawn, defaults to []
    :type transition_times: list, optional
    :param markers: time values at which transitions are marked somehow, I forget, defaults to []
    :type markers: list, optional
    :param interval_labels: list of labels of intervals defined by markers, defaults to []
    :type interval_labels: list, optional
    :param y2names: names of quantities to be plotted on a secondary y axis, defaults to []
    :type y2names: list, optional
    """
    # disable debug-level logging and above since matplotlib has a lot of debug statements
    default_units = {
        'Temperature': 'K',
        'Pressure': 'bar',
        'Density': 'kg/m^3',
        'Potential': 'kJ/mol'
    }
    units = kwargs.get('units', default_units)
    logging.disable(logging.DEBUG)
    size = kwargs.get('size', (16, 4 * len(names)))
    legend = kwargs.get('legend', False)
    fig, ax = plt.subplots(len(names), 1, figsize=size)
    plt.xlabel('time(ps)')
    cmapname = kwargs.get('colormap', 'plasma')
    # yunits=kwargs.get('yunits',None)
    cmap = cm.get_cmap(cmapname)
    # print(f'in global_trace:\n{df.head().to_string()}')

    interval_times = []
    if interval_labels:
        for i in range(1, len(transition_times)):
            interval_times.append(
                (transition_times[i] + transition_times[i - 1]) / 2
            )
    for l, t in zip(interval_labels, interval_times):
        logger.info(f'{t} {l}')
    assert len(interval_labels) == len(interval_times)
    L, R = -1, -1
    if len(markers) > 1:
        L, R = markers[0], markers[-1]
        in_tt = [x for x in transition_times if L < x < R]
        marked_df = df[(df['time(ps)'] > L) & (df['time(ps)'] < R)]
        fig, ax = plt.subplots(len(names) * 2, 1, figsize=size)
        for i, colname in enumerate(names):
            out_ax = ax[0] if len(names) == 1 else ax[i * 2]
            in_ax = ax[1] if len(names) == 1 else ax[i * 2 + 1]
            out_ax.plot(df['time(ps)'], df[colname], label=colname)
            ylabel = colname
            if ylabel in units:
                ylabel += f' ({units[ylabel]})'
            out_ax.set_ylabel(ylabel)
            out_ax.set_xlabel('time(ps)')
            if len(y2names) > i:
                out_ax2 = out_ax.twinx()
                out_ax2.plot(
                    df['time(ps)'],
                    df[y2names[i]],
                    label=y2names[i],
                    color='black'
                )
                ylabel = y2names[i]
                if ylabel in units:
                    ylabel += f' ({units[ylabel]})'
                out_ax2.set_ylabel(ylabel)
                out_ax2.set_xlabel('time(ps)')
            if len(transition_times) > 0:
                colors = [
                    cmap(i / len(transition_times))
                    for i in range(len(transition_times))
                ]
                ylim = out_ax.get_ylim()
                out_ax.vlines(
                    transition_times,
                    ylim[0],
                    ylim[1],
                    color=colors,
                    linewidth=0.75,
                    alpha=0.5
                )
                # for x,l in zip(interval_times,interval_labels):
                #     out_ax.text(x,0.9*ylim[1],l,fontsize=8)
            in_ax.plot(marked_df['time(ps)'], marked_df[colname], label=colname)
            ylabel = colname
            if ylabel in units:
                ylabel += f' ({units[ylabel]})'
            in_ax.set_xlabel('time(ps)')
            in_ax.set_yabel(ylabel)
            if len(y2names) > i:
                in_ax2 = in_ax.twinx()
                in_ax2.plot(
                    marked_df['time(ps)'],
                    marked_df[y2names[i]],
                    label=y2names[i],
                    color='black'
                )
                ylabel = y2names[i]
                if ylabel in units:
                    ylabel += f' ({units[ylabel]})'
                in_ax2.set_ylabel(ylabel)
                in_ax2.set_xlabel('time(ps)')
            if len(transition_times) > 0:
                colors = [
                    cmap(i / len(transition_times))
                    for i in range(len(transition_times))
                ]
                ylim = in_ax.get_ylim()
                in_ax.vlines(
                    in_tt,
                    ylim[0],
                    ylim[1],
                    color=colors,
                    linewidth=0.75,
                    alpha=0.5
                )
                # for x,l in zip(interval_times,interval_labels):
                #     if L<x<R:
                #         out_ax.text(x,0.9*ylim[1],l,fontsize=8)
    else:
        fig, ax = plt.subplots(len(names), 1, figsize=size)
        for i, colname in enumerate(names):
            the_ax = ax if len(names) == 1 else ax[i]
            the_ax.plot(df['time(ps)'], df[colname], label=colname)
            ylabel = colname
            if ylabel in units:
                ylabel += f' ({units[ylabel]})'
            the_ax.set_ylabel(ylabel)
            the_ax.set_xlabel('time(ps)')
            if len(y2names) > i:
                the_ax2 = the_ax.twinx()
                the_ax2.plot(
                    df['time(ps)'],
                    df[y2names[i]],
                    label=y2names[i],
                    color='black'
                )
                ylabel = y2names[i]
                if ylabel in units:
                    ylabel += f' ({units[ylabel]})'
                the_ax2.set_ylabel(ylabel)
            if len(transition_times) > 0:
                colors = [
                    cmap(i / len(transition_times))
                    for i in range(len(transition_times))
                ]
                ylim = the_ax.get_ylim()
                the_ax.vlines(
                    transition_times,
                    ylim[0],
                    ylim[1],
                    color=colors,
                    linewidth=0.5,
                    alpha=0.5
                )

    # plt.xlabel('time(ps)')
    if legend:
        plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)


def network_graph(G, filename, **kwargs):
    """network_graph draws a custom formatted network plot from graph G

    :param G: a graph from networkx
    :type G: nx.Graph or nx.DiGraph
    :param filename: name of output image filename
    :type filename: str
    """
    logging.disable(logging.DEBUG)
    arrows = kwargs.get('arrows', False)
    figsize = kwargs.get('figsize', (32, 32))
    node_size = kwargs.get('node_size', 200)
    with_labels = kwargs.get('with_labels', False)
    cmap = cm.get_cmap('seismic')
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.axis('off')
    molnames = list(
        set([n.get('molecule_name', 'anonymous') for k, n in G.nodes.items()])
    )
    nmolname = len(molnames)
    cx = []
    for n in G.nodes.values():
        idx = molnames.index(n.get('molecule_name', 'anonymous'))
        cx.append((float(idx) + 0.25) / (nmolname + 1))
    nx.draw_networkx(
        G,
        ax=ax,
        arrows=arrows,
        node_size=node_size,
        node_color=cx,
        cmap=cmap,
        with_labels=with_labels
    )
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)


# below are representive diagnostic output lines to establish extraction patterns
_template_1 = '2022-08-11 17:40:36,969 HTPolyNet.runtime.my_logger INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********'
_template_1_token_idx = [2, 3, 5, 7]
_template_2 = '2022-09-03 19:32:46,830 HTPolyNet.curecontroller.do_iter INFO> Iteration 1 current conversion 0.283 or 1082 bonds'
_template_2_token_idx = [2, 3, 6, 7]
_template_2_data_idx = {
    'iter': (int, 5),
    'conv': (float, 8),
    'nbonds': (int, 10)
}


def _token_match(l, template, pat_idx):
    """_token_match returns True if tokens indexed by pat_idx in the space-split l and template match

    :param l: probe string
    :type l: str
    :param template: template string
    :type template: str
    :param pat_idx: list of token indices
    :type pat_idx: list
    :return: True if tokens in l match those in template
    :rtype: bool
    """
    if len(l.split()) != len(template.split()):
        return
    return all([l.split()[t] == template.split()[t] for t in pat_idx])


def _parse_data(dat, l, idx_dict):
    """_parse_data parses data in a matching line

    :param dat: pandas dataframe
    :type dat: pd.DataFame
    :param l: line
    :type l: str
    :param idx_dict: dictionary of column-name:(type-converter,token-index)
    :type idx_dict: dict
    """
    tok = l.split()
    for k, v in idx_dict.items():
        conv, s = v
        dat[k].append(conv(tok[s]))


def diagnostics_graphs(logfiles, filename, **kwargs):
    """diagnostics_graphs extracts selected data from the diagnostic output and generates plots

    :param logfiles: list of names of diagnostic log files to treat in parallel
    :type logfiles: list
    :param filename: name of output image file
    :type filename: str
    """
    xmax = kwargs.get('xmax', -1)
    figsize = kwargs.get('figsize', (12, 6))
    logging.disable(logging.DEBUG)
    df: dict[pd.DataFrame] = {}
    for logfile in logfiles:
        bn, ex = os.path.splitext(logfile)
        with open(logfile, 'r') as f:
            lines = f.read().split('\n')
        logger.info(f'read {len(lines)} lines from {logfile}')
        data = {}
        data['time'] = []
        data['iter'] = []
        data['conv'] = []
        data['nbonds'] = []
        counter = 0
        for l in lines:
            if _token_match(l, _template_1, _template_1_token_idx):
                # print('you should only see this once')
                counter += 1
                assert not counter > 1
                data['time'].append(
                    datetime.strptime(
                        ' '.join(l.split()[0:2]), '%Y-%m-%d %H:%M:%S,%f'
                    )
                )
                data['iter'].append(0)
                data['conv'].append(0.0)
                data['nbonds'].append(0)
            elif _token_match(l, _template_2, _template_2_token_idx):
                data['time'].append(
                    datetime.strptime(
                        ' '.join(l.split()[0:2]), '%Y-%m-%d %H:%M:%S,%f'
                    )
                )
                # print('data tok',f'{l.split()}')
                _parse_data(data, l, _template_2_data_idx)
        # print('data',f'{data}')
        df[logfile] = pd.DataFrame(data)
        time_idx = list(df[logfile].columns).index('time')
        df[logfile]['elapsed'] = (
            df[logfile]['time'] - df[logfile].iloc[0, time_idx]
        ).astype(int) / 1.e9 / 3600.0
        df[logfile].to_csv(f'{bn}.csv', index=False, sep=' ', header=True)
    fig, ax = plt.subplots(1, 2, sharex=True, figsize=figsize)
    ax[0].set_ylim([0, 1])
    ax[0].set_xlabel('runtime (h)')
    ax[0].set_ylabel('conversion')
    ax[1].set_xlabel('runtime (h)')
    ax[1].set_ylabel('iteration')
    if xmax > -1:
        ax[0].set_xlim([0, xmax])
    for logfile in logfiles:
        ax[0].plot(df[logfile]['elapsed'], df[logfile]['conv'], label=logfile)
        ax[1].plot(df[logfile]['elapsed'], df[logfile].index + 1, label=logfile)
    plt.legend()
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)


def init_molecule_graph(proj_dir):
    """init_molecule_graph creates and initializes a inter-molecular graph to show network connectivity

    :param proj_dir: name of project directory
    :type proj_dir: str
    :return: a nodes-only Graph enumerating all molecules; this will be further processed elsewhere to add connectivity information
    :rtype: networkx.Graph
    """
    gro = os.path.join(proj_dir, 'systems/init/init.gro')
    top = os.path.join(proj_dir, 'systems/init/init.top')
    grx = os.path.join(proj_dir, 'systems/init/init.grx')
    TC = TopoCoord(grofilename=gro, topfilename=top, grxfilename=grx)
    G = nx.Graph()
    adf = TC.Coordinates.A
    mm = set(zip(adf['molecule'], adf['molecule_name']))
    for mx, mn in mm:
        G.add_node(mx, molecule_name=mn)
    return G


def plots(args):
    """plots handles the plots subcommand

    :param args: command-line arguments
    :type args: argparse.Namespace
    """
    loglevel_numeric = getattr(logging, args.loglevel.upper())
    logging.basicConfig(
        format='%(levelname)s> %(message)s', level=loglevel_numeric
    )
    if not args.no_banner:
        banner(logger.info)

    if args.source == 'build':
        build_plots(args)
    elif args.source == 'diag':
        diag_plots(args)
    elif args.source == 'post':
        post_plots(args)
    else:
        logger.error(f'Source {args.source} is not recognized.')


def diag_plots(args):
    diags = args.diags
    plotfile = args.plotfile
    if not plotfile:
        plotfile = 'cure_info.png'
    if len(diags) > 0:
        diagnostics_graphs(diags, plotfile)


def build_plots(args):
    GmxNames = {'t': 'Temperature', 'd': 'Density', 'p': 'Potential'}
    plot_types = args.buildplot
    trace_types = []
    if 't' in plot_types:
        trace_types = [GmxNames[i] for i in args.traces]
    for p in args.proj:
        if trace_types:
            df, transition_times, cure_markers, interval_labels = density_evolution(
                p
            )
            global_trace(
                df,
                trace_types,
                os.path.join(p, 'buildtraces.png'),
                transition_times=transition_times,
                markers=[],
                interval_labels=interval_labels,
                y2names=['nbonds', 'nbonds'],
                legend=True
            )
            df.to_csv(
                os.path.join(p, 'buildtraces.csv'),
                index=False,
                header=True,
                float_format='{:.3f}'.format
            )
        if any([x in plot_types for x in 'gnc']):
            G = init_molecule_graph(p)
            n = 1
            while os.path.exists(
                os.path.join(p, f'systems/iter-{n}/2-cure_update-bonds.csv')
            ):
                logger.info(f'iter-{n}/2-cure_update-bonds.csv')
                g = graph_from_bondsfile(
                    os.path.join(p, f'systems/iter-{n}/2-cure_update-bonds.csv')
                )
                G = nx.compose(G, g)
                if 'g' in plot_types:
                    network_graph(
                        G, os.path.join(p, f'plots/iter-{n}-graph.png')
                    )
                n += 1
        if 'g' in plot_types:
            network_graph(G, os.path.join(p, 'graph.png'))
        if 'c' in plot_types:
            clu = clusters(G)
            clu.to_csv(
                os.path.join(p, 'clusters.csv'),
                sep=' ',
                header=True,
                index=False
            )
            logger.info(f'{os.path.join(p,"clusters.csv")} created.')
            # cluster_plot(clu,os.path.join(p,"clusters.png"))
        if 'n' in plot_types:
            am = mwbxl(G)
            logger.info(
                f'Avg homo-N between xlinks: {np.average(am["n"],weights=am["counts"]):.2f}'
            )
            am.to_csv(
                os.path.join(p, 'dist_bw_xlinks.csv'),
                sep=' ',
                index=False,
                header=True
            )
            logger.info(f'{os.path.join(p,"dist_bw_xlinks.csv")} created.')
            # dist_bw_xlinks_plot(am,os.path.join(p,"dist_bw_xlinks.png"))


def do_tg_plots(
    phases, projdirs, outfile='tg.png', save_data='data.csv', n_points=[10, 20]
):
    res = {}
    means = {}
    stds = {}
    rate = []
    Tgs = []
    nproj = len(projdirs)
    for i, phase in enumerate(phases):
        p = phase['ladder']
        rate.append(p['deltaT'] / (p['ps_per_rise'] + p['ps_per_run']))
        res[i] = []
        for d in projdirs:
            df = pd.read_csv(
                os.path.join(d, p['subdir'], 'ladder.csv'),
                index_col=None,
                header=0
            )
            t0 = p['warmup_ps']
            ps_per_step = p['ps_per_rise'] + p['ps_per_run']
            final_ps = df['time(ps)'].iloc[-1]
            curr_ps = t0
            T = []
            rho = []
            V = []
            Tstd = []
            rhostd = []
            Vstd = []
            while curr_ps < final_ps:
                ll = curr_ps + p['ps_per_rise'] + 0.5 * p['ps_per_run']
                ul = ll + 0.5 * p['ps_per_run']
                tdf = df[(df['time(ps)'] >= ll) & (df['time(ps)'] <= ul)]
                T.append(tdf['Temperature'].mean())
                rho.append(tdf['Density'].mean())
                V.append(tdf['Volume'].mean())
                Tstd.append(tdf['Temperature'].std())
                rhostd.append(tdf['Density'].std())
                Vstd.append(tdf['Volume'].std())
                curr_ps += ps_per_step

            res[i].append(
                pd.DataFrame(
                    {
                        'Temperature': T,
                        'Volume': V,
                        'Density': rho,
                        'Temperature-std': Tstd,
                        'Volume-std': Vstd,
                        'Density-std': rhostd
                    }
                )
            )

        df_concat0 = pd.concat(res[i], axis=1)
        means[i] = df_concat0.stack().groupby(level=[0, 1]).mean().unstack()
        stds[i] = df_concat0.stack().groupby(level=[0, 1]).std().unstack()
        means[i]['Volume-std'] = stds[i]['Volume']
        means[i]['Density-std'] = stds[i]['Density']

    fig, ax = plt.subplots(1, 2, figsize=(10, 6), sharex=True, sharey=True)
    for i, v in enumerate(means.keys()):
        m = means[v]
        std = stds[v]
        if np.isnan(np.sum(std['Density'])):
            std['Density'] = np.zeros(len(std['Density']))
        m.sort_values('Temperature', axis=0, inplace=True)
        ax[i].set_xlabel('Temperature [K]')
        ax[i].set_ylabel('Density [kg/m$^3$]')
        if nproj < 2:
            ax[i].scatter(m['Temperature'], m['Density'])
        else:
            ax[i].errorbar(m['Temperature'], m['Density'], std['Density'])
        Tg, c, h = compute_tg(m['Temperature'], m['Density'], n_points=n_points)
        Tgs.append(Tg)
        if Tg != -1:
            ax[i].plot(
                m['Temperature'],
                c[0] * m['Temperature'] + c[1],
                color='blue',
                alpha=0.5
            )
            ax[i].plot(
                m['Temperature'],
                h[0] * m['Temperature'] + h[1],
                color='red',
                alpha=0.5
            )
            ax[i].scatter([Tg], [c[0] * Tg + c[1]], marker='o', color='black')
            ax[i].text(
                Tg, c[0] * Tg + c[1], f'   {Tg:.2f} K', verticalalignment='top'
            )
            means[v]['glassy-line-Density'] = c[0] * m['Temperature'] + c[1]
            means[v]['rubbery-line-Density'] = h[0] * m['Temperature'] + h[1]
        means[v].to_csv(
            f'ladder{v}-{save_data}', index=False, header=True, sep=' '
        )
        logger.info(f'ladder{v}-{save_data} created.')
    plt.savefig(outfile, bbox_inches='tight')
    plt.close(fig)
    hTg = Tgs[0]
    cTg = Tgs[1]
    logger.info(f'{outfile} created.')
    logger.info(
        f'heating Tg = {hTg:.2f} K ({(hTg-273.15):.2f} C) at {rate[0]:.5f} K/ps ({rate[0]*1.e12:.3e} K/s)'
    )
    logger.info(
        f'cooling Tg = {cTg:.2f} K ({(cTg-273.15):.2f} C) at {rate[1]:.5f} K/ps ({rate[1]*1.e12:.3e} K/s)'
    )


def do_E_plots(
    phases, projdirs, outfile='e.png', fit_domain=[10, 200], save_data='E.csv'
):
    MPa_per_bar = 1.e-1
    # average over replicas and directions (here, phases)
    all_stress_strains = []
    for p in phases:
        params = p['deform']
        dir = params['direction']
        # Box-X-strain,Pres-XX-stress
        strain_name = f'Box-{dir.upper()}-strain'
        stress_name = f'Pres-{dir.upper()}{dir.upper()}-stress'
        for d in projdirs:
            df = pd.read_csv(
                os.path.join(d, params['subdir'], f'deform-{dir}.csv'),
                index_col=None,
                header=0
            )
            all_stress_strains.append(
                pd.DataFrame(
                    {
                        'strain': df[strain_name],
                        'stress': (df[stress_name] * MPa_per_bar)
                    }
                )
            )
    df_concat = pd.concat(all_stress_strains, axis=1)
    mean_stress_strains = df_concat.stack().groupby(level=[0, 1]
                                                   ).mean().unstack()
    stds_stress_strains = df_concat.stack().groupby(level=[0, 1]).std().unstack()
    mean_stress_strains['stress-std'] = stds_stress_strains['stress']
    mean_stress_strains.to_csv(save_data, header=True, index=False, sep=' ')
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.set_ylabel('Stress, [MPa]')
    ax.set_xlabel('Strain [*]')
    ax.errorbar(
        mean_stress_strains['strain'],
        mean_stress_strains['stress'],
        stds_stress_strains['stress'],
        alpha=0.2
    )
    ax.plot(mean_stress_strains['strain'], mean_stress_strains['stress'])
    E, R2 = compute_E(
        mean_stress_strains['strain'],
        mean_stress_strains['stress'],
        fit_domain=fit_domain
    )
    fitline = E * mean_stress_strains['strain']
    half_domain = int(mean_stress_strains['strain'].shape[0] / 2)
    line_domain = [
        0, fit_domain[1] if fit_domain[1] > half_domain else half_domain
    ]
    X = np.array(mean_stress_strains['strain'])[line_domain[0]:line_domain[1]]
    Y = np.array(fitline)[line_domain[0]:line_domain[1]]
    ax.plot(X, Y, 'k--', alpha=0.7)
    plt.savefig(outfile, bbox_inches='tight')
    plt.close(fig)
    logger.info(
        f'{outfile} and {save_data} created. E = {E/1000.0:.3f} GPa (R^2 {R2:.3f})'
    )


def post_plots(args):
    n_points = args.n_points
    phases = []
    # print(args.cfg)
    for c in args.cfg:
        with open(c, 'r') as f:
            phases += yaml.safe_load(f)
    phasenames = [list(x.keys())[0] for x in phases]
    # print(phasenames)
    mdf = []
    for p in args.proj:
        if 'anneal' in phasenames and 'equilibrate' in phasenames:
            mdf.append(postsim_density_evolution(p))
    if mdf:
        multi_trace(
            mdf,
            xnames=['time(ps)'] * len(mdf),
            ynames=['Density'] * len(mdf),
            labels=args.proj,
            ylabel='Density [kg/m$^3$]',
            outfile='anneal-equil-density.png'
        )
        logger.info('-'.join(phasenames) + '-density.png created.')
        m = []
        for d in mdf:
            seg = d.iloc[int(0.9 * d.shape[0]):]
            logger.debug(
                f'averaging {seg.shape[0]} final density values out of {d.shape[0]}'
            )
            m.append(seg['Density'].mean())
        m = np.array(m)
        logger.info(
            f'mean density {m.mean():.0f} kg/m^3 ({m.mean()/1000.0:.3f} g/cc)'
        )

    ladder_phases = [
        idx for idx, value in enumerate(phasenames) if value == 'ladder'
    ]
    # print(ladder_phases)
    if len(ladder_phases) > 0:    # perform Tg calculations on each phase
        do_tg_plots(
            [phases[i] for i in ladder_phases], args.proj, n_points=n_points
        )

    deform_phases = [
        idx for idx, value in enumerate(phasenames) if value == 'deform'
    ]
    if len(deform_phases) > 0:
        do_E_plots([phases[i] for i in deform_phases], args.proj)
