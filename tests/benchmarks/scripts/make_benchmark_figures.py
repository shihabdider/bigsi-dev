import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch
plt.style.use('seaborn-whitegrid')

# Globals
ERROR_RATES = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
# SEQ_LENGTHS = [1000, 2000, 3000, 4000, 5000, 7500, 10000]
SEQ_LENGTHS = [1000, 2000, 3000, 4000, 5000, 7500, 10000, 12500, 15000, 17500,
               20000]


def make_simulation_df(file, parameter_name, metric_name):
    data = []
    with open(file, 'r') as handle:
        for line in handle:
            split_line = line.rstrip().split()
            experiment = int(split_line[0].split('/')[-2].split('_')[1])
            parameter = int(split_line[0].split('/')[-1].split('.')[0])
            metric = float(split_line[1])
            data_row = {
                'experiment number': experiment,
                parameter_name: parameter,
                metric_name: metric
            }

            data.append(data_row)

    df = pd.DataFrame(data)
    df = df.pivot(index='experiment number', columns=parameter_name)
    df.columns = df.columns.droplevel(0)
    return df


def compute_errorbars(means, stds):
    lower_errors = []
    upper_errors = []
    for i, mean in enumerate(means):
        lower_error = stds.iloc[i]*2
        if (mean - lower_error) < 0:
            lower_error = mean
        lower_errors.append(lower_error)

        upper_error = stds.iloc[i]*2
        if upper_error + mean > 1.0:
            upper_error = 1.0 - mean
        upper_errors.append(upper_error)

    errors = pd.DataFrame(
        {
            'lower': lower_errors,
            'upper': upper_errors
        }
    )

    return errors


def get_simulation_stats(benchmark_name, parameter_name):

    metrics = ['sensitivity', 'specificity']

    benchmark_stats = {}
    for metric in metrics:
        metrics = make_simulation_df(
            'metrics/{0}_{1}.txt'.format(benchmark_name, metric),
            parameter_name,
            metric
        )

        means = metrics.mean(axis=0)
        stds = metrics.std(axis=0)
        errors = compute_errorbars(means, stds)

        stat_names = ['means', 'stds', 'errors']
        stats = dict(zip(stat_names, [means, stds, errors]))
        for name in stat_names:
            benchmark_stats['{0}_{1}'.format(metric, name)] = stats[name]

    return benchmark_stats


def make_experiments_figure(xs, ys, metrics):
    '''
    Generates template figure for experiments

    args:
        xs: list of parameters (e.g substitution rates and query lengths)
        ys: list of benchmark statistics (e.g pan trog and gorilla)

    returns: figure and axis objects
    '''

    fig, axs = plt.subplots(1, len(ys), sharey='row', figsize=(12, 6))

    for i, ax in enumerate(axs):
        for metric in metrics:
            # xs = ys[i]['{0}_means'.format(metric['name'])].index
            ax.plot(
                # xs,
                xs[i] if i < len(xs) else xs[0],
                ys[i]['{0}_means'.format(metric['name'])],
                label=metric['name'] if i == 0 else "",
                marker='o'
            )
            ax.errorbar(
                # xs,
                xs[i] if i < len(xs) else xs[0],
                ys[i]['{0}_means'.format(metric['name'])],
                yerr=[
                    ys[i]['{0}_errors'.format(metric['name'])].lower,
                    ys[i]['{0}_errors'.format(metric['name'])].upper
                ],
                fmt='-',
                color=metric['color'],
            )

    return fig, axs


def make_simulation_trials_figure(figure_output=False):
    sub_stats = get_simulation_stats('sub_rate_95_32M', 'sub_rate')
    length_stats = get_simulation_stats('query_length_95_32M', 'query_length')

    sub_stats_5kb = get_simulation_stats('sub_rate_5000', 'sub_rate')
    length_stats_006 = get_simulation_stats('query_length_006', 'query_length')

    xs = [ERROR_RATES, SEQ_LENGTHS, ERROR_RATES, SEQ_LENGTHS]
    ys = [sub_stats, length_stats, sub_stats_5kb, length_stats_006]

    metrics = [
        {'name': 'sensitivity', 'color': 'blue'},
        {'name': 'specificity', 'color': 'orange'}
    ]

    fig, axs = make_experiments_figure(xs, ys, metrics)

    # Theoretical Accuracies

    theory_error_rate_df_10000 = pd.read_csv(
        './metrics/theory_error_rate_acc_10000.txt'
    )
    theory_error_rate_df_5000 = pd.read_csv(
        './metrics/theory_error_rate_acc_5000.txt'
    )

    theory_query_size_df_000 = pd.read_csv(
        './metrics/theory_query_size_acc_00.txt'
    )
    theory_query_size_df_006 = pd.read_csv(
        './metrics/theory_query_size_acc_006.txt'
    )

    # Theoretical Curves

    metrics = ['sensitivities', 'specificities']
    for metric in metrics:
        axs[0].plot(
            ERROR_RATES,
            theory_error_rate_df_10000[metric],
            label='theoretical {0}'.format(metric)
        )

        axs[1].plot(
            SEQ_LENGTHS,
            theory_query_size_df_000[metric],
        )

        axs[2].plot(
            ERROR_RATES,
            theory_error_rate_df_5000[metric],
        )

        axs[3].plot(
            SEQ_LENGTHS,
            theory_query_size_df_006[metric],
        )

    fig.suptitle('Flashmap accuracy on simulated human genome sequences')

    # Substitution Rate
    axs[0].set_title('10Kb query')
    axs[0].axvline(x=0.06, linestyle='--', color='grey')
    axs[0].text(0.045, 0.2, '0.06 substitution rate threshold', rotation=90)
    axs[0].set_xlabel('Substitutions per site')

    axs[2].set_title('5Kb query')
    axs[2].axvline(x=0.06, linestyle='--', color='grey')
    axs[2].text(0.045, 0.2, '0.06 substitution rate threshold', rotation=90)
    axs[2].set_xlabel('Substitutions per site')

    # Query Length
    axs[1].set_title('Exact matching query')
    axs[1].axvline(x=5000, linestyle='--', color='grey')
    axs[1].text(5500, 0.3, '5kb query threshold', rotation=90)
    # axs[1].set_xscale('log')
    axs[1].set_xlabel('Query length (kb)')
    axs[1].set_xticks([1000, 5000, 10000, 15000, 20000])

    axs[3].set_title('0.06 substitution rate query')
    axs[3].axvline(x=5000, linestyle='--', color='grey')
    axs[3].text(5500, 0.3, '5kb query threshold', rotation=90)
    # axs[3].set_xscale('log')
    axs[3].set_xlabel('Query length (kb)')
    axs[3].set_xticks([1000, 5000, 10000, 15000, 20000])

    fig.legend(loc='center right', frameon=False)

    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()


def get_read_metrics(file):
    metrics = []
    with open(file, 'r') as handle:
        for line in handle:
            metric = float(line.split(' ')[1].rstrip())
            metrics.append(metric)

    return metrics

def make_mammal_figure(figure_output=False):
    ys = []
    species = ['pan_trog', 'gorilla', 'bonobo']
    for specie in species:
        y = get_simulation_stats(specie, 'query_length')
        ys.append(y)

    xs = [SEQ_LENGTHS]
    metrics = [
        {'name': 'sensitivity', 'color': 'blue'},
        {'name': 'specificity', 'color': 'orange'}
    ]

    fig, axs = make_experiments_figure(xs, ys, metrics)
    fig.suptitle(
        ('Flashmap accuracy for mammalian genomes '
         'mapped to human reference hg38')
    )

    for ax in axs:
        ax.axvline(x=5000, linestyle='--', color='grey')
        ax.text(5500, 0.3, '5kb query threshold', rotation=90)
        ax.set_xscale('log')
        ax.set_xlabel('Query length (kb)')
        ax.set_xticks([1000, 2500, 5000, 10000, 20000])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    # Set titles
    axs[0].set_title('Chimpanzee')
    axs[1].set_title('Gorilla')
    axs[2].set_title('Bonobo')
    # axs[3].set_title('Basset')
    fig.legend()

    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()



def make_synth_figure():
    error_rates, error_sensitivities = get_error_metrics(
        'synth_seq_error_sensitivity.txt')

    _, error_specificities = get_error_metrics(
        'synth_seq_error_specificity.txt')

    seq_lengths, seq_length_sensitivities = get_length_metrics(
        'synth_seq_query_len_sensitivity.txt')

    _, seq_length_specificities = get_length_metrics(
        'synth_seq_query_len_specificity.txt')

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on synthetic data')

    # Sub errors
    ax1.plot(error_rates, error_sensitivities, label='sensitivity',
             marker='o')
    #ax1.errorbar(error_rates, subs_sensitivity_means,
    #             yerr=subs_sensitivity_errors, fmt='-', color='blue')
    ax1.plot(error_rates, error_specificities, label='specificity',
             marker='o')
    #ax1.errorbar(error_rates, subs_specificity_means,
    #             yerr=[subs_specificity_errors_lower,
    #                   subs_specificity_errors_upper], fmt='-', color='orange')
    ax1.axvline(x=0.06, linestyle='--', color='grey')
    ax1.text(0.045, 0.3, '0.06 substitution rate threshold', rotation=90)
    ax1.set_xlabel('Substitutions per site')
    fig.legend()

    ax2.plot(seq_lengths, seq_length_sensitivities, label='sensitivity',
             marker='o')
    #ax2.errorbar(seq_lengths, length_sensitivity_means,
    #             yerr=length_sensitivity_errors, fmt='-', color='blue')

    ax2.plot(seq_lengths, seq_length_specificities, label='specificity',
             marker='o')
    #ax2.errorbar(seq_lengths, length_specificity_means,
    #             yerr=length_specificity_errors, fmt='-', color='orange')
    ax2.axvline(x=5000, linestyle='--', color='grey')
    ax2.text(5500, 0.3, '5kb query threshold', rotation=90)
    ax2.set_xscale('log')
    ax2.set_xlabel('Query length (kb)')
    #ax2.legend()
    plt.show()
    #plt.savefig('flashmap_accuracy_synthetic_data.png')


def make_read_figure():
    # Nanopore reads
    nanopore_sensitivities = get_read_metrics(
        'metrics/nanopore_read_sensitivities_003.txt')

    nanopore_specificities = get_read_metrics(
        'metrics/nanopore_read_specificities_003.txt')

    
    # Pacbio reads
    pacbio_sensitivities = get_read_metrics(
        'metrics/pacbio_read_sensitivities_003.txt')

    pacbio_specificities = get_read_metrics(
        'metrics/pacbio_read_specificities_003.txt')

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on read datasets')

    read_dataset_names = ['Nanopore', 'Pacbio']
    # Sensitivity
    ax1.scatter([read_dataset_names[0]]*len(nanopore_sensitivities), nanopore_sensitivities, marker='o')
    ax1.scatter([read_dataset_names[1]]*len(pacbio_sensitivities), pacbio_sensitivities, marker='o')
    ax1.set_xlabel('Sensitivity')

    # Pacbio
    ax2.scatter([read_dataset_names[0]]*len(nanopore_specificities), nanopore_specificities, marker='o')
    ax2.scatter([read_dataset_names[1]]*len(pacbio_specificities), pacbio_specificities, marker='o')
    ax2.set_xlabel('Specificity')
    plt.show()
    #plt.savefig('figures/flashmap_reads_accuracy_003.png')


def make_runtime_figure():
    # Varying BIGSI/target size
    bigsi_sizes = [100, 400, 800, 1200, 1600, 2000, 2500, 3000]
    all_bigsi_runtimes = []
    num_trials = 100
    with open('metrics/synth_bigsi_size_performance_times.txt', 'r') as handle:
        for lst in chunked(handle.readlines(), num_trials):
            runtimes = []
            for line in lst:
                runtimes.append(float(line.rstrip()))

            all_bigsi_runtimes.append(runtimes)

    bigsi_runtime_means = [np.mean(runtimes) for runtimes in 
                           all_bigsi_runtimes]
    bigsi_runtime_stds = [np.std(runtimes) for runtimes in all_bigsi_runtimes]
    bigsi_runtime_errors = [2*std for std in bigsi_runtime_stds]

    # Varying query size
    query_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 300]
    all_query_runtimes = []
    num_trials = 100
    with open('metrics/synth_query_size_performance_times.txt', 'r') as handle:
        for lst in chunked(handle.readlines(), num_trials):
            runtimes = []
            for line in lst:
                runtimes.append(float(line.rstrip()))

            all_query_runtimes.append(runtimes)

    query_runtime_means = [np.mean(runtimes) for runtimes in 
                           all_query_runtimes]
    query_runtime_stds = [np.std(runtimes) for runtimes in all_query_runtimes]
    query_runtime_errors = [2*std for std in query_runtime_stds]


    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap Runtimes')

    # BIGSI Size
    ax1.plot(bigsi_sizes, bigsi_runtime_means, marker='o')
    ax1.errorbar(bigsi_sizes, bigsi_runtime_means, yerr=bigsi_runtime_errors,
                 fmt='-', color='blue')
    ax1.set_ylabel('Wall-clock Runtimes (s)')
    ax1.set_xlabel('BIGSI size (MB)')

    # Query size
    ax2.plot(query_sizes, query_runtime_means, marker='o')
    ax2.errorbar(query_sizes, query_runtime_means, yerr=query_runtime_errors,
                 fmt='-', color='orange')
    ax2.set_xlabel('Query size (kb)')
    #plt.show()
    plt.savefig('figures/flashmap_runtimes.png')


def make_jaccard_test_figure(window_size, savefig=False):
    #query_sizes = [5000, 10000, 20000, 40000, 80000, 160000, 300000]
    query_sizes = [5000, 20000, 80000, 300000]
    sub_rates = ['000', '001', '002', '003', '004', '005', '006',
                 '007', '008', '009', '010']
    ref_size = 16

    fig, axs = plt.subplots(len(query_sizes), sharey='row')
    plt.subplots_adjust(hspace=0.8)
    fig.suptitle('Deviation From True Containment Score')
    x_axis_text = 'Substitution Rate'
    y_axis_text = 'True Jaccard Containment - Winnowed Jaccard Containment'
    fig.text(0.5, 0.04, x_axis_text, ha='center')
    fig.text(0.04, 0.5, y_axis_text, va='center', rotation='vertical')

    for i, query_size in enumerate(query_sizes):
        jaccard_sizes_diffs = []
        for rate_index, rate in enumerate(sub_rates):
            filename = 'metrics/jaccard/{0}_{1}_{2}_w{3}.txt'.format(
                ref_size, query_size, rate, window_size)
            jaccard_diffs = []
            with open(filename, 'r') as handle:
                for line in handle:
                    jaccard_diffs.append(float(line.strip()))
            jaccard_sizes_diffs.append(jaccard_diffs)

        sub_rates_floats = [int(rate)/100 for rate in sub_rates]
        means = [np.mean(diffs) for diffs in jaccard_sizes_diffs]
        #mean = (np.mean(means))
        stds = [np.std(diffs) for diffs in jaccard_sizes_diffs]
        #upper_95 = mean + 2*np.mean(stds)
        #lower_95 = mean - 2*np.mean(stds)

        axs[i].set_title('{}kb Query'.format(int(query_size)/1000))
        #axs[i].set_ylim(ymin=-0.3, ymax=0.15)
        axs[i].plot(sub_rates_floats, means, marker='o')
        axs[i].errorbar(sub_rates_floats, means,
                        yerr=[2*std for std in stds],
                        fmt='-')
        axs[i].axhline(y=0)

    figure_output = 'figures/jaccard_deviation_{0}_{1}_w{2}.txt'.format(
                ref_size, query_size, window_size)
    if savefig:
        plt.savefig(figure_output)
    else:
        plt.show()

def make_window_size_boundary_points(genome_name):
    error_rates = []
    query_lengths = []
    window_sizes = []
    with open(
        './metrics/window_size_estimates_{0}.txt'.format(genome_name), 'r'
    ) as handle:
        for line in handle:
            split_line = line.rstrip().split(',')
            query_lengths.append(float(split_line[0]))
            error_rates.append(float(split_line[1]))
            window_sizes.append(float(split_line[2]))

    sketch_window_sizes = {
        'worm': 66,
        'hg38': 2000,
        'plant': 8000
    }

    use_bigsi = ['bigsi filter' 
                 if window_size <= sketch_window_sizes[genome_name] 
                 else 'mashmap'
                 for window_size in window_sizes]

    df = pd.DataFrame(
        {'x': query_lengths, 'y': error_rates, 'z': use_bigsi}
    )

    # To plot the line use plt.plot(x, y) with the right xs and the right ys
    boundary_points = []
    for i, length in enumerate(set(query_lengths)):
        row = df.loc[(df['x'] == length) & (df['z'] == 'bigsi filter')]
        if not row.empty:
            row = row[row['y'] == max(row['y'])]
            boundary_points.append(
                {
                    'x': float(row['x']), 
                    'y': float(row['y'])
                }
            )

    boundary_points = pd.DataFrame(boundary_points).sort_values(by='x')
    return boundary_points


def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p




def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax2 : the big main axes
    ax1 : the zoomed axes
    The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=2, loc2a=3, loc1b=1, loc2b=4, 
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def test():
    # file = './metrics/sub_rate_995_32M_sensitivity.txt'
    # make_simulation_trials_figure()
    t = np.arange(1000)
    x1 = np.sin(np.arange(1000)/10.)
    x2 = np.cos(np.arange(1000)/10.)

    plt.figure()
    ax1 = plt.subplot(2, 2, 1)
    plt.plot(t, x1)
    ax2 = plt.subplot(2, 2, 2)
    plt.plot(t, x2)

    ax3 = plt.subplot(2, 2, 3)
    plt.plot(t[100:200], x1[100:200])
    zoom_effect02(ax3, ax1)

    ax4 = plt.subplot(2, 2, 4)
    plt.plot(t[600:700], x2[600:700])
    zoom_effect02(ax4, ax2)
    plt.show()


# test()

def make_window_size_zoomed(main_axis, position, data, x_lim):
    ax = plt.subplot(position)
    plt.plot(
        data.x[:x_lim], 
        data.y[:x_lim]
    )
    ax.fill_between(
        data.x[:x_lim], 
        94, data.y[:x_lim], 
        where=(
            (data.y[:x_lim] >= 94) & 
            (data.x[:x_lim] >= 5000)
        ),
    )
    ax.fill_between(
        data.x[:x_lim], 
        100, data.y[:x_lim],
    )
    ax.axhline(y=94, linestyle='--', color='grey')
    ax.axvline(x=5000, linestyle='--', color='grey')
    zoom_effect02(ax, main_axis)


def window_size_figure(figure_output=False):
    genomes = ['worm', 'hg38', 'plant']
    window_size_boundary_points = {
    }
    for genome in genomes:
        boundary_points = make_window_size_boundary_points(genome)
        window_size_boundary_points[genome] = boundary_points

    plt.figure(figsize=(10, 6))
    axs = []
    for i, genome in enumerate(genomes):
        ax = plt.subplot(2, 3, i + 1)
        plt.plot(
                window_size_boundary_points[genome].x, 
                window_size_boundary_points[genome].y
        )
        ax.fill_between(
            window_size_boundary_points[genome].x, 
            94, window_size_boundary_points[genome].y, 
            where=(
                (window_size_boundary_points[genome].y >= 94) & 
                (window_size_boundary_points[genome].x >= 5000)
            ),
            label='Bigsi filter is required' if i == 0 else ""
        )
        ax.fill_between(
            window_size_boundary_points[genome].x, 
            100, window_size_boundary_points[genome].y,
            label='MashMap only' if i == 0 else "",
        )
        ax.axhline(y=94, linestyle='--', color='grey')
        ax.axvline(x=5000, linestyle='--', color='grey')
        axs.append(ax)

    make_window_size_zoomed(
        axs[1], 
        223, 
        window_size_boundary_points['hg38'], 
        20
    )

    make_window_size_zoomed(
        axs[2], 
        224, 
        window_size_boundary_points['plant'], 
        100
    )

    plt.suptitle(
        'Decision boundary for using BIGSI filter for different genome sizes'
    )

    plt.figtext(0.5, 0.04, 'Query Length', ha='center')
    plt.figtext(0.04, 0.5, 'Percent Identity', va='center',
                rotation='vertical')
    plt.figlegend()

    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()


# window_size_figure(figure_output='figures/window_size_figure')

# make_read_figure()
# make_runtime_figure()
make_simulation_trials_figure() #figure_output='figures/simulation_trials_95_32M')
# make_mammal_figure(figure_output='figures/mammal_accuracy')
# make_synth_figure()
# make_jaccard_test_figure(50)
