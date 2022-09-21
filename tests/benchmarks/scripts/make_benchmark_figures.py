import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

# Globals
ERROR_RATES = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
SEQ_LENGTHS = [1000, 2000, 3000, 4000, 5000, 7500, 10000, 12500, 15000, 17500,
               20000, 40000, 80000]


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

    fig, axs = plt.subplots(1, 2, sharey='row', figsize=(8, 6))

    for i, ax in enumerate(axs):
        for metric in metrics:
            ax.plot(
                xs[i] if i < len(xs) else xs[0],
                ys[i]['{0}_means'.format(metric['name'])],
                label=metric['name'] if i == 0 else "",
                marker='o'
            )
            ax.errorbar(
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
    xs = [ERROR_RATES, SEQ_LENGTHS]
    ys = [sub_stats, length_stats]
    metrics = [
        {'name': 'sensitivity', 'color': 'blue'},
        {'name': 'specificity', 'color': 'orange'}
    ]

    fig, axs = make_experiments_figure(xs, ys, metrics)

    # Theoritical Curve
    error_rate_theory_sensitivities = [0.9999999999999992, 0.9999999973133248,
                                       0.999999035776888, 0.9999633461915861,
                                       0.9996118046476143, 0.9985398391125027,
                                       0.9951159587019267, 0.9887248080755489,
                                       0.9786179973528009, 0.9625853664976465]
    error_rate_theory_specificities = [1.0, 1.0, 1.0, 0.9999999999630376,
                                       0.8999965025321311, 0.6881394058910561,
                                       0.2188627586201921,
                                       5.96189764223709e-14, 0.0, 0.0]
    query_size_theory_sensitivities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                       1.0, 1.0, 1.0, 1.0, 1.0]
    query_size_theory_specificities = [0.999999585864842, 0.9999999999999574,
                                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                       1.0, 1.0, 1.0]


    # Theoritical Curves
    axs[0].plot(
        ERROR_RATES,
        error_rate_theory_sensitivities,
        label='theoretical sensitivity'
    )
    axs[0].plot(
        ERROR_RATES,
        error_rate_theory_specificities,
        label='theoretical specificity'
    )

    axs[1].plot(
        SEQ_LENGTHS,
        query_size_theory_sensitivities[:11],
    )
    axs[1].plot(
        SEQ_LENGTHS,
        query_size_theory_specificities[:11],
    )

    # Substitution Rate
    # ax1.axhline(y=0.95, linestyle='--', color='grey')
    fig.suptitle('Flashmap accuracy on simulated human genome sequences')
    axs[0].axvline(x=0.05, linestyle='--', color='grey')
    axs[0].text(0.045, 0.1, '0.05 substitution rate threshold', rotation=90)
    axs[0].set_xlabel('Substitutions per site')
    # axs[1].set_ylim(ymin=0.7)

    # Query Length
    axs[1].axvline(x=5000, linestyle='--', color='grey')
    axs[1].text(5500, 0.3, '5kb query threshold', rotation=90)
    # axs[2].set_xscale('log')
    axs[1].set_xlabel('Query length (kb)')
    fig.legend(frameon=True)

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
    pan_trog = get_simulation_stats('pan_trog', 'query_length')
    gorilla = get_simulation_stats('gorilla', 'query_length')

    xs = [SEQ_LENGTHS]
    ys = [pan_trog, gorilla]
    metrics = [
        {'name': 'sensitivity', 'color': 'blue'},
        {'name': 'specificity', 'color': 'orange'}
    ]

    fig, axs = make_experiments_figure(xs, ys, metrics)
    fig.suptitle(
        ('Flashmap accuracy for mammalian genomes '
         'mapped to human reference hg38')
    )

    # Pan Trog
    axs[0].set_title('Chimpanzee')
    axs[0].axvline(x=5000, linestyle='--', color='grey')
    axs[0].text(5500, 0.2, '5kb query threshold', rotation=90)
    axs[0].set_xscale('log')
    axs[0].set_xlabel('Query length (kb)')

    # Gorilla
    # axs[1].set_ylim(ymin=0)
    axs[1].set_title('Gorilla')
    axs[1].text(5500, 0.2, '5kb query threshold', rotation=90)
    axs[1].axvline(x=5000, linestyle='--', color='grey')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Query length (kb)')

    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()


def test():
    # file = './metrics/sub_rate_995_32M_sensitivity.txt'
    # make_simulation_trials_figure()
    make_mammal_figure()

# test()


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
    ax1.axvline(x=0.05, linestyle='--', color='grey')
    ax1.text(0.045, 0.3, '0.05 substitution rate threshold', rotation=90)
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

def window_size_figure(figure_output=False):
    error_rates = []
    query_lengths = []
    window_sizes = []
    with open('./metrics/window_size_estimates.txt', 'r') as handle:
        for line in handle:
            split_line = line.rstrip().split(',')
            query_lengths.append(int(split_line[0]))
            error_rates.append(int(split_line[1]))
            window_sizes.append(float(split_line[2]))

    use_bigsi = ['bigsi filter' if window_size <= 2000 else 'mashmap'
                 for window_size in window_sizes]


    import pandas as pd

    #create DataFrame
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
    print(boundary_points)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(boundary_points.x, boundary_points.y)

    ax.fill_between(
        boundary_points.x, 94, boundary_points.y, 
        where=(
            (boundary_points.y >= 94) & (boundary_points.x >= 5000)
        ),
        label='Bigsi filter is required'
    )
    ax.fill_between(
        boundary_points.x, 100, boundary_points.y,
        label='MashMap only',
    )
    # groups = df.groupby('z')
    # for name, group in groups:
    #     ax.plot(group.x, group.y, marker='.', linestyle='', markersize=12,
    #              label=name)

    plt.title('Decision boundary for using BIGSI filter')
    plt.xlabel('Query Length')
    plt.ylabel('Percent Identity')
    plt.axhline(y=94, linestyle='--', color='grey')
    plt.axvline(x=5000, linestyle='--', color='grey')
    plt.xscale('log')
    plt.legend(loc='upper right', frameon=True)
    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()


window_size_figure()  # figure_output='figures/window_size_figure')
# make_read_figure()
# make_runtime_figure()
# make_simulation_trials_figure(10,
#                               figure_output='figures/simulation_trials_95_32M')
# make_mammal_figure(100)
# make_synth_figure()
# make_jaccard_test_figure(50)
