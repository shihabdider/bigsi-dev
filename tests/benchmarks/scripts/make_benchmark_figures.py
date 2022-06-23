import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
seq_lengths = [1000, 2000, 3000, 4000, 5000, 10000, 20000, 40000, 80000,
               160000, 200000, 250000, 300000]


def chunked(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_parameter_metrics(file, num_trials, parameter):
    parameter_metrics = []
    with open(file, 'r') as handle:
        for lst in chunked(handle.readlines(), num_trials):
            metrics = []
            for i, line in enumerate(lst):
                split_line = line.split(' ')
                metric = float(split_line[1])
                metrics.append(metric)
            parameter_metrics.append(metrics)

    return parameter_metrics


def get_read_metrics(file):
    metrics = []
    with open(file, 'r') as handle:
        for line in handle:
            metric = float(line.split(' ')[1].rstrip())
            metrics.append(metric)

    return metrics


def get_simulation_stats(benchmark_name, num_trials, benchmark_parameters):
    # Sensitivity
    benchmark_sensitivities = get_parameter_metrics(
        'metrics/{0}_sensitivities.txt'.format(benchmark_name), num_trials)

    benchmark_sensitivity_means = [np.mean(sensitivities) for sensitivities in
                                   benchmark_sensitivities]
    benchmark_sensitivity_stds = [np.std(sensitivities) for sensitivities in
                                  benchmark_sensitivities]

    benchmark_sensitivity_errors_upper = []
    for i, mean in enumerate(benchmark_sensitivity_means):
        if (mean + 2*benchmark_sensitivity_stds[i] > 1):
            error = 1 - mean
            benchmark_sensitivity_errors_upper.append(error)
        else:
            benchmark_sensitivity_errors_upper.append(benchmark_sensitivity_stds[i]*2)

    benchmark_sensitivity_errors_lower = [max((2*std, 0)) for std in
            benchmark_sensitivity_stds]

    # Specificity
    benchmark_specificities = get_parameter_metrics(
        'metrics/{0}_specificities.txt'.format(benchmark_name), num_trials)

    benchmark_specificity_means = [np.mean(specificities) for specificities in
                                   benchmark_specificities]
    benchmark_specificity_stds = [np.std(specificities) for specificities in
                                  benchmark_specificities]
    benchmark_specificity_errors = [2*std for std in 
                                    benchmark_specificity_stds]

    benchmark_stats = {
            'specificity_means': benchmark_specificity_means,
            'specificity_stds': benchmark_specificity_stds,
            'specificity_errors': benchmark_specificity_errors,
            'sensitivity_means': benchmark_sensitivity_means,
            'sensitivity_stds': benchmark_sensitivity_stds,
            'sensitivity_errors_upper': benchmark_sensitivity_errors_upper,
            'sensitivity_errors_lower': benchmark_sensitivity_errors_lower,
            }

    return benchmark_stats


def make_mammal_figure(num_trials, figure_output=False):
    pan_trog = get_simulation_stats('pan_trog', num_trials, seq_lengths)
    gorilla = get_simulation_stats('gorilla', num_trials, seq_lengths)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on mammalian genomes')

    # Pan Trog
    ax1.set_title('Chimpanzee')
    ax1.plot(seq_lengths, pan_trog['sensitivity_means'], label='sensitivity',
             marker='o')
    ax1.errorbar(seq_lengths, pan_trog['sensitivity_means'],
                 yerr=[pan_trog['sensitivity_errors_lower'],
                       pan_trog['sensitivity_errors_upper']],
                 fmt='-', color='blue')

    ax1.plot(seq_lengths, pan_trog['specificity_means'], label='specificity',
             marker='o')
    ax1.errorbar(seq_lengths, pan_trog['specificity_means'],
                 yerr=pan_trog['specificity_errors'], fmt='-', color='orange')
    ax1.axvline(x=5000, linestyle='--', color='grey')
    ax1.text(5500, 0.2, '5kb query threshold', rotation=90)
    ax1.set_xscale('log')
    ax1.set_xlabel('Query length (kb)')
    #ax1.set_ylim(ymin=0)

    # Gorilla
    ax2.set_title('Gorilla')
    ax2.plot(seq_lengths, gorilla['sensitivity_means'], label='sensitivity',
             marker='o')
    ax2.errorbar(seq_lengths, gorilla['sensitivity_means'],
                 yerr=[gorilla['sensitivity_errors_lower'],
                       gorilla['sensitivity_errors_upper']],
                 fmt='-', color='blue')

    ax2.plot(seq_lengths, gorilla['specificity_means'], label='specificity',
             marker='o')
    ax2.errorbar(seq_lengths, gorilla['specificity_means'],
                 yerr=gorilla['specificity_errors'], fmt='-', color='orange')
    ax2.text(5500, 0.2, '5kb query threshold', rotation=90)
    ax2.axvline(x=5000, linestyle='--', color='grey')
    ax2.set_xscale('log')
    ax2.set_xlabel('Query length (kb)')
    if figure_output:
        plt.savefig(figure_output)
    else:
        plt.show()


def make_simulation_trials_figure(num_trials, figure_output):
    subs = get_simulation_stats('sub_rate', num_trials, error_rates)
    length = get_simulation_stats('query_length', num_trials, seq_lengths)

    # Theoritical Curve
    # error_rate_theory_sensitivities = [0.9999999999999992, 0.9999999973133248, 
    #                                    0.999999035776888, 0.9999633461915861, 
    #                                    0.9996118046476143, 0.9985398391125027, 
    #                                    0.9951159587019267, 0.9887248080755489, 
    #                                    0.9786179973528009, 0.9625853664976465]
    # error_rate_theory_specificities = [1.0, 1.0, 1.0, 0.9999999999630376, 
    #                                    0.9999965025321311, 0.9881394058910561, 
    #                                    0.2188627586201921, 
    #                                    5.96189764223709e-14, 0.0, 0.0]
    # query_size_theory_sensitivities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
    #                                    1.0, 1.0, 1.0, 1.0, 1.0]
    # query_size_theory_specificities = [0.999999585864842, 0.9999999999999574, 
    #                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
    #                                    1.0, 1.0, 1.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on simulated data')

    # Sub errors
    ax1.plot(error_rates, subs['sensitivity_means'], label='sensitivity', 
             marker='o')
    # ax1.plot(error_rates, error_rate_theory_sensitivities,
    #          label='theoretical sensitivity')
    # ax1.plot(error_rates, error_rate_theory_specificities,
    #         label='theoretical specificity')
    ax1.errorbar(error_rates, subs['sensitivity_means'],
                 #yerr=subs_sensitivity_errors, 
                 yerr=[subs['sensitivity_errors_lower'],
                       subs['sensitivity_errors_upper']],
                 fmt='-', color='blue',)
    ax1.plot(error_rates, subs['specificity_means'], label='specificity',
             marker='o')
    ax1.errorbar(error_rates, subs['specificity_means'],
                 #yerr=[subs_specificity_errors_lower, subs_specificity_errors_upper], 
                 yerr=subs['specificity_errors'],
                 fmt='-', color='orange')
    #ax1.axhline(y=0.95, linestyle='--', color='grey')
    ax1.axvline(x=0.05, linestyle='--', color='grey')
    ax1.text(0.045, 0.1, '0.05 substitution rate threshold', rotation=90)
    ax1.set_xlabel('Substitutions per site')
    #ax1.set_ylim(ymin=0.7)
    fig.legend(frameon=True)

    # Query Length
    ax2.plot(seq_lengths, length['sensitivity_means'], label='sensitivity', 
             marker='o')
    # ax2.plot(seq_lengths, query_size_theory_sensitivities, 
    #          label='theoretical sensitivity')
    # ax2.plot(seq_lengths, query_size_theory_specificities,
    #          label='theoretical specificity')
    ax2.errorbar(seq_lengths, length['sensitivity_means'],
                 yerr=[length['sensitivity_errors_lower'],
                 length['sensitivity_errors_upper']],
                 #yerr=length_sensitivity_errors, 
                 fmt='-', color='blue')

    ax2.plot(seq_lengths, length['specificity_means'], label='specificity',
             marker='o')
    ax2.errorbar(seq_lengths, length['specificity_means'],
                 yerr=length['specificity_errors'], fmt='-', color='orange')
    ax2.axvline(x=5000, linestyle='--', color='grey')
    ax2.text(5500, 0.3, '5kb query threshold', rotation=90)
    ax2.set_xscale('log')
    ax2.set_xlabel('Query length (kb)')
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


def make_jaccard_test_figure(window_size, query_size, savefig=False):
    query_sizes = [5000, 10000, 20000, 40000, 80000, 160000, 300000]
    sub_rates = ['000', '001', '002', '003', '004', '005', '006', 
                 '007', '008', '009', '010']
    ref_size = 16

    fig, ax = plt.subplots(sharey='row')
    fig.suptitle(
        'Containment Score Deviation For {0}kb Query'.format(query_size/1000))
    x_axis_text = 'Substitution Rate'
    y_axis_text = 'True Jaccard Containment - Winnowed Jaccard Containment'
    fig.text(0.5, 0.04, x_axis_text, ha='center')
    fig.text(0.04, 0.5, y_axis_text, va='center', rotation='vertical')

    for query_size in query_sizes:
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

        ax.plot(sub_rates_floats, means, marker='o')
        ax.errorbar(sub_rates_floats, means,
                    yerr=[2*std for std in stds], 
                    fmt='-')
        ax.axhline(y=0)

    figure_output = 'figures/jaccard_deviation_{0}_{1}_w{2}.txt'.format(
                ref_size, query_size, window_size)
    if savefig:
        plt.savefig(figure_output)
    else:
        plt.show()

#make_read_figure()
#make_runtime_figure()
#make_simulation_trials_figure(100)
#make_mammal_figure(100)
#make_synth_figure()
make_jaccard_test_figure(100, 10000)
