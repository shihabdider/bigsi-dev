import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

#error_rates = []
#with open('error_rate_sensitivities.txt', 'r') as handle:
#    for line in handle:
#        split_line = line.split(' ')
#        error_rate = int(split_line[0][13:15])/100
#        error_rates.append(error_rate)

#seq_lengths = []
#seq_len_sensitivities = []
#with open('seq_length_sensitivities.txt', 'r') as handle:
#    for line in handle:
#        split_line = line.split(' ')
#        seq_length = int(split_line[0][12:-13])
#        print(seq_length)
#        seq_lengths.append(seq_length)
#
#        sensitivity = float(split_line[1])
#        seq_len_sensitivities.append(sensitivity)

#seq_len_specificities = []
#with open('seq_length_specificities.txt', 'r') as handle:
#    for line in handle:
#        split_line = line.split(' ')
#        specificity = float(split_line[1])
#        seq_len_specificities.append(specificity)

def chunked(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


#seq_lengths = [1000, 2000, 3000, 4000, 5000, 7000, 10000, 20000, 40000, 80000, 
#               160000, 200000, 250000, 300000]
#trial_sensitivities = []
#with open('seq_length_trials_sensitivity.txt', 'r') as handle:
#    for lst in chunked(handle.readlines(), 20):
#        seq_len_sensitivities = []
#        for i, line in enumerate(lst):
#            ignore_trial = (i >= 2 and i <= 5)
#            #ignore_trial = False
#            if not ignore_trial:
#                split_line = line.split(' ')
#                #seq_length = int(split_line[0][26:-11])
#                #print(seq_length)
#                #if seq_length not in seq_lengths:
#                #    seq_lengths.append(seq_length)
#
#                sensitivity = float(split_line[1])
#                seq_len_sensitivities.append(sensitivity)
#        trial_sensitivities.append(seq_len_sensitivities)
#
#trial_specificities = []
#with open('seq_length_trials_specificity.txt', 'r') as handle:
#    for lst in chunked(handle.readlines(), 20):
#        seq_len_specificities = []
#        for i, line in enumerate(lst):
#            ignore_trial = i >= 2 and i <= 5
#            #ignore_trial = False
#            if not ignore_trial:
#                split_line = line.split(' ')
#                specificity = float(split_line[1])
#                seq_len_specificities.append(specificity)
#        trial_specificities.append(seq_len_specificities)
#
#length_sensitivity_means = [np.mean(sensitivities) for sensitivities in 
#                            trial_sensitivities]
#length_sensitivity_stds = [np.std(sensitivities) for sensitivities in 
#                           trial_sensitivities]
##length_sensitivity_errors = [2*std for std in length_sensitivity_stds]
#length_sensitivity_errors_upper = []
#for i, mean in enumerate(length_sensitivity_means):
#    if (mean + 2*length_sensitivity_stds[i] > 1):
#        error = 1 - mean
#        length_sensitivity_errors_upper.append(error)
#    else:
#        length_sensitivity_errors_upper.append(length_sensitivity_stds[i]*2)
#
#length_sensitivity_errors_lower = [max((2*std, 0)) for std in 
#                                   length_sensitivity_stds]
#
# 
#
#length_specificity_means = [np.mean(specificities) for specificities in 
#                            trial_specificities]
#length_specificity_stds = [np.std(specificities) for specificities in 
#                           trial_specificities]
#length_specificity_errors = [2*std for std in length_specificity_stds]
#
#subs_sensitivity_means = [0.9672668236, 0.9978902053, 0.9892268481, 
#                          0.9612782473, 0.9686281217, 0.7647099155, 
#                          0.4649572138, 0.1941832427, 0.07847325855, 
#                          0.03383950359]
#subs_sensitivity_stds = [0.01820500844, 0.003751457372, 0.007053695757, 
#                         0.01663856124, 0.01458290775, 0.04421938972, 
#                         0.06239228437, 0.03726337989, 0.02210031014, 
#                         0.01430279011]
#
#subs_sensitivity_errors = [2*std for std in subs_sensitivity_stds]
#
#subs_specificity_means = [0.9722387901, 0.9781456601, 0.9865311132, 
#                          0.9939467173, 0.9666281629, 0.9973872353,	
#                          0.9968382777, 0.9973861561, 0.9973865158, 
#                          0.9973865139]
#subs_specificity_stds = [0.04951011097, 0.0567989133, 0.02330796448,	
#                         0.005852751996, 0.11714118, 0.001232437158, 
#                         0.002693404477,	0.000003467098229, 
#                         0.000003387384166, 0.000004087748032]
#
#subs_specificity_errors_upper = []
#for i, mean in enumerate(subs_specificity_means):
#    if (mean + 2*subs_specificity_stds[i] > 1):
#        error = 1 - mean
#        subs_specificity_errors_upper.append(error)
#    else:
#        subs_specificity_errors_upper.append(subs_specificity_stds[i]*2)
#
#subs_specificity_errors_lower = [max((2*std, 0)) for std in 
#                                 subs_specificity_stds]

def get_error_metrics(file, num_trials):
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    error_rate_metrics = []
    with open(file, 'r') as handle:
        for lst in chunked(handle.readlines(), num_trials):
            metrics = []
            for i, line in enumerate(lst):
                split_line = line.split(' ')
                #error_rate = int(
                #    split_line[0].split('/')[-1].split('.')[0]
                #)/100
                #if error_rate not in error_rates:
                #    error_rates.append(error_rate)

                metric = float(split_line[1])
                metrics.append(metric)
            error_rate_metrics.append(metrics)

    return error_rates, error_rate_metrics


def get_length_metrics(file, num_trials):

    seq_lengths = [1000, 2000, 3000, 4000, 5000, 10000, 20000, 40000, 80000, 
                   160000, 200000, 250000, 300000]
    seq_length_metrics = []
    with open(file, 'r') as handle:
        for lst in chunked(handle.readlines(), num_trials):
            metrics = []
            for i, line in enumerate(lst):
                split_line = line.split(' ')
                #seq_length = int(
                #    split_line[0].split('/')[-1].split('.')[0]
                #)
                #if seq_length not in seq_lengths:
                #    seq_lengths.append(seq_length)

                metric = float(split_line[1])
                metrics.append(metric)
            seq_length_metrics.append(metrics)

    return seq_lengths, seq_length_metrics


def get_read_metrics(file):
    metrics = []
    with open(file, 'r') as handle:
        for line in handle:
            metric = float(line.split(' ')[1].rstrip())
            metrics.append(metric)

    return metrics


def make_mammal_figure():
    pan_trog_lengths, pan_trog_sensitivities = get_length_metrics(
        'metrics/pan_trog_sensitivities_90CI.txt', num_trials=10)

    _, pan_trog_specificities = get_error_metrics(
        'metrics/pan_trog_specificities_90CI.txt', num_trials=10)

    gorilla_lengths, gorilla_sensitivities = get_length_metrics(
        'metrics/gorilla_sensitivities_03_unbound.txt', num_trials=100)

    _, gorilla_specificities = get_error_metrics(
        'metrics/gorilla_specificities_03_unbound.txt', num_trials=100)

    # Pan Trog (Chimp)
    pan_trog_sensitivity_means = [np.mean(sensitivities) for sensitivities in 
                                  pan_trog_sensitivities]
    pan_trog_sensitivity_stds = [np.std(sensitivities) for sensitivities in 
                                 pan_trog_sensitivities]

    pan_trog_sensitivity_errors_upper = []
    for i, mean in enumerate(pan_trog_sensitivity_means):
        if (mean + 2*pan_trog_sensitivity_stds[i] > 1):
            error = 1 - mean
            pan_trog_sensitivity_errors_upper.append(error)
        else:
            pan_trog_sensitivity_errors_upper.append(pan_trog_sensitivity_stds[i]*2)
    
    pan_trog_sensitivity_errors_lower = [max((2*std, 0)) for std in 
                                        pan_trog_sensitivity_stds]

    pan_trog_specificity_means = [np.mean(specificities) for specificities in 
                                  pan_trog_specificities]
    pan_trog_specificity_stds = [np.std(specificities) for specificities in 
                                 pan_trog_specificities]
    pan_trog_specificity_errors = [2*std for std in pan_trog_specificity_stds]

    # Gorilla
    gorilla_sensitivity_means = [np.mean(sensitivities) for sensitivities in 
                                  gorilla_sensitivities]
    gorilla_sensitivity_stds = [np.std(sensitivities) for sensitivities in 
                                 gorilla_sensitivities]

    gorilla_sensitivity_errors_upper = []
    for i, mean in enumerate(gorilla_sensitivity_means):
        if (mean + 2*gorilla_sensitivity_stds[i] > 1):
            error = 1 - mean
            gorilla_sensitivity_errors_upper.append(error)
        else:
            gorilla_sensitivity_errors_upper.append(gorilla_sensitivity_stds[i]*2)
    
    gorilla_sensitivity_errors_lower = [max((2*std, 0)) for std in 
                                        gorilla_sensitivity_stds]

    gorilla_specificity_means = [np.mean(specificities) for specificities in 
                                  gorilla_specificities]
    gorilla_specificity_stds = [np.std(specificities) for specificities in 
                                 gorilla_specificities]
    gorilla_specificity_errors = [2*std for std in gorilla_specificity_stds]


    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on mammalian genomes')

    # Pan Trog
    ax1.set_title('Chimpanzee')
    ax1.plot(pan_trog_lengths, pan_trog_sensitivity_means, label='sensitivity', 
             marker='o')
    ax1.errorbar(pan_trog_lengths, pan_trog_sensitivity_means, 
                 yerr=[pan_trog_sensitivity_errors_lower, 
                       pan_trog_sensitivity_errors_upper], 
                 fmt='-', color='blue')

    ax1.plot(pan_trog_lengths, pan_trog_specificity_means, label='specificity', 
             marker='o')
    ax1.errorbar(pan_trog_lengths, pan_trog_specificity_means, 
                 yerr=pan_trog_specificity_errors, fmt='-', color='orange')
    ax1.axvline(x=5000, linestyle='--', color='grey')
    ax1.text(5500, 0.2, '5kb query threshold', rotation=90)
    ax1.set_xscale('log')
    ax1.set_xlabel('Query length (kb)')
    #ax1.set_ylim(ymin=0)

    # Gorilla
    ax2.set_title('Gorilla')
    ax2.plot(gorilla_lengths, gorilla_sensitivity_means, label='sensitivity', 
             marker='o')
    ax2.errorbar(gorilla_lengths, gorilla_sensitivity_means, 
                 yerr=[gorilla_sensitivity_errors_lower, 
                       gorilla_sensitivity_errors_upper], 
                 fmt='-', color='blue')

    ax2.plot(gorilla_lengths, gorilla_specificity_means, label='specificity', 
             marker='o')
    ax2.errorbar(gorilla_lengths, gorilla_specificity_means, 
                 yerr=gorilla_specificity_errors, fmt='-', color='orange')
    ax2.text(5500, 0.2, '5kb query threshold', rotation=90)
    ax2.axvline(x=5000, linestyle='--', color='grey')
    ax2.set_xscale('log')
    ax2.set_xlabel('Query length (kb)')
    #ax2.legend()
    plt.show()
    #plt.savefig('figures/flashmap_mammal_accuracy_03.png')


def make_simulation_trials_figure():
    error_rates, error_sensitivities = get_error_metrics(
        'metrics/adaptive_error_error_sensitivity_90CI.txt', num_trials=100)

    _, error_specificities = get_error_metrics(
        'metrics/adaptive_error_error_specificity_90CI.txt', num_trials=100)

    seq_lengths, seq_length_sensitivities = get_length_metrics(
        'metrics/adaptive_error_length_sensitivity_90CI.txt', num_trials=100)

    _, seq_length_specificities = get_length_metrics(
        'metrics/adaptive_error_length_specificity_90CI.txt', num_trials=100)

    # Sub Rate
    subs_sensitivity_means = [np.mean(sensitivities) for sensitivities in 
                              error_sensitivities]
    subs_sensitivity_stds = [np.std(sensitivities) for sensitivities in 
                             error_sensitivities]
    subs_sensitivity_errors = [2*std for std in subs_sensitivity_stds]

    subs_sensitivity_errors_upper = []
    for i, mean in enumerate(subs_sensitivity_means):
        if (mean + 2*subs_sensitivity_stds[i] > 1):
            error = 1 - mean
            subs_sensitivity_errors_upper.append(error)
        else:
            subs_sensitivity_errors_upper.append(subs_sensitivity_stds[i]*2)
    
    subs_sensitivity_errors_lower = [max((2*std, 0)) for std in 
                                     subs_sensitivity_stds]


    subs_specificity_means = [np.mean(specificities) for specificities in 
                              error_specificities]
    subs_specificity_stds = [np.std(specificities) for specificities in 
                             error_specificities]
    subs_specificity_errors = [2*std for std in subs_specificity_stds]

    # Query Length
    length_sensitivity_means = [np.mean(sensitivities) for sensitivities in 
                                seq_length_sensitivities]
    length_sensitivity_stds = [np.std(sensitivities) for sensitivities in 
                               seq_length_sensitivities]
    length_sensitivity_errors = [2*std for std in length_sensitivity_stds]

    length_sensitivity_errors_upper = []
    for i, mean in enumerate(length_sensitivity_means):
        if (mean + 2*length_sensitivity_stds[i] > 1):
            error = 1 - mean
            length_sensitivity_errors_upper.append(error)
        else:
            length_sensitivity_errors_upper.append(length_sensitivity_stds[i]*2)
    
    length_sensitivity_errors_lower = [max((2*std, 0)) for std in 
                                     length_sensitivity_stds]

    length_specificity_means = [np.mean(specificities) for specificities in 
                                seq_length_specificities]
    length_specificity_stds = [np.std(specificities) for specificities in 
                               seq_length_specificities]
    length_specificity_errors = [2*std for std in length_specificity_stds]


    # Theoritical Curve
    error_rate_theory_sensitivities = [0.9999999999999992, 0.9999999973133248, 
                                       0.999999035776888, 0.9999633461915861, 
                                       0.9996118046476143, 0.9985398391125027, 
                                       0.9951159587019267, 0.9887248080755489, 
                                       0.9786179973528009, 0.9625853664976465]
    error_rate_theory_specificities = [1.0, 1.0, 1.0, 0.9999999999630376, 
                                       0.9999965025321311, 0.9881394058910561, 
                                       0.2188627586201921, 
                                       5.96189764223709e-14, 0.0, 0.0]
    query_size_theory_sensitivities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                       1.0, 1.0, 1.0, 1.0, 1.0]
    query_size_theory_specificities = [0.999999585864842, 0.9999999999999574, 
                                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                       1.0, 1.0, 1.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')
    fig.suptitle('Flashmap accuracy on simulated data')


    # Sub errors
    ax1.plot(error_rates, subs_sensitivity_means, label='sensitivity', 
             marker='o')
    # ax1.plot(error_rates, error_rate_theory_sensitivities,
    #          label='theoretical sensitivity')
    # ax1.plot(error_rates, error_rate_theory_specificities,
    #         label='theoretical specificity')
    ax1.errorbar(error_rates, subs_sensitivity_means, 
                 #yerr=subs_sensitivity_errors, 
                 yerr=[subs_sensitivity_errors_lower, 
                       subs_sensitivity_errors_upper], 
                 fmt='-', color='blue',)
    ax1.plot(error_rates, subs_specificity_means, label='specificity', 
             marker='o')
    ax1.errorbar(error_rates, subs_specificity_means, 
                 #yerr=[subs_specificity_errors_lower, subs_specificity_errors_upper], 
                 yerr=subs_specificity_errors,
                 fmt='-', color='orange')
    #ax1.axhline(y=0.95, linestyle='--', color='grey')
    ax1.axvline(x=0.05, linestyle='--', color='grey')
    ax1.text(0.045, 0.1, '0.05 substitution rate threshold', rotation=90)
    ax1.set_xlabel('Substitutions per site')
    #ax1.set_ylim(ymin=0.7)
    fig.legend(frameon=True)

    # Query Length
    ax2.plot(seq_lengths, length_sensitivity_means, label='sensitivity', 
             marker='o')
    # ax2.plot(seq_lengths, query_size_theory_sensitivities, 
    #          label='theoretical sensitivity')
    # ax2.plot(seq_lengths, query_size_theory_specificities,
    #          label='theoretical specificity')
    ax2.errorbar(seq_lengths, length_sensitivity_means, 
                 yerr=[length_sensitivity_errors_lower, length_sensitivity_errors_upper], 
                 #yerr=length_sensitivity_errors, 
                 fmt='-', color='blue')

    ax2.plot(seq_lengths, length_specificity_means, label='specificity', 
             marker='o')
    ax2.errorbar(seq_lengths, length_specificity_means, 
                 yerr=length_specificity_errors, fmt='-', color='orange')
    ax2.axvline(x=5000, linestyle='--', color='grey')
    ax2.text(5500, 0.3, '5kb query threshold', rotation=90)
    ax2.set_xscale('log')
    ax2.set_xlabel('Query length (kb)')
    plt.show()
    #plt.savefig('figures/flashmap_adaptive_error_accuracy_03_unbound.png')


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


#make_read_figure()
#make_runtime_figure()
make_simulation_trials_figure()
#make_mammal_figure()
#make_synth_figure()
