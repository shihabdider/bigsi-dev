import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

error_rates = []
sensitivities = []
with open('error_rate_sensitivities.txt', 'r') as handle:
    for line in handle:
        split_line = line.split(' ')
        error_rate = int(split_line[0][13:15])/100
        error_rates.append(error_rate)

        sensitivity = float(split_line[1])
        sensitivities.append(sensitivity)

specificities = []
with open('error_rate_specificities.txt', 'r') as handle:
    for line in handle:
        split_line = line.split(' ')
        specificity = float(split_line[1])
        specificities.append(specificity)

#plt.plot(error_rates, sensitivities, label='sensitivity', marker='o')
#plt.plot(error_rates, specificities, label='specificity', marker='o')
#plt.title('Flashmap Performance at Varying Error Rates')
#plt.xlabel('Error rate')
#plt.ylabel('Performance')
#plt.legend()
#plt.savefig('flashmap_performance_error_rate.png')


seq_lengths = []
sensitivities = []
with open('seq_length_sensitivities.txt', 'r') as handle:
    for line in handle:
        split_line = line.split(' ')
        seq_length = int(split_line[0][12:-13])
        print(seq_length)
        seq_lengths.append(seq_length)

        sensitivity = float(split_line[1])
        sensitivities.append(sensitivity)

specificities = []
with open('seq_length_specificities.txt', 'r') as handle:
    for line in handle:
        split_line = line.split(' ')
        specificity = float(split_line[1])
        specificities.append(specificity)

plt.plot(seq_lengths, sensitivities, label='sensitivity', marker='o')
plt.plot(seq_lengths, specificities, label='specificity', marker='o')
plt.title('Flashmap Performance at Varying Sequence Lengths')
plt.xlabel('Sequence Length (Kbp)')
plt.ylabel('Performance')
plt.legend()
plt.savefig('flashmap_performance_seq_length.png')
