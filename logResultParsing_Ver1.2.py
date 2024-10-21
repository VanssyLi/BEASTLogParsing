import numpy as np
import arviz as az
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import argparse


def add_space(template, output):
    inputfile = open(template, 'r').readlines()
    write_file = open(output, 'w')

    for line in inputfile:
        write_file.write(line)
        if '# keywords: skygrid' in line:
            new_line = str('\n')
            write_file.write(new_line)
        else:
            continue
    write_file.close()
    return output


def read_log(filename):
    global log, parameter_names, taxa_name, burnin_percent
    add_space(filename, filename)
    # Read the log file and transpose the data
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                log = pd.read_table(f, index_col=0)

    # burn in
    burn_in = int(log.shape[0] * burnin_percent)
    log = log[burn_in:log.shape[0]]
    parameter_names = list(log.columns)
    # print(list(log.columns))

    for param in parameter_names:  # Select the date columns
        if 'age' and 'ND' in param:
            log = log[param]
            log.columns = param
            taxa_name = param
            # print(taxa_name)

    return log


def stat_summary(log, taxa_name):
    # Summarize the statistic data
    log_np = log.to_numpy()

    ess = az.ess(log_np, method='mean')
    estimated_means = np.mean(log_np, axis=0)
    stdev = np.std(log_np, axis=0)
    hpd_intervals = az.hdi(log_np, hdi_prob=0.95)

    res_dict = {}
    res_dict[taxa_name] = [estimated_means, stdev, hpd_intervals[0], hpd_intervals[1], ess]

    res_df = pd.DataFrame.from_dict(res_dict, orient='index')
    res_df.columns = ['mean', 'stdev', '95%lower', '95%upper', 'ess']

    return res_df


def plot_res(df, filename):
    fig, ax = plt.subplots(figsize=(10, 25))
    ax.errorbar(df['mean'],
                df.index, xerr=[df['mean'] - df['95%lower'], df['95%upper'] - df['mean']],
                fmt='|', color='black',
                ecolor='#939393', elinewidth=10)
    ax.set_xlim(0, 2E6)
    ax.invert_xaxis()
    ax.grid(True)
    ax.grid(axis='y', alpha=0)
    ax.set_xlabel('Estimated age (Ma)')

    plt.yticks(size=16)
    plt.yticks(size=15)

    plt.savefig(filename + '.png', dpi=300)
    # plt.show()


def load_save(path, filename):
    # read all file in a directory
    res_all = pd.DataFrame([])
    for p in Path(path).glob('*.log'):
        res_all = res_all.append(stat_summary(read_log(p), taxa_name))
    res_all.to_csv(filename + '.csv')

    # deal with index(taxa) name
    for index in res_all.index.values.tolist():
        taxa = str(index).split('_')[0].lstrip('age(').rstrip(')')
        res_all.rename(index={index: taxa}, inplace=True)
    res_all = res_all.sort_index(ascending=False)

    print(res_all)
    plot_res(res_all, filename)


'''path = '/Users/vanssyli/Master/SU/CPG/BEASTdata/dating_comparisons/DPVetal2021/date_all/log/'
filename = 'allResults.csv'
burnin_percent = 0.1


load_save(path, filename)'''



parser = argparse.ArgumentParser(prog='logResultParsing')

parser.add_argument('-p', '--path', type=str, help='The path for log files')
parser.add_argument('-f', '--filename', type=str, help='The output file name')
parser.add_argument('-b', '--burnin_percentage', type=int, help='Burn in first 10 percent (0.1 default) or 20 percent (0.2)')

args = parser.parse_args()

load_save(path=args.path,
          filename=args.filename)
