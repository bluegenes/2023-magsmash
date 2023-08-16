import os
import sys
import argparse
import pandas as pd

# List of file paths
#file_paths = ['file1.csv', 'file2.csv', 'file3.csv']  # Replace with your file paths


def main(args):

    gather_info = []
    # Iterate through each file and read into pandas dataframe
    for file_path in args.gather_files:
        data = pd.read_csv(file_path)
        gather_info.append(data)

#    import pdb;pdb.set_trace()
    # now combine the dataframes
    combined_data = pd.concat(gather_info)#, ignore_index=True)

    # Pivot the data to create the summary matrix
    summary_matrix = combined_data.pivot(index='query_name', columns='name', values='f_unique_weighted')
    summary_matrix['total'] = summary_matrix.sum(axis=1)

    # Save the summary matrix to a file
    summary_matrix.reset_index(inplace=True)
    summary_matrix[['query_name', 'total']].to_csv(args.output_summary, index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument('--gather-files', nargs='+', required=True,
                        help='List of file paths to gather')
    p.add_argument('--output-summary', required=True,
                        help='Output file path for the summary csv')
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

# # write argparse to input files and parameter
# def cmdline(sys_args):
#     "Command line entry point w/argparse action."
#     p = argparse.ArgumentParser()
#     p.add_argument('--pocp-table', default='brady_pocp_table.tab', help='Original POCP table')
#     p.add_argument('--sourmash-compare-csv', default='output.pocp/brady.protein.sc5.compare.csv', help='Sourmash compare csv')
#     p.add_argument('--output-comparison-csv')
#     p.add_argument('--output-plot', default='plots/POCP-vs-cANI.sc5.png', help='Output plot')
#     args = p.parse_args()
#     return main(args)

# if __name__ == '__main__':
#     returncode = cmdline(sys.argv[1:])
#     sys.exit(returncode)
