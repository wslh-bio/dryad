import argparse
import pandas as pd
import os

#this gets us the root dir of the project
base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('dryad_matrix_1',
    help='Path to validated snp_dists_matrix'
    )
parser.add_argument('dryad_matrix_2',
    help='Path to test spriggan_report.csv'
    )
args = parser.parse_args()

validated_df = pd.read_csv(os.path.abspath(args.dryad_matrix_1),sep=',')
to_test_df = pd.read_csv(os.path.abspath(args.dryad_matrix_2),sep=',')

if validated_df.equals(to_test_df):
    print("These dataframes match and pass validation.")
else:
    
    print("WARNING: Failed validation. These snp matrices do not match.")