import argparse
import pandas as pd
import os
import logging

### Setting up loggin structure
logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

### Setting up argparse arguments
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('dryad_matrix_1',
    help='Path to validated snp_dists_matrix'
    )
parser.add_argument('dryad_matrix_2',
    help='Path to test spriggan_report.csv'
    )
args = parser.parse_args()

### Getting df from matrix input
validated_df = pd.read_csv(os.path.abspath(args.dryad_matrix_1),sep='\t')
to_test_df = pd.read_csv(os.path.abspath(args.dryad_matrix_2),sep='\t')

### Determiniing if matrices are the same, if not, why?
if validated_df.equals(to_test_df):
    logging.info("These dataframes match and pass validation.")
else:
    logging.info("Failed validation. These snp matrices do not match!")
    try:
        diff = validated_df.compare(to_test_df)
        logging.info(f'The discrepancies are noted below: \n{diff}')
    except ValueError:
        logging.info("The snp matrices have different sizes! \nOnly compare matrices that have the same amount of row and column vectors.")
