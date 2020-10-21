#! /usr/bin/env python

# Script to:
# - detect the build of the target and base datasets
# - update the base build to match the target build if need be



######################
# Importing packages #
######################

import argparse
from snps import SNPs
import io



#####################
# Parsing arguments #
#####################

parser = argparse.ArgumentParser(description='Detect target and base datasets build and update the base build to match the target if need be')
parser.add_argument('-t', '--input_target', help='Input BIM file (a combination of all BIM files, transformed into a 23andme-like format')
parser.add_argument('-b', '--input_base', help='Input base file, transformed into a 23andme-like format')

args = vars(parser.parse_args())

# Args to variable
input_target = args['input_target']
input_base = args['input_base']



###############################################################################
# Detect builds and update the base's build if it does not match the target's #
###############################################################################

target = SNPs(input_target, output_dir='.')
base = SNPs(input_base, output_dir='.')

if base.build != target.build:
    base.remap_snps(target.build)
    updated_base = base.save_snps("new_base_coordinates.txt", sep="\t", header=True)


