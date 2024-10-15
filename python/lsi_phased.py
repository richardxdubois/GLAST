from fermipy.gtanalysis import GTAnalysis
import argparse

# Command line arguments
parser = argparse.ArgumentParser(description='run phased analysis')

parser.add_argument('--config', default='config_0.yaml', help="config file")
parser.add_argument('--evfile', default='LSI_phase_0.txt', help="events file list")
parser.add_argument('--output', default='pickle_0', help="output pickle file")
parser.add_argument('--source', default='4FGL J1826.2-1450', help="source name")
parser.add_argument('--overwrite', action='store_true', help="overwrite files?")


args = parser.parse_args()

print("Entering lsi_phased.py using", args.config)

gta = GTAnalysis(args.config, logging={'verbosity': 3})

print(args.overwrite)
gta.setup(overwrite=False)

# Free Normalization of all Sources within 3 deg of ROI center
gta.free_sources(distance=3.0, pars='norm')

# Free all parameters of isotropic and galactic diffuse components
gta.free_source('galdiff')
gta.free_source('isodiff')
gta.free_source(args.source)

fit_results = gta.fit()
print('Fit Quality: ', fit_results['fit_quality'])
print(gta.roi[args.source])

sed = gta.sed(args.source)

gta.write_roi(args.output, make_plots=True)
