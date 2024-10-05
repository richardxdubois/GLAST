from fermipy.gtanalysis import GTAnalysis
import argparse

# Command line arguments
parser = argparse.ArgumentParser(description='run phased analysis')

parser.add_argument('--config', default='config_0.yaml', help="config file")
parser.add_argument('--evfile', default='LSI_phase_0.txt', help="events file list")
parser.add_argument('--output', default='pickle_0', help="output pickle file")


args = parser.parse_args()

print("Entering lsi_phased.py using", args.config)

gta = GTAnalysis(args.config,logging={'verbosity': 3})

gta.setup(overwrite=True)
#gta.config["data"]["evfile"] = args.evfile

# Free Normalization of all Sources within 3 deg of ROI center
gta.free_sources(distance=3.0,pars='norm')

# Free all parameters of isotropic and galactic diffuse components
gta.free_source('galdiff')
gta.free_source('isodiff')
gta.free_source('4FGL J0240.5+6113')

fit_results = gta.fit()
print('Fit Quality: ',fit_results['fit_quality'])
print(gta.roi['4FGL J0240.5+6113'])

gta.write_roi(args.output,make_plots=True)
