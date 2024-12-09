from fermipy.gtanalysis import GTAnalysis
import argparse

# Command line arguments
parser = argparse.ArgumentParser(description='run phased analysis')

parser.add_argument('--config', default='config_0.yaml', help="config file")
parser.add_argument('--evfile', default='LSI_phase_0.txt', help="events file list")
parser.add_argument('--output', default='pickle_0', help="output pickle file")
parser.add_argument('--source', default='4FGL J1826.2-1450', help="source name")
parser.add_argument('--overwrite', action='store_true', help="overwrite files?")
parser.add_argument('--freeze', default=3., type=float, help="source freeze radius")
parser.add_argument('--lowE', action='store_true', help="low energy fit (PL)")
parser.add_argument('--gated', action='store_true', help="replace gated pulsar with log parabola")

args = parser.parse_args()

print("Entering lsi_phased.py using", args.config)

gta = GTAnalysis(args.config, logging={'verbosity': 3})

print(args.overwrite)
gta.setup(overwrite=args.overwrite)

gta.optimize()

model = gta.roi.get_source_by_name(args.source)
source_glon = model['glon']
source_glat = model['glat']

if args.lowE:
    print("switching ", args.source, "to PL model")

    gta.delete_source(args.source)

    # Add SourceA to the model
    gta.add_source(args.source, { 'glon': source_glon, 'glat': source_glat,
                    'SpectrumType' : 'PowerLaw', 'Index': 2.0,
                    'Scale': 1000, 'Prefactor': 1e-11,
                    'SpatialModel': 'PointSource'})

if args.gated:
    print("switching pulsar", args.source, "to LogParabola model")

    gta.delete_source(args.source)

    # Add SourceA to the model
    gta.add_source(args.source, { 'glon': source_glon, 'glat': source_glat,
                    'SpectrumType' : 'LogParabola', 'alpha': 2.0, 'beta': 1.e-4,
                    'Scale': 1000, 'norm': 1e-11,
                    'SpatialModel': 'PointSource'})

# Free Normalization of all Sources within 3 deg of ROI center
if args.freeze > 0.:
    gta.free_sources(distance=args.freeze, pars='norm')

# Free all parameters of isotropic and galactic diffuse components
gta.free_source('galdiff')
gta.free_source('isodiff')
gta.free_source(args.source)

fit_results = gta.fit()
print('Fit Quality: ', fit_results['fit_quality'])
print(gta.roi[args.source])

sed = gta.sed(args.source, make_plots=True, prefix=args.output)
print(sed["eflux"])

gta.write_roi(args.output, make_plots=True)
