from fermipy.gtanalysis import GTAnalysis
import argparse
import re

# Command line arguments
parser = argparse.ArgumentParser(description='run phased analysis')

parser.add_argument('--config', default='config_0.yaml', help="config file")
parser.add_argument('--evfile', default='LSI_phase_0.txt', help="events file list")
parser.add_argument('--output', default='pickle_0', help="output pickle file")
parser.add_argument('--source', default='4FGL J1826.2-1450', help="source name")
parser.add_argument('--overwrite', action='store_true', help="overwrite files?")
parser.add_argument('--freeze', default=3., type=float, help="source freeze radius")
parser.add_argument('--lowE', action='store_true', help="low energy fit (PL)")
parser.add_argument('--cutoffPL', action='store_true', help="cutoff power law")
parser.add_argument('--gated', default='', help="replace named gated pulsar")
parser.add_argument('--pulsar_logP', action='store_true', help="set pulsar to log parabola")
parser.add_argument('--add_fixed', action='store_true', help="add fixed source at target location")

args = parser.parse_args()

print("Entering lsi_phased.py using", args.config)

gta = GTAnalysis(args.config, logging={'verbosity': 3})

print(args.overwrite)
gta.setup(overwrite=args.overwrite)

#gta.optimize()

model = gta.roi.get_source_by_name(args.source)
print(model)
source_glon = model['glon']
source_glat = model['glat']

if args.lowE:
    print("switching ", args.source, "to PL model")

    gta.delete_source(args.source)

    # Add Source back to the model
    gta.add_source(args.source, { 'glon': source_glon, 'glat': source_glat,
                    'SpectrumType' : 'PowerLaw', 'Index': 2.0,
                    'Scale': 1000, 'Prefactor': 1e-11,
                    'SpatialModel': 'PointSource'})

if args.cutoffPL:
    print("switching ", args.source, "to cutoff PL model")

    gta.delete_source(args.source)

    # Add Source back to the model

    gta.add_source(args.source, { 'glon': source_glon, 'glat': source_glat,
                    'SpectrumType' : 'PLSuperExpCutoff4', 'IndexS': -2.0,
                    'Scale': 1000, 'Prefactor': 1e-11, 'ExpfactorS': 0.4, 'Index2': 0.5,
                    'SpatialModel': 'PointSource'})


# Free Normalization of all Sources within specified deg of ROI center

gta.free_sources(distance=args.freeze, pars='norm')
print("freeing sources to ", args.freeze, "deg")

if args.add_fixed:
    fixed_source_name = args.source + '_fixed'
    gta.add_source(fixed_source_name, { 'glon': source_glon, 'glat': source_glat,
                    'SpectrumType' : 'LogParabola', 'alpha': 2.164,
                    'Scale': 1000, 'norm': 3.696e-11, 'beta': 0.1728, 'Eb':1178.,
                    'SpatialModel': 'PointSource'}, free=False)
    print("adding fixed source: ", fixed_source_name)



# Free norm parameters of isotropic and galactic diffuse components
gta.free_source('galdiff', pars='norm')
gta.free_source('isodiff', pars='norm')
gta.free_source(args.source)

if args.gated != '':

    model_p = gta.roi.get_source_by_name(args.gated)
    p_glon = model_p['glon']
    p_glat = model_p['glat']

    gta.delete_source(args.gated)

    if args.pulsar_logP:
        # Add Source back to the model
        gta.add_source(args.gated, {'glon': p_glon, 'glat': p_glat,
                                     'SpectrumType': 'LogParabola', 'alpha': 2.0, 'beta': 1.e-4,
                                     'Scale': 1000, 'norm': 1e-11,
                                     'SpatialModel': 'PointSource'})
        print("switching pulsar", args.gated, "to LogParabola model, and freeing")

        gta.free_source(args.gated)

free_sources = []
source_model = gta.roi.get_sources()
for m in source_model:
    # Check if the source has any free parameters
    if m.is_free:
        free_sources.append(m.name)

print("Currently free sources in the model:")
for source in free_sources:
    print(source)

fit_results = gta.fit()
print('Fit Quality: ', fit_results['fit_quality'])
print(gta.roi[args.source])

if args.pulsar_logP:
    print(gta.roi[args.gated])

if args.add_fixed:
    print(gta.roi[fixed_source_name])


sed = gta.sed(args.source, make_plots=True, prefix=args.output)
print(sed["eflux"])

map_name = args.output + "_" + re.sub(r'\s+', '_', args.source) + "_"
r_model = {'Index' : 2.0, 'SpatialModel' : 'Gaussian', 'SpatialWidth' : 0.3 }
resid_map = gta.residmap(map_name, make_plots=True, model=r_model)
gta.plotter.make_residmap_plots(resid_map, roi=None)

ts_map = gta.tsmap(map_name, make_plots=True, model=r_model)
gta.plotter.make_tsmap_plots(ts_map, roi=None)

gta.write_roi(args.output, make_plots=True)
