import numpy as np
from datetime import datetime
import subprocess
import argparse
import pickle
import os
import shutil
import yaml
from pathlib import Path

from prepare_gti import prepare_gti


class prepare_super_orbital():
    def __init__(self, input_yaml):

        self.input_yaml = input_yaml

        self.daysecs = 86400.

        with open(self.input_yaml, "r") as f:
            data = yaml.safe_load(f)

        self.do_orb_only = data['do_orb_only']
        self.do_super_orbital = data['do_super_orbital']

        self.orbit_period = data['orbit_period'] * self.daysecs
        self.super_orbital_period = data['super_orbit_period'] * self.daysecs
        self.zero_phase = data['zero_phase']

        self.config = data["config"]
        self.batch = data["batch"]

        self.num_orb_bins = data["num_orb_bins"]
        self.num_super_bins = data["num_super_bins"]
        self.top_dir = data["top_dir"]
        self.base_orb_dir = data["base_orb_phase"]
        self.base_super_dir = data["base_super_phase"]
        self.input_file = data["input_file"]
        self.current_dir = None
        self.update_configs_only = False
        try:
            self.update_configs_only = data['update_configs_only']
        except KeyError:
            pass

        print("Finished loading input yaml file")

        self.prep_gti = prepare_gti()

        with open(self.config, "r") as f_config:
            save_config = f_config.read()

        with open(self.batch, "r") as f_batch:
            save_batch = f_batch.read()
            for pickle_version in save_batch.split():
                if pickle_version.startswith("pickle"):
                    break

        yaml_name = Path(input_yaml).stem
        pickle_config = "prepare_" + yaml_name + "_" + pickle_version + "_config.pkl"
        pickle_config_data = [data, save_config, save_batch]

        with open(pickle_config, 'wb') as pickle_file:
            pickle.dump(pickle_config_data, pickle_file)

    def prepare_file(self, template_file_path, tgt_file_path, bin=""):

        with open(template_file_path, 'r') as file:
            filedata = file.read()

        filedata = filedata.replace("BIN", bin)
        filedata = filedata.replace("CWD", self.current_dir)

        with open(tgt_file_path, 'w') as file:
            file.write(filedata)

        return 0

    def prepare_fits_file(self, input_file, tgt_file, bin):

        shutil.copyfile(input_file, tgt_file)

        # update GTIs in tgt_file
        orb = self.orbit_period
        if not self.do_orb_only:
            orb = self.super_orbital_period

        rc = self.prep_gti.get_infile_stuff(infile=tgt_file)
        rc = self.prep_gti.extrapolate_t0(t0_MJD=self.zero_phase, orb=orb)

        if self.do_orb_only:
            num_bins = self.num_orb_bins
        else:
            num_bins = self.num_super_bins
        rc = self.prep_gti.make_new_gti(phase_bin=bin, num_bins=num_bins)
        rc = self.prep_gti.merge_gtis()

        return 0

    def create_directory(self, path=None):
        try:
            os.makedirs(path)
            print("Directory ", path, " created")
        except OSError:
            print("Directory ", path, " already exists")
            pass

        return

    def fill_dir(self, dir, dir_pattern, bin):
        sub_dir = os.path.join(dir, dir_pattern + bin)
        self.current_dir = sub_dir

        rc = self.create_directory(sub_dir)

        config_in = os.path.join(self.top_dir, self.config)
        config_out = os.path.join(sub_dir, "config_" + bin + ".yaml")
        rc = self.prepare_file(template_file_path=config_in, tgt_file_path=config_out, bin=bin)

        batch_in = os.path.join(self.top_dir, self.batch)
        batch_out = os.path.join(sub_dir, "fpy_slurm.txt")
        rc = self.prepare_file(template_file_path=batch_in, tgt_file_path=batch_out, bin=bin)

        if self.update_configs_only is False:
            input_file = os.path.join(dir, self.input_file + ".fits")
            tgt_file = os.path.join(sub_dir, "ft1_gti_updated-" + bin + ".fits")
            rc = self.prepare_fits_file(input_file=input_file, tgt_file=tgt_file, bin=bin)

        return 0

    def populate(self):

        for i in range(self.num_orb_bins):

            orb_dir = os.path.join(self.top_dir, self.base_orb_dir + str(i))
            print("working on", orb_dir)

            if self.do_orb_only:
                rc = self.fill_dir(self.top_dir, self.base_orb_dir, str(i))
                self.current_dir = orb_dir

            if self.do_super_orbital:

                for j in range(self.num_super_bins):

                    str_j = str(j)

                    super_dir = os.path.join(orb_dir, self.base_super_dir + str_j)
                    rc = self.fill_dir(orb_dir, self.base_super_dir, str(j))
                    self.current_dir = super_dir
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--app_config',
                        default="/Users/richarddubois/Code/GLAST/tmp/dbg/super_only_test/app_config.yaml",
                        help="overall app config file")

    args = parser.parse_args()

    super = prepare_super_orbital(input_yaml=args.app_config)

    super.populate()
