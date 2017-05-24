import yaml
import os.path
import shutil
import subprocess

if __name__ == '__main__':
    config = yaml.load(open('config.yml', 'r'))
    directories = config['directories']

    for d in os.listdir(directories['analyzed_mri']):
        dirname = os.path.relpath(d)

        # iterate through subjects
        # all subject folders start with s0
        if dirname.startswith('s0'):

            print("Current subject is " + dirname)

            transitions_file = os.path.join(
                    directories['raw_behavioral'],
                    dirname,
                    (dirname +'_hc_transitions.yml'))

            # if there's a transitions file for the subject, read it in
            if os.path.isfile(transitions_file):

                hc_transitions = yaml.load(open(transitions_file, 'r'))

                # loop across hemispheres (keys)
                for k in hc_transitions:
                    print("Current hemisphere is " + k)
                    output_dir = os.path.join(
                            directories['analyzed_mri'],
                            dirname,
                            'ROIs',
                            ('ashs_' + k))

                    if os.path.exists(output_dir):

                        # make combined ROIs: whole ROI, anterior, posterior

                        # combined 35/36
                        subprocess.call(["fslmaths", os.path.join(output_dir,'35.nii'), "-add", os.path.join(output_dir,'36.nii'), os.path.join(output_dir,'35_36.nii')])

                        # combined CA2, CA3, DG
                        subprocess.call(["fslmaths", os.path.join(output_dir,'CA2.nii'), "-add", os.path.join(output_dir,'CA3.nii'), "-add", os.path.join(output_dir,'DG.nii'), os.path.join(output_dir,'CA2_3_DG.nii')])
                        subprocess.call(["fslmaths", os.path.join(output_dir,'CA3.nii'), "-add", os.path.join(output_dir,'DG.nii'), os.path.join(output_dir,'CA3_DG.nii')])

                        # whole hippo
                        subprocess.call(["fslmaths", os.path.join(output_dir,'CA1.nii'), "-add", os.path.join(output_dir,'CA2.nii'), "-add", os.path.join(output_dir,'CA3.nii'), "-add", os.path.join(output_dir,'DG.nii'), "-add", os.path.join(output_dir,'subiculum.nii'), os.path.join(output_dir, 'whole_hippo.nii')])
