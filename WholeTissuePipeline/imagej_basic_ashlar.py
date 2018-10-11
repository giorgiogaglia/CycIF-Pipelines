# @File(label="Select a slide to process") filename
# @File(label="Select the output location", style="directory") output_dir
# @String(label="Experiment name (base name for output files)") experiment_name
# @Float(label="Flat field smoothing parameter (0 for automatic)", value=0.1) lambda_flat
# @Float(label="Dark field smoothing parameter (0 for automatic)", value=0.01) lambda_dark

# Takes a slide (or other multi-series BioFormats-compatible file set) and
# generates flat- and dark-field correction profile images with BaSiC. The
# output format is two multi-series TIFF files (one for flat and one for dark)
# which is the input format used by Ashlar.

# Invocation for running from the commandline:
#
# ImageJ --ij2 --headless --run imagej_basic_ashlar.py "filename='input.ext',output_dir='output',experiment_name='my_experiment'"

import sys
from ij import IJ, WindowManager
from ij.macro import Interpreter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import BaSiC_ as Basic


def main():

    Interpreter.batchMode = True

    if (lambda_flat == 0) ^ (lambda_dark == 0):
        print ("ERROR: Both of lambda_flat and lambda_dark must be zero,"
               " or both non-zero.")
        return
    lambda_estimate = "Automatic" if lambda_flat == 0 else "Manual"

    print "Loading images..."

    options = ImporterOptions()
    options.setId(str(filename))
    options.setOpenAllSeries(True)
    options.setConcatenate(True)
    options.setSplitChannels(True)
    imps = BF.openImagePlus(options)

    num_channels = len(imps)
    w = imps[0].getWidth()
    h = imps[0].getHeight()
    ff_imp = IJ.createImage("Flat-field", w, h, num_channels, 32);
    df_imp = IJ.createImage("Dark-field", w, h, num_channels, 32);

    basic = Basic()
    Basic_noOfSlices = Basic.getDeclaredField('noOfSlices')
    Basic_noOfSlices.setAccessible(True)

    for channel, imp in enumerate(imps):
        title = imp.getTitle()
        print "Processing:", title
        x, y, c, z, t = imp.getDimensions()
        assert z == 1 and c == 1
        imp.setDimensions(1, t, 1)

        WindowManager.setTempCurrentImage(imp)
        Basic_noOfSlices.setInt(basic, t)
        basic.exec(
            imp, None, None,
            "Estimate shading profiles", "Estimate both flat-field and dark-field",
            lambda_estimate, lambda_flat, lambda_dark,
            "Ignore", "Compute shading only"
        )
        ff_channel = WindowManager.getImage('Flat-field:' + title)
        ff_channel.copy()
        ff_imp.setSlice(channel + 1)
        ff_imp.paste()
        ff_channel.close()

        df_channel = WindowManager.getImage('Dark-field:' + title)
        df_channel.copy()
        df_imp.setSlice(channel + 1)
        df_imp.paste()
        df_channel.close()

        imp.close()

    # Setting the active slice back to 1 seems to fix an issue where
    # the last slice was empty in the saved TIFFs. Not sure why.
    ff_imp.setSlice(1)
    df_imp.setSlice(1)

    ff_filename = '%s/%s-ffp-basic.tif' % (output_dir, experiment_name)
    IJ.saveAsTiff(ff_imp, ff_filename)
    ff_imp.show()
    ff_imp.close()

    df_filename = '%s/%s-dfp-basic.tif' % (output_dir, experiment_name)
    IJ.saveAsTiff(df_imp, df_filename)
    df_imp.show()
    df_imp.close()

    print "Done!"


main()
