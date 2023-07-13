import os
import shutil

def get_data_filename(filename):
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(test_dir, "data")
    fname = os.path.join(data_dir, filename)
    if not os.path.exists(fname):
        raise OSError("File %s not found!" % fname)
    return fname

def copy_to_tmpdir(filename, tmpdir):
    basename = os.path.basename(filename)
    new_filename = os.path.join(tmpdir, basename)
    shutil.copyfile(filename, new_filename)
    return new_filename