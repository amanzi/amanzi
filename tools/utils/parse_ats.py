import sys,os
import shutil
import h5py
import numpy as np
import argparse

def get_keys(dat, time_range=None):
    """Collect all time index keys, optionally only those in a given 2-tuple of start and end time."""
    keys = sorted(dat[next(iter(dat.keys()))].keys(), key=int)

    if time_range is not None:
        all_times = get_times(dat, keys)
        keys = [key for key,time in zip(keys,all_times) if time_range[0] <= time <= time_range[1]]
    return keys

def get_times(dat, keys=None):
    """Get the times in the file, optionally at a given list of time index keys."""
    a_field = next(iter(dat.keys()))

    if keys is None:
        keys = get_keys(dat)

    times = [dat[a_field][key].attrs['Time'] for key in keys]
    return times

def get_keys_and_times(dat, time_range=None):
    """Get the time index keys and times in a file, optionally only within a range of (start,end) times."""
    if time_range is not None:
        all_keys, all_times = get_keys_and_times(dat, None)
        keys = [key for key,time in zip(all_keys,all_times) if time_range[0] <= time <= time_range[1]]
        times = [time for time in all_times if time_range[0] <= time <= time_range[1]]        
    else:
        keys = get_keys(dat)
        times = get_times(dat, keys)

    return keys, times

    

def readATS(directory='.', base="visdump_data.h5", inds=None, time_range=None,
            timeunits='yr'):
    """Main access point to data."""
    dat = h5py.File(os.path.join(directory,base),'r')

    if inds is None:
        keys, times = get_keys_and_times(dat, time_range)
    else:
        keys = get_keys(dat, time_range)
        keys = [keys[ind] for ind in inds]
        times = get_times(dat, keys)

    times = np.array(times)
    if timeunits == 'd':
        times = np.array([t*365.25 for t in times])
    elif timeunits == 's':
        times = np.array([t*365.25*86400 for t in times])
    elif timeunits != 'yr':
        raise RuntimeError("Invalid time unit: %s"%timeunits)
    return keys, times, dat

def getSurfaceData(keys, dat, name):
    if not len(name.split(".")) == 3:
        name = name + ".cell.0"
    res = np.array([dat[name][key][0] for key in keys])
    if len(res.shape) == 2 and res.shape[1] == 1:
        res = res[:,0]
    return res

def subsetFile(directory=".", base="visdump_data.h5", outfile="my_visdump_data.h5", inds=None, interval=1, time_range=None, names=None, dry_run=False):
    """Read one file, write another"""
    null_time_range = False
    if (time_range[0] == 0.0) and (time_range[1] == 1.e99):
        null_time_range = True
    
    keys, times, dat = readATS(directory,base,inds,time_range)

    if interval > 1:
        keys = keys[::interval]
        times = times[::interval]

    if (interval == 1) and (null_time_range) and (inds == None) and dry_run:
        print("Available times (count = %d):"%len(times), times)
        return (None,None)
    elif dry_run:
        print("Matched times (count = %d):"%len(times), times)
        return (None, None)
        
    print("Transfering %d times to %s"%(len(times),outfile))
    
    if names is None:
        names = dat.keys()
    
    out = h5py.File(outfile)
    for name in names:
        grp = out.create_group(name)
        for key in keys:
            grp.create_dataset(key, data=dat[name][key][:])

    out.create_dataset('Times', data=np.array(times))
    out.close()

    return keys, times

    
def subsetXMFFile(out_directory, base="visdump_data.VisIt.xmf", inds=None, interval=1, time_range=None, names=None, dry_run=False):
    """Read one simulation set, write another."""
    import xml.etree.ElementTree as etree

    null_time_range = False
    if (time_range[0] == 0.0) and (time_range[1] == 1.e99):
        null_time_range = True

    xmf_in = etree.parse(base)
    files = list(xmf_in.getroot()[0][0])
    if interval > 1:
        files = files[::interval]

    if (interval == 1) and (null_time_range) and (inds == None) and dry_run:
        print("Available times (count = %d):"%len(times), times)
        return (None,None)
    elif dry_run:
        print("Matched times (count = %d):"%len(times), times)
        return (None, None)

    to_remove = [f for f in xmf_in.getroot()[0][0] if f not in files]
    for f in to_remove:
        xmf_in.getroot()[0][0].remove(f)

    # write the VisIt.xmf file
    xmf_in.write(os.path.join(out_directory,os.path.split(base)[-1]))

    # write the .N.xmf files
    in_dir_list = os.path.split(base)
    if len(in_dir_list) > 1:
        in_dir = os.path.join(*in_dir_list[:-1])
    else:
        in_dir = "."

    if names is None:
        for f in xmf_in.getroot()[0][0]:
            # simply copy over the filenames
            fname = f.get("href")
            in_f = os.path.join(in_dir, f.get("href"))

            fname_strip = os.path.split(fname)[-1]
            out_f = os.path.join(out_directory, fname_strip)
            assert not os.path.exists(out_f)
            shutil.copyfile(in_f, out_f)

    else:
        for f in xmf_in.getroot()[0][0]:
            # write new files with a subset of the data
            fname = f.get("href")
            in_f = os.path.join(in_dir, f.get("href"))

            fname_strip = os.path.split(fname)[-1]
            out_f = os.path.join(out_directory, fname_strip)
            try:
                assert not os.path.exists(out_f)
            except AssertionError:
                print("Exists:", out_f)
                raise AssertionError

            step_xmf = etree.parse(in_f)
            to_remove = [e for e in step_xmf.getroot()[0][0] if e.tag == "Attribute" and e.get("Name") not in names]
            for e in to_remove:
                step_xmf.getroot()[0][0].remove(e)

            step_xmf.write(out_f)
    
    # find the h5 file, and assume that all use the same h5
    an_xmf_in = xmf_in.getroot()[0][0][0].get("href")
    in_f = os.path.join(in_dir, f.get("href"))
    an_step_xmf = etree.parse(in_f)
    in_data = an_step_xmf.getroot()[0][0][-1]
    assert in_data.tag == "Attribute"
    fname_h5 = in_data[-1].text.strip("\n").strip().split(":")[0]
    in_h5 = os.path.join(in_dir, fname_h5)

    fname_h5_strip = os.path.split(fname_h5)[-1]
    out_h5 = os.path.join(out_directory, fname_h5_strip)
    assert not os.path.exists(out_h5)

    subsetFile(base=in_h5, outfile=out_h5, inds=inds, interval=interval, time_range=time_range, names=names, dry_run=False)

    # find the mesh file, copy it over
    mesh_name = an_step_xmf.getroot()[0][0][0][-1].text.strip("\n").strip().split(":")[0]
    in_mesh = os.path.join(in_dir, mesh_name)
    mesh_name_strip = os.path.split(mesh_name)[-1]
    out_mesh = os.path.join(out_directory, mesh_name_strip)
    assert not os.path.exists(out_mesh)
    shutil.copyfile(in_mesh, out_mesh)
    
    

        
