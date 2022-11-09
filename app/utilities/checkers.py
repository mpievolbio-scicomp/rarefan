import glob
import sys
import os
import logging
import pandas
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def count_lines(fname):
    """ Count number of lines in file.

    :param fname: Name (path) of the file to be counted.
    :type  fname: str (path)

    """

    with open(fname, 'r') as fh:
        try:
            number_of_lines = len(fh.readlines())
        except IOError:
            logger.warning("%s is empty.", fname)
            return 0

    return number_of_lines


def parse_results(outdir,
                  reference_strain,
                  ):
    """ Parse output files from a given job and infer how the job has ended.

    :param outdir: The rarefan run output directory.
    :type  run_id: str (path)

    :param reference_strain: The reference strain
    :type  reference_strain: str

    :return: A dictionary of counts for number of rayts, number of overrepresented reps and number of repins.
    :rtype: dict

    """

    # Strip extension from reference strain filename.
    reference_strain = ".".join(reference_strain.split(".")[:-1])
    assoc_fname = os.path.join(outdir, "repin_rayt_association.txt.fas")
    overrep_fname = os.path.join(outdir, reference_strain + ".overrep")

    filenames = {"rayts": assoc_fname,
                 "nmers": overrep_fname
                 }
    results = {"counts": {"rayts": 0,
                          "nmers": 0,
                          "repins": {},
                          },
               "status": {"rayts": 0,
                          "nmers": 0,
                          "repins": 0,
                          }
               }

    for fname in "rayts", "nmers":
        try:
            results["counts"][fname] = count_lines(filenames[fname])
        except:
            results['status'][fname] = 1
            logger.warning("Could not count lines in RAYT file", fname)

    # RAYTs must be divided by 2 to get the correct number.
    results['counts']['rayts'] = results['counts']['rayts'] // 2

    # Parse presAbs files into pandas dataframes
    pres_abs_frames = dict()

    pres_abs_fnames = glob.glob(os.path.join(outdir, "presAbs_*.txt"))
    pres_abs_fnames = sorted(pres_abs_fnames)
    repin_checks = [0]*len(pres_abs_fnames)
    for i,pres_abs_fname in enumerate(pres_abs_fnames):
        try:
            pres_abs_frames[i] = pandas.read_csv(pres_abs_fname, sep="\t", index_col=0).to_dict(orient='index')
            repin_checks[i] = 0
        except:
            pres_abs_frames[i] = pandas.DataFrame()
            repin_checks[i] = 1
            logger.warning("Could not count lines in RAYT file", fname)

    results['counts']['repins'] = pres_abs_frames

    return results

if __name__ == "__main__":
    parsed = parse_results(sys.argv[1], 'Nmen_2594.fas')

    print(parsed)

