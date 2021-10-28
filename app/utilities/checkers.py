import sys
import os
import logging
logger = logging.getLogger('routes')


def count_lines(fname):
    """ Count number of lines in file.

    :param fname: Name (path) of the file to be counted.
    :type  fname: str (path)

    """

    with open(fname, 'r') as fh:
        try:
            number_of_lines = len(fh.readlines())
        except IOError:
            logging.warning("%s is empty.", fname)
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

    assoc_fname = os.path.join(outdir, "repin_rayt_association.txt.fas")
    overrep_fname = os.path.join(outdir, reference_strain + ".overrep")

    filenames = {"rayts": assoc_fname,
                 "seeds": overrep_fname
                 }
    results = {"counts":
               {"rayts": 0,
                "seeds": 0,
                "repins": {},
                },
               "status": {"rayts": 0,
                          "seeds": 0,
                          "repins": 0 
                          }
               }

    for fname in "rayts", "seeds":
        try:
            results["counts"][fname] = count_lines(filenames[fname])
        except:
            results['status'][fname] = 1

    # RAYTs must be divided by 2 to get the correct number.
    results['counts']['rayts'] = results['counts']['rayts'] // 2

    # Aggregate "*_largestCluster" filenames.
    cluster_files = [os.path.join(outdir,
                                  f"{reference_strain}_{i}",
                                  f"{reference_strain}_{i}_largestCluster.nodes")
                     for i in range(5)]

    # Count repins for each type 0..5
    for i, cluster_file in enumerate(cluster_files):
        try:
            results['counts']['repins'][i] = count_lines(cluster_file)
        except:
            results['status']['repins'] = 1

    return results

